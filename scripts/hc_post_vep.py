"""
hc_post_vep.py
Merges VEP annotations into the HC observations, performs a full join
against the Mutect2 cache, attaches QC metrics, and writes the master
telomere call table.

Usage:
    python hc_post_vep.py \
        --known-parquet  hc_obs_known.parquet \
        --novel-parquet  hc_obs_novel_waiting.parquet \
        --vep-output     hc_novel_variants_annotated.txt \
        --m2-cache       variant_classification_table_v3.csv \
        --metrics-csv    hs_metrics_summary.csv \
        --output-parquet master_telomere_calls_final.parquet \
        --output-csv     master_telomere_calls_final.csv
"""

import argparse
import io
import os
import subprocess

import polars as pl


VEP_TO_PRIMARY = {
    "SYMBOL":             "SYMBOL.PRIMARY",
    "Consequence":        "Consequence.PRIMARY",
    "HGVSp":              "HGVSp.PRIMARY",
    "HGVSc":              "HGVSc.PRIMARY",
    "Feature":            "Feature.PRIMARY",
    "CDS_position":       "CDS_position.PRIMARY",
    "Protein_position":   "Protein_position.PRIMARY",
    "Codons":             "Codons.PRIMARY",
    "IMPACT":             "IMPACT.PRIMARY",
    "STRAND":             "STRAND.PRIMARY",
    "HGNC_ID":            "HGNC_ID.PRIMARY",
    "EXON":               "EXON.PRIMARY",
    "Existing_variation": "Existing_variation.PRIMARY",
}

TARGET_GENES = ["TERT", "TERC", "POT1", "DKC1", "VEP_REJECTED"]

METRICS_TO_DROP = [
    "MEAN_TARGET_COVERAGE", "PCT_TARGET_BASES_20X", "PCT_EXC_DUPE",
    "PCT_TARGET_BASES_10X", "PCT_TARGET_BASES_30X",
]


def standardize_ids(df: pl.DataFrame) -> pl.DataFrame:
    if "CHROM" in df.columns:
        df = df.with_columns(pl.col("CHROM").cast(pl.Utf8))
    if "POS" in df.columns:
        df = df.with_columns(pl.col("POS").cast(pl.Int64))
    if "SAMPLE_ID" in df.columns:
        df = df.with_columns(
            pl.col("SAMPLE_ID").cast(pl.Utf8).str.strip_chars()
              .str.replace(r"_full_targets.*", "")
              .str.replace(r"_v[0-9].*", "")
              .alias("SAMPLE_ID")
        )
    return df


def parse_vep(vep_output_path: str) -> pl.DataFrame:
    cmd = f"grep -v '##' {vep_output_path} | sed 's/^#//'"
    cleaned = subprocess.check_output(cmd, shell=True).decode("utf-8")
    df_vep = pl.read_csv(
        io.StringIO(cleaned),
        separator="\t",
        null_values=["-"],
        infer_schema_length=0,
    )
    rename_map = {k: v for k, v in VEP_TO_PRIMARY.items() if k in df_vep.columns}
    return (
        df_vep
          .rename(rename_map)
          .rename({"Uploaded_variation": "VARIANT_ID"})
          .unique(subset=["VARIANT_ID"])
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--known-parquet",  required=True)
    ap.add_argument("--novel-parquet",  required=True)
    ap.add_argument("--vep-output",     required=True)
    ap.add_argument("--m2-cache",       required=True)
    ap.add_argument("--metrics-csv",    required=True)
    ap.add_argument("--output-parquet", required=True)
    ap.add_argument("--output-csv",     required=True)
    args = ap.parse_args()

    # ── Load ──────────────────────────────────────────────────────────────────
    print("Loading HC observations...")
    df_known = pl.read_parquet(args.known_parquet)
    df_novel = pl.read_parquet(args.novel_parquet)

    print("Parsing VEP output...")
    df_vep = parse_vep(args.vep_output)

    print("Loading Mutect2 cache...")
    df_m2 = (
        pl.read_csv(args.m2_cache, null_values=["NA", ""], infer_schema_length=10000)
          .rename({"sample_id": "SAMPLE_ID"} if "sample_id" in
                  pl.read_csv(args.m2_cache, n_rows=1).columns else {})
          .drop([c for c in METRICS_TO_DROP
                 if c in pl.read_csv(args.m2_cache, n_rows=1).columns])
    )

    # ── Merge VEP into novel observations ────────────────────────────────────
    df_novel_annotated = df_novel.join(df_vep, on="VARIANT_ID", how="left")
    df_hc_full = pl.concat([df_known, df_novel_annotated], how="diagonal")
    df_hc_full = standardize_ids(df_hc_full)
    df_m2      = standardize_ids(df_m2)

    df_hc_full = df_hc_full.with_columns(pl.lit(True).alias("HC_PRESENT"))

    # ── Full join against M2 cache ───────────────────────────────────────────
    print("Performing master join...")
    join_keys = ["CHROM", "POS", "REF", "ALT", "SAMPLE_ID"]
    df_final = df_hc_full.join(df_m2, on=join_keys, how="full", suffix="_m2")

    # Coalesce overlapping columns
    all_base_cols = {c.removesuffix("_m2") for c in df_final.columns}
    merge_exprs = []
    for col in all_base_cols:
        m2_col = f"{col}_m2"
        if col in df_final.columns and m2_col in df_final.columns:
            merge_exprs.append(pl.coalesce([pl.col(col), pl.col(m2_col)]).alias(col))

    if merge_exprs:
        df_final = df_final.with_columns(merge_exprs)

    df_final = df_final.drop([c for c in df_final.columns if c.endswith("_m2")])

    # ── Metric alignment ─────────────────────────────────────────────────────
    df_final = df_final.with_columns([
        pl.coalesce([
            pl.col("AD"),
            pl.col("T_ALT_COUNT").cast(pl.Float64, strict=False)
              if "T_ALT_COUNT" in df_final.columns else pl.lit(None, dtype=pl.Float64),
        ]).alias("AD"),
        pl.coalesce([
            pl.col("DP"),
            pl.col("T_TOTAL_COUNT").cast(pl.Float64, strict=False)
              if "T_TOTAL_COUNT" in df_final.columns else pl.lit(None, dtype=pl.Float64),
        ]).alias("DP"),
    ])

    # ── Attach HS metrics ────────────────────────────────────────────────────
    if os.path.exists(args.metrics_csv):
        df_hs = standardize_ids(
            pl.read_csv(args.metrics_csv).rename(
                {"Sample_ID": "SAMPLE_ID"}
                if "Sample_ID" in pl.read_csv(args.metrics_csv, n_rows=1).columns
                else {}
            )
        )
        df_final = df_final.join(df_hs, on="SAMPLE_ID", how="left")

    # ── Final labelling and gene filter ─────────────────────────────────────
    if "briefcase" in df_final.columns:
        df_final = df_final.rename({"briefcase": "DISEASE_GROUP"})

    df_final = df_final.with_columns([
        pl.col("SYMBOL.PRIMARY").fill_null("VEP_REJECTED"),
        pl.when(
            pl.col("HC_PRESENT").is_not_null() &
            pl.col("T_ALT_COUNT").is_not_null()
              if "T_ALT_COUNT" in df_final.columns
              else pl.col("HC_PRESENT").is_not_null()
        ).then(pl.lit("Both"))
          .when(pl.col("HC_PRESENT").is_not_null()).then(pl.lit("HaplotypeCaller"))
          .otherwise(pl.lit("Mutect2"))
          .alias("CALLER"),
    ]).filter(pl.col("SYMBOL.PRIMARY").is_in(TARGET_GENES))

    df_out = df_final.drop(
        [c for c in ["HC_PRESENT"] if c in df_final.columns]
    )

    # ── Save ─────────────────────────────────────────────────────────────────
    for path in [args.output_parquet, args.output_csv]:
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    df_out.write_parquet(args.output_parquet)
    df_out.write_csv(args.output_csv)
    print(f"Done — {df_out.height:,} variants written to {args.output_parquet}")


if __name__ == "__main__":
    main()
