"""
get_cn_info.py
Fetches per-sample gene copy number (from RobustCNV) and the germline
heterozygous SNP anchor VAF (from Mutect2), then applies the ploidy-corrected
VAF classifier to label each variant as Germline / Somatic / Ambiguous.

Usage:
    python get_cn_info.py \
        --input-parquet  master_telomere_calls_final.parquet \
        --paths-csv      file_paths_mapping.csv \
        --output-parquet master_telomere_calls_classified.parquet \
        --output-csv     master_telomere_calls_classified.csv
"""

import argparse
import os

import polars as pl


CLASSIFICATION_THRESHOLDS = [0.20, 0.25, 0.30]


# ── Per-file helpers ──────────────────────────────────────────────────────────

def fetch_ct_from_path(path: str, gene_symbol: str) -> float:
    """Return copy number (Ct) for a gene from RobustCNV output."""
    if path == "NOT_FOUND" or not path or not os.path.exists(path):
        return 2.0
    try:
        cnv = pl.read_csv(path, separator="\t", quote_char='"', infer_schema_length=0)
        cnv.columns = [c.replace('"', "").strip() for c in cnv.columns]
        gene_info = (
            cnv.with_columns(
                pl.col("Gene").str.replace_all('"', "").str.strip_chars().alias("Gene_Clean")
            ).filter(pl.col("Gene_Clean") == gene_symbol)
        )
        if gene_info.is_empty():
            return 2.0
        return 2 * (2 ** float(gene_info["GeneSegmentMedian"][0]))
    except Exception:
        return 2.0


def fetch_anchor_from_path(path: str) -> float:
    """
    Estimate the tumour-fraction anchor VAF from heterozygous germline SNPs
    in the Mutect2 somatic file.  Falls back to 0.48 (diploid germline).
    """
    if path == "NOT_FOUND" or not path or not os.path.exists(path):
        return 0.48
    try:
        mutect = pl.read_csv(path, separator="\t", infer_schema_length=1000)
        mutect.columns = [c.upper().replace('"', "").strip() for c in mutect.columns]

        if "VAF" not in mutect.columns:
            return 0.48

        mutect = mutect.with_columns([
            pl.col("VAF").cast(pl.Float64, strict=False),
            pl.col("GNOMAD_NON_CANCER_AF_POPMAX").cast(pl.Float64, strict=False)
              if "GNOMAD_NON_CANCER_AF_POPMAX" in mutect.columns
              else pl.lit(None, dtype=pl.Float64).alias("GNOMAD_NON_CANCER_AF_POPMAX"),
        ])

        anchors = mutect.filter(
            (pl.col("GNOMAD_NON_CANCER_AF_POPMAX") > 0.0001) &
            (pl.col("VAF").is_between(0.30, 0.70))
        )
        if anchors.height > 0:
            return anchors["VAF"].median()
        return 0.48
    except Exception:
        return 0.48


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-parquet",  required=True)
    ap.add_argument("--paths-csv",      required=True)
    ap.add_argument("--output-parquet", required=True)
    ap.add_argument("--output-csv",     required=True)
    args = ap.parse_args()

    print("Loading data and path mapping...")
    df       = pl.read_parquet(args.input_parquet)
    path_map = pl.read_csv(args.paths_csv)

    # Join paths into the main dataframe
    df = df.join(
        path_map.select(["cbio_sample_id", "cnv_path", "mutect_path"]),
        left_on="SAMPLE_ID",
        right_on="cbio_sample_id",
        how="left",
    )
    for col in ["cnv_path", "mutect_path"]:
        if col not in df.columns:
            df = df.with_columns(pl.lit("NOT_FOUND").alias(col))

    print(f"Extracting Ct and anchor VAF for {df.height:,} observations...")
    df = df.with_columns([
        pl.struct(["cnv_path", "SYMBOL.PRIMARY"]).map_elements(
            lambda x: fetch_ct_from_path(x["cnv_path"], x["SYMBOL.PRIMARY"]),
            return_dtype=pl.Float64,
        ).alias("Ct"),
        pl.col("mutect_path").map_elements(
            fetch_anchor_from_path,
            return_dtype=pl.Float64,
        ).alias("sample_anchor"),
    ])

    # p_c_vaf: expected VAF for a diploid copy at the tumour purity
    df = df.with_columns(
        (pl.col("sample_anchor") / (pl.col("Ct") / 2)).alias("p_c_vaf")
    )

    print("Applying classification scenarios...")
    for t in CLASSIFICATION_THRESHOLDS:
        df = df.with_columns(
            pl.when(pl.col("VAF") > 0.85)
              .then(pl.lit("Germline"))
              .when((pl.col("VAF") - pl.col("p_c_vaf")).abs() <= 0.12)
              .then(pl.lit("Germline"))
              .when(
                  ((pl.col("VAF") - pl.col("p_c_vaf")).abs() > t) &
                  (pl.col("VAF") < pl.col("sample_anchor"))
              ).then(pl.lit("Somatic"))
              .otherwise(pl.lit("Ambiguous"))
              .alias(f"Classification_{t}")
        )

    df = df.with_columns(
        ((pl.col("VAF") * pl.col("Ct")) / (2 * pl.col("sample_anchor")))
        .alias("Multiplicity")
    )

    # Remove helper columns before saving
    df_out = df.drop([c for c in ["cnv_path", "mutect_path", "p_c_vaf"] if c in df.columns])

    for path in [args.output_parquet, args.output_csv]:
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    df_out.write_parquet(args.output_parquet)
    df_out.write_csv(args.output_csv)
    print(f"Done — {df_out.height:,} variants written to {args.output_parquet}")


if __name__ == "__main__":
    main()
