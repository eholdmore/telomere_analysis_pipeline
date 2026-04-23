"""
hc_pre_vep.py
Collects filtered HaplotypeCaller TSVs, checks each variant against the
Mutect2 annotation cache, and emits:
  - hc_obs_known.parquet   — variants already annotated in the M2 cache
  - hc_obs_novel_waiting.parquet — variants that need VEP annotation
  - hc_novel_variants_to_annotate.vcf — VCF for VEP input

Usage:
    python hc_pre_vep.py \
        --tsv-glob "filtered_hc/*.tsv" \
        --m2-cache variant_classification_table_v3.csv \
        --output-known  hc_obs_known.parquet \
        --output-novel  hc_obs_novel_waiting.parquet \
        --output-vcf    hc_novel_variants_to_annotate.vcf
"""

import argparse
import glob
import os

import polars as pl


HC_SCHEMA = {
    "CHROM":  pl.Utf8,
    "POS":    pl.Int64,
    "REF":    pl.Utf8,
    "ALT":    pl.Utf8,
    "QD":     pl.Float64,
    "MQRS":   pl.Float64,
    "AD":     pl.Float64,
    "DP":     pl.Float64,
}

HC_COLUMNS = [
    "CHROM", "POS", "REF", "ALT", "RAW_FILTER", "QUAL",
    "QD", "MQRS", "RPRS", "AD", "DP", "SAMPLE_ID", "briefcase",
]


def load_hc_tsvs(tsv_glob: str) -> pl.DataFrame:
    files = sorted(glob.glob(tsv_glob))
    if not files:
        raise FileNotFoundError(f"No TSV files matched: {tsv_glob}")
    print(f"  Loading {len(files)} HC TSV files...")
    frames = []
    for fp in files:
        try:
            df = pl.read_csv(
                fp,
                separator="\t",
                has_header=False,
                null_values=["."],
                new_columns=HC_COLUMNS,
                schema_overrides=HC_SCHEMA,
            )
            frames.append(df)
        except Exception as e:
            print(f"  WARNING: could not read {fp}: {e}")
    return pl.concat(frames, how="diagonal")


def apply_gatk_filters(df: pl.DataFrame) -> pl.DataFrame:
    """Reproduce the HC quality-flagging logic from hc_pre-vep.py."""
    df = df.with_columns([
        pl.col("SAMPLE_ID").str.strip_chars().alias("SAMPLE_ID"),
        (
            pl.col("CHROM") + ":" +
            pl.col("POS").cast(pl.Utf8) + ":" +
            pl.col("REF") + ":" +
            pl.col("ALT")
        ).alias("VARIANT_ID"),
        (pl.col("AD") / pl.col("DP")).fill_nan(0.0).alias("VAF"),
    ])
    is_snv = (
        pl.col("REF").str.len_chars() == 1
    ) & (
        pl.col("ALT").str.len_chars() == 1
    )
    df = df.with_columns(
        pl.when(is_snv)
          .then(
              pl.when(
                  (pl.col("QD").fill_null(10) < 3) |
                  (pl.col("MQRS").fill_null(0) < -4)
              ).then(pl.lit("lowQual_SNV"))
               .otherwise(pl.lit("PASS"))
          )
          .otherwise(
              pl.when(
                  (pl.col("QD").fill_null(10) < 3) |
                  (pl.col("MQRS").fill_null(0) < -5)
              ).then(pl.lit("lowQual_Indel"))
               .otherwise(pl.lit("PASS"))
          )
          .alias("HC_FILTER")
    )
    return df


def load_m2_cache(m2_cache_path: str) -> pl.DataFrame:
    cols = ["CHROM", "POS", "REF", "ALT", "SYMBOL.PRIMARY",
            "Consequence.PRIMARY", "HGVSp.PRIMARY"]
    return (
        pl.read_csv(m2_cache_path, null_values=["NA", ""], infer_schema_length=10000)
          .select([c for c in cols if c in
                   pl.read_csv(m2_cache_path, n_rows=1, null_values=["NA", ""]).columns])
          .unique(subset=["CHROM", "POS", "REF", "ALT"])
          .with_columns([
              pl.col("CHROM").cast(pl.Utf8),
              pl.col("POS").cast(pl.Int64),
          ])
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv-glob",      required=True)
    ap.add_argument("--m2-cache",      required=True)
    ap.add_argument("--output-known",  required=True)
    ap.add_argument("--output-novel",  required=True)
    ap.add_argument("--output-vcf",    required=True)
    args = ap.parse_args()

    print("Loading HC TSVs...")
    df_hc = load_hc_tsvs(args.tsv_glob)
    df_hc = apply_gatk_filters(df_hc)
    print(f"  Total HC observations: {df_hc.height:,}")

    print("Loading Mutect2 cache...")
    m2_annotated = load_m2_cache(args.m2_cache)

    join_keys = ["CHROM", "POS", "REF", "ALT"]
    df_hc_known = df_hc.join(m2_annotated, on=join_keys, how="inner")
    df_hc_novel = df_hc.join(m2_annotated, on=join_keys, how="anti")
    print(f"  Known (in M2 cache):  {df_hc_known.height:,}")
    print(f"  Novel (needs VEP):    {df_hc_novel.height:,}")

    # Write parquets
    for df, path in [(df_hc_known, args.output_known), (df_hc_novel, args.output_novel)]:
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        df.write_parquet(path)

    # Write VCF
    vcf_prep = (
        df_hc_novel
          .select(["CHROM", "POS", "VARIANT_ID", "REF", "ALT"])
          .unique()
          .sort(["CHROM", "POS"])
          .with_columns([
              pl.lit(".").alias("QUAL"),
              pl.lit("PASS").alias("FILTER"),
              pl.lit(".").alias("INFO"),
          ])
    )
    print(f"Writing VCF for {vcf_prep.height:,} unique novel sites -> {args.output_vcf}")
    os.makedirs(os.path.dirname(os.path.abspath(args.output_vcf)), exist_ok=True)
    with open(args.output_vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        vcf_prep.select(
            ["CHROM", "POS", "VARIANT_ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        ).write_csv(f, separator="\t", include_header=False)

    print("hc_pre_vep done.")


if __name__ == "__main__":
    main()
