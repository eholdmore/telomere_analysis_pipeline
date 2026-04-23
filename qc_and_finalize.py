"""
qc_and_finalize.py
Applies sample- and variant-level QC flags, hard filters, generates the
final output table, five publication-quality figures, and two summary CSVs.

Outputs (all written to --output-dir):
    final_variant_table.csv / .parquet
    qc_summary.csv
    variant_summary.csv
    figures/vaf_by_classification.png
    figures/cn_distribution.png
    figures/anchor_vaf.png
    figures/coverage_vs_20x.png
    figures/qd_vs_qual.png

Usage:
    python qc_and_finalize.py \
        --input-parquet master_telomere_calls_classified.parquet \
        --output-dir    output/
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import polars as pl
import pandas as pd


# ── QC thresholds (mirrors qc_metrics.R) ─────────────────────────────────────
MIN_MEAN_COVG  = 100.0
MIN_20X_PCT    = 0.95
MAX_DUPE_PCT   = 0.80
MIN_HC_QD      = 3.0
MIN_HC_MQRS    = -4.0
MIN_HC_RPRS    = -3.0
MIN_ALT_READS  = 9          # hard filter: ALT count > 9


CLASSIFICATION_COL = "Classification_0.25"   # primary classification threshold


# ── Colour palette ────────────────────────────────────────────────────────────
PALETTE = {
    "Germline":   "#2c7bb6",
    "Somatic":    "#d7191c",
    "Ambiguous":  "#fdae61",
    "PASS":       "#1a9641",
    "WARN":       "#fdae61",
    "FAIL":       "#d7191c",
}


# ── Helpers ───────────────────────────────────────────────────────────────────

def to_float(series: pl.Series) -> pl.Series:
    return series.cast(pl.Float64, strict=False)


def safe_col(df: pl.DataFrame, col: str, default=None) -> pl.Series:
    if col in df.columns:
        return df[col]
    if default is not None:
        return pl.Series([default] * df.height)
    raise KeyError(f"Column not found: {col}")


# ── Sample-level QC ───────────────────────────────────────────────────────────

def add_sample_qc_flags(df: pl.DataFrame) -> pl.DataFrame:
    has_covg  = "MEAN_TARGET_COVERAGE" in df.columns
    has_20x   = "PCT_TARGET_BASES_20X"  in df.columns
    has_dupe  = "PCT_EXC_DUPE"          in df.columns

    exprs = []
    flag_cols = []

    if has_covg:
        exprs.append(
            (pl.col("MEAN_TARGET_COVERAGE").cast(pl.Float64, strict=False) < MIN_MEAN_COVG)
            .alias("FLAG_LOW_COVG")
        )
        flag_cols.append("FLAG_LOW_COVG")
    if has_20x:
        exprs.append(
            (pl.col("PCT_TARGET_BASES_20X").cast(pl.Float64, strict=False) < MIN_20X_PCT)
            .alias("FLAG_LOW_20X")
        )
        flag_cols.append("FLAG_LOW_20X")
    if has_dupe:
        exprs.append(
            (pl.col("PCT_EXC_DUPE").cast(pl.Float64, strict=False) > MAX_DUPE_PCT)
            .alias("FLAG_HIGH_DUPE")
        )
        flag_cols.append("FLAG_HIGH_DUPE")

    if exprs:
        df = df.with_columns(exprs)

    # NUM_SAMPLE_FLAGS and SAMPLE_QC_STATUS
    if flag_cols:
        df = df.with_columns(
            sum(pl.col(c).cast(pl.Int32) for c in flag_cols).alias("NUM_SAMPLE_FLAGS")
        )
        df = df.with_columns(
            pl.when(pl.col("NUM_SAMPLE_FLAGS") >= 2).then(pl.lit("FAIL"))
              .when(pl.col("NUM_SAMPLE_FLAGS") == 1).then(pl.lit("WARN"))
              .otherwise(pl.lit("PASS"))
              .alias("SAMPLE_QC_STATUS")
        )
    else:
        df = df.with_columns([
            pl.lit(0).alias("NUM_SAMPLE_FLAGS"),
            pl.lit("PASS").alias("SAMPLE_QC_STATUS"),
        ])

    return df


# ── Variant-level QC ─────────────────────────────────────────────────────────

def add_variant_qc_flags(df: pl.DataFrame) -> pl.DataFrame:
    exprs = []
    flag_cols = []

    if "QD" in df.columns:
        exprs.append(
            (pl.col("QD").cast(pl.Float64, strict=False) < MIN_HC_QD)
            .fill_null(False).alias("FLAG_HC_LOW_QD")
        )
        flag_cols.append("FLAG_HC_LOW_QD")

    if "QUAL" in df.columns:
        exprs.append(
            (pl.col("QUAL").cast(pl.Float64, strict=False) < 30)
            .fill_null(False).alias("FLAG_HC_LOW_QUAL")
        )
        flag_cols.append("FLAG_HC_LOW_QUAL")

    if "MQRS" in df.columns:
        exprs.append(
            (pl.col("MQRS").cast(pl.Float64, strict=False) < MIN_HC_MQRS)
            .fill_null(False).alias("FLAG_HC_MQRS")
        )
        flag_cols.append("FLAG_HC_MQRS")

    if "RPRS" in df.columns:
        exprs.append(
            (pl.col("RPRS").cast(pl.Float64, strict=False) < MIN_HC_RPRS)
            .fill_null(False).alias("FLAG_HC_RPRS")
        )
        flag_cols.append("FLAG_HC_RPRS")

    if exprs:
        df = df.with_columns(exprs)

    if flag_cols:
        df = df.with_columns(
            sum(pl.col(c).cast(pl.Int32) for c in flag_cols).alias("NUM_FLAGS")
        )
        df = df.with_columns(
            pl.when(pl.col("NUM_FLAGS") == 0)
              .then(pl.lit("HighConf"))
              .otherwise(pl.lit("LowConf"))
              .alias("HC_VARIANT_QUALITY")
        )
    else:
        df = df.with_columns([
            pl.lit(0).alias("NUM_FLAGS"),
            pl.lit("HighConf").alias("HC_VARIANT_QUALITY"),
        ])

    # Germline PON / rescue flags (placeholder logic — adapt to your PON)
    if "GERMLINE_PON" not in df.columns:
        df = df.with_columns(pl.lit("FALSE").alias("GERMLINE_PON"))
    if "GERMLINE_RESCUE" not in df.columns:
        df = df.with_columns(pl.lit("FALSE").alias("GERMLINE_RESCUE"))

    return df


# ── Hard filter ───────────────────────────────────────────────────────────────

def apply_hard_filters(df: pl.DataFrame) -> pl.DataFrame:
    """
    Keep rows where:
      SAMPLE_QC_STATUS == PASS
      HC_VARIANT_QUALITY == HighConf
      GERMLINE_PON == FALSE
      AD > MIN_ALT_READS
    """
    mask = pl.lit(True)
    if "SAMPLE_QC_STATUS" in df.columns:
        mask = mask & (pl.col("SAMPLE_QC_STATUS") == "PASS")
    if "HC_VARIANT_QUALITY" in df.columns:
        mask = mask & (pl.col("HC_VARIANT_QUALITY") == "HighConf")
    if "GERMLINE_PON" in df.columns:
        mask = mask & (pl.col("GERMLINE_PON") == "FALSE")
    if "AD" in df.columns:
        mask = mask & (pl.col("AD").cast(pl.Float64, strict=False) > MIN_ALT_READS)
    return df.filter(mask)


# ── Figures ───────────────────────────────────────────────────────────────────

def fig_vaf_by_classification(df: pd.DataFrame, out_path: str):
    if CLASSIFICATION_COL not in df.columns or "VAF" not in df.columns:
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    groups = sorted(df[CLASSIFICATION_COL].dropna().unique())
    for label in groups:
        sub = df.loc[df[CLASSIFICATION_COL] == label, "VAF"].dropna()
        ax.hist(sub, bins=50, alpha=0.6, label=label,
                color=PALETTE.get(label, "grey"))
    ax.set_xlabel("Variant Allele Fraction (VAF)")
    ax.set_ylabel("Count")
    ax.set_title("VAF Distribution by Classification")
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def fig_cn_distribution(df: pd.DataFrame, out_path: str):
    if "Ct" not in df.columns:
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    ct = df["Ct"].dropna()
    ax.hist(ct, bins=40, color="#4393c3", edgecolor="white", linewidth=0.5)
    ax.axvline(2.0, color="red", linestyle="--", linewidth=1.2, label="Diploid (Ct=2)")
    ax.set_xlabel("Gene Copy Number (Ct)")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of Gene Copy Number")
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def fig_anchor_vaf(df: pd.DataFrame, out_path: str):
    if "sample_anchor" not in df.columns:
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    vals = df["sample_anchor"].dropna()
    ax.hist(vals, bins=40, color="#74add1", edgecolor="white", linewidth=0.5)
    ax.axvline(0.48, color="red", linestyle="--", linewidth=1.2, label="Default (0.48)")
    ax.set_xlabel("Sample Anchor VAF (median het-SNP VAF)")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of Sample Anchor VAF")
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def fig_coverage_vs_20x(df: pd.DataFrame, out_path: str):
    if "MEAN_TARGET_COVERAGE" not in df.columns or "PCT_TARGET_BASES_20X" not in df.columns:
        return
    sub = (
        df[["SAMPLE_ID", "MEAN_TARGET_COVERAGE", "PCT_TARGET_BASES_20X", "SAMPLE_QC_STATUS"]]
        .drop_duplicates("SAMPLE_ID")
        .dropna(subset=["MEAN_TARGET_COVERAGE", "PCT_TARGET_BASES_20X"])
    )
    status_col = "SAMPLE_QC_STATUS" if "SAMPLE_QC_STATUS" in sub.columns else None

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axhspan(MIN_20X_PCT, 1.0, xmin=MIN_MEAN_COVG / (sub["MEAN_TARGET_COVERAGE"].max() * 1.1),
               alpha=0.08, color="green", label="Safe zone")
    ax.axvline(MIN_MEAN_COVG, linestyle=":", color="grey", linewidth=1)
    ax.axhline(MIN_20X_PCT,   linestyle=":", color="grey", linewidth=1)

    if status_col:
        for status, color in PALETTE.items():
            mask = sub[status_col] == status
            ax.scatter(sub.loc[mask, "MEAN_TARGET_COVERAGE"],
                       sub.loc[mask, "PCT_TARGET_BASES_20X"],
                       c=color, label=status, alpha=0.7, s=20)
    else:
        ax.scatter(sub["MEAN_TARGET_COVERAGE"], sub["PCT_TARGET_BASES_20X"],
                   c="#4393c3", alpha=0.5, s=20)

    ax.set_xlabel("Mean Target Coverage (×)")
    ax.set_ylabel("% Bases at 20×")
    ax.set_title("Sample QC: Coverage vs. 20× Breadth")
    ax.legend(markerscale=2)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def fig_qd_vs_qual(df: pd.DataFrame, out_path: str):
    if "QD" not in df.columns or "HC_VARIANT_QUALITY" not in df.columns:
        return
    sub = df[["QD", "HC_VARIANT_QUALITY"]].dropna()
    groups = ["HighConf", "LowConf"]
    fig, ax = plt.subplots(figsize=(7, 5))
    for label, color in [("HighConf", PALETTE["PASS"]), ("LowConf", PALETTE["FAIL"])]:
        vals = sub.loc[sub["HC_VARIANT_QUALITY"] == label, "QD"]
        ax.hist(vals, bins=40, alpha=0.65, label=label, color=color)
    ax.axvline(MIN_HC_QD, linestyle="--", color="black", linewidth=1.2, label=f"QD={MIN_HC_QD}")
    ax.set_xlabel("Quality by Depth (QD)")
    ax.set_ylabel("Count")
    ax.set_title("QD Distribution by HC Variant Quality")
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


# ── Summaries ─────────────────────────────────────────────────────────────────

def write_qc_summary(df: pd.DataFrame, out_path: str):
    qc_cols = [c for c in [
        "SAMPLE_ID", "SAMPLE_QC_STATUS", "NUM_SAMPLE_FLAGS",
        "MEAN_TARGET_COVERAGE", "PCT_TARGET_BASES_20X", "PCT_TARGET_BASES_30X",
        "PCT_EXC_DUPE", "MEDIAN_TARGET_COVERAGE",
        "FLAG_LOW_COVG", "FLAG_LOW_20X", "FLAG_HIGH_DUPE",
    ] if c in df.columns]

    (df[qc_cols]
     .drop_duplicates("SAMPLE_ID")
     .to_csv(out_path, index=False))


def write_variant_summary(df: pd.DataFrame, out_path: str):
    rows = []
    for gene in df.get("SYMBOL.PRIMARY", pd.Series(dtype=str)).dropna().unique():
        sub = df[df.get("SYMBOL.PRIMARY") == gene] if "SYMBOL.PRIMARY" in df.columns else df
        for cls in (df.get(CLASSIFICATION_COL, pd.Series(dtype=str)).dropna().unique()
                    if CLASSIFICATION_COL in df.columns else []):
            s2 = sub[sub.get(CLASSIFICATION_COL) == cls] if CLASSIFICATION_COL in sub.columns else sub
            for csq in (s2.get("Consequence.PRIMARY", pd.Series(dtype=str)).dropna().unique()
                        if "Consequence.PRIMARY" in s2.columns else []):
                s3 = s2[s2.get("Consequence.PRIMARY") == csq] if "Consequence.PRIMARY" in s2.columns else s2
                rows.append({
                    "Gene":        gene,
                    "Classification": cls,
                    "Consequence": csq,
                    "N_variants":  len(s3),
                })
    pd.DataFrame(rows).to_csv(out_path, index=False)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-parquet",    required=True)
    ap.add_argument("--output-dir",       required=True)
    ap.add_argument("--prefilter-parquet", required=True,
                    help="Path for the pre-hard-filter parquet (all QC flags attached)")
    args = ap.parse_args()

    fig_dir = os.path.join(args.output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    print("Loading classified variants...")
    df_pl = pl.read_parquet(args.input_parquet)

    df_pl = add_sample_qc_flags(df_pl)
    df_pl = add_variant_qc_flags(df_pl)

    # ── Pre-filter snapshot (all variants, QC flags attached) ─────────────────
    os.makedirs(os.path.dirname(os.path.abspath(args.prefilter_parquet)), exist_ok=True)
    df_pl.write_parquet(args.prefilter_parquet)
    print(f"  Pre-filter:  {df_pl.height:,} variants → {args.prefilter_parquet}")

    df_filtered = apply_hard_filters(df_pl)
    print(f"  Post-filter: {df_filtered.height:,} variants")

    # Save final tables
    df_filtered.write_parquet(os.path.join(args.output_dir, "final_variant_table.parquet"))
    df_filtered.write_csv(os.path.join(args.output_dir, "final_variant_table.csv"))

    # Convert to pandas for plotting / summaries (full pre-filter df for QC)
    df_pd   = df_pl.to_pandas()
    df_f_pd = df_filtered.to_pandas()

    print("Writing summaries...")
    write_qc_summary(df_pd,   os.path.join(args.output_dir, "qc_summary.csv"))
    write_variant_summary(df_f_pd, os.path.join(args.output_dir, "variant_summary.csv"))

    print("Generating figures...")
    fig_vaf_by_classification(df_f_pd, os.path.join(fig_dir, "vaf_by_classification.png"))
    fig_cn_distribution(df_f_pd,       os.path.join(fig_dir, "cn_distribution.png"))
    fig_anchor_vaf(df_pd,              os.path.join(fig_dir, "anchor_vaf.png"))
    fig_coverage_vs_20x(df_pd,         os.path.join(fig_dir, "coverage_vs_20x.png"))
    fig_qd_vs_qual(df_pd,              os.path.join(fig_dir, "qd_vs_qual.png"))

    print(f"All outputs written to {args.output_dir}")


if __name__ == "__main__":
    main()
