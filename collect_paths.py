"""
collect_paths.py
Builds file_paths_mapping.csv from the sample manifest.

Usage:
    python collect_paths.py \
        --manifest transferred_samples.csv \
        --data-root /data \
        --output file_paths_mapping.csv
"""

import argparse
import csv
import glob
import os
import sys


# ── Priority selectors (ported from create_Mu2_RobustCNV_paths.py) ────────────

def select_cnv_file(files):
    """Prefer somatic CNV files over germline."""
    if not files:
        return None
    somatic = [f for f in files if "somatic" in f.lower()]
    if somatic:
        return somatic[0]
    standard = [f for f in files if "germline" not in f.lower()]
    if standard:
        return standard[0]
    return files[0]


def select_mutect_file(files):
    """Prefer somatic_output / final_snvindels over MNV or paired files."""
    if not files:
        return None
    priority_keywords = ["somatic_output.tsv", "final_snvindels.tsv", "snv_indel.tsv"]
    for kw in priority_keywords:
        match = [f for f in files if kw in f.lower()]
        if match:
            return match[0]
    non_mnv = [f for f in files if "mnv" not in f.lower()]
    return non_mnv[0] if non_mnv else files[0]


def find_bam(data_root, briefcase, cbio_id):
    pattern = os.path.join(data_root, briefcase, cbio_id, "aligned", cbio_id + "*.bam")
    matches = glob.glob(pattern)
    return matches[0] if matches else "NOT_FOUND"


def find_metrics(data_root, briefcase, cbio_id):
    """Locate HS metrics file; tries several common naming patterns."""
    patterns = [
        os.path.join(data_root, briefcase, cbio_id, "metrics", cbio_id + "*hsMetrics.txt"),
        os.path.join(data_root, briefcase, cbio_id, "metrics", cbio_id + "*hs_metrics.txt"),
        os.path.join(data_root, briefcase, cbio_id, "metrics", cbio_id + ".txt"),
    ]
    for p in patterns:
        matches = glob.glob(p)
        if matches:
            return matches[0]
    return "NOT_FOUND"


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest",  required=True, help="Path to transferred_samples.csv")
    ap.add_argument("--data-root", required=True, help="Root data directory (e.g. /data)")
    ap.add_argument("--output",    required=True, help="Output CSV path")
    args = ap.parse_args()

    found_data = []
    print(f"Searching for files under {args.data_root} ...")

    with open(args.manifest) as f:
        reader = csv.DictReader(f)
        for row in reader:
            briefcase = row["briefcase"]
            cbio_id   = row["cbio_sample_id"]

            # CNV
            cnv_pattern = os.path.join(
                args.data_root, briefcase, cbio_id,
                "analysis", "RobustCNV_*", "*Gene_level_calls_ALL.tsv"
            )
            selected_cnv = select_cnv_file(glob.glob(cnv_pattern))

            # Mutect2
            mutect_pattern = os.path.join(
                args.data_root, briefcase, cbio_id,
                "analysis", "MuTect_2", "*.tsv"
            )
            selected_mutect = select_mutect_file(glob.glob(mutect_pattern))

            # BAM
            bam = find_bam(args.data_root, briefcase, cbio_id)

            # HS Metrics
            metrics = find_metrics(args.data_root, briefcase, cbio_id)

            found_data.append({
                "cbio_sample_id": cbio_id,
                "briefcase":      briefcase,
                "cnv_path":       selected_cnv    or "NOT_FOUND",
                "mutect_path":    selected_mutect or "NOT_FOUND",
                "bam_path":       bam,
                "metrics_path":   metrics,
            })

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["cbio_sample_id", "briefcase", "cnv_path", "mutect_path", "bam_path", "metrics_path"]
        )
        writer.writeheader()
        writer.writerows(found_data)

    missing_bam    = sum(1 for x in found_data if x["bam_path"]    == "NOT_FOUND")
    missing_cnv    = sum(1 for x in found_data if x["cnv_path"]    == "NOT_FOUND")
    missing_mutect = sum(1 for x in found_data if x["mutect_path"] == "NOT_FOUND")
    missing_met    = sum(1 for x in found_data if x["metrics_path"]== "NOT_FOUND")

    print(f"Done — {len(found_data)} samples written to {args.output}")
    print(f"  Missing BAM:     {missing_bam}")
    print(f"  Missing CNV:     {missing_cnv}")
    print(f"  Missing MuTect2: {missing_mutect}")
    print(f"  Missing Metrics: {missing_met}")

    if missing_bam > 0:
        print("WARNING: Some samples are missing BAM files and will be skipped in HaplotypeCaller.", file=sys.stderr)


if __name__ == "__main__":
    main()
