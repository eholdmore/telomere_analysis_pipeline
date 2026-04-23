"""
hs_metrics_summary.py
Reads per-sample Picard HS metrics files (paths from file_paths_mapping.csv)
and writes a consolidated QC summary CSV.

Usage:
    python hs_metrics_summary.py \
        --paths-csv file_paths_mapping.csv \
        --output hs_metrics_summary.csv
"""

import argparse
import csv
import os


METRIC_COLS = [
    "MEAN_TARGET_COVERAGE",
    "MEDIAN_TARGET_COVERAGE",
    "PCT_TARGET_BASES_10X",
    "PCT_TARGET_BASES_20X",
    "PCT_TARGET_BASES_30X",
    "PCT_EXC_DUPE",
]


def parse_hs_metrics(filepath, sample_id):
    """Parse a Picard HS metrics file and return a dict of the metrics we need."""
    try:
        with open(filepath) as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith("## METRICS CLASS"):
                headers = lines[i + 1].strip().split("\t")
                values  = lines[i + 2].strip().split("\t")
                row = dict(zip(headers, values))
                row["SAMPLE_ID"] = sample_id
                return row
    except Exception as e:
        print(f"  WARNING: could not parse {filepath}: {e}")
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paths-csv", required=True, help="file_paths_mapping.csv")
    ap.add_argument("--output",    required=True, help="Output CSV path")
    args = ap.parse_args()

    metrics_data = []
    missing = []

    with open(args.paths_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id    = row["cbio_sample_id"]
            metrics_path = row.get("metrics_path", "NOT_FOUND")

            if metrics_path == "NOT_FOUND" or not os.path.exists(metrics_path):
                missing.append(sample_id)
                continue

            parsed = parse_hs_metrics(metrics_path, sample_id)
            if parsed:
                metrics_data.append(parsed)

    if not metrics_data:
        print("ERROR: No HS metrics data found. Check paths in file_paths_mapping.csv.")
        raise SystemExit(1)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    fieldnames = ["SAMPLE_ID"] + METRIC_COLS
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(metrics_data)

    print(f"HS metrics summary written to {args.output} ({len(metrics_data)} samples)")
    if missing:
        print(f"  WARNING: {len(missing)} samples had no metrics file: {missing[:5]}{'...' if len(missing)>5 else ''}")


if __name__ == "__main__":
    main()
