#!/usr/bin/env bash
# setup_venv.sh
# Creates (or refreshes) the Python virtual environment needed by the pipeline.
# Run once on the login node before submitting the Snakemake workflow.
#
# Usage:
#   bash setup_venv.sh [optional: path/to/venv]
#
# Default venv location: ~/telomere_numpy_venv

set -euo pipefail

VENV_DIR="${1:-${HOME}/telomere_numpy_venv}"

echo "=== Telomere Pipeline — venv setup ==="
echo "Target:  ${VENV_DIR}"

# ── Create or reuse ──────────────────────────────────────────────────────────
if [ ! -d "${VENV_DIR}" ]; then
    echo "Creating new virtual environment..."
    python3 -m venv "${VENV_DIR}"
else
    echo "Virtual environment already exists — upgrading packages..."
fi

source "${VENV_DIR}/bin/activate"

# ── Core packages ────────────────────────────────────────────────────────────
pip install --upgrade pip wheel

pip install \
    polars \
    pandas \
    pyarrow \
    matplotlib \
    snakemake \
    snakemake-executor-plugin-slurm

echo ""
echo "=== Setup complete ==="
echo "Activate with:  source ${VENV_DIR}/bin/activate"
echo ""
echo "Test with:"
echo "  python -c 'import polars, pandas, matplotlib, snakemake; print(\"OK\")'"
