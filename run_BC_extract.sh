#!/usr/bin/env bash
#SBATCH --job-name=BC_extract
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=BC_extract_%j.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.email@university.ch

set -euo pipefail

# Find the directory containing this script
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Activate Conda if CONDA_BASE is specified
if [[ -n "${CONDA_BASE:-}" ]]; then
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "${BC_CONDA_ENV:-BC_extraction}"
fi

# Run barcode extraction
python BC_extract.py --parameters BC_extract_params.json

echo "Barcode extraction completed successfully."
