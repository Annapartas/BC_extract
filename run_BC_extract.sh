#!/usr/bin/env bash
#SBATCH --job-name=BC_extract
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=BC_extract_%j.log
set -euo pipefail

# Run with: bash run_BC_extract.sh OR sbatch run_BC_extract.sh
# Install dependencies first (see README). The environment must be active
# or its executables must be on PATH in the compute job.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Optional: activate Conda automatically on compute nodes by setting
# CONDA_BASE to your installation directory, e.g. /home/user/miniforge3.
if [[ -n "${CONDA_BASE:-}" ]]; then
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "${BC_CONDA_ENV:-BC_extraction}"
fi

# Avoid Cutadapt oversubscribing CPUs on a SLURM allocation.
# Configure processing.cutadapt_cores and pigz_threads in the JSON to
# fit --cpus-per-task (default: 4).
python "${SCRIPT_DIR}/BC_extract.py" --parameters "${SCRIPT_DIR}/BC_extract_params.json" "$@"
