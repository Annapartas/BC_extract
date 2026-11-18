# BC_extraction

Extract and count barcode fragments from FASTQ reads using **Cutadapt** and exact IUPAC sequence validation. Both read orientations are supported.

## Files

- `BC_extract.py` — extraction and validation pipeline
- `BC_extract_params.json` — input, sequence, processing, and output settings
- `run_BC_extract.sh` — local or SLURM launcher
- `environment.yml` — Conda dependencies

## Setup

```bash
git clone https://github.com/Annapartas/BC_extract.git
cd BC_extract
conda env create -f environment.yml
conda activate BC_extraction
```

Replace `YOUR_USERNAME` with your GitHub username. The cluster must have Conda (or compatible software) installed.

## Configure

Edit `BC_extract_params.json` before running:

- `input.files`: FASTQ glob(s), e.g. `"/path/to/reads/*.fastq.gz"`
- `sequences`: upstream flank, expected IUPAC fragment pattern, downstream flank
- `processing`: Cutadapt cores, pigz threads, compression, optional intermediate files
- `output.directory`: destination for results

**Paths in the JSON are relative to the directory from which the job is launched**, not the repository. Absolute input/output paths are recommended on a cluster. Set `cutadapt_cores` and `pigz_threads` to no more than the CPUs allocated to the job.

## Run

From a login/workstation shell (with the environment active):

```bash
bash run_BC_extract.sh
```

On a **SLURM cluster**, adjust the `#SBATCH` resource lines in `run_BC_extract.sh` to your cluster's requirements, then submit:

```bash
sbatch run_BC_extract.sh
```

SLURM jobs may not inherit your activated Conda environment. To activate it within the job, set `CONDA_BASE` in your submission environment to your Conda installation directory, e.g. `export CONDA_BASE="$HOME/miniforge3"`; the script then activates the `BC_extraction` environment. Alternatively, configure the cluster's environment/modules so that `python`, `cutadapt`, and `pigz` are on `PATH` in the job.

You can also run the Python script directly:

```bash
python BC_extract.py --parameters BC_extract_params.json
```

## Outputs

Results are written to `BC_extraction_output/` by default:

```text
BC_extraction_output/
├── pipeline_summary.tsv
└── SAMPLE/
    ├── SAMPLE_validated.fastq.gz
    ├── SAMPLE_fragment_counts.tsv
    ├── SAMPLE_cutadapt_report.txt
    ├── SAMPLE_cutadapt_report.json
    ├── SAMPLE_pipeline_summary.json
    ├── SAMPLE_flank_unmatched.fastq.gz
    └── SAMPLE_pattern_rejected.fastq.gz
```

The last two FASTQ files are optional and controlled in the JSON. The pipeline requires exact 15-nt flanks, extracts the intervening fragment, and retains only fragments matching the complete IUPAC pattern. Each input FASTQ is processed separately.
