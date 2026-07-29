#!/bin/bash
#SBATCH --job-name=t2_strand
#SBATCH --account=bio250374
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/anvil/projects/x-bio250374/daniel/TESorter2_dev/_regression/logs/strand_%j.out
#SBATCH --error=/anvil/projects/x-bio250374/daniel/TESorter2_dev/_regression/logs/strand_%j.err
set -euo pipefail
module load python/3.9.5
REG=/anvil/projects/x-bio250374/daniel/TESorter2_dev/_regression
DB=/anvil/projects/x-bio250374/daniel/TESorter2_dev/tesorter2/database/AnnoSINE_core.hmm
"$REG/venv/bin/python" "$REG/strand_probe.py" "$DB" \
    /anvil/projects/x-bio250374/daniel/minirice.fa 100
