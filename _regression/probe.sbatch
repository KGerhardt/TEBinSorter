#!/bin/bash
#SBATCH --job-name=t2_strand
#SBATCH --account=bio250374
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/strand_%j.out
#SBATCH --error=logs/strand_%j.err
set -euo pipefail
module load python/3.9.5

# See run_regression.sbatch for why the checkout is resolved this way.
REG=""
for _c in "${TESORTER2_DEV:+$TESORTER2_DEV/_regression}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)" \
          "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$_c" && -f "$_c/run_regression.sbatch" ]] && { REG="$_c"; break; }
done
[[ -n "$REG" ]] || { echo "cannot locate _regression/; set TESORTER2_DEV" >&2; exit 1; }
DEV="$(dirname "$REG")"
mkdir -p "$REG/logs"

DB="$DEV/tesorter2/database/AnnoSINE_core.hmm"
GENOME="${TESORTER2_TEST_GENOME:-/anvil/projects/x-bio250374/daniel/minirice.fa}"
"$REG/venv/bin/python" "$REG/strand_probe.py" "$DB" "$GENOME" 100
