#!/bin/bash
#SBATCH --job-name=t2_bundled
#SBATCH --account=bio250374
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/anvil/projects/x-bio250374/daniel/TESorter2_dev/_regression/logs/bundled_%j.out
#SBATCH --error=/anvil/projects/x-bio250374/daniel/TESorter2_dev/_regression/logs/bundled_%j.err

# Zero-config check: a clean venv with ONLY the wheel installed (no checkout,
# no --db-dir, no TESORTER2_DB) must find the packaged databases and reproduce
# the pre-refactor baseline in old_out/.
set -euo pipefail

module load python/3.9.5
REG=/anvil/projects/x-bio250374/daniel/TESorter2_dev/_regression
INPUT=/anvil/projects/x-bio250374/daniel/bath_vs_hmmer/phase3_data/edta_ltr_100.fa
export PATH="/home/x-dmahechalope/.conda/envs/EDTA-RepGen/bin:$PATH"

unset TESORTER2_DB || true
cd /tmp                       # deliberately outside the checkout
rm -rf "$REG/bundled_out"

echo "=== tesorter2 with NO --db-dir, from a wheel-only install ==="
"$REG/venv_clean/bin/tesorter2" "$INPUT" -d rexdb -p 8 \
    -o "$REG/bundled_out" --prefix t

echo
echo "=== compare against pre-refactor baseline (old_out) ==="
status=0
shopt -s nullglob
for n in "$REG"/bundled_out/*.tsv "$REG"/bundled_out/*.aa; do
    f=$(basename "$n"); o="$REG/old_out/$f"
    if [[ ! -f "$o" ]]; then echo "MISSING in baseline: $f"; status=1; continue; fi
    # old_out predates the Sequence-Ontology columns (SO_name/SO_ID, cols 9-10),
    # so for .cls.tsv compare only the deterministic original columns 1-8; the
    # .aa translation is compared in full.
    case "$f" in
        *.tsv) ocmp() { cut -f1-8 "$1"; };;
        *)     ocmp() { cat "$1"; };;
    esac
    if diff -q <(ocmp "$o") <(ocmp "$n") >/dev/null; then
        echo "IDENTICAL  $f  ($(wc -l < "$n") lines)"
    else
        echo "DIFFERS    $f"; diff <(ocmp "$o") <(ocmp "$n") | head -15; status=1
    fi
done

echo
if [[ $status -eq 0 ]]; then
    echo "RESULT: PASS - zero-config install reproduces the baseline"
else
    echo "RESULT: FAIL"
fi
exit $status
