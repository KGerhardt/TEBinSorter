#!/bin/bash
#SBATCH --job-name=t2_multidb
#SBATCH --account=bio250374
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/multidb_%j.out
#SBATCH --error=logs/multidb_%j.err

# End-to-end check of the multi-database path (-d rexdb,gydb), which no other
# regression covers: every other case runs a single database, where clade
# harmonization is a deliberate no-op. Asserts that the cross-database vote runs
# on harmonized labels (see tesorter2/clade_harmonize.py) and that the per-
# database SecondaryHits audit trail keeps native names.
set -euo pipefail

module load python/3.9.5

# Locate the suite rather than hardcoding it; see run_regression.sbatch for why
# each candidate is probed for a marker file.
REG=""
for _c in "${TESORTER2_DEV:+$TESORTER2_DEV/_regression}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)" \
          "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$_c" && -f "$_c/multidb_e2e.sbatch" ]] && { REG="$_c"; break; }
done
[[ -n "$REG" ]] || { echo "cannot locate _regression/; set TESORTER2_DEV" >&2; exit 1; }
DEV="$(dirname "$REG")"
VENV=$REG/venv
mkdir -p "$REG/logs"

export PATH="${TESORTER2_TOOLS:-/home/x-dmahechalope/.conda/envs/EDTA-RepGen/bin}:$PATH"
INPUT="${TESORTER2_TEST_LIB:-/anvil/projects/x-bio250374/daniel/bath_vs_hmmer/phase3_data/edta_ltr_100.fa}"

cd "$REG"
rm -rf multidb_out single_out

# The venv holds an installed COPY of the package, not an editable install, so
# $VENV/bin/tesorter2 would run whatever was last pip-installed. Run the
# checkout through PYTHONPATH instead: same interpreter and dependencies, but
# the code under test is the working tree, and the venv is left untouched.
export PYTHONPATH="$DEV${PYTHONPATH:+:$PYTHONPATH}"
T2=("$VENV/bin/python" -m tesorter2.pipeline)

echo "=== [1/3] Multi-database run (-d rexdb,gydb) ==="
"${T2[@]}" "$INPUT" \
    -d rexdb,gydb -p 8 --db-dir "$DEV/tesorter2/database" \
    -o "$REG/multidb_out" --prefix t

echo
echo "=== [2/3] Per-database outputs must all be present ==="
status=0
for f in t.cls.tsv t.rexdb.cls.tsv t.gydb.cls.tsv; do
    if [[ -f "$REG/multidb_out/$f" ]]; then
        echo "OK      $f  ($(( $(wc -l < "$REG/multidb_out/$f") - 1 )) elements)"
    else
        echo "MISSING $f"; status=1
    fi
done

echo
echo "=== Compatibility contract: per-DB files stay at TEsorter's 7 columns ==="
# EDTA reads *.{db}db.cls.tsv positionally as (id, Order, Superfamily) and must
# never see harmonized names or extra columns. Harmonization rewrites Order and
# Superfamily too, so a leak here silently changes EDTA's classification.
for db in rexdb gydb; do
    f="$REG/multidb_out/t.$db.cls.tsv"
    ncol=$(head -1 "$f" | awk -F'\t' '{print NF}')
    if [[ "$ncol" -ne 7 ]]; then
        echo "FAIL: t.$db.cls.tsv has $ncol columns, expected 7"; status=1
    elif grep -qE $'\t'"(Bel-Pao|Caulimoviridae|Retrovirus)"$'\t' "$f"; then
        echo "FAIL: harmonized superfamily leaked into t.$db.cls.tsv"; status=1
    else
        echo "OK      t.$db.cls.tsv  7 columns, native names only"
    fi
done

echo
echo "=== [3/3] Harmonization assertions on the combined call ==="
"$VENV/bin/python" "$REG/check_multidb.py" \
    "$REG/multidb_out/t.cls.tsv" || status=1

echo
echo "=== Single-database no-op: -d rexdb must keep native names ==="
# Harmonization is gated on >= 2 databases voting, so a single-database run has
# to come out on native REXdb spellings. run_regression.sbatch compares this
# same path against the pre-refactor baseline; here we only assert the gate.
"${T2[@]}" "$INPUT" \
    -d rexdb -p 8 --db-dir "$DEV/tesorter2/database" \
    -o "$REG/single_out" --prefix t >/dev/null 2>&1
if diff -q <(cut -f1-7 "$REG/single_out/t.cls.tsv") \
           <(cut -f1-7 "$REG/single_out/t.rexdb.cls.tsv") >/dev/null; then
    echo "OK: single-db combined call is identical to the per-database call"
else
    echo "NOTE: single-db combined differs from per-database (BLAST pass-2 may"
    echo "      reclassify; checking only that no unified name appeared)"
    if cut -f4 "$REG/single_out/t.cls.tsv" | grep -qE '^(Ale/Retrofit|Tekay/Del|Ivana/Oryco)$'; then
        echo "FAIL: unified clade name in a single-database run"
        status=1
    else
        echo "OK: no unified clade names in the single-database run"
    fi
fi

echo
if [[ $status -eq 0 ]]; then
    echo "RESULT: PASS - multi-database harmonization is active and scoped correctly"
else
    echo "RESULT: FAIL - see above"
fi
exit $status
