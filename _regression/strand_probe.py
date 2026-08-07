"""
Does the current DNA path miss minus-strand hits?

Takes real sequences, searches the SINE (DNA) database with the engine the
pipeline uses today (hmmsearch, DNA alphabet) and with nhmmer, on both the
sequences as given and their reverse complement. If hmmsearch is strand-blind,
its hit count collapses on the reverse-complemented input while nhmmer's does
not.
"""
import sys
import pyhmmer
from pyhmmer import easel, plan7

DB = sys.argv[1]
FASTA = sys.argv[2]
LIMIT = int(sys.argv[3]) if len(sys.argv) > 3 else 200

dna = easel.Alphabet.dna()

with plan7.HMMFile(DB) as fh:
    hmms = list(fh)
print(f"database: {DB}")
print(f"  models: {len(hmms)}  alphabet: {'DNA' if hmms[0].alphabet.is_dna() else hmms[0].alphabet}")

# Load sequences, chunking anything long: hmmsearch has a hard 100k-residue
# per-sequence limit, so long targets must be windowed to compare the engines
# on the same input at all.
WIN = 50_000
seqs = []
with easel.SequenceFile(FASTA, digital=False) as sf:
    for s in sf:
        txt = s.sequence
        name = s.name.decode() if isinstance(s.name, bytes) else s.name
        for i in range(0, len(txt), WIN):
            chunk = txt[i:i + WIN]
            ts = easel.TextSequence(name=f"{name}_{i}".encode(), sequence=chunk)
            seqs.append(ts.digitize(dna))
            if len(seqs) >= LIMIT:
                break
        if len(seqs) >= LIMIT:
            break
seqs = seqs[:LIMIT]
print(f"  sequences: {len(seqs)} windows (<= {WIN} bp) from {FASTA}")
print(f"  total bp: {sum(len(s) for s in seqs):,}")

fwd = easel.DigitalSequenceBlock(dna, seqs)
rc_seqs = []
for s in seqs:
    t = s.copy()
    t.reverse_complement(inplace=True)
    rc_seqs.append(t)
rc = easel.DigitalSequenceBlock(dna, rc_seqs)


def count_hmmsearch(block):
    """Exactly how pipeline.run_database_legacy searches DNA databases today."""
    n = 0
    for top in pyhmmer.hmmsearch(hmms, block, bias_filter=False,
                                 Z=len(hmms), domZ=len(hmms), E=1e10, domE=1e10):
        for hit in top:
            if hit.evalue < 1e-3:
                n += 1
    return n


def count_nhmmer(block):
    n = 0
    for top in pyhmmer.nhmmer(hmms, block, E=1e10):
        for hit in top:
            if hit.evalue < 1e-3:
                n += 1
    return n


def best_evalues(fn_name, block, k=3):
    """Show the strongest hits so we can tell signal from noise."""
    ev = []
    if fn_name == "hmmsearch":
        it = pyhmmer.hmmsearch(hmms, block, bias_filter=False,
                               Z=len(hmms), domZ=len(hmms), E=1e10, domE=1e10)
    else:
        it = pyhmmer.nhmmer(hmms, block, E=1e10)
    for top in it:
        for hit in top:
            ev.append(hit.evalue)
    ev.sort()
    return ev[:k]


print()
print(f"{'engine':<12} {'forward':>10} {'revcomp':>10}   verdict")
print("-" * 56)
hs_f, hs_r = count_hmmsearch(fwd), count_hmmsearch(rc)
# A both-strand engine returns the SAME count on an input and its revcomp.
verdict = "single-strand only" if hs_f and hs_f != hs_r else "both strands"
print(f"{'hmmsearch':<12} {hs_f:>10} {hs_r:>10}   {verdict}  <- current DNA path")

nh_f, nh_r = count_nhmmer(fwd), count_nhmmer(rc)
verdict = "symmetric" if nh_f and abs(nh_f - nh_r) <= max(1, 0.1 * nh_f) else "asymmetric?"
print(f"{'nhmmer':<12} {nh_f:>10} {nh_r:>10}   {verdict}")

print()
print("strongest E-values (sanity: is there real SINE signal at all?)")
for eng in ("hmmsearch", "nhmmer"):
    print(f"  {eng:<10} fwd={best_evalues(eng, fwd)}")
    print(f"  {eng:<10} rc ={best_evalues(eng, rc)}")
