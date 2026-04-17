"""
Hot N-mer pre-filter for jewelry_hammer.

Extracts emission-score-derived k-mers ("hot N-mers") from HMM models.
A hot N-mer is a k-length sequence where every base has a positive
log-odds emission score at its corresponding model position. These are
sequences the model considers better than random at every position.

Hot N-mers are ranked by summed log-odds score.
"""

import math
import os
from collections import defaultdict
from itertools import product

import numpy as np

from hmm import load_hmms

BASES = "ACGT"
BASE_IDX = {b: i for i, b in enumerate(BASES)}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
BACKGROUND = {"A": 0.2869, "C": 0.1899, "G": 0.2161, "T": 0.3072}


def revcomp(seq):
    return "".join(COMP.get(b, "N") for b in reversed(seq))


def compute_log_odds(hmm):
    """Compute log-odds emission scores for all positions. Returns (M+1) x 4 array."""
    M = hmm.M
    me = hmm.match_emissions
    log_odds = np.zeros((M + 1, 4))
    for i in range(1, M + 1):
        for j in range(4):
            p = me[i][j]
            bg = BACKGROUND[BASES[j]]
            log_odds[i, j] = math.log2(p / bg) if p > 0 else -100
    return log_odds


def extract_hot_nmers(hmm, k=6):
    """
    Extract hot N-mers: k-length sequences where every base scores
    positive (better than background) at its model position.

    Returns:
        list of (kmer_string, model_start_pos, total_score) sorted by
        total_score descending
    """
    M = hmm.M
    log_odds = compute_log_odds(hmm)

    # At each position, which bases are positive?
    positive_bases = []  # pos index 0 = model pos 1
    for i in range(1, M + 1):
        pos_bases = []
        for j in range(4):
            if log_odds[i, j] > 0:
                pos_bases.append((j, log_odds[i, j]))
        positive_bases.append(pos_bases)

    hot_nmers = []

    for start in range(M - k + 1):
        window = [positive_bases[start + p] for p in range(k)]

        # Skip if any position has zero positive bases
        if any(len(wb) == 0 for wb in window):
            continue

        # Enumerate all-positive combinations
        for combo in product(*window):
            kmer = "".join(BASES[base_idx] for base_idx, _ in combo)
            total_score = sum(score for _, score in combo)
            model_pos = start + 1  # 1-based

            hot_nmers.append((kmer, model_pos, total_score))

    # Sort by score descending
    hot_nmers.sort(key=lambda x: -x[2])

    return hot_nmers


# ---------------------------------------------------------------------------
# CLI: inspect hot N-mers for models
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    hmm_dir = sys.argv[1] if len(sys.argv) > 1 else "hmm_plants"
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 6

    hmm_paths = sorted(
        os.path.join(hmm_dir, f) for f in os.listdir(hmm_dir)
        if f.endswith(".hmm"))

    print(f"k={k}, {len(hmm_paths)} models")
    print(f"{'Model':<35} {'M':>4} {'HotNmers':>9} {'Unique':>7} "
          f"{'TopScore':>8} {'MedScore':>8} {'PosCoverage':>11}")
    print("-" * 95)

    for path in hmm_paths:
        hmm = load_hmms(path)[0]
        name = hmm.name if isinstance(hmm.name, str) else hmm.name.decode()
        M = hmm.M

        hot = extract_hot_nmers(hmm, k)

        unique_kmers = len(set(kmer for kmer, _, _ in hot))
        positions = set(pos for _, pos, _ in hot)
        pos_coverage = len(positions) / M if M > 0 else 0

        if hot:
            scores = [s for _, _, s in hot]
            top = scores[0]
            med = scores[len(scores) // 2]
        else:
            top = med = 0

        print(f"{name:<35} {M:>4} {len(hot):>9} {unique_kmers:>7} "
              f"{top:>8.2f} {med:>8.2f} {pos_coverage:>10.1%}")

    # Show top 20 hot N-mers for first 3 models
    print("\n--- Top 20 hot N-mers per model (first 3) ---")
    for path in hmm_paths[:3]:
        hmm = load_hmms(path)[0]
        name = hmm.name if isinstance(hmm.name, str) else hmm.name.decode()
        hot = extract_hot_nmers(hmm, k)

        print(f"\n{name} (M={hmm.M}):")
        for kmer, pos, score in hot[:20]:
            rc = revcomp(kmer)
            print(f"  {kmer} (rc={rc})  pos={pos:>3}  score={score:.2f}")
