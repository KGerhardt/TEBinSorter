"""
decompose_hmm.py — Decompose HMM models into high-scoring sub-model domains.

For each HMM, finds non-overlapping windows of N positions that maximize
emission log-odds score. These sub-models are the most informative
regions of the model — ideal for fast pre-screening.

Adapted from jewelry_hammer for TEBinSorter.
"""

import argparse
import sys
import numpy as np

from pyhmmer.plan7 import HMMFile
import pyhmmer.easel as easel

# Background frequencies
DNA_BG = np.array([0.2869, 0.1899, 0.2161, 0.3072])  # A C G T

# Robinson-Robinson amino acid background (HMMER default)
AA_BG = np.array([
    0.0787945, 0.0151600, 0.0535222, 0.0668298,  # A C D E
    0.0397062, 0.0695071, 0.0229198, 0.0590092,  # F G H I
    0.0594422, 0.0963728, 0.0237718, 0.0414386,  # K L M N
    0.0482904, 0.0395639, 0.0540978, 0.0683364,  # P Q R S
    0.0540687, 0.0673417, 0.0114135, 0.0304133,  # T V W Y
])


def compute_log_odds(hmm):
    """Compute log-odds emission scores for all positions.
    Auto-detects DNA vs amino acid alphabet.
    Returns (M+1, K) array where K is alphabet size.
    """
    M = hmm.M
    me = hmm.match_emissions

    if hmm.alphabet == easel.Alphabet.dna():
        bg = DNA_BG
        K = 4
    else:
        bg = AA_BG
        K = 20

    log_odds = np.zeros((M + 1, K))
    for i in range(1, M + 1):
        for j in range(K):
            p = me[i][j]
            log_odds[i, j] = np.log2(p / bg[j]) if p > 0 else -100
    return log_odds


def score_positions(hmm):
    """Compute per-position informativeness score.
    Returns array of length M with the max log-odds at each position.
    """
    lo = compute_log_odds(hmm)
    return np.max(lo[1:], axis=1)  # skip row 0


def find_domains(position_scores, window_size, max_domains=None,
                 min_score_per_base=0.0, max_overlap_frac=0.0):
    """Greedy window selection by total score with configurable overlap.

    Args:
        position_scores: array of per-position scores, length M
        window_size: sub-model size N
        max_domains: optional cap on number of domains
        min_score_per_base: minimum average score per base to accept a window
        max_overlap_frac: maximum fraction of window that can overlap with
                          already-selected domains (0.0 = no overlap,
                          0.33 = up to 33% overlap)

    Returns:
        list of (start_pos, end_pos, total_score) sorted by model position
    """
    M = len(position_scores)
    if M < window_size:
        return [(0, M, float(np.sum(position_scores)))]

    max_overlap = int(window_size * max_overlap_frac)

    # Score all windows
    n_windows = M - window_size + 1
    window_scores = np.zeros(n_windows)
    # Sliding sum
    window_scores[0] = np.sum(position_scores[:window_size])
    for i in range(1, n_windows):
        window_scores[i] = (window_scores[i-1]
                            - position_scores[i-1]
                            + position_scores[i + window_size - 1])

    # Greedy selection: pick best window, allow partial overlap
    used = np.zeros(M, dtype=int)  # count how many domains cover each position
    domains = []

    while True:
        if max_domains is not None and len(domains) >= max_domains:
            break

        # Mask windows that exceed overlap budget
        valid_scores = window_scores.copy()
        for i in range(n_windows):
            overlap_count = int(np.sum(used[i:i + window_size] > 0))
            if overlap_count > max_overlap:
                valid_scores[i] = -np.inf

        best_idx = np.argmax(valid_scores)
        best_score = valid_scores[best_idx]

        if best_score == -np.inf:
            break
        if best_score / window_size < min_score_per_base:
            break

        start = best_idx
        end = best_idx + window_size
        domains.append((start, end, float(best_score)))
        used[start:end] += 1

    # Sort by position
    domains.sort(key=lambda d: d[0])
    return domains


def decompose_model(hmm, window_size, max_domains=None,
                    min_score_per_base=0.0, max_overlap_frac=0.0):
    """Decompose one HMM into sub-model domains.

    Returns:
        name: model name
        M: model length
        domains: list of (start, end, score) tuples (0-based model positions)
        position_scores: per-position score array
    """
    name = hmm.name if isinstance(hmm.name, str) else hmm.name.decode()
    M = hmm.M

    position_scores = score_positions(hmm)
    domains = find_domains(position_scores, window_size, max_domains,
                           min_score_per_base, max_overlap_frac)

    return name, M, domains, position_scores


def decompose_database(hmm_path, window_size=64, max_domains=None,
                       min_score_per_base=0.0, max_overlap_frac=0.0):
    """Decompose all models in an HMM database file.

    Args:
        hmm_path: path to HMM database file (multi-model)
        window_size: sub-model window size
        max_domains: optional cap per model
        min_score_per_base: minimum avg score to accept a window
        max_overlap_frac: maximum overlap fraction between windows

    Returns:
        dict of {model_name: [(start, end, score), ...]}
    """
    from hmm import load_hmms

    hmms = load_hmms(hmm_path)
    results = {}

    for hmm in hmms:
        name, M, domains, _ = decompose_model(
            hmm, window_size, max_domains, min_score_per_base,
            max_overlap_frac)
        results[name] = domains

    return results


def splice_sub_hmm(src_hmm, start, end):
    """
    Splice a sub-HMM from a parent model by copying emission and
    transition probabilities for positions start:end.

    Args:
        src_hmm: parent pyhmmer HMM
        start: 0-based start position in parent model
        end: 0-based end position (exclusive)

    Returns:
        pyhmmer HMM object ready for searching
    """
    import pyhmmer.plan7 as plan7
    import io

    sub_M = end - start
    alphabet = src_hmm.alphabet
    K = alphabet.K  # 4 for DNA, 20 for amino

    name = src_hmm.name if isinstance(src_hmm.name, str) else src_hmm.name
    sub_name = f"{name}__sub_{start}-{end}"

    sub = plan7.HMM(alphabet, M=sub_M, name=sub_name)

    # Copy emissions and transitions from parent
    for i in range(1, sub_M + 1):
        src_i = int(start + i)
        for j in range(K):
            sub.match_emissions[i][j] = float(src_hmm.match_emissions[src_i][j])
            sub.insert_emissions[i][j] = float(src_hmm.insert_emissions[src_i][j])
        for j in range(7):
            sub.transition_probabilities[i][j] = float(
                src_hmm.transition_probabilities[src_i][j])

    # Entry state transitions
    for j in range(7):
        sub.transition_probabilities[0][j] = float(
            src_hmm.transition_probabilities[0][j])

    sub.set_composition()
    sub.consensus = src_hmm.consensus[start:end]

    # Write, inject STATS lines, reload for a fully valid HMM
    buf = io.BytesIO()
    sub.write(buf)
    text = buf.getvalue().decode()

    if "STATS LOCAL MSV" not in text:
        lines = text.split("\n")
        new_lines = []
        for line in lines:
            new_lines.append(line)
            if line.startswith("MAP"):
                new_lines.append("STATS LOCAL MSV       -9.9014  0.70957")
                new_lines.append("STATS LOCAL VITERBI  -10.7224  0.70957")
                new_lines.append("STATS LOCAL FORWARD   -4.1637  0.70957")
        text = "\n".join(new_lines)

    buf2 = io.BytesIO(text.encode())
    with plan7.HMMFile(buf2) as hf:
        return next(hf)


def build_sub_hmms(hmm_path, window_size=25, max_domains=None,
                   min_score_per_base=0.0, max_overlap_frac=0.0):
    """
    Decompose all models in a database and build searchable sub-HMMs.

    Args:
        hmm_path: path to HMM database file
        window_size: sub-model window size (default 25 for amino, 64 for DNA)
        max_domains: optional cap per model
        min_score_per_base: minimum avg score to accept a window
        max_overlap_frac: maximum overlap fraction between windows

    Returns:
        list of (sub_hmm, parent_name, start, end) tuples
    """
    from hmm import load_hmms

    hmms = load_hmms(hmm_path)
    sub_hmms = []

    for hmm in hmms:
        name = hmm.name
        M = hmm.M

        if M <= window_size:
            sub_hmms.append((hmm, name, 0, M))
            continue

        _, _, domains, _ = decompose_model(
            hmm, window_size, max_domains, min_score_per_base,
            max_overlap_frac)

        for start, end, score in domains:
            sub = splice_sub_hmm(hmm, start, end)
            sub_hmms.append((sub, name, start, end))

    return sub_hmms


def main():
    parser = argparse.ArgumentParser(
        description="Decompose HMM models into high-scoring sub-model domains")
    parser.add_argument("hmm_file",
                        help="HMM database file (single or multi-model)")
    parser.add_argument("--window", type=int, default=64,
                        help="Sub-model window size (default 64)")
    parser.add_argument("--max-domains", type=int, default=None,
                        help="Max domains per model (default: no limit)")
    parser.add_argument("--min-score-per-base", type=float, default=0.0,
                        help="Minimum avg score/base to accept domain")
    parser.add_argument("--output", default=None,
                        help="Output TSV (default: stdout)")
    args = parser.parse_args()

    from hmm import load_hmms
    hmms = load_hmms(args.hmm_file)

    results = []
    for hmm in hmms:
        name, M, domains, pos_scores = decompose_model(
            hmm, args.window, args.max_domains, args.min_score_per_base)
        results.append((name, M, domains, pos_scores))

    out = open(args.output, "w") if args.output else sys.stdout

    out.write(f"model\tM\tn_domains\tdomain_bp\tcoverage\t"
              f"domains\tavg_score_per_base\n")

    total_M = 0
    total_domain_bp = 0

    for name, M, domains, pos_scores in sorted(results, key=lambda r: -r[1]):
        n_dom = len(domains)
        domain_bp = sum(e - s for s, e, _ in domains)
        coverage = domain_bp / M if M > 0 else 0
        total_M += M
        total_domain_bp += domain_bp

        dom_strs = []
        for s, e, sc in domains:
            avg = sc / (e - s)
            dom_strs.append(f"{s}-{e}({sc:.1f})")

        avg_spb = sum(sc for _, _, sc in domains) / domain_bp if domain_bp else 0

        out.write(f"{name}\t{M}\t{n_dom}\t{domain_bp}\t{coverage:.2f}\t"
                  f"{','.join(dom_strs)}\t{avg_spb:.2f}\n")

    if out is not sys.stdout:
        out.close()

    n_models = len(results)
    n_domains = sum(len(d) for _, _, d, _ in results)
    print(f"\n{n_models} models -> {n_domains} domains at window={args.window}",
          file=sys.stderr)
    print(f"  Total model positions: {total_M:,}", file=sys.stderr)
    print(f"  Total domain positions: {total_domain_bp:,} "
          f"({100*total_domain_bp/total_M:.1f}% of model space)",
          file=sys.stderr)
    print(f"  Avg domains/model: {n_domains/n_models:.1f}", file=sys.stderr)


if __name__ == "__main__":
    main()
