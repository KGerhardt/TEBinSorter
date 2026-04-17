"""
deconflict.py — Fast in-memory deconfliction and export of HMM search results.

Loads all hits from SQLite into numpy arrays, computes best-per-family
assignments in one pass, exports clean flat files.
"""

import sqlite3
import numpy as np


def _get_family(model_name):
    """Extract domain family from model name."""
    if ":" in model_name:
        return model_name.split(":")[1].split("-")[-1]
    return model_name.split("_")[0]


def _get_base_seq(target_name):
    """Strip frame suffix to get original sequence name."""
    return target_name.rsplit("|", 1)[0]


def load_hits(db_path, table="legacy_hits", database=None):
    """Load hits from SQLite into structured numpy arrays.

    Args:
        db_path: path to results .db file
        table: table name
        database: filter to specific database name (or None for all)

    Returns:
        dict with numpy arrays: target, model, score, evalue, acc,
        hmm_from, hmm_to, model_len, env_from, env_to,
        plus derived: base_seq, family, hmm_cov, norm_score
    """
    conn = sqlite3.connect(db_path)

    query = f"""
        SELECT target_name, query_name, dom_score, i_evalue, acc,
               hmm_from, hmm_to, query_len, env_from, env_to
        FROM {table}
    """
    params = ()
    if database:
        query += " WHERE database = ?"
        params = (database,)

    rows = conn.execute(query, params).fetchall()
    conn.close()

    n = len(rows)
    if n == 0:
        return None

    target = np.array([r[0] for r in rows])
    model = np.array([r[1] for r in rows])
    score = np.array([r[2] for r in rows], dtype=np.float64)
    evalue = np.array([r[3] for r in rows], dtype=np.float64)
    acc = np.array([r[4] for r in rows], dtype=np.float64)
    hmm_from = np.array([r[5] for r in rows], dtype=np.int32)
    hmm_to = np.array([r[6] for r in rows], dtype=np.int32)
    model_len = np.array([r[7] for r in rows], dtype=np.int32)
    env_from = np.array([r[8] for r in rows], dtype=np.int32)
    env_to = np.array([r[9] for r in rows], dtype=np.int32)

    base_seq = np.array([_get_base_seq(t) for t in target])
    family = np.array([_get_family(m) for m in model])
    hmm_cov = 100.0 * (hmm_to - hmm_from + 1) / model_len
    norm_score = score / model_len

    return {
        "target": target, "model": model, "score": score,
        "evalue": evalue, "acc": acc,
        "hmm_from": hmm_from, "hmm_to": hmm_to, "model_len": model_len,
        "env_from": env_from, "env_to": env_to,
        "base_seq": base_seq, "family": family,
        "hmm_cov": hmm_cov, "norm_score": norm_score,
    }


def best_per_family(hits):
    """Find the best-scoring model per (base_sequence, domain_family).

    Args:
        hits: dict from load_hits()

    Returns:
        numpy index array into hits for the best entries
    """
    keys = np.char.add(np.char.add(hits["base_seq"], "|"), hits["family"])
    sort_idx = np.argsort(-hits["score"])
    _, first_idx = np.unique(keys[sort_idx], return_index=True)
    return sort_idx[first_idx]


def best_per_frame(hits):
    """Find the single best-scoring hit per translated frame.

    Args:
        hits: dict from load_hits()

    Returns:
        numpy index array into hits for the best entries
    """
    sort_idx = np.argsort(-hits["score"])
    _, first_idx = np.unique(hits["target"][sort_idx], return_index=True)
    return sort_idx[first_idx]


def best_per_seq_model(hits):
    """Find the best-scoring domain per (target_frame, model) pair.

    Args:
        hits: dict from load_hits()

    Returns:
        numpy index array into hits for the best entries
    """
    keys = np.char.add(np.char.add(hits["target"], "|"), hits["model"])
    sort_idx = np.argsort(-hits["score"])
    _, first_idx = np.unique(keys[sort_idx], return_index=True)
    return sort_idx[first_idx]


def export_best_tsv(hits, indices, out_path, nucl_lengths=None):
    """Export selected hits as a clean TSV.

    Args:
        hits: dict from load_hits()
        indices: index array from best_per_* functions
        out_path: output file path
        nucl_lengths: optional {seq_name: length} for coordinate conversion
    """
    from sequence import parse_frame_suffix, aa_to_nucl_coords

    columns = [
        "seq_id", "model", "family", "strand", "frame",
        "nuc_start", "nuc_end", "env_from_aa", "env_to_aa",
        "hmm_from", "hmm_to", "hmm_cov",
        "score", "norm_score", "evalue", "accuracy",
    ]

    with open(out_path, "w") as f:
        f.write("\t".join(columns) + "\n")

        for i in indices:
            target = hits["target"][i]
            model = hits["model"][i]
            family = hits["family"][i]
            seq_id, strand, frame = parse_frame_suffix(target)

            env_from = int(hits["env_from"][i])
            env_to = int(hits["env_to"][i])

            if (strand in ("+", "-") and nucl_lengths
                    and seq_id in nucl_lengths):
                nuc_start, nuc_end = aa_to_nucl_coords(
                    env_from, env_to, strand, frame, nucl_lengths[seq_id])
            else:
                nuc_start, nuc_end = env_from, env_to
                if strand == ".":
                    frame = "."

            vals = [
                seq_id, model, family, strand, frame,
                nuc_start, nuc_end, env_from, env_to,
                int(hits["hmm_from"][i]), int(hits["hmm_to"][i]),
                f"{hits['hmm_cov'][i]:.1f}",
                f"{hits['score'][i]:.1f}",
                f"{hits['norm_score'][i]:.4f}",
                f"{hits['evalue'][i]:.2e}",
                f"{hits['acc'][i]:.3f}",
            ]
            f.write("\t".join(str(v) for v in vals) + "\n")
