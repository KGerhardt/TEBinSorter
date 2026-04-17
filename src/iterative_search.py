"""
iterative_search.py — Facet screen → graph-guided confirmation → legacy fallback.

Orchestrates the iterative search process:
  1. Facet screen: bulk sub-HMM search, produces per-frame ranked candidate list
  2. Confirmation rounds: batched full-model searches guided by similarity graph,
     frames drop out as they confirm or exhaust candidates
  3. Legacy fallback: remaining unresolved frames get full nobias search

The work planner balances batches by estimated cost (M^2 * total_seqlen)
and is reused across all rounds.
"""

import json
import logging
import time
from collections import defaultdict

import pyhmmer
import pyhmmer.easel as easel

from decompose_hmm import build_sub_hmms_from_file, _load_hmms
from model_graph import get_or_build_graph
from search import (tophits_to_domtbl, parse_domtbl_text, _collect_hits,
                    _partition_hmms_by_size)

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Work planner: balance (model, frame_set) jobs by cost
# ---------------------------------------------------------------------------

def plan_batches(jobs, hmms_dict, seq_lens, n_batches):
    """Plan cost-balanced batches from a list of (model_name, [frame_names]).

    Each batch is a list of (model_name, [frame_names]) tuples whose
    total estimated cost (M^2 * sum_seqlen) is roughly equal.

    Args:
        jobs: list of (model_name, [frame_names])
        hmms_dict: {name: HMM} for M values
        seq_lens: {frame_name: length}
        n_batches: target number of batches

    Returns:
        list of batches, each a list of (model_name, [frame_names])
    """
    # Compute cost per job
    costed = []
    for model, frames in jobs:
        if model not in hmms_dict:
            continue
        m_sq = hmms_dict[model].M ** 2
        total_len = sum(seq_lens.get(f, 0) for f in frames)
        cost = m_sq * total_len
        costed.append((model, frames, cost))

    # Sort by cost descending
    costed.sort(key=lambda x: -x[2])

    # Greedy bin packing
    batches = [[] for _ in range(min(n_batches, len(costed)))]
    batch_costs = [0.0] * len(batches)

    for model, frames, cost in costed:
        # Put in lightest batch
        min_idx = batch_costs.index(min(batch_costs))
        batches[min_idx].append((model, frames))
        batch_costs[min_idx] += cost

    return [b for b in batches if b]


def run_batch(batch, hmms_dict, seq_block, name_to_idx, alphabet, Z):
    """Execute one batch: search each model against its assigned frames.

    Args:
        batch: list of (model_name, [frame_names])
        hmms_dict: {name: HMM}
        seq_block: full DigitalSequenceBlock
        name_to_idx: {frame_name: index in seq_block}
        alphabet: easel.Alphabet
        Z: database size for E-value calculation

    Returns:
        list of hit dicts
    """
    all_hits = []

    for model_name, frame_names in batch:
        if model_name not in hmms_dict:
            continue

        # Build subset block for these frames
        indices = [name_to_idx[f] for f in frame_names if f in name_to_idx]
        if not indices:
            continue
        subset = [seq_block[i] for i in indices]
        sub_block = easel.DigitalSequenceBlock(alphabet, subset)

        hmm = hmms_dict[model_name]
        for top_hits in pyhmmer.hmmsearch(
            [hmm], sub_block,
            bias_filter=False,
            Z=Z, domZ=Z, E=1e10, domE=1e10,
        ):
            if len(top_hits) == 0:
                continue
            text = tophits_to_domtbl(top_hits, header=False)
            all_hits.extend(parse_domtbl_text(text))

    return all_hits


# ---------------------------------------------------------------------------
# Facet screen
# ---------------------------------------------------------------------------

def facet_screen(hmm_path, seq_block, window_size=64, max_overlap_frac=0.33):
    """Run sub-HMM screen against all frames.

    Returns:
        frame_candidates: {frame_name: [(parent_model, score), ...]}
            sorted by score descending, deduplicated by parent model
    """
    log.info("  Facet screen: building sub-HMMs")
    t0 = time.time()
    sub_hmms = build_sub_hmms_from_file(
        hmm_path, window_size=window_size,
        max_overlap_frac=max_overlap_frac)
    parent_map = {s[0].name: s[1] for s in sub_hmms}
    just_subs = [s[0] for s in sub_hmms]
    t1 = time.time()
    log.info(f"    {len(just_subs)} sub-HMMs built in {t1 - t0:.1f}s")

    log.info("  Facet screen: searching")
    normal, outliers = _partition_hmms_by_size(just_subs)
    sub_Z = len(just_subs)

    raw_hits = []
    if normal:
        raw_hits.extend(_collect_hits(pyhmmer.hmmsearch(
            normal, seq_block,
            bias_filter=False,
            Z=sub_Z, domZ=sub_Z, E=1e10, domE=1e10,
            parallel="queries",
        )))
    if outliers:
        raw_hits.extend(_collect_hits(pyhmmer.hmmsearch(
            outliers, seq_block,
            bias_filter=False,
            Z=sub_Z, domZ=sub_Z, E=1e10, domE=1e10,
            parallel="targets",
        )))

    t2 = time.time()
    log.info(f"    {len(raw_hits)} raw facet hits in {t2 - t1:.1f}s")

    # Build per-frame candidate list: best score per parent model
    frame_raw = defaultdict(dict)
    for hit in raw_hits:
        frame = hit["target_name"]
        parent = parent_map.get(hit["query_name"], hit["query_name"])
        score = hit["dom_score"]
        if parent not in frame_raw[frame] or score > frame_raw[frame][parent]:
            frame_raw[frame][parent] = score

    frame_candidates = {}
    for frame, model_scores in frame_raw.items():
        ranked = sorted(model_scores.items(), key=lambda x: -x[1])
        frame_candidates[frame] = ranked

    log.info(f"    {len(frame_candidates)} frames with candidates")
    return frame_candidates


# ---------------------------------------------------------------------------
# Iterative confirmation with graph-guided expansion
# ---------------------------------------------------------------------------

def iterative_confirm(frame_candidates, graph, hmms_dict, seq_block,
                      seq_lens, alphabet, n_workers=4,
                      max_rounds=5, expand_k=10):
    """Iterative confirmation: facet hits → graph expansion → convergence.

    Args:
        frame_candidates: from facet_screen()
        graph: ModelGraph for neighbor traversal
        hmms_dict: {name: HMM}
        seq_block: full DigitalSequenceBlock
        seq_lens: {frame_name: length}
        alphabet: easel.Alphabet
        n_workers: controls batch granularity
        max_rounds: maximum confirmation rounds
        expand_k: how many graph neighbors to expand per confirmed model

    Returns:
        confirmed_hits: list of all hit dicts from confirmation searches
        unresolved_frames: set of frame names that need legacy fallback
    """
    Z = len(hmms_dict)
    name_to_idx = {seq_block[i].name: i for i in range(len(seq_block))}

    all_frame_names = set(name_to_idx.keys())
    frames_with_candidates = set(frame_candidates.keys())
    resolved_frames = set()
    all_confirmed_hits = []

    # Per-frame state:
    #   searched: models we've run full search on
    #   excluded: models ruled out by graph termination (scores only declining)
    #   queue: remaining candidates to try, ordered by expected promise
    #   best_per_family: best score seen per domain family
    #   status: 'active', 'resolved' (all families converged), 'exhausted' (no more candidates)
    frame_state = {}
    for frame, candidates in frame_candidates.items():
        frame_state[frame] = {
            "searched": set(),
            "excluded": set(),
            "queue": list(candidates),
            "best_per_family": defaultdict(float),
            "confirmed_families": set(),
            "declining_count": defaultdict(int),  # consecutive declines per family
            "status": "active",
        }

    DECLINE_LIMIT = 3  # consecutive score declines before terminating a family

    for round_num in range(1, max_rounds + 1):
        # Build jobs: for each active frame, pick its next model to try
        jobs = defaultdict(list)  # model -> [frames]
        frames_this_round = 0

        for frame, state in frame_state.items():
            if state["status"] != "active":
                continue

            # Find next untried, non-excluded model
            queue = state["queue"]
            next_model = None
            while queue:
                model, _ = queue[0]
                if (model not in state["searched"]
                        and model not in state["excluded"]):
                    next_model = model
                    break
                queue.pop(0)

            if next_model is None:
                state["status"] = "exhausted"
                continue

            jobs[next_model].append(frame)
            frames_this_round += 1

        if not jobs:
            log.info(f"    Round {round_num}: no jobs, converged")
            break

        # Plan and execute batches
        job_list = list(jobs.items())
        n_batches = max(n_workers * 8, 1)
        batches = plan_batches(job_list, hmms_dict, seq_lens, n_batches)

        t0 = time.time()
        round_hits = []
        for batch in batches:
            round_hits.extend(run_batch(
                batch, hmms_dict, seq_block, name_to_idx, alphabet, Z))
        t1 = time.time()

        # Index hits by frame
        hit_by_frame = defaultdict(list)
        for hit in round_hits:
            hit_by_frame[hit["target_name"]].append(hit)

        # Update frame state
        n_confirmed = 0
        for model, frames in jobs.items():
            for frame in frames:
                state = frame_state[frame]
                state["searched"].add(model)

                family = graph.family_of(model)
                frame_hits = [h for h in hit_by_frame.get(frame, [])
                              if h["query_name"] == model]

                if frame_hits:
                    n_confirmed += 1
                    all_confirmed_hits.extend(frame_hits)
                    best_hit_score = max(h["dom_score"] for h in frame_hits)

                    if best_hit_score > state["best_per_family"][family]:
                        # Score improved -> reset decline counter, expand graph
                        state["best_per_family"][family] = best_hit_score
                        state["declining_count"][family] = 0
                        state["confirmed_families"].add(family)

                        for neighbor, _ in graph.neighbors(model, k=expand_k):
                            if (neighbor not in state["searched"]
                                    and neighbor not in state["excluded"]):
                                state["queue"].append(
                                    (neighbor, best_hit_score * 0.9))
                    else:
                        # Score declined for this family
                        state["declining_count"][family] += 1
                        if state["declining_count"][family] >= DECLINE_LIMIT:
                            # Terminate this family's search path
                            # Exclude remaining graph neighbors for this family
                            for neighbor, _ in graph.neighbors(model, k=50):
                                if graph.family_of(neighbor) == family:
                                    state["excluded"].add(neighbor)
                else:
                    # No hit at all -> this model is a dead end
                    # Don't expand its neighbors
                    pass

                # Pop consumed entry from queue
                queue = state["queue"]
                while queue and queue[0][0] in state["searched"]:
                    queue.pop(0)

        # Check for newly resolved frames (all queue exhausted or excluded)
        for frame, state in frame_state.items():
            if state["status"] != "active":
                continue
            queue = state["queue"]
            has_work = any(
                m not in state["searched"] and m not in state["excluded"]
                for m, _ in queue
            )
            if not has_work:
                state["status"] = "resolved"

        n_active = sum(1 for s in frame_state.values()
                       if s["status"] == "active")
        n_resolved = sum(1 for s in frame_state.values()
                         if s["status"] == "resolved")
        n_exhausted = sum(1 for s in frame_state.values()
                          if s["status"] == "exhausted")
        n_models = len(jobs)

        log.info(f"    Round {round_num}: {n_models} models, "
                 f"{frames_this_round} frames, "
                 f"{len(round_hits)} hits, "
                 f"{n_confirmed} confirmed in {t1 - t0:.1f}s "
                 f"[active={n_active} resolved={n_resolved} "
                 f"exhausted={n_exhausted}]")

    # Unresolved = frames with no candidates + active/exhausted frames
    unresolved = all_frame_names - frames_with_candidates
    for frame, state in frame_state.items():
        if state["status"] in ("exhausted",) and not state["confirmed_families"]:
            # Had candidates but never confirmed anything -> legacy
            unresolved.add(frame)

    return all_confirmed_hits, unresolved


# ---------------------------------------------------------------------------
# Full iterative pipeline
# ---------------------------------------------------------------------------

def iterative_search(hmm_path, seq_block, seq_fasta, alphabet,
                     n_workers=4, window_size=64, max_overlap_frac=0.33,
                     max_rounds=5, expand_k=10):
    """Full iterative search: facets → confirm → legacy.

    Args:
        hmm_path: path to HMM database
        seq_block: DigitalSequenceBlock of all frames
        seq_fasta: path to FASTA for seq_lens
        alphabet: easel.Alphabet
        n_workers: worker count for batching
        window_size: facet sub-HMM window
        max_overlap_frac: facet overlap
        max_rounds: max confirmation rounds
        expand_k: graph expansion depth

    Returns:
        list of all hit dicts
    """
    import pyfastx

    hmms = _load_hmms(hmm_path)
    hmms_dict = {h.name: h for h in hmms}
    Z = len(hmms)

    fa = pyfastx.Fasta(seq_fasta, build_index=True)
    seq_lens = {rec.name: len(rec) for rec in fa}

    graph = get_or_build_graph(hmm_path)
    log.info(f"  Graph: {graph}")

    # Step 1: Facet screen
    frame_candidates = facet_screen(
        hmm_path, seq_block, window_size, max_overlap_frac)

    # Step 2: Iterative confirmation
    log.info("  Iterative confirmation:")
    confirmed, unresolved = iterative_confirm(
        frame_candidates, graph, hmms_dict, seq_block,
        seq_lens, alphabet, n_workers, max_rounds, expand_k)

    log.info(f"  Confirmed: {len(confirmed)} hits, "
             f"unresolved: {len(unresolved)} frames")

    # Step 3: Legacy fallback on unresolved
    all_hits = list(confirmed)

    if unresolved:
        log.info(f"  Legacy fallback on {len(unresolved)} frames")
        name_to_idx = {seq_block[i].name: i for i in range(len(seq_block))}
        leftover = [seq_block[name_to_idx[f]] for f in unresolved
                    if f in name_to_idx]

        if leftover:
            leftover_block = easel.DigitalSequenceBlock(alphabet, leftover)
            t0 = time.time()

            normal, outliers = _partition_hmms_by_size(hmms)
            if normal:
                all_hits.extend(_collect_hits(pyhmmer.hmmsearch(
                    normal, leftover_block,
                    bias_filter=False, Z=Z, domZ=Z, E=1e10, domE=1e10,
                    parallel="queries",
                )))
            if outliers:
                all_hits.extend(_collect_hits(pyhmmer.hmmsearch(
                    outliers, leftover_block,
                    bias_filter=False, Z=Z, domZ=Z, E=1e10, domE=1e10,
                    parallel="targets",
                )))
            t1 = time.time()
            log.info(f"    {len(all_hits) - len(confirmed)} legacy hits "
                     f"in {t1 - t0:.1f}s")

    return all_hits
