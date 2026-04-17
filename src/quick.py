"""
Quick mode: sub-HMM triage + single-model confirmation + legacy fallback.

1. Screen all frames with spliced sub-HMMs → assign each frame to its
   highest-scoring sub-HMM's parent model
2. Confirm each assignment with one nobias full-model search
3. Unassigned leftovers get full legacy search

Trades a small amount of recall for major speed improvement on the
confirmation step: most frames only need 1 full-model search instead
of 266.
"""

import logging
import time
from collections import defaultdict

import pyhmmer
import pyhmmer.easel as easel

from decompose_hmm import build_sub_hmms_from_file, _load_hmms
from search import (tophits_to_domtbl, parse_domtbl_text, _collect_hits,
                    _partition_hmms_by_size, build_sequence_block)

log = logging.getLogger(__name__)


def _search_models_against_block(hmms, seq_block, Z, bias_filter=False):
    """Run hmmsearch with IQR partitioning, return hit list."""
    normal, outliers = _partition_hmms_by_size(hmms)
    all_hits = []

    if normal:
        all_hits.extend(_collect_hits(pyhmmer.hmmsearch(
            normal, seq_block,
            bias_filter=bias_filter,
            Z=Z, domZ=Z, E=1e10, domE=1e10,
            parallel="queries",
        )))

    if outliers:
        all_hits.extend(_collect_hits(pyhmmer.hmmsearch(
            outliers, seq_block,
            bias_filter=bias_filter,
            Z=Z, domZ=Z, E=1e10, domE=1e10,
            parallel="targets",
        )))

    return all_hits


def quick_search(hmm_path, seq_block, seq_fasta, alphabet,
                 window_size=64, max_overlap_frac=0.33):
    """
    Quick mode search: sub-HMM triage → confirm → legacy fallback.

    Args:
        hmm_path: path to HMM database file
        seq_block: DigitalSequenceBlock of all frames
        seq_fasta: path to the FASTA (for building subset blocks)
        alphabet: easel.Alphabet
        window_size: sub-HMM window size
        max_overlap_frac: sub-HMM overlap fraction

    Returns:
        all_hits: list of hit dicts (combined from all tiers)
    """
    hmms = _load_hmms(hmm_path)
    hmms_dict = {h.name: h for h in hmms}
    Z = len(hmms)

    # --- Tier 1: Sub-HMM screen ---
    log.info("    Tier 1: sub-HMM screen")
    t0 = time.time()
    sub_hmms = build_sub_hmms_from_file(hmm_path, window_size=window_size,
                                         max_overlap_frac=max_overlap_frac)
    parent_map = {s[0].name: s[1] for s in sub_hmms}
    just_subs = [s[0] for s in sub_hmms]
    t1 = time.time()
    log.info(f"      Built {len(just_subs)} sub-HMMs in {t1 - t0:.1f}s")

    # Search sub-HMMs against all frames
    sub_Z = len(just_subs)
    normal_subs, outlier_subs = _partition_hmms_by_size(just_subs)

    sub_hits = []
    if normal_subs:
        sub_hits.extend(_collect_hits(pyhmmer.hmmsearch(
            normal_subs, seq_block,
            bias_filter=False,
            Z=sub_Z, domZ=sub_Z, E=1e10, domE=1e10,
            parallel="queries",
        )))
    if outlier_subs:
        sub_hits.extend(_collect_hits(pyhmmer.hmmsearch(
            outlier_subs, seq_block,
            bias_filter=False,
            Z=sub_Z, domZ=sub_Z, E=1e10, domE=1e10,
            parallel="targets",
        )))

    t2 = time.time()
    log.info(f"      {len(sub_hits)} sub-HMM hits in {t2 - t1:.1f}s")

    # Rank sub-HMM hits per frame by score descending
    frame_candidates = defaultdict(list)  # frame -> [(parent_model, score)]
    for hit in sub_hits:
        frame = hit["target_name"]
        parent = parent_map.get(hit["query_name"], hit["query_name"])
        score = hit["dom_score"]
        frame_candidates[frame].append((parent, score))

    # Deduplicate and sort: unique models per frame, best score first
    for frame in frame_candidates:
        seen = {}
        for model, score in frame_candidates[frame]:
            if model not in seen or score > seen[model]:
                seen[model] = score
        frame_candidates[frame] = sorted(seen.items(), key=lambda x: -x[1])

    all_frame_names = set()
    for i in range(len(seq_block)):
        all_frame_names.add(seq_block[i].name)
    frames_with_candidates = set(frame_candidates.keys())
    unassigned_frames = all_frame_names - frames_with_candidates

    log.info(f"      {len(frames_with_candidates)} frames with candidates, "
             f"{len(unassigned_frames)} with no sub-HMM signal")

    # --- Tier 2: Confirm by searching candidates in score order ---
    # For each model, collect the frames that want to try it.
    # Process models from most-requested to least. A frame drops out
    # as soon as any model confirms it.
    log.info("    Tier 2: confirming assignments (score-ordered)")
    t3 = time.time()

    name_to_idx = {seq_block[i].name: i for i in range(len(seq_block))}

    # Build priority queue: which model each frame wants to try next
    frame_queue_idx = {f: 0 for f in frames_with_candidates}
    confirmed_frames = set()
    confirmed_hits = []
    models_searched = set()

    while True:
        # Collect: for each unconfirmed frame, what model does it want next?
        model_frames = defaultdict(list)
        for frame in frames_with_candidates - confirmed_frames:
            idx = frame_queue_idx[frame]
            candidates = frame_candidates[frame]
            if idx >= len(candidates):
                # Exhausted all candidates -> goes to legacy
                unassigned_frames.add(frame)
                confirmed_frames.add(frame)  # remove from loop
                continue
            model, _ = candidates[idx]
            model_frames[model].append(frame)

        if not model_frames:
            break

        # Search each needed model against its requesting frames
        for model_name, frames in model_frames.items():
            if model_name not in hmms_dict:
                for f in frames:
                    frame_queue_idx[f] += 1
                continue

            # Only search frames not yet confirmed
            active_frames = [f for f in frames if f not in confirmed_frames]
            if not active_frames:
                continue

            subset = [seq_block[name_to_idx[f]] for f in active_frames
                      if f in name_to_idx]
            if not subset:
                continue
            sub_block = easel.DigitalSequenceBlock(alphabet, subset)

            hmm = hmms_dict[model_name]
            models_searched.add(model_name)

            for top_hits in pyhmmer.hmmsearch(
                [hmm], sub_block,
                bias_filter=False,
                Z=Z, domZ=Z, E=1e10, domE=1e10,
            ):
                if len(top_hits) == 0:
                    continue
                text = tophits_to_domtbl(top_hits, header=False)
                hits = parse_domtbl_text(text)
                hit_names = set(h["target_name"] for h in hits)
                confirmed_hits.extend(hits)

                for f in active_frames:
                    if f in hit_names:
                        confirmed_frames.add(f)
                    else:
                        frame_queue_idx[f] += 1

    t4 = time.time()
    # Frames that exhausted candidates without confirmation
    for frame in frames_with_candidates - confirmed_frames:
        unassigned_frames.add(frame)

    log.info(f"      {len(confirmed_hits)} confirmed hits, "
             f"{len(models_searched)} models searched in {t4 - t3:.1f}s")
    log.info(f"      {len(unassigned_frames)} frames to legacy fallback")

    # --- Tier 3: Legacy search on unassigned frames ---
    all_hits = list(confirmed_hits)

    if unassigned_frames:
        log.info("    Tier 3: legacy search on leftovers")
        t5 = time.time()

        # Build subset block for unassigned frames
        leftover = [seq_block[name_to_idx[f]] for f in unassigned_frames
                    if f in name_to_idx]
        if leftover:
            leftover_block = easel.DigitalSequenceBlock(alphabet, leftover)
            legacy_hits = _search_models_against_block(
                hmms, leftover_block, Z, bias_filter=False)
            all_hits.extend(legacy_hits)
            t6 = time.time()
            log.info(f"      {len(legacy_hits)} legacy hits in {t6 - t5:.1f}s")
    else:
        log.info("    Tier 3: no leftovers, skipped")

    return all_hits
