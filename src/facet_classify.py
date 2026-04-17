"""
facet_classify.py — Facet-based TE classification pipeline.

1. Facet screen: sub-HMM search, all frames, all models
2. Filter: keep facet hits with score >= threshold (default 10)
3. Verify: full-model nobias search for the top facet per (frame, family)
4. Report: verified hits + probable hits with confidence tiers

Confidence tiers based on facet score (empirically determined):
    score >= 20: 100% verification rate
    score >= 15: 99.9% verification rate
    score >= 10: 86.5% verification rate
    score < 10:  unreliable, excluded by default
"""

import logging
import time
from collections import defaultdict

import pyhmmer
import pyhmmer.easel as easel

from decompose_hmm import build_sub_hmms_from_file, build_sub_hmms_tiered, _load_hmms
from model_graph import get_or_build_graph
from search import (tophits_to_domtbl, parse_domtbl_text, _collect_hits,
                    _partition_hmms_by_size)

log = logging.getLogger(__name__)

# Confidence tiers per facet size (empirically calibrated)
# AA models (high information per position)
AA_CONF_TIERS = {
    96: [(30.0, "confirmed_100"), (20.0, "probable_99.9"), (15.0, "probable_99"), (10.0, "probable_91")],
    64: [(15.0, "confirmed_100"), (10.0, "probable_97")],
    48: [(15.0, "confirmed_100"), (10.0, "probable_97")],  # limited data, use 64's
    32: [(10.0, "probable_97")],  # limited data
}

# DNA models (low information per position, scores are much lower)
DNA_CONF_TIERS = {
    192: [(20.0, "confirmed_100"), (0.0, "probable_90"), (-20.0, "probable_87")],
    128: [(50.0, "confirmed_100"), (-5.0, "probable_87")],
    64:  [(10.0, "confirmed_100"), (5.0, "probable_70")],
}

CONF_TIERS_DEFAULT = [(20.0, "confirmed_100"), (15.0, "probable_99"), (10.0, "probable_90")]


def _get_confidence(score, facet_M, is_dna=False):
    """Get confidence tier for a facet hit based on score and facet size."""
    tier_table = DNA_CONF_TIERS if is_dna else AA_CONF_TIERS

    tiers = CONF_TIERS_DEFAULT
    for size in sorted(tier_table.keys(), reverse=True):
        if facet_M >= size:
            tiers = tier_table[size]
            break

    for threshold, tier_name in tiers:
        if score >= threshold:
            return tier_name
    return "low"


def _get_family(model_name):
    if ":" in model_name:
        return model_name.split(":")[1].split("-")[-1]
    return model_name.split("_")[0]


def facet_screen(hmm_path, seq_block, window_size=64, max_overlap_frac=0.33,
                 checkpoint_path=None):
    """Run sub-HMM screen. Returns {frame: {parent_model: best_score}}."""
    import json
    import os

    if checkpoint_path and os.path.exists(checkpoint_path):
        log.info(f"  Loading facet checkpoint: {checkpoint_path}")
        t0 = time.time()
        with open(checkpoint_path) as f:
            raw = json.load(f)
        # Convert back from JSON lists
        result = {k: dict(v) for k, v in raw.items()}
        log.info(f"    {len(result)} frames in {time.time() - t0:.1f}s")
        return result

    log.info("  Building sub-HMMs (tiered)")
    t0 = time.time()
    sub_hmms = build_sub_hmms_tiered(
        hmm_path, max_overlap_frac=max_overlap_frac)
    parent_map = {s[0].name: s[1] for s in sub_hmms}
    # Sort largest first for optimal thread utilization
    sub_hmms.sort(key=lambda s: -s[0].M)
    just_subs = [s[0] for s in sub_hmms]
    t1 = time.time()
    log.info(f"    {len(just_subs)} sub-HMMs in {t1 - t0:.1f}s")

    log.info("  Searching facets")
    raw_hits = _collect_hits(pyhmmer.hmmsearch(
        just_subs, seq_block,
        bias_filter=False,
        Z=len(just_subs), domZ=len(just_subs),
        E=1e10, domE=1e10,
        parallel="queries",
    ))
    t2 = time.time()
    log.info(f"    {len(raw_hits)} raw hits in {t2 - t1:.1f}s")

    # Build facet size lookup
    facet_sizes = {s[0].name: s[0].M for s in sub_hmms}

    # Best facet score per (frame, parent_model), with facet size
    result = defaultdict(dict)  # frame -> {parent: (score, facet_M)}
    for hit in raw_hits:
        frame = hit["target_name"]
        sub_name = hit["query_name"]
        parent = parent_map.get(sub_name, sub_name)
        score = hit["dom_score"]
        fM = facet_sizes.get(sub_name, hit["query_len"])
        if parent not in result[frame] or score > result[frame][parent][0]:
            result[frame][parent] = (score, fM)

    result = dict(result)
    log.info(f"    {len(result)} frames with signal")

    if checkpoint_path:
        with open(checkpoint_path, "w") as f:
            json.dump({k: list(v.items()) for k, v in result.items()}, f)
        log.info(f"    Saved checkpoint: {checkpoint_path}")

    return result


def classify_frames(facet_hits, hmms_dict, seq_block, alphabet,
                    min_facet_score=10.0, n_workers=4, is_dna=False):
    """
    Classify frames using facet scores + targeted full-model verification.

    Args:
        facet_hits: {frame: {model: score}} from facet_screen
        hmms_dict: {model_name: HMM}
        seq_block: DigitalSequenceBlock
        alphabet: easel.Alphabet
        min_facet_score: minimum facet score to consider
        n_workers: for future parallel batching

    Returns:
        classifications: list of dicts, one per (frame, family) assignment:
            {frame, family, model, facet_score, confidence,
             verified, full_score (if verified)}
        unclassified_frames: set of frame names with no facet signal
    """
    Z = len(hmms_dict)
    name_to_idx = {seq_block[i].name: i for i in range(len(seq_block))}

    all_frame_names = set(name_to_idx.keys())
    frames_with_signal = set(facet_hits.keys())
    unclassified = all_frame_names - frames_with_signal

    # Step 1: Find best facet per (frame, family), filtered by min score
    top_per_family = {}  # (frame, family) -> (model, score, facet_M)
    other_hits = defaultdict(list)  # (frame, family) -> [(model, score, facet_M), ...]

    for frame, model_scores in facet_hits.items():
        by_family = defaultdict(list)
        for model, (score, fM) in model_scores.items():
            if score < min_facet_score:
                continue
            family = _get_family(model)
            by_family[family].append((model, score, fM))

        for family, candidates in by_family.items():
            candidates.sort(key=lambda x: -x[1])
            top_per_family[(frame, family)] = candidates[0]
            if len(candidates) > 1:
                other_hits[(frame, family)] = candidates[1:]

    log.info(f"  {len(top_per_family)} (frame, family) assignments to verify")

    # Step 2: Verify top picks with full-model nobias search
    # Group frames by their top model for batched searching
    model_frames = defaultdict(list)
    for (frame, family), (model, score, fM) in top_per_family.items():
        model_frames[model].append((frame, family))

    log.info(f"  Verifying against {len(model_frames)} unique models")
    t0 = time.time()

    verified_results = {}  # (frame, model) -> full_score

    # Batch: collect all models and their frames, search in bulk
    for model_name, frame_families in model_frames.items():
        if model_name not in hmms_dict:
            continue

        frames = [ff[0] for ff in frame_families]
        indices = [name_to_idx[f] for f in frames if f in name_to_idx]
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
            for hit in parse_domtbl_text(text):
                verified_results[(hit["target_name"], model_name)] = hit

    t1 = time.time()
    log.info(f"  Verification: {len(verified_results)} hits in {t1 - t0:.1f}s")

    # Step 3: Build classification output
    classifications = []

    for (frame, family), (model, facet_score, facet_M) in top_per_family.items():
        confidence = _get_confidence(facet_score, facet_M, is_dna)

        entry = {
            "frame": frame,
            "family": family,
            "model": model,
            "facet_score": facet_score,
            "facet_M": facet_M,
            "confidence": confidence,
            "verified": False,
            "full_score": None,
            "full_hit": None,
        }

        # Check verification
        vkey = (frame, model)
        if vkey in verified_results:
            entry["verified"] = True
            entry["full_score"] = verified_results[vkey]["dom_score"]
            entry["full_hit"] = verified_results[vkey]

        classifications.append(entry)

        # Add secondary hits as unverified probable assignments
        for alt_model, alt_score, alt_fM in other_hits.get((frame, family), []):
            classifications.append({
                "frame": frame,
                "family": family,
                "model": alt_model,
                "facet_score": alt_score,
                "facet_M": alt_fM,
                "confidence": _get_confidence(alt_score, alt_fM),
                "verified": False,
                "full_score": None,
                "full_hit": None,
                "is_secondary": True,
            })

    n_verified = sum(1 for c in classifications
                     if c.get("verified") and not c.get("is_secondary"))
    n_primary = sum(1 for c in classifications if not c.get("is_secondary"))
    log.info(f"  {n_primary} primary assignments, "
             f"{n_verified} verified ({100*n_verified/max(n_primary,1):.1f}%)")

    return classifications, unclassified


def facet_classify(hmm_path, seq_block, seq_fasta, alphabet,
                   n_workers=4, window_size=64, max_overlap_frac=0.33,
                   min_facet_score=10.0, checkpoint_dir=None):
    """
    Full facet classification pipeline.

    Returns:
        classifications: list of classification dicts
        legacy_hits: list of hit dicts from legacy fallback on unclassified frames
    """
    import os
    import pyfastx

    hmms = _load_hmms(hmm_path)
    hmms_dict = {h.name: h for h in hmms}
    Z = len(hmms)

    # Detect DNA vs AA
    is_dna = hmms[0].alphabet == easel.Alphabet.dna() if hmms else False

    # DNA models use lower score thresholds
    if is_dna and min_facet_score > 0:
        min_facet_score = -20.0
        log.info(f"  DNA database detected, using min_facet_score={min_facet_score}")

    # Step 1: Facet screen
    ckpt = None
    if checkpoint_dir:
        db_base = os.path.basename(hmm_path)
        ckpt = os.path.join(checkpoint_dir, f"{db_base}.facets.json")

    log.info("  --- Facet classification ---")
    facet_hits = facet_screen(
        hmm_path, seq_block, window_size, max_overlap_frac,
        checkpoint_path=ckpt)

    # Step 2: Classify
    classifications, unclassified = classify_frames(
        facet_hits, hmms_dict, seq_block, alphabet,
        min_facet_score=min_facet_score, n_workers=n_workers,
        is_dna=is_dna)

    # Step 3: Legacy fallback on unclassified
    legacy_hits = []
    if unclassified:
        log.info(f"  Legacy fallback on {len(unclassified)} unclassified frames")
        name_to_idx = {seq_block[i].name: i for i in range(len(seq_block))}
        leftover = [seq_block[name_to_idx[f]] for f in unclassified
                    if f in name_to_idx]

        if leftover:
            t0 = time.time()
            leftover_block = easel.DigitalSequenceBlock(alphabet, leftover)
            normal, outliers = _partition_hmms_by_size(hmms)

            if normal:
                legacy_hits.extend(_collect_hits(pyhmmer.hmmsearch(
                    normal, leftover_block,
                    bias_filter=False, Z=Z, domZ=Z, E=1e10, domE=1e10,
                    parallel="queries",
                )))
            if outliers:
                legacy_hits.extend(_collect_hits(pyhmmer.hmmsearch(
                    outliers, leftover_block,
                    bias_filter=False, Z=Z, domZ=Z, E=1e10, domE=1e10,
                    parallel="targets",
                )))
            t1 = time.time()
            log.info(f"    {len(legacy_hits)} legacy hits in {t1 - t0:.1f}s")

    return classifications, legacy_hits


def export_classifications_tsv(classifications, out_path):
    """Export classifications to TSV."""
    columns = [
        "frame", "family", "model", "facet_score", "facet_M",
        "confidence", "verified", "full_score", "is_secondary",
    ]

    with open(out_path, "w") as f:
        f.write("\t".join(columns) + "\n")
        for c in classifications:
            vals = [
                c["frame"],
                c["family"],
                c["model"],
                f"{c['facet_score']:.1f}",
                c.get("facet_M", ""),
                c["confidence"],
                "yes" if c.get("verified") else "no",
                f"{c['full_score']:.1f}" if c.get("full_score") is not None else "",
                "yes" if c.get("is_secondary") else "no",
            ]
            f.write("\t".join(str(v) for v in vals) + "\n")
