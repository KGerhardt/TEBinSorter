"""
cross_family.py — Complete missing family assignments for facet-classified frames.

After facet classification, some frames have hits for some families (e.g. RT, INT)
but not others (e.g. GAG, TPase). This module searches the missing families'
models against those frames to close the gap.

Only searches the specific (model, frame) pairs that are needed -- no redundant work.
"""

import logging
import time
from collections import defaultdict

import pyhmmer
import pyhmmer.easel as easel

from search import tophits_to_domtbl, parse_domtbl_text

log = logging.getLogger(__name__)


def _get_family(model_name):
    if ":" in model_name:
        return model_name.split(":")[1].split("-")[-1]
    return model_name.split("_")[0]


def find_missing_families(classifications, hmms_dict):
    """Identify which families each classified frame is missing.

    Args:
        classifications: list of classification dicts from facet_classify
        hmms_dict: {model_name: HMM} for all models

    Returns:
        missing: dict of {frame_name: set(family_names)} that need searching
        all_families: set of all known family names
    """
    # All families in the database
    all_families = set(_get_family(name) for name in hmms_dict)

    # Which families did each frame find?
    frame_found = defaultdict(set)
    classified_frames = set()
    for c in classifications:
        if c.get("is_secondary"):
            continue
        frame_found[c["frame"]].add(c["family"])
        classified_frames.add(c["frame"])

    # Missing = all_families - found for each classified frame
    missing = {}
    for frame in classified_frames:
        missed = all_families - frame_found[frame]
        if missed:
            missing[frame] = missed

    return missing, all_families


def search_missing(missing, hmms_dict, seq_block, alphabet, Z=None):
    """Search missing families for classified frames.

    Groups by model for efficient batching: each model searches only
    the frames that need it.

    Args:
        missing: {frame: set(families)} from find_missing_families
        hmms_dict: {model_name: HMM}
        seq_block: DigitalSequenceBlock
        alphabet: easel.Alphabet
        Z: database size for E-values (default: len(hmms_dict))

    Returns:
        list of hit dicts
    """
    if not missing:
        return []

    if Z is None:
        Z = len(hmms_dict)

    name_to_idx = {seq_block[i].name: i for i in range(len(seq_block))}

    # Invert: for each model, which frames need it?
    # A frame needs a model if the model's family is in the frame's missing set
    model_frames = defaultdict(set)
    for frame, missed_families in missing.items():
        for model_name, hmm in hmms_dict.items():
            if _get_family(model_name) in missed_families:
                model_frames[model_name].add(frame)

    # Batch: group models that share similar frame sets
    # Simple approach: just search each model against its frames
    # But batch multiple models into one hmmsearch call for efficiency
    log.info(f"    {len(missing)} frames missing families, "
             f"{len(model_frames)} models to search")

    t0 = time.time()

    # Group models by their frame set size for rough batching
    # Sort largest frame sets first (most work first)
    sorted_models = sorted(model_frames.items(), key=lambda x: -len(x[1]))

    # Batch: collect models into groups, each group searches a combined frame set
    BATCH_SIZE = 20  # models per hmmsearch call
    all_hits = []
    batch_models = []
    batch_frames = set()

    for model_name, frames in sorted_models:
        batch_models.append(model_name)
        batch_frames |= frames

        if len(batch_models) >= BATCH_SIZE:
            hits = _run_batch(batch_models, batch_frames, hmms_dict,
                              seq_block, name_to_idx, alphabet, Z,
                              model_frames)
            all_hits.extend(hits)
            batch_models = []
            batch_frames = set()

    # Flush remaining
    if batch_models:
        hits = _run_batch(batch_models, batch_frames, hmms_dict,
                          seq_block, name_to_idx, alphabet, Z,
                          model_frames)
        all_hits.extend(hits)

    t1 = time.time()
    log.info(f"    {len(all_hits)} cross-family hits in {t1 - t0:.1f}s")

    return all_hits


def _run_batch(model_names, frame_names, hmms_dict, seq_block,
               name_to_idx, alphabet, Z, model_frames):
    """Run one batch of models against their needed frames."""
    hmms = [hmms_dict[m] for m in model_names if m in hmms_dict]
    if not hmms:
        return []

    # Build subset block
    indices = [name_to_idx[f] for f in frame_names if f in name_to_idx]
    if not indices:
        return []

    subset = [seq_block[i] for i in indices]
    sub_block = easel.DigitalSequenceBlock(alphabet, subset)

    raw_hits = []
    for top_hits in pyhmmer.hmmsearch(
        hmms, sub_block,
        bias_filter=False,
        Z=Z, domZ=Z, E=1e10, domE=1e10,
    ):
        if len(top_hits) == 0:
            continue
        text = tophits_to_domtbl(top_hits, header=False)
        for hit in parse_domtbl_text(text):
            # Only keep hits for frames that actually need this model's family
            model = hit["query_name"]
            frame = hit["target_name"]
            if frame in model_frames.get(model, set()):
                raw_hits.append(hit)

    return raw_hits
