"""
classifier.py — Config-driven TE classification from HMM domain hits.

Replicates TEsorter's classification logic: domain deconfliction (hmm2best),
order/superfamily/clade assignment, completeness checking. Parameterized
by per-database configs rather than hardcoded per-database classes.

Two stages:
  1. hmm2best: select best domain hit per (sequence, domain_type),
     with database-specific remapping and overlap rules
  2. classify: assign order/superfamily/clade from domain architecture
"""

import logging
import numpy as np
from collections import Counter, defaultdict

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Database configs
# ---------------------------------------------------------------------------

REXDB_CONFIG = {
    "name": "rexdb",
    "domain_remap": {"aRH": "RH", "TPase": "INT"},
    "overlap_aware": True,
    "structures": {
        ("LTR", "Copia"):   ["GAG", "PROT", "INT", "RT", "RH"],
        ("LTR", "Gypsy"):   ["GAG", "PROT", "RT", "RH", "INT"],
        ("LTR", "Bel-Pao"): ["GAG", "PROT", "RT", "RH", "INT"],
    },
    "clade_parser": "rexdb",
    "clade_restrict": {"Copia", "Gypsy"},  # only these get clade names
}

GYDB_CONFIG = {
    "name": "gydb",
    "domain_remap": {},
    "overlap_aware": False,
    "structures": {
        ("LTR", "Copia"):          ["GAG", "AP", "INT", "RT", "RNaseH"],
        ("LTR", "Gypsy"):          ["GAG", "AP", "RT", "RNaseH", "INT"],
        ("LTR", "Pao"):            ["GAG", "AP", "RT", "RNaseH", "INT"],
        ("LTR", "Retroviridae"):   ["GAG", "AP", "RT", "RNaseH", "INT", "ENV"],
        ("LTR", "Caulimoviridae"): ["GAG", "AP", "RT", "RNaseH"],
    },
    "clade_parser": "gydb",
    "clade_restrict": None,
}

LINE_CONFIG = {
    "name": "line",
    "domain_remap": {},
    "overlap_aware": False,
    "structures": {},
    "clade_parser": "rexdb",
    "clade_restrict": None,
}

TIR_CONFIG = {
    "name": "tir",
    "domain_remap": {},
    "overlap_aware": False,
    "structures": {},
    "clade_parser": "rexdb",
    "clade_restrict": None,
}

SINE_CONFIG = {
    "name": "sine",
    "domain_remap": {},
    "overlap_aware": False,
    "structures": {},
    "clade_parser": "sine",
    "clade_restrict": None,
}

DB_CONFIGS = {
    "rexdb": REXDB_CONFIG,
    "gydb": GYDB_CONFIG,
    "line": LINE_CONFIG,
    "tir": TIR_CONFIG,
    "sine": SINE_CONFIG,
    "sine-so": SINE_CONFIG,
}


# ---------------------------------------------------------------------------
# Clade parsing
# ---------------------------------------------------------------------------

def parse_clade_rexdb(model_name):
    """Parse REXdb model name -> (order, superfamily, clade, gene).

    Format: Class_I/LTR/Ty1_copia/SIRE:Ty1-RT
    """
    if ":" not in model_name:
        return "Unknown", "unknown", model_name, model_name

    clade_path, domain = model_name.split(":", 1)
    gene = domain.split("-")[-1]

    if clade_path.startswith("Class_I/LTR/Ty1_copia"):
        order, superfamily = "LTR", "Copia"
    elif clade_path.startswith("Class_I/LTR/Ty3_gypsy"):
        order, superfamily = "LTR", "Gypsy"
    elif clade_path.startswith("Class_I/LTR/"):
        parts = clade_path.split("/")
        order, superfamily = parts[1], parts[2] if len(parts) > 2 else "unknown"
    elif clade_path.startswith("Class_I/"):
        parts = clade_path.split("/")
        order = parts[1]
        superfamily = parts[2] if len(parts) > 2 else "unknown"
    elif clade_path.startswith("Class_II/"):
        parts = clade_path.split("/")
        order = parts[2] if len(parts) > 2 else "unknown"
        superfamily = parts[3] if len(parts) > 3 else "unknown"
    elif clade_path.startswith("NA"):
        order, superfamily = "LTR", "Retrovirus"
    else:
        order, superfamily = "Unknown", "unknown"

    clade = clade_path.split("/")[-1] if "/" in clade_path else clade_path

    return order, superfamily, clade, gene


def parse_clade_gydb(model_name):
    """Parse GyDB model name -> (gene, clade).

    Format: AP_copia, RT_gypsy, GAG_lentiviridae
    """
    parts = model_name.split("_", 1)
    gene = parts[0]
    clade = parts[1] if len(parts) > 1 else model_name
    return gene, clade


# ---------------------------------------------------------------------------
# hmm2best: domain deconfliction
# ---------------------------------------------------------------------------

def hmm2best(hits, config, compat_rounding=False):
    """Select best domain hits per (base_seq, domain_type).

    Args:
        hits: dict from deconflict.load_hits()
        config: database config dict
        compat_rounding: if True, round norm_score to 2 decimal places
                         (replicates TEsorter rounding bug)

    Returns:
        numpy index array of selected hits
    """
    n = len(hits["score"])
    if n == 0:
        return np.array([], dtype=int)

    remap = config["domain_remap"]
    overlap_aware = config["overlap_aware"]

    # Apply domain remapping
    domain_types = hits["family"].copy()
    for old, new in remap.items():
        mask = domain_types == old
        domain_types[mask] = new

    # Normalized score
    norm_score = hits["score"] / hits["model_len"]
    if compat_rounding:
        norm_score = np.round(norm_score, 2)

    # Group by (base_seq, domain_type)
    # Process in score-descending order within each group
    best = {}  # (base_seq, domain_type) -> index

    sort_idx = np.argsort(-norm_score)

    for i in sort_idx:
        base = hits["base_seq"][i]
        dtype = domain_types[i]
        key = (base, dtype)

        if key not in best:
            best[key] = i
            continue

        if not overlap_aware:
            # Simple: already have a better score, skip
            continue

        # Overlap-aware logic (REXdb)
        curr = best[key]
        curr_score = norm_score[curr]
        new_score = norm_score[i]

        if new_score <= curr_score:
            continue

        # Only replace if same gene OR overlapping envelopes
        curr_model = hits["model"][curr]
        new_model = hits["model"][i]

        # Same gene = same domain prefix (e.g. both Ty1-RT)
        if ":" in curr_model and ":" in new_model:
            curr_gene = curr_model.split(":")[1]
            new_gene = new_model.split(":")[1]
            same_gene = curr_gene == new_gene
        else:
            same_gene = curr_model.split("_")[0] == new_model.split("_")[0]

        if same_gene:
            best[key] = i
            continue

        # Check envelope overlap
        curr_start, curr_end = hits["env_from"][curr], hits["env_to"][curr]
        new_start, new_end = hits["env_from"][i], hits["env_to"][i]
        if new_start <= curr_end and new_end >= curr_start:
            best[key] = i

    return np.array(list(best.values()), dtype=int)


def apply_filters(hits, indices, min_cov=20.0, max_evalue=1e-3,
                  min_acc=0.5, min_norm_score=0.1, compat_rounding=False):
    """Apply filters to selected hits.

    Uses full precision by default. --compat-tesorter-rounding rounds
    norm_score to 2 decimal places before threshold comparison.
    """
    idx = np.array(indices)
    cov = hits["hmm_cov"][idx]
    evalue = hits["evalue"][idx]
    acc = hits["acc"][idx]
    nscore = hits["score"][idx] / hits["model_len"][idx]
    if compat_rounding:
        nscore = np.round(nscore, 2)

    mask = (cov >= min_cov) & (evalue <= max_evalue) & (acc >= min_acc) & (nscore >= min_norm_score)
    return idx[mask]


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def classify_element(genes, clades, models, config):
    """Classify a single TE element from its domain hits.

    Args:
        genes: list of gene names (domain types) in positional order
        clades: list of clade names corresponding to each gene
        models: list of full model names
        config: database config dict

    Returns:
        (order, superfamily, clade, complete)
    """
    parser = config["clade_parser"]

    if parser == "rexdb":
        return _classify_rexdb(genes, clades, models, config)
    elif parser == "gydb":
        return _classify_gydb(genes, clades, models, config)
    elif parser == "sine":
        return "SINE", "unknown", "unknown", "unknown"
    else:
        return "Unknown", "unknown", "unknown", "unknown"


def _classify_rexdb(genes, clades, models, config):
    """REXdb classification logic."""
    clade_count = Counter(clades)
    max_clade = max(clade_count, key=lambda x: clade_count[x])

    order, superfamily, _, _ = parse_clade_rexdb(
        [m for m, c in zip(models, clades) if c == max_clade][0])

    counts = list(clade_count.values())
    if len(clade_count) == 1 or (clade_count[max_clade] > 1 and
                                  len(counts) > 1 and counts[0] > counts[1]):
        display_clade = max_clade.split("/")[-1]
    elif len(clade_count) > 1:
        display_clade = "mixture"
        superfamilies = [parse_clade_rexdb(m)[1] for m in models]
        if len(Counter(superfamilies)) > 1:
            superfamily = "mixture"
            orders = [parse_clade_rexdb(m)[0] for m in models]
            if len(Counter(orders)) > 1:
                order = "mixture"
    else:
        display_clade = max_clade.split("/")[-1]

    # Check completeness
    structures = config["structures"]
    try:
        expected = structures[(order, superfamily)]
        present = [g for g in genes if g in set(expected)]
        complete = "yes" if expected == present else "no"
    except KeyError:
        complete = "unknown"

    # Restrict clade names
    restrict = config.get("clade_restrict")
    if restrict and superfamily not in restrict:
        display_clade = "unknown"
    if display_clade.startswith("Ty"):
        display_clade = "unknown"

    return order, superfamily, display_clade, complete


def _classify_gydb(genes, clades, models, config):
    """GyDB classification logic. Requires clade_map."""
    clade_count = Counter(clades)
    max_clade = max(clade_count, key=lambda x: clade_count[x])

    # Look up order/superfamily from clade map
    clade_map = config.get("_clade_map", {})
    order, superfamily = clade_map.get(max_clade, ("Unknown", "unknown"))

    if len(clade_count) == 1 or clade_count[max_clade] > 1:
        display_clade = max_clade
    elif len(clade_count) > 1:
        display_clade = "mixture"
        superfamilies = [clade_map.get(c, [None, None])[1] for c in clades]
        if len(Counter(superfamilies)) > 1:
            superfamily = "mixture"
            orders = [clade_map.get(c, [None, None])[0] for c in clades]
            if len(Counter(orders)) > 1:
                order = "mixture"
    else:
        display_clade = max_clade

    structures = config["structures"]
    try:
        expected = structures[(order, superfamily)]
        present = [g for g in genes if g in set(expected)]
        complete = "yes" if expected == present else "no"
    except KeyError:
        complete = "unknown"

    return order, superfamily, display_clade, complete


# ---------------------------------------------------------------------------
# Full classification pipeline
# ---------------------------------------------------------------------------

def classify_sequences(hits, config, gydb_clade_map=None, compat_rounding=False):
    """Full classification: hmm2best -> filter -> classify per sequence.

    Args:
        hits: dict from deconflict.load_hits()
        config: database config dict
        gydb_clade_map: optional {clade: (order, superfamily)} for GyDB

    Returns:
        list of dicts with keys: id, order, superfamily, clade, complete,
        strand, domains
    """
    if gydb_clade_map and config["clade_parser"] == "gydb":
        config = dict(config)
        config["_clade_map"] = gydb_clade_map

    # Step 1: hmm2best
    best_idx = hmm2best(hits, config, compat_rounding=compat_rounding)
    log.info(f"  hmm2best: {len(best_idx)} domain assignments")

    # Step 2: filter
    filtered_idx = apply_filters(hits, best_idx, compat_rounding=compat_rounding)
    log.info(f"  After filter: {len(filtered_idx)} assignments")

    # Step 3: group by base_seq and classify
    seq_domains = defaultdict(list)
    for i in filtered_idx:
        base = hits["base_seq"][i]
        seq_domains[base].append(i)

    results = []
    for base_seq, indices in seq_domains.items():
        # Determine strand
        # Use frame info from target names
        strands = []
        for i in indices:
            target = hits["target"][i]
            if "|fwd" in target:
                strands.append("+")
            elif "|rev" in target:
                strands.append("-")
            else:
                strands.append(".")

        unique_strands = set(strands)
        if len(unique_strands) > 1:
            strand = "?"
        elif len(unique_strands) == 1:
            strand = strands[0]
        else:
            continue

        # Sort by position (env_from)
        sorted_indices = sorted(indices, key=lambda i: hits["env_from"][i])
        if strand == "-":
            sorted_indices.reverse()

        # Extract genes, clades, models
        remap = config["domain_remap"]
        genes = []
        clades = []
        models = []
        domain_strs = []

        for i in sorted_indices:
            model = hits["model"][i]
            if config["clade_parser"] == "rexdb":
                _, _, clade, gene = parse_clade_rexdb(model)
            elif config["clade_parser"] == "gydb":
                gene, clade = parse_clade_gydb(model)
            else:
                gene, clade = "SINE", "SINE"

            # Apply remapping for display
            display_gene = remap.get(gene, gene)
            genes.append(display_gene)
            clades.append(clade)
            models.append(model)
            domain_strs.append(f"{display_gene}|{clade}")

        order, superfamily, max_clade, complete = classify_element(
            genes, clades, models, config)

        results.append({
            "id": base_seq,
            "order": order,
            "superfamily": superfamily,
            "clade": max_clade,
            "complete": complete,
            "strand": strand,
            "domains": " ".join(domain_strs),
        })

    log.info(f"  Classified: {len(results)} sequences")
    return results


def export_classification_tsv(results, out_path):
    """Export classification results as TSV (TEsorter cls.tsv format)."""
    columns = ["#TE", "Order", "Superfamily", "Clade", "Complete",
               "Strand", "Domains"]
    with open(out_path, "w") as f:
        f.write("\t".join(columns) + "\n")
        for r in results:
            line = [r["id"], r["order"], r["superfamily"], r["clade"],
                    r["complete"], r["strand"], r["domains"]]
            f.write("\t".join(line) + "\n")
