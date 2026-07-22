"""
so_map.py — Map TEsorter2 classifications onto Sequence Ontology terms.

Classifications are emitted as Order/Superfamily/Clade, which is TEsorter's
vocabulary rather than a standard one. This resolves those labels to Sequence
Ontology names and accessions so the GFF3 is ontology-compliant and the
classifications are comparable with other annotation tools.

The authority is EDTA's TE_Sequence_Ontology.txt (bundled in data/), whose
third column lists the aliases in use across TE databases. Those aliases
already include "Order/Superfamily" spellings such as LTR/Copia and TIR/hAT,
so most labels resolve directly against it.

Lookup is a ladder, most specific first:
    1. "Order/Superfamily"   e.g. LTR/Copia   -> Copia_LTR_retrotransposon
    2. "Superfamily"         e.g. Helitron    -> helitron
    3. "Order"               e.g. LINE        -> LINE_element
    4. repeat_region (SO:0000657), the generic fallback

Matching is case- and separator-insensitive, since databases disagree on both
(REXdb carries both Tc1_Mariner and Tc1_mariner).
"""

import os

_DATA = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")
SO_FILE = os.path.join(_DATA, "TE_Sequence_Ontology.txt")

# Generic fallback for labels with no established SO term.
UNKNOWN_SO = ("repeat_region", "SO:0000657")

# Labels TEsorter2 emits that EDTA's aliases do not cover. Each maps onto the
# closest established SO term; none of these have a term of their own.
_EXTRA_ALIASES = {
    # Order-level: no SO term for these orders, fall back to their class.
    "dirs":           ("LTR_retrotransposon", "SO:0000186"),
    "penelope":       ("non_LTR_retrotransposon", "SO:0000189"),
    "pararetrovirus": ("LTR_retrotransposon", "SO:0000186"),
    "maverick":       ("DNA_transposon", "SO:0000182"),
    # LTR/Retrovirus is a retroviral lineage, not a TE superfamily with a term.
    "ltr/retrovirus": ("LTR_retrotransposon", "SO:0000186"),
    # "TIR" is ambiguous in the EDTA file: it is an alias of both
    # terminal_inverted_repeat_element (SO:0000208, the element) and
    # terminal_inverted_repeat (SO:0000481, the repeat feature itself). As an
    # order label it means the element.
    "tir":            ("terminal_inverted_repeat_element", "SO:0000208"),
}

# Superfamily alias prefixes to try per order. Databases spell the same
# superfamily under different prefixes -- EDTA lists Zator as DNA/Zator
# (RepeatMasker convention), not TIR/Zator -- so try each order's known
# spellings before falling back to the order itself.
_ORDER_PREFIXES = {
    "tir":  ("TIR", "DNA", "MITE"),
    "line": ("LINE", "nonLTR"),
    "ltr":  ("LTR",),
    "sine": ("SINE",),
}

# For lineages with no SO term of their own, EDTA invents a descriptive name
# and files it under a generic accession (CR1_LINE_retrotransposon is not a
# real SO term; its SO:0000194 is LINE_element's). Those names must not be
# used as a GFF3 feature type, so map the generic accessions back to the real
# term. IDs absent here are ones whose EDTA name IS the canonical SO name.
_CANONICAL_BY_ID = {
    "SO:0000182": "DNA_transposon",
    "SO:0000186": "LTR_retrotransposon",
    "SO:0000189": "non_LTR_retrotransposon",
    "SO:0000194": "LINE_element",
    "SO:0000206": "SINE_element",
    "SO:0000208": "terminal_inverted_repeat_element",
    "SO:0000544": "helitron",
    "SO:0000657": "repeat_region",
}


def _norm(s):
    """Normalize an alias for matching: case- and separator-insensitive."""
    return s.strip().lower().replace("-", "_").replace(" ", "")


def load_so_index(path=None):
    """Parse the EDTA ontology file into {normalized_alias: (so_name, so_id)}.

    Column 1 is the exact SO name, column 2 the SO accession, column 3 a
    comma-separated alias list. Commented and blank lines are ignored.
    """
    path = path or SO_FILE
    index = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            so_name, so_id = cols[0].strip(), cols[1].strip()
            aliases = [so_name]
            if len(cols) > 2:
                aliases += [a for a in cols[2].split(",")]
            for a in aliases:
                a = _norm(a)
                if a:
                    index.setdefault(a, (so_name, so_id))
    index.update(_EXTRA_ALIASES)
    return index


_INDEX = None


def _index():
    global _INDEX
    if _INDEX is None:
        _INDEX = load_so_index()
    return _INDEX


def _is_set(v):
    return bool(v) and v.strip().lower() not in ("unknown", "na", "")


def lookup(order, superfamily=None):
    """Resolve an Order/Superfamily label to (so_name, so_id).

    The name is EDTA's, which for lineages without their own SO term is a
    descriptive name filed under a generic accession. Use so_gff_type() for a
    name that is guaranteed to be a real SO term.

    Returns the generic repeat_region term when nothing matches.
    """
    idx = _index()
    order = (order or "").strip()
    superfamily = (superfamily or "").strip()

    candidates = []
    if _is_set(order) and _is_set(superfamily):
        candidates.append(f"{order}/{superfamily}")
        for prefix in _ORDER_PREFIXES.get(_norm(order), ()):
            candidates.append(f"{prefix}/{superfamily}")
        candidates.append(superfamily)
    if _is_set(order):
        candidates.append(order)

    for c in candidates:
        hit = idx.get(_norm(c))
        if hit:
            return hit
    return UNKNOWN_SO


def so_gff_type(order, superfamily=None):
    """Resolve to (so_term, so_id) where so_term is a real SO name.

    Suitable for the GFF3 type column, which must hold a valid ontology term.
    """
    name, so_id = lookup(order, superfamily)
    return _CANONICAL_BY_ID.get(so_id, name), so_id
