"""Cross-database clade/superfamily name harmonization.

REXdb and GyDB classify TEs with different nomenclatures (GyDB emits lowercase
clade tokens like ``athila`` and superfamily ``Pao``; REXdb emits ``Athila`` and
``Bel-Pao``; REXdb's plant pararetrovirus branch is order ``pararetrovirus``
where GyDB has superfamily ``Caulimoviridae``). When several databases classify
the same element, ``reconcile_classifications`` votes on these labels by string
identity -- so without harmonization the only level at which evidence can pool is
Order. This module maps each database's runtime ``(order, superfamily, clade)``
onto a unified taxonomy, so the cross-database vote pools at every level where
the two databases genuinely agree (e.g. GyDB ``athila`` and REXdb ``Athila``
both become ``Gypsy/Athila``).

The mapping table ``database/clade_harmonization.tsv`` contains ONLY the rows
that change something; any triple absent from the table is already canonical and
returned unchanged (identity fallback). The table is keyed on the runtime triple
that ``classify_sequences`` actually emits, not on raw HMM names.

Both tables are looked up in the resolved database directory first, then in the
packaged one, so a custom ``TESORTER2_DB`` collection that does not carry them
still harmonizes using the shipped tables.
"""
import logging
import os

from .paths import get_db_dir, packaged_db_dir

log = logging.getLogger(__name__)

_HARMONIZATION_TSV = "clade_harmonization.tsv"
_SCOPE_TSV = "clade_scope.tsv"

# Module-level cache: {path: {(db, order, sf, clade): (order_h, sf_h, clade_h)}}
_CACHE = {}
# Scope cache: {path: {db: {harmonized_superfamily: resolving_bool}}}
_SCOPE_CACHE = {}


def _resolve_table(filename, db_dir=None):
    """Locate a harmonization table: resolved database dir, else packaged one.

    Returns the packaged path even when it does not exist, so the caller emits
    one warning naming a stable location rather than a transient override.
    """
    candidate = os.path.join(get_db_dir(db_dir), filename)
    if os.path.exists(candidate):
        return candidate
    return os.path.join(packaged_db_dir(), filename)


def load_harmonization(path=None):
    """Load the harmonization table into a dict keyed by the runtime triple.

    Returns {(db, order_native, superfamily_native, clade_native):
             (order_h, superfamily_h, clade_h, lineage_h)}. Cached per path.
    Missing file yields an empty table (harmonization then becomes a pure
    identity no-op). ``lineage_h`` is "" for every row except the handful where
    one database resolves below the shared clade level (the REXdb Tat group);
    7-column rows are accepted and read as an empty lineage.
    """
    if path is None:
        path = _resolve_table(_HARMONIZATION_TSV)
    if path in _CACHE:
        return _CACHE[path]

    table = {}
    if not os.path.exists(path):
        log.warning(f"Clade harmonization table not found: {path} "
                    f"(combined classification will use native names)")
        _CACHE[path] = table
        return table

    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if parts[0] == "db":  # header
                continue
            if len(parts) == 7:
                parts = parts + [""]
            if len(parts) != 8:
                continue
            db, o_n, sf_n, c_n, o_h, sf_h, c_h, l_h = parts
            table[(db, o_n, sf_n, c_n)] = (o_h, sf_h, c_h, l_h)

    _CACHE[path] = table
    return table


def harmonize(db, order, superfamily, clade, table=None):
    """Map one database's native ``(order, superfamily, clade)`` to the unified
    taxonomy. Triples absent from the table are already canonical and returned
    unchanged. Returns the triple only -- use ``harmonize_lineage`` for the
    sub-clade level.
    """
    if table is None:
        table = load_harmonization()
    return table.get((db, order, superfamily, clade),
                     (order, superfamily, clade, ""))[:3]


def harmonize_lineage(db, order, superfamily, clade, table=None):
    """The sub-clade lineage for a native triple, or "" when the database does
    not resolve below the shared clade level.

    Only populated where one database splits a clade the other can only name as
    a whole (REXdb's Ogre/Retand/TatI-III under GyDB's single ``tat``). Keeping
    it out of ``clade`` is what stops the clade name from depending on which
    database happened to win the vote.
    """
    if table is None:
        table = load_harmonization()
    return table.get((db, order, superfamily, clade),
                     (order, superfamily, clade, ""))[3]


def load_scope(path=None):
    """Load the a-priori scope mask: {db: {harmonized_superfamily: resolving}}.

    A database "resolves" a superfamily when it models >= 2 distinct clades
    there; otherwise that branch is a catch-all and the database cannot
    adjudicate the superfamily at lineage level. Missing file yields an empty
    mask (scope-aware fusion then degrades to the plain weighted vote).
    """
    if path is None:
        path = _resolve_table(_SCOPE_TSV)
    if path in _SCOPE_CACHE:
        return _SCOPE_CACHE[path]

    scope = {}
    if not os.path.exists(path):
        log.warning(f"Clade scope mask not found: {path} "
                    f"(combined classification will not be scope-aware)")
        _SCOPE_CACHE[path] = scope
        return scope

    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if parts[0] == "db":  # header
                continue
            if len(parts) != 4:
                continue
            db, sf_h, _n, resolving = parts
            scope.setdefault(db, {})[sf_h] = (resolving == "yes")

    _SCOPE_CACHE[path] = scope
    return scope


def resolves(db, superfamily, scope=None):
    """True iff ``db`` can adjudicate ``superfamily`` at lineage level (models
    >= 2 clades there). Unknown (db, superfamily) pairs default to False, so a
    database with no scope entry is treated as non-resolving (catch-all)."""
    if scope is None:
        scope = load_scope()
    return scope.get(db, {}).get(superfamily, False)
