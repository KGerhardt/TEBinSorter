"""
paths.py — Centralized resolution of the TEsorter2 database directory.

The HMM databases ship inside the package (tesorter2/database/), as they do in
TEsorter, so an install needs no download step. Resolution order:

    1. explicit override (e.g. a --db-dir CLI argument)
    2. the TESORTER2_DB environment variable
    3. the packaged tesorter2/database/

Only set 1 or 2 to point at a custom database collection.
"""

import os

_PKG_DIR = os.path.dirname(os.path.realpath(__file__))


def packaged_db_dir():
    """The database/ directory shipped inside the package."""
    return os.path.join(_PKG_DIR, "database")


def get_db_dir(override=None):
    """Resolve the database directory. See module docstring for precedence."""
    if override:
        return os.path.abspath(override)
    env = os.environ.get("TESORTER2_DB")
    if env:
        return os.path.abspath(env)
    return packaged_db_dir()
