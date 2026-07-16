"""
paths.py — Centralized resolution of the TESorter2 database directory.

The HMM databases are distributed separately from the code (they are hundreds
of MB) and fetched with `tesorter2-download-db`. This module resolves where
they live, with the following precedence:

    1. explicit override (e.g. a --db-dir CLI argument)
    2. the TESORTER2_DB environment variable
    3. a repo-bundled ``database/`` next to the installed package (dev checkout)
    4. the user data dir populated by download-db:
       ``$XDG_DATA_HOME/tesorter2/database`` (default ``~/.local/share/...``)
"""

import os

_PKG_DIR = os.path.dirname(os.path.realpath(__file__))


def user_data_dir():
    """Base data directory where download-db installs the databases."""
    xdg = os.environ.get("XDG_DATA_HOME") or os.path.join(
        os.path.expanduser("~"), ".local", "share")
    return os.path.join(xdg, "tesorter2")


def _bundled_db_dir():
    """database/ shipped alongside the source tree (development checkout)."""
    return os.path.normpath(os.path.join(_PKG_DIR, "..", "database"))


def _has_hmms(d):
    try:
        return any(f.endswith(".hmm") for f in os.listdir(d))
    except OSError:
        return False


def get_db_dir(override=None):
    """Resolve the database directory. See module docstring for precedence.

    Never raises: returns the best candidate even if it does not yet exist, so
    callers can produce a helpful "run tesorter2-download-db" message.
    """
    if override:
        return os.path.abspath(override)
    env = os.environ.get("TESORTER2_DB")
    if env:
        return os.path.abspath(env)
    bundled = _bundled_db_dir()
    if _has_hmms(bundled):
        return bundled
    return os.path.join(user_data_dir(), "database")
