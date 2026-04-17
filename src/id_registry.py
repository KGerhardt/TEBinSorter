"""
id_registry.py — Deterministic integer ID registry for string identifiers.

Maps string names (model names, sequence names) to stable 64-bit integer
IDs via hashing. The mapping is deterministic: same string always produces
the same ID regardless of run order or timing.

Lookup tables are persisted in SQLite for human-readable recovery.
The hash function is the authority; the table is just a cache.

Usage:
    reg = IDRegistry(conn)
    model_id = reg.model_id("Class_I/LTR/Ty1_copia/SIRE:Ty1-RT")
    seq_id = reg.seq_id("Os0008_INT#LTR/Copia")
    frame_id = reg.frame_id("Os0008_INT#LTR/Copia|fwd2")
"""

import hashlib
import sqlite3
from functools import lru_cache


REGISTRY_SCHEMA = """
CREATE TABLE IF NOT EXISTS id_models (
    id      INTEGER PRIMARY KEY,
    name    TEXT NOT NULL UNIQUE,
    family  TEXT,
    M       INTEGER
);

CREATE TABLE IF NOT EXISTS id_sequences (
    id      INTEGER PRIMARY KEY,
    name    TEXT NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS id_frames (
    id          INTEGER PRIMARY KEY,
    name        TEXT NOT NULL UNIQUE,
    seq_id      INTEGER,
    direction   INTEGER,
    frame       INTEGER
);

CREATE TABLE IF NOT EXISTS id_families (
    id      INTEGER PRIMARY KEY,
    name    TEXT NOT NULL UNIQUE
);
"""


def stable_hash(name):
    """Deterministic 63-bit positive integer hash of a string.

    Uses SHA-256 truncated to 63 bits (positive signed int64).
    Collision probability is negligible for any practical database size.
    """
    h = hashlib.sha256(name.encode("utf-8")).digest()
    return int.from_bytes(h[:8], "little") & 0x7FFFFFFFFFFFFFFF


def _parse_family(model_name):
    """Extract domain family from model name."""
    if ":" in model_name:
        return model_name.split(":")[1].split("-")[-1]
    return model_name.split("_")[0]


def _parse_frame(frame_name):
    """Parse frame name into (seq_name, direction, frame_num).

    Direction: 1 = forward, -1 = reverse, 0 = unknown.
    Frame: 0-indexed (0, 1, 2).
    """
    parts = frame_name.rsplit("|", 1)
    if len(parts) != 2:
        return frame_name, 0, 0

    seq_name, suffix = parts
    if suffix.startswith("fwd"):
        direction = 1
        frame_num = int(suffix[-1]) - 1
    elif suffix.startswith("rev"):
        direction = -1
        frame_num = int(suffix[-1]) - 1
    else:
        direction = 0
        frame_num = 0

    return seq_name, direction, frame_num


class IDRegistry:
    """Manages string-to-integer ID mappings with SQLite persistence."""

    def __init__(self, conn):
        """Initialize with a SQLite connection. Creates tables if needed."""
        self.conn = conn
        conn.executescript(REGISTRY_SCHEMA)
        conn.commit()

    @lru_cache(maxsize=4096)
    def model_id(self, name, M=None):
        """Get or register a model ID."""
        mid = stable_hash(name)
        family = _parse_family(name)
        fid = self.family_id(family)
        self.conn.execute(
            "INSERT OR IGNORE INTO id_models (id, name, family, M) "
            "VALUES (?, ?, ?, ?)",
            (mid, name, family, M),
        )
        return mid

    @lru_cache(maxsize=4096)
    def family_id(self, name):
        """Get or register a family ID."""
        fid = stable_hash(f"family:{name}")
        self.conn.execute(
            "INSERT OR IGNORE INTO id_families (id, name) VALUES (?, ?)",
            (fid, name),
        )
        return fid

    @lru_cache(maxsize=65536)
    def seq_id(self, name):
        """Get or register a sequence ID."""
        sid = stable_hash(name)
        self.conn.execute(
            "INSERT OR IGNORE INTO id_sequences (id, name) VALUES (?, ?)",
            (sid, name),
        )
        return sid

    @lru_cache(maxsize=65536)
    def frame_id(self, name):
        """Get or register a frame ID."""
        fid = stable_hash(name)
        seq_name, direction, frame_num = _parse_frame(name)
        sid = self.seq_id(seq_name)
        self.conn.execute(
            "INSERT OR IGNORE INTO id_frames "
            "(id, name, seq_id, direction, frame) VALUES (?, ?, ?, ?, ?)",
            (fid, name, sid, direction, frame_num),
        )
        return fid

    def commit(self):
        """Flush pending registrations to disk."""
        self.conn.commit()

    # --- Reverse lookups ---

    def model_name(self, mid):
        """Look up model name from ID."""
        row = self.conn.execute(
            "SELECT name FROM id_models WHERE id = ?", (mid,)
        ).fetchone()
        return row[0] if row else None

    def seq_name(self, sid):
        """Look up sequence name from ID."""
        row = self.conn.execute(
            "SELECT name FROM id_sequences WHERE id = ?", (sid,)
        ).fetchone()
        return row[0] if row else None

    def frame_name(self, fid):
        """Look up frame name from ID."""
        row = self.conn.execute(
            "SELECT name FROM id_frames WHERE id = ?", (fid,)
        ).fetchone()
        return row[0] if row else None

    def family_name(self, fid):
        """Look up family name from ID."""
        row = self.conn.execute(
            "SELECT name FROM id_families WHERE id = ?", (fid,)
        ).fetchone()
        return row[0] if row else None

    # --- Batch operations ---

    def register_models(self, hmms):
        """Register all models from a list of pyhmmer HMM objects."""
        for h in hmms:
            name = h.name if isinstance(h.name, str) else h.name.decode()
            self.model_id(name, M=h.M)
        self.commit()

    def register_frames(self, frame_names):
        """Register a list of frame names."""
        for name in frame_names:
            self.frame_id(name)
        self.commit()
