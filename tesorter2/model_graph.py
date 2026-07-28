"""
model_graph.py — Cross-model similarity graph for HMM databases.

Searches every model's consensus sequence against every model to build
a directed similarity graph. The graph encodes which models recognize
similar sequences, enabling efficient traversal during classification.

The graph is persisted as a SQLite database alongside the HMM file,
precomputed once and reused across runs.

Standalone module: depends only on pyhmmer, numpy, sqlite3.

API:
    build_graph(hmms)                  → ModelGraph from HMM list
    build_graph_from_file(hmm_path)    → ModelGraph from file
    ModelGraph.save(db_path)           → persist to SQLite
    ModelGraph.load(db_path)           → load from SQLite
    ModelGraph.neighbors(model, k)     → top-k similar models
    ModelGraph.traversal_order(entry)  → graph-search order from entry point
"""

import io
import os
import sqlite3
from collections import defaultdict

import numpy as np
import pyhmmer
import pyhmmer.easel as easel
import pyhmmer.plan7 as plan7


GRAPH_SCHEMA = """
CREATE TABLE IF NOT EXISTS models (
    name    TEXT PRIMARY KEY,
    M       INTEGER NOT NULL,
    family  TEXT
);

CREATE TABLE IF NOT EXISTS edges (
    query   TEXT NOT NULL,
    target  TEXT NOT NULL,
    score   REAL NOT NULL,
    PRIMARY KEY (query, target)
);

CREATE TABLE IF NOT EXISTS metadata (
    key     TEXT PRIMARY KEY,
    value   TEXT
);

CREATE INDEX IF NOT EXISTS idx_edges_query ON edges(query);
CREATE INDEX IF NOT EXISTS idx_edges_score ON edges(score);
"""


def _infer_family(model_name):
    """Infer domain family from model name.

    Handles both REXdb format (Class_I/LTR/Ty1_copia/SIRE:Ty1-RT → RT)
    and GyDB format (RT_copia → RT).
    """
    if ":" in model_name:
        # REXdb: domain is after the colon, family is after the last dash
        domain = model_name.split(":")[1]
        return domain.split("-")[-1]
    else:
        # GyDB: family is the prefix before the first underscore
        return model_name.split("_")[0]


class ModelGraph:
    """Directed similarity graph between HMM models."""

    def __init__(self):
        self.models = {}       # name → {M, family}
        self.edges = {}        # (query, target) → score
        self._neighbors = {}   # query → [(target, score), ...] sorted desc

    def add_model(self, name, M, family=None):
        if family is None:
            family = _infer_family(name)
        self.models[name] = {"M": M, "family": family}

    def add_edge(self, query, target, score):
        if query == target:
            return
        self.edges[(query, target)] = score

    def _build_neighbor_cache(self):
        """Build sorted neighbor lists from edges."""
        adj = defaultdict(list)
        for (q, t), score in self.edges.items():
            adj[q].append((t, score))
        self._neighbors = {
            q: sorted(pairs, key=lambda x: -x[1])
            for q, pairs in adj.items()
        }

    def neighbors(self, model, k=None):
        """Get top-k neighbors of a model, sorted by score descending.

        Args:
            model: model name
            k: max neighbors to return (None = all)

        Returns:
            list of (neighbor_name, score) tuples
        """
        if not self._neighbors:
            self._build_neighbor_cache()
        nbrs = self._neighbors.get(model, [])
        return nbrs[:k] if k else nbrs

    def family_of(self, model):
        """Get the domain family of a model."""
        if model in self.models:
            return self.models[model]["family"]
        return _infer_family(model)

    def traversal_order(self, entry_model, max_depth=10):
        """Generate graph traversal order from an entry point.

        Breadth-first by score: at each step, expand the best-scoring
        unvisited neighbor. Used for guided search after a facet hit.

        Args:
            entry_model: starting model name
            max_depth: max neighbors to expand per visited node

        Yields:
            model names in traversal order (entry_model first)
        """
        visited = set()
        # Priority queue: (model, score_that_led_here)
        queue = [(entry_model, float("inf"))]

        while queue:
            # Pop highest-scoring candidate
            queue.sort(key=lambda x: -x[1])
            model, _ = queue.pop(0)

            if model in visited:
                continue
            visited.add(model)
            yield model

            # Expand neighbors
            for neighbor, score in self.neighbors(model, k=max_depth):
                if neighbor not in visited:
                    queue.append((neighbor, score))

    def save(self, db_path):
        """Persist graph to SQLite."""
        conn = sqlite3.connect(db_path)
        conn.executescript(GRAPH_SCHEMA)

        conn.executemany(
            "INSERT OR REPLACE INTO models (name, M, family) VALUES (?, ?, ?)",
            [(name, info["M"], info["family"])
             for name, info in self.models.items()],
        )

        conn.executemany(
            "INSERT OR REPLACE INTO edges (query, target, score) VALUES (?, ?, ?)",
            [(q, t, s) for (q, t), s in self.edges.items()],
        )

        conn.executemany(
            "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
            [("n_models", str(len(self.models))),
             ("n_edges", str(len(self.edges)))],
        )

        conn.commit()
        conn.close()

    @classmethod
    def load(cls, db_path):
        """Load graph from SQLite."""
        g = cls()
        conn = sqlite3.connect(db_path)

        for name, M, family in conn.execute(
                "SELECT name, M, family FROM models"):
            g.models[name] = {"M": M, "family": family}

        for query, target, score in conn.execute(
                "SELECT query, target, score FROM edges"):
            g.edges[(query, target)] = score

        conn.close()
        g._build_neighbor_cache()
        return g

    def __repr__(self):
        return (f"ModelGraph({len(self.models)} models, "
                f"{len(self.edges)} edges)")


def build_graph(hmms):
    """Build a similarity graph from a list of pyhmmer HMM objects.

    Searches every model's consensus sequence against every model
    using hmmsearch. Edges represent cross-model recognition.

    Args:
        hmms: list of pyhmmer.plan7.HMM objects

    Returns:
        ModelGraph
    """
    g = ModelGraph()
    abc = hmms[0].alphabet

    # Register models
    for h in hmms:
        name = h.name if isinstance(h.name, str) else h.name.decode()
        g.add_model(name, h.M)

    # Build consensus sequence block
    cons_seqs = []
    for h in hmms:
        name = h.name if isinstance(h.name, str) else h.name.decode()
        ts = easel.TextSequence(
            name=name.encode("ascii"),
            sequence=h.consensus,
        )
        cons_seqs.append(ts.digitize(abc))
    cons_block = easel.DigitalSequenceBlock(abc, cons_seqs)

    # Search all models against all consensuses
    Z = len(hmms)
    for top_hits in pyhmmer.hmmsearch(
        hmms, cons_block,
        bias_filter=False,
        Z=Z, domZ=Z, E=1e10, domE=1e10,
    ):
        if len(top_hits) == 0:
            continue
        buf = io.BytesIO()
        top_hits.write(buf, format="targets", header=False)
        text = buf.getvalue().decode().strip()
        if not text:
            continue
        for line in text.split("\n"):
            cols = line.split()
            if len(cols) < 6:
                continue
            target, query, score = cols[0], cols[2], float(cols[5])
            g.add_edge(query, target, score)

    g._build_neighbor_cache()
    return g


def build_graph_from_file(hmm_path):
    """Build a similarity graph from an HMM database file.

    Args:
        hmm_path: path to HMM file (single or multi-model)

    Returns:
        ModelGraph
    """
    hmms = []
    with plan7.HMMFile(hmm_path) as hf:
        for hmm in hf:
            hmms.append(hmm)
    return build_graph(hmms)


def get_or_build_graph(hmm_path, graph_dir=None):
    """Load a cached graph or build and cache one.

    Graph is stored as {hmm_basename}.graph.db alongside the HMM file
    (or in graph_dir if specified).

    Args:
        hmm_path: path to HMM database file
        graph_dir: optional directory for graph files

    Returns:
        ModelGraph
    """
    if graph_dir is None:
        graph_dir = os.path.dirname(os.path.abspath(hmm_path))

    graph_path = f"{hmm_path}.similarity_graph.db"
    if graph_dir:
        basename = os.path.basename(hmm_path)
        graph_path = os.path.join(graph_dir,
                                   f"{basename}.similarity_graph.db")

    if os.path.exists(graph_path):
        return ModelGraph.load(graph_path)

    graph = build_graph_from_file(hmm_path)
    graph.save(graph_path)
    return graph


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    import argparse
    import time

    parser = argparse.ArgumentParser(
        description="Build or inspect cross-model similarity graph")
    parser.add_argument("hmm_file",
                        help="HMM database file")
    parser.add_argument("--rebuild", action="store_true",
                        help="Force rebuild even if cached graph exists")
    parser.add_argument("--show-neighbors", type=str, default=None,
                        help="Show top neighbors for a specific model")
    parser.add_argument("--top-k", type=int, default=10,
                        help="Number of neighbors to show [default: 10]")
    args = parser.parse_args()

    graph_path = f"{args.hmm_file}.similarity_graph.db"

    if args.rebuild or not os.path.exists(graph_path):
        print(f"Building graph for {args.hmm_file}...", flush=True)
        t0 = time.time()
        graph = build_graph_from_file(args.hmm_file)
        graph.save(graph_path)
        t1 = time.time()
        print(f"  {graph} in {t1 - t0:.1f}s -> {graph_path}")
    else:
        graph = ModelGraph.load(graph_path)
        print(f"Loaded {graph} from {graph_path}")

    # Summary
    families = defaultdict(int)
    for info in graph.models.values():
        families[info["family"]] += 1
    print(f"\nFamilies ({len(families)}):")
    for fam, count in sorted(families.items(), key=lambda x: -x[1]):
        print(f"  {fam}: {count} models")

    # Cross-family edges
    cross = 0
    within = 0
    for (q, t), score in graph.edges.items():
        if graph.family_of(q) == graph.family_of(t):
            within += 1
        else:
            cross += 1
    print(f"\nEdges: {within} within-family, {cross} cross-family")

    if args.show_neighbors:
        model = args.show_neighbors
        print(f"\nTop {args.top_k} neighbors of {model}:")
        for nbr, score in graph.neighbors(model, k=args.top_k):
            fam = graph.family_of(nbr)
            print(f"  {score:7.1f}  {fam:8s}  {nbr}")


if __name__ == "__main__":
    main()
