#!/usr/bin/env python3
"""Summarize a self-alignment PAF (with cg:Z: CIGAR).

Designed for `WFA_TEsorter -q X.fa -t X.fa --cigar > X.paf` where the same
FASTA is both query and target (so every sequence appears with a self-hit
at qname == tname). Self-hits are excluded from all three outputs.

Produces three files at <prefix>.*:

  <prefix>.best_per_query.tsv
      One row per query: the single best non-self target.
      Same columns as WFA_best_per_query.py:
          qname  pass  pid  qcov  tcov  best_tname

  <prefix>.rbh.tsv
      Reciprocal best-hit pairs: A's best non-self target is B AND B's
      best non-self target is A. Pairs are emitted once with seqA<seqB
      lexically. Metrics are reported from both directions:
          seqA  seqB  pid_AB  qcov_AB  tcov_AB
                      pid_BA  qcov_BA  tcov_BA  both_pass

  <prefix>.clusters.tsv
      Single-linkage clusters from non-self alignments that pass all
      three thresholds (--min-pid AND --min-qcov AND --min-tcov).
      Edges are added with union-find; isolates become singleton clusters.
      Representative = longest member (ties broken by name asc). Members
      are listed by length desc then name asc, so the representative is
      the first row of each cluster. Clusters ranked by size (desc), then
      representative name. One row per member:
          cluster_id  representative  member  cluster_size

Per-row metrics (from cg:Z: CIGAR, NM:i: when present):
  pid  = matches / (matches + mismatches + gap_opens)   # gap-compressed
  qcov = M_bases / qlen
  tcov = M_bases / tlen
  pass = pid >= --min-pid AND qcov >= --min-qcov AND tcov >= --min-tcov

For "best" target per query the sort key is (passed, pid*qcov*tcov), so a
passing alignment always beats a failing one and ties break by the raw
product. RBH reports the reciprocity regardless of pass/fail and adds a
both_pass column so downstream filtering is one awk away.

Streams the PAF; only one record per unique query is held in RAM.
Requires `cg:Z:` on every row (run WFA_TEsorter / minimap2 with --cigar).
"""

import argparse
import re
import sys


CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def parse_cigar(cigar):
    """Return (m_bases, i_bases, d_bases, gap_opens). M counts =/X too."""
    m = i = d = gap_opens = 0
    for n_str, op in CIGAR_RE.findall(cigar):
        n = int(n_str)
        if op in 'M=X':
            m += n
        elif op == 'I':
            i += n
            gap_opens += 1
        elif op == 'D':
            d += n
            gap_opens += 1
    return m, i, d, gap_opens


def row_metrics(fields):
    """Compute metrics for one PAF row, or None to skip.

    Returns (qname, qlen, tname, tlen, pid, qcov, tcov).
    """
    qname = fields[0]
    qlen = int(fields[1])
    tname = fields[5]
    tlen = int(fields[6])
    paf_matches = int(fields[9])

    cigar = None
    nm = None
    for tag in fields[12:]:
        if tag.startswith('cg:Z:'):
            cigar = tag[5:]
        elif tag.startswith('NM:i:'):
            nm = int(tag[5:])

    if cigar is None or qlen <= 0 or tlen <= 0:
        return None

    m_bases, i_bases, d_bases, gap_opens = parse_cigar(cigar)

    if nm is not None:
        # NM = mismatches + indel_bases
        mismatches = nm - i_bases - d_bases
        matches = m_bases - mismatches
    else:
        # fall back to PAF column 10 (residue matches)
        matches = paf_matches
        mismatches = m_bases - matches

    if matches < 0 or mismatches < 0:
        return None

    denom = matches + mismatches + gap_opens
    pid = matches / denom if denom > 0 else 0.0
    qcov = m_bases / qlen
    tcov = m_bases / tlen

    return qname, qlen, tname, tlen, pid, qcov, tcov


class UnionFind:
    """Union-find with path compression and union-by-size."""

    def __init__(self):
        self.parent = {}
        self.size = {}

    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.size[x] = 1

    def find(self, x):
        self.add(x)
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        # iterative path compression
        while self.parent[x] != root:
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('paf', nargs='?', default='-',
                    help='Input PAF (default: stdin). Use "-" for stdin.')
    ap.add_argument('-o', '--output', required=True,
                    help='Output prefix; writes <prefix>.best_per_query.tsv, '
                         '<prefix>.rbh.tsv, <prefix>.clusters.tsv')
    ap.add_argument('--min-pid', type=float, default=0.70,
                    help='Min gap-compressed identity (default: 0.70).')
    ap.add_argument('--min-qcov', type=float, default=0.70,
                    help='Min query coverage (default: 0.70).')
    ap.add_argument('--min-tcov', type=float, default=0.70,
                    help='Min target coverage (default: 0.70).')
    ap.add_argument('-v', '--verbose', action='store_true',
                    help='Report progress and summary stats to stderr.')
    args = ap.parse_args()

    fin = sys.stdin if args.paf == '-' else open(args.paf)

    # qname -> (sort_key, passed, pid, qcov, tcov, tname)  (non-self only)
    best = {}
    seq_len = {}
    uf = UnionFind()
    n_rows = n_skipped = n_self = n_edges_passing = 0

    for line in fin:
        if not line or line[0] == '#':
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 12:
            n_skipped += 1
            continue
        n_rows += 1
        out = row_metrics(fields)
        if out is None:
            n_skipped += 1
            continue
        qname, qlen, tname, tlen, pid, qcov, tcov = out

        # register every sequence we see (incl. via self-hit) so isolates
        # still appear as singleton clusters in the output
        if qname not in seq_len:
            seq_len[qname] = qlen
            uf.add(qname)
        if tname not in seq_len:
            seq_len[tname] = tlen
            uf.add(tname)

        if qname == tname:
            n_self += 1
            continue

        passed = (pid >= args.min_pid and
                  qcov >= args.min_qcov and
                  tcov >= args.min_tcov)
        key = (1 if passed else 0, pid * qcov * tcov)
        prev = best.get(qname)
        if prev is None or key > prev[0]:
            best[qname] = (key, passed, pid, qcov, tcov, tname)

        if passed:
            uf.union(qname, tname)
            n_edges_passing += 1

        if args.verbose and n_rows % 1_000_000 == 0:
            print(f'  ... {n_rows:,} rows  {len(seq_len):,} seqs  '
                  f'{n_edges_passing:,} passing pairs',
                  file=sys.stderr)

    if fin is not sys.stdin:
        fin.close()

    # ---------- single-best per query (non-self) ----------
    best_path = f'{args.output}.best_per_query.tsv'
    with open(best_path, 'w') as fb:
        fb.write('qname\tpass\tpid\tqcov\ttcov\tbest_tname\n')
        for qname, (_, passed, pid, qcov, tcov, tname) in best.items():
            fb.write(f'{qname}\t{"pass" if passed else "fail"}\t'
                     f'{pid:.4f}\t{qcov:.4f}\t{tcov:.4f}\t{tname}\n')

    # ---------- reciprocal best hits ----------
    rbh = []
    seen_pairs = set()
    for q, qrec in best.items():
        t = qrec[5]
        trec = best.get(t)
        if trec is None or trec[5] != q:
            continue
        a, b = (q, t) if q < t else (t, q)
        if (a, b) in seen_pairs:
            continue
        seen_pairs.add((a, b))
        ab = best[a]   # a -> b: (key, passed, pid, qcov, tcov, _)
        ba = best[b]   # b -> a
        both_pass = ab[1] and ba[1]
        rbh.append((a, b,
                    ab[2], ab[3], ab[4],
                    ba[2], ba[3], ba[4],
                    both_pass))
    rbh.sort()

    rbh_path = f'{args.output}.rbh.tsv'
    with open(rbh_path, 'w') as fr:
        fr.write('seqA\tseqB\tpid_AB\tqcov_AB\ttcov_AB\t'
                 'pid_BA\tqcov_BA\ttcov_BA\tboth_pass\n')
        for a, b, p_ab, q_ab, t_ab, p_ba, q_ba, t_ba, bp in rbh:
            fr.write(f'{a}\t{b}\t'
                     f'{p_ab:.4f}\t{q_ab:.4f}\t{t_ab:.4f}\t'
                     f'{p_ba:.4f}\t{q_ba:.4f}\t{t_ba:.4f}\t'
                     f'{"yes" if bp else "no"}\n')

    # ---------- clusters (single-linkage on passing non-self edges) ----------
    members = {}
    for s in seq_len:
        members.setdefault(uf.find(s), []).append(s)

    clusters = []  # list of (size, representative, sorted_members)
    for r, mems in members.items():
        # sort by length desc, then name asc; representative = first
        mems_sorted = sorted(mems, key=lambda s: (-seq_len.get(s, 0), s))
        clusters.append((len(mems), mems_sorted[0], mems_sorted))
    # rank: largest first, then by representative name
    clusters.sort(key=lambda c: (-c[0], c[1]))

    cl_path = f'{args.output}.clusters.tsv'
    with open(cl_path, 'w') as fc:
        fc.write('cluster_id\trepresentative\tmember\tcluster_size\n')
        for idx, (size, rep, mems) in enumerate(clusters, start=1):
            cid = f'C{idx:06d}'
            for m in mems:
                fc.write(f'{cid}\t{rep}\t{m}\t{size}\n')

    if args.verbose:
        n_pass_q = sum(1 for v in best.values() if v[1])
        n_singletons = sum(1 for c in clusters if c[0] == 1)
        n_nontriv = len(clusters) - n_singletons
        n_rbh_pass = sum(1 for r in rbh if r[8])
        print(f'rows: {n_rows:,}  skipped: {n_skipped:,}  '
              f'self-hits: {n_self:,}', file=sys.stderr)
        print(f'sequences seen: {len(seq_len):,}', file=sys.stderr)
        print(f'queries with non-self hit: {len(best):,}  '
              f'(pass: {n_pass_q:,}, fail: {len(best) - n_pass_q:,})',
              file=sys.stderr)
        print(f'rbh pairs: {len(rbh):,}  (both-pass: {n_rbh_pass:,})',
              file=sys.stderr)
        print(f'clusters: {len(clusters):,}  '
              f'(non-singleton: {n_nontriv:,}, singletons: {n_singletons:,})',
              file=sys.stderr)
        print(f'wrote {best_path}', file=sys.stderr)
        print(f'wrote {rbh_path}', file=sys.stderr)
        print(f'wrote {cl_path}', file=sys.stderr)


if __name__ == '__main__':
    main()
