#!/usr/bin/env python3
"""Pick the single best target per query from a PAF (with cg:Z: CIGAR).

For each PAF row we compute, from the CIGAR:
  pid  = matches / (matches + mismatches + gap_opens)   # gap-compressed
  qcov = M_bases / qlen   (M_bases = bases not in I or D ops)
  tcov = M_bases / tlen
  pass = pid >= --min-pid AND qcov >= --min-qcov AND tcov >= --min-tcov

Per query, we keep the alignment that maximizes pass*pid*qcov*tcov;
ties (e.g. all failing -> product 0) are broken by raw pid*qcov*tcov.

Output TSV:  qname  pass  pid  qcov  tcov  best_tname

Notes:
- Streams the PAF line by line; only one row per query is held in memory.
- Requires `cg:Z:` CIGAR on every row (run minimap2/WFA_TEsorter with --cigar).
- mismatches are derived as NM - I_bases - D_bases (NM = mismatches + indels).
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
    """Compute (qname, tname, pid, qcov, tcov) for one PAF row, or None to skip."""
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

    return qname, tname, pid, qcov, tcov


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('paf', nargs='?', default='-',
                    help='Input PAF (default: stdin). Use "-" for stdin.')
    ap.add_argument('-o', '--output', default='-',
                    help='Output TSV (default: stdout).')
    ap.add_argument('--min-pid', type=float, default=0.70,
                    help='Min gap-compressed identity to pass (default: 0.70).')
    ap.add_argument('--min-qcov', type=float, default=0.70,
                    help='Min query coverage to pass (default: 0.70).')
    ap.add_argument('--min-tcov', type=float, default=0.70,
                    help='Min target coverage to pass (default: 0.70).')
    ap.add_argument('-v', '--verbose', action='store_true',
                    help='Report progress to stderr.')
    args = ap.parse_args()

    fin = sys.stdin if args.paf == '-' else open(args.paf)
    fout = sys.stdout if args.output == '-' else open(args.output, 'w')

    # qname -> (sort_key, passed, pid, qcov, tcov, tname)
    best = {}
    n_rows = n_skipped = 0

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
        qname, tname, pid, qcov, tcov = out
        passed = (pid >= args.min_pid and
                  qcov >= args.min_qcov and
                  tcov >= args.min_tcov)
        # primary: pass beats fail; secondary: raw product
        key = (1 if passed else 0, pid * qcov * tcov)
        prev = best.get(qname)
        if prev is None or key > prev[0]:
            best[qname] = (key, passed, pid, qcov, tcov, tname)
        if args.verbose and n_rows % 1_000_000 == 0:
            print(f'  ... {n_rows:,} rows, {len(best):,} unique queries',
                  file=sys.stderr)

    fout.write('qname\tpass\tpid\tqcov\ttcov\tbest_tname\n')
    for qname, (_, passed, pid, qcov, tcov, tname) in best.items():
        fout.write(f'{qname}\t{"pass" if passed else "fail"}\t'
                   f'{pid:.4f}\t{qcov:.4f}\t{tcov:.4f}\t{tname}\n')

    if args.verbose:
        n_pass = sum(1 for v in best.values() if v[1])
        print(f'rows: {n_rows:,}  skipped: {n_skipped:,}  '
              f'queries: {len(best):,}  pass: {n_pass:,}  '
              f'fail: {len(best) - n_pass:,}', file=sys.stderr)

    if fin is not sys.stdin:
        fin.close()
    if fout is not sys.stdout:
        fout.close()


if __name__ == '__main__':
    main()
