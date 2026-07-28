#!/usr/bin/env python3
"""Compare a minimap2-derived classified TSV to a WFA goldstandard TSV.

Both files share the schema:
    qname  pass  pid  qcov  tcov  best_tname

Universe = WFA qnames. Queries missing from the prediction file are treated as
predicted_fail (matches downstream behavior of classify_ltr_paf_fast.py, which
drops queries with no PAF alignment carrying dv:f / de:f).

Pair agreement (TP only): best_tname stems compared after stripping any
"#..." suffix. Records counts of agree / disagree.

Output: JSON to --out (default stdout), single-line summary to stderr.
"""

import argparse
import json
import math
import sys


def strip_class(tname):
    if tname is None:
        return None
    i = tname.find("#")
    return tname if i < 0 else tname[:i]


def load_tsv(path, has_header):
    """Return {qname: (pass_bool, best_tname)}."""
    out = {}
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            if has_header and lineno == 1:
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            if len(f) < 6:
                continue
            qname = f[0]
            passing = (f[1].lower() == "pass")
            tname = f[5]
            out[qname] = (passing, tname)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gold", required=True, help="WFA goldstandard TSV (with header)")
    ap.add_argument("--pred", required=True, help="prediction TSV (no header)")
    ap.add_argument("--gold-has-header", action="store_true", default=True)
    ap.add_argument("--pred-has-header", action="store_true", default=False)
    ap.add_argument("--out", default="-", help="output JSON path (default stdout)")
    ap.add_argument("--label", default="", help="optional label echoed into JSON")
    args = ap.parse_args()

    gold = load_tsv(args.gold, has_header=args.gold_has_header)
    pred = load_tsv(args.pred, has_header=args.pred_has_header)

    universe = set(gold.keys())
    extras = set(pred.keys()) - universe

    TP = TN = FP = FN = 0
    tp_agree = tp_disagree = 0
    missing_from_pred = 0

    for q in universe:
        g_pass, g_t = gold[q]
        if q in pred:
            p_pass, p_t = pred[q]
        else:
            p_pass, p_t = False, None
            missing_from_pred += 1

        if g_pass and p_pass:
            TP += 1
            if strip_class(g_t) == strip_class(p_t):
                tp_agree += 1
            else:
                tp_disagree += 1
        elif (not g_pass) and (not p_pass):
            TN += 1
        elif (not g_pass) and p_pass:
            FP += 1
        else:
            FN += 1

    total = TP + TN + FP + FN
    accuracy = (TP + TN) / total if total else 0.0
    precision = TP / (TP + FP) if (TP + FP) else 0.0
    recall = TP / (TP + FN) if (TP + FN) else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0

    denom = math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc = ((TP * TN) - (FP * FN)) / denom if denom else 0.0

    out = {
        "label": args.label,
        "total": total,
        "TP": TP, "TN": TN, "FP": FP, "FN": FN,
        "accuracy": round(accuracy, 6),
        "precision": round(precision, 6),
        "recall": round(recall, 6),
        "f1": round(f1, 6),
        "mcc": round(mcc, 6),
        "tp_agree": tp_agree,
        "tp_disagree": tp_disagree,
        "tp_agree_frac": round(tp_agree / TP, 6) if TP else 0.0,
        "missing_from_pred": missing_from_pred,
        "extras": len(extras),
    }

    if args.out == "-":
        json.dump(out, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        with open(args.out, "w") as fh:
            json.dump(out, fh, indent=2)
            fh.write("\n")

    sys.stderr.write(
        "[eval] %s F1=%.4f acc=%.4f P=%.4f R=%.4f MCC=%.4f "
        "TP=%d TN=%d FP=%d FN=%d agree=%d/%d miss=%d extra=%d\n"
        % (args.label or args.pred, f1, accuracy, precision, recall, mcc,
           TP, TN, FP, FN, tp_agree, TP, missing_from_pred, len(extras))
    )


if __name__ == "__main__":
    main()
