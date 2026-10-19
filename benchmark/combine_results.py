#!/usr/bin/env python3
"""Collect all evaluate.py JSON outputs under a directory into one summary table."""
import json, glob, os, sys

def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "results"
    files = [f for f in sorted(glob.glob(os.path.join(root, "**", "*.json"), recursive=True))
             if not any(part.startswith("_") for part in f.split(os.sep))]
    if not files:
        print(f"No *.json under {root}/ (run the pipelines first)"); return
    cols = ["tool", "library", "strandedness", "calls", "TP", "FN", "FP", "precision", "recall",
            "f1", "specificity", "accuracy", "quant_spearman", "strand_accuracy", "recall_strand"]
    rows = []
    for f in files:
        d = json.load(open(f))
        variant = os.path.basename(f).replace("eval_", "").replace(".json", "")
        d["tool"] = f"{d['tool']}"
        d["library"] = variant
        rows.append(d)
    w = {c: max(len(c), *(len(f"{r.get(c,'')}" if not isinstance(r.get(c), float)
                          else f"{r.get(c):.3f}") for r in rows)) for c in cols}
    line = "  ".join(c.ljust(w[c]) for c in cols)
    print(line); print("-" * len(line))
    for r in sorted(rows, key=lambda x: (x["tool"], x["library"])):
        cells = []
        for c in cols:
            v = r.get(c, "")
            v = f"{v:.3f}" if isinstance(v, float) else f"{v}"
            cells.append(v.ljust(w[c]))
        print("  ".join(cells))
    out = os.path.join(root, "summary.tsv")
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda x: (x["tool"], x["library"])):
            fh.write("\t".join(f"{r.get(c,'')}" for c in cols) + "\n")
    print(f"\nwritten: {out}")

if __name__ == "__main__":
    main()
