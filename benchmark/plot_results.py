#!/usr/bin/env python3
"""
Plot benchmark results produced by evaluate.py (one JSON per tool x library x variant).

Outputs (PNG, into --out):
  recall_vs_coverage.png        recall as a function of junction-read coverage bin, one line per run
  per_category_recall.png       heatmap of recall per test-case category (A-M) x run
  precision_recall_f1.png       grouped bars of precision / recall / F1 / specificity per run
  quantification.png            Spearman correlation of detected vs expected junction reads per run
  strand_accuracy.png           strand of matched calls: correct / wrong / none, per run

Usage:
  python3 plot_results.py --results results --out results/plots
"""
import argparse, glob, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COV_BINS = ["1-2", "3-9", "10-49", "50+"]

def frac(s):
    a, b = s.split("/"); a, b = int(a), int(b)
    return (a / b if b else float("nan")), a, b

def load(results):
    runs = []
    for f in sorted(glob.glob(os.path.join(results, "**", "*.json"), recursive=True)):
        rel = os.path.relpath(f, results)
        if any(p.startswith("_") for p in rel.split(os.sep)):   # skip _selftest/_demo unless --results points inside it
            continue
        d = json.load(open(f))
        label = os.path.basename(f)[:-5].replace("eval_", "")
        d["label"] = f"{d.get('tool','?')}:{label}"
        runs.append(d)
    return runs

def plot_recall_vs_coverage(runs, out):
    plt.figure(figsize=(7, 4.5))
    x = np.arange(len(COV_BINS))
    for r in runs:
        cr = r.get("coverage_recall", {})
        y = [frac(cr[b])[0] if b in cr else np.nan for b in COV_BINS]
        plt.plot(x, y, marker="o", label=r["label"])
    plt.xticks(x, [f"{b}\nreads" for b in COV_BINS])
    plt.ylim(-0.03, 1.03); plt.ylabel("recall"); plt.xlabel("junction-read coverage")
    plt.title("Sensitivity vs coverage")
    plt.grid(alpha=0.3); plt.legend(fontsize=7, ncol=2)
    plt.tight_layout(); plt.savefig(os.path.join(out, "recall_vs_coverage.png"), dpi=140); plt.close()

def plot_category_heatmap(runs, out):
    cats = set()
    for r in runs: cats |= set(r.get("per_category_recall", {}))
    cats = sorted(cats)
    M = np.full((len(cats), len(runs)), np.nan)
    for j, r in enumerate(runs):
        pc = r.get("per_category_recall", {})
        for i, c in enumerate(cats):
            if c in pc: M[i, j] = frac(pc[c])[0]
    fig, ax = plt.subplots(figsize=(max(6, 1.1*len(runs)+3), max(4, 0.42*len(cats)+1)))
    im = ax.imshow(M, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(runs))); ax.set_xticklabels([r["label"] for r in runs], rotation=40, ha="right", fontsize=7)
    ax.set_yticks(range(len(cats))); ax.set_yticklabels(cats, fontsize=8)
    for i in range(len(cats)):
        for j in range(len(runs)):
            if not np.isnan(M[i, j]):
                ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=6.5,
                        color="black")
    fig.colorbar(im, ax=ax, label="recall", shrink=0.8)
    ax.set_title("Per-category recall")
    plt.tight_layout(); plt.savefig(os.path.join(out, "per_category_recall.png"), dpi=140); plt.close()

def plot_prf(runs, out):
    metrics = ["precision", "recall", "f1", "specificity"]
    x = np.arange(len(runs)); w = 0.2
    plt.figure(figsize=(max(6, 1.2*len(runs)+2), 4.5))
    for k, m in enumerate(metrics):
        vals = [r.get(m, np.nan) for r in runs]
        plt.bar(x + (k-1.5)*w, vals, w, label=m)
    plt.xticks(x, [r["label"] for r in runs], rotation=40, ha="right", fontsize=7)
    plt.ylim(0, 1.05); plt.ylabel("score"); plt.title("Precision / Recall / F1 / Specificity")
    plt.legend(fontsize=8); plt.grid(axis="y", alpha=0.3)
    plt.tight_layout(); plt.savefig(os.path.join(out, "precision_recall_f1.png"), dpi=140); plt.close()

def plot_quant(runs, out):
    x = np.arange(len(runs))
    vals = [r.get("quant_spearman", np.nan) for r in runs]
    plt.figure(figsize=(max(6, 1.2*len(runs)+2), 4))
    plt.bar(x, vals, color="#4C72B0")
    plt.xticks(x, [r["label"] for r in runs], rotation=40, ha="right", fontsize=7)
    plt.ylim(0, 1.05); plt.ylabel("Spearman r"); plt.title("Quantification: detected vs expected junction reads")
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout(); plt.savefig(os.path.join(out, "quantification.png"), dpi=140); plt.close()

def plot_strand(runs, out):
    x = np.arange(len(runs))
    ok = np.array([r.get("strand_correct", 0) for r in runs], float)
    bad = np.array([r.get("strand_wrong", 0) for r in runs], float)
    unk = np.array([r.get("strand_unknown", 0) for r in runs], float)
    tp = np.where(ok + bad + unk > 0, ok + bad + unk, np.nan)
    plt.figure(figsize=(max(6, 1.2*len(runs)+2), 4.5))
    plt.bar(x, ok/tp, color="#55A868", label="correct strand")
    plt.bar(x, bad/tp, bottom=ok/tp, color="#C44E52", label="wrong strand")
    plt.bar(x, unk/tp, bottom=(ok+bad)/tp, color="#BBBBBB", label="no strand")
    plt.xticks(x, [f"{r['label']} [{r.get('strandedness', 'NA')}]" for r in runs],
               rotation=40, ha="right", fontsize=7)
    plt.ylim(0, 1.05); plt.ylabel("fraction of matched BSJs")
    plt.title("Strand of detected BSJs vs truth")
    plt.legend(fontsize=8); plt.grid(axis="y", alpha=0.3)
    plt.tight_layout(); plt.savefig(os.path.join(out, "strand_accuracy.png"), dpi=140); plt.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = a.out or os.path.join(a.results, "plots")
    os.makedirs(out, exist_ok=True)
    runs = load(a.results)
    if not runs:
        print(f"No JSON results under {a.results}/ — run the pipelines first."); return
    plot_recall_vs_coverage(runs, out)
    plot_category_heatmap(runs, out)
    plot_prf(runs, out)
    plot_quant(runs, out)
    plot_strand(runs, out)
    print(f"Wrote 5 plots for {len(runs)} run(s) -> {out}/")

if __name__ == "__main__":
    main()
