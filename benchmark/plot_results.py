#!/usr/bin/env python3
"""
Plot benchmark results produced by evaluate.py (one JSON per tool x library x variant).

Every figure is small multiples with a shared library axis, so runs stay readable as
tools and libraries are added. Strandedness variants (*_rev, *_unstr) sit next to their
base library. Tools keep a fixed color + marker in every figure.

Outputs (PNG, into --out):
  recall_vs_coverage.png        recall per junction-read coverage bin: library x bin heatmap, one panel per tool
  per_category_recall.png       recall per test-case category (A-M): category x library heatmap, one panel per tool
  precision_recall_f1.png       precision / recall / F1 / specificity: one panel per metric, one dot per tool
  quantification.png            Spearman of detected vs expected junction reads, one dot per tool
  strand_accuracy.png           every true BSJ: correct / wrong / no strand / not detected, one panel per tool
  strand_heatmap.png            the same outcome per individual true BSJ: BSJ x library, one panel per tool
                                (needs --dataset for the truth list)

Usage:
  python3 plot_results.py --results results --out results/plots --dataset ../dataset
"""
import argparse, glob, json, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Patch

COV_BINS = ["1-2", "3-9", "10-49", "50+"]

# Fixed tool identity: color (categorical slots, in order) + marker as secondary encoding.
TOOL_ORDER = ["bsj_finder", "ciri2", "circexplorer2_bed", "circexplorer2_known"]
TOOL_LABEL = {"bsj_finder": "bsj_finder", "ciri2": "CIRI2",
              "circexplorer2_bed": "CIRCexplorer2 (de novo)",
              "circexplorer2_known": "CIRCexplorer2 (known)"}
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
MARKERS = ["o", "s", "D", "^", "v", "P", "X", "h"]

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"
GOOD, CRITICAL, NONE = "#0ca30c", "#d03b3b", "#c3c2b7"
SEQ = LinearSegmentedColormap.from_list(
    "seq_blue", ["#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#104281"])

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8, "text.color": INK,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 9, "axes.facecolor": SURFACE, "figure.facecolor": SURFACE,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
})

LIB_RE = re.compile(r"^(pe|se)(\d+)_(lowerr|higherr)(?:_(rev|unstr))?$")

def lib_key(lib):
    """pe before se, long reads first, lowerr before higherr, then fwd/rev/unstr."""
    m = LIB_RE.match(lib)
    if not m:
        return (1, lib)
    layout, rlen, err, strand = m.groups()
    return (0, layout != "pe", -int(rlen), err != "lowerr", ["", "rev", "unstr"].index(strand or ""))

def frac(s):
    a, b = s.split("/"); a, b = int(a), int(b)
    return a / b if b else float("nan")

def load(results):
    runs = []
    for f in sorted(glob.glob(os.path.join(results, "**", "*.json"), recursive=True)):
        rel = os.path.relpath(f, results)
        if any(p.startswith("_") for p in rel.split(os.sep)):   # skip _selftest/_demo/_old unless --results points inside it
            continue
        d = json.load(open(f))
        d.setdefault("tool", "?")
        if d.get("library", "NA") == "NA":
            d["library"] = os.path.basename(f)[:-5].replace("eval_", "")
        runs.append(d)
    return runs

def axes_of(runs):
    tools = sorted({r["tool"] for r in runs},
                   key=lambda t: (TOOL_ORDER.index(t) if t in TOOL_ORDER else len(TOOL_ORDER), t))
    libs = sorted({r["library"] for r in runs}, key=lib_key)
    by = {(r["tool"], r["library"]): r for r in runs}
    style = {t: (SLOTS[i % len(SLOTS)], MARKERS[i % len(MARKERS)]) for i, t in enumerate(tools)}
    return tools, libs, by, style

def lib_ticks(ax, libs, axis="x"):
    """library tick labels; strand variants indented under their base so the grouping reads."""
    labels = [l for l in libs]
    if axis == "x":
        ax.set_xticks(range(len(libs))); ax.set_xticklabels(labels, rotation=55, ha="right")
    else:
        ax.set_yticks(range(len(libs))); ax.set_yticklabels(labels)

def group_separators(ax, libs, axis="x"):
    """hairlines between base-library groups (a base and its _rev/_unstr variants)."""
    base = [re.sub(r"_(rev|unstr)$", "", l) for l in libs]
    for i in range(1, len(libs)):
        if base[i] != base[i - 1]:
            (ax.axvline if axis == "x" else ax.axhline)(i - 0.5, color=GRID, lw=0.8, zorder=0)

def tool_legend(fig, tools, style, ncol=None):
    handles = [plt.Line2D([], [], ls="", marker=style[t][1], ms=7, mfc=style[t][0], mec=SURFACE,
                          mew=1.2, label=TOOL_LABEL.get(t, t)) for t in tools]
    fig.legend(handles=handles, loc="upper center", ncol=ncol or len(tools), frameon=False,
               bbox_to_anchor=(0.5, 1.0), fontsize=8)

def heat(ax, M, ann, title):
    im = ax.imshow(M, cmap=SEQ, vmin=0, vmax=1, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if ann[i][j]:
                v = M[i, j]
                ax.text(j, i, ann[i][j], ha="center", va="center", fontsize=6,
                        color=SURFACE if (not np.isnan(v) and v > 0.6) else INK)
    ax.set_title(title, loc="left")
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    return im

def save(fig, out, name):
    fig.savefig(os.path.join(out, name), dpi=150, bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- figures
def plot_recall_vs_coverage(runs, out):
    tools, libs, by, _ = axes_of(runs)
    fig, axs = plt.subplots(1, len(tools), figsize=(2.2 * len(tools) + 1.5, 0.28 * len(libs) + 1.6),
                            sharey=True, squeeze=False)
    for ax, t in zip(axs[0], tools):
        M = np.full((len(libs), len(COV_BINS)), np.nan); ann = [[""] * len(COV_BINS) for _ in libs]
        for i, l in enumerate(libs):
            cr = by.get((t, l), {}).get("coverage_recall", {})
            for j, b in enumerate(COV_BINS):
                if b in cr: M[i, j] = frac(cr[b]); ann[i][j] = cr[b]
        im = heat(ax, M, ann, TOOL_LABEL.get(t, t))
        ax.set_xticks(range(len(COV_BINS))); ax.set_xticklabels(COV_BINS)
        ax.set_xlabel("junction reads")
        group_separators(ax, libs, axis="y")
    lib_ticks(axs[0][0], libs, axis="y")
    fig.colorbar(im, ax=axs[0].tolist(), label="recall", shrink=0.6, pad=0.02)
    fig.suptitle("Sensitivity vs junction-read coverage (detected / true BSJs per bin)",
                 x=0.02, ha="left", fontsize=10)
    save(fig, out, "recall_vs_coverage.png")

def plot_category_heatmap(runs, out):
    tools, libs, by, _ = axes_of(runs)
    cats = sorted({c for r in runs for c in r.get("per_category_recall", {})})
    ncol = 2 if len(tools) > 2 else len(tools)
    nrow = int(np.ceil(len(tools) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(0.42 * len(libs) * ncol + 3, (0.3 * len(cats) + 0.6) * nrow + 1.4),
                            sharex=True, sharey=True, squeeze=False)
    for k, t in enumerate(tools):
        ax = axs[k // ncol][k % ncol]
        M = np.full((len(cats), len(libs)), np.nan); ann = [[""] * len(libs) for _ in cats]
        for j, l in enumerate(libs):
            pc = by.get((t, l), {}).get("per_category_recall", {})
            for i, c in enumerate(cats):
                if c in pc: M[i, j] = frac(pc[c]); ann[i][j] = pc[c]
        im = heat(ax, M, ann, TOOL_LABEL.get(t, t))
        group_separators(ax, libs)
        ax.set_yticks(range(len(cats))); ax.set_yticklabels(cats)
        if k + ncol >= len(tools):          # library labels on the bottom panel of each column only
            lib_ticks(ax, libs); ax.tick_params(labelbottom=True)
    for k in range(len(tools), nrow * ncol): axs[k // ncol][k % ncol].axis("off")
    fig.colorbar(im, ax=axs.ravel().tolist(), label="recall", shrink=0.5, pad=0.02)
    fig.suptitle("Per-category recall (detected / true BSJs)", x=0.02, ha="left", fontsize=10)
    save(fig, out, "per_category_recall.png")

def dot_panel(ax, tools, libs, by, style, metric):
    n = len(tools); w = 0.6
    for k, t in enumerate(tools):
        xs, ys = [], []
        for j, l in enumerate(libs):
            v = by.get((t, l), {}).get(metric, np.nan)
            if v is not None and not (isinstance(v, float) and np.isnan(v)):
                xs.append(j + (k - (n - 1) / 2) * w / max(n - 1, 1)); ys.append(v)
        c, m = style[t]
        ax.scatter(xs, ys, s=34, marker=m, color=c, edgecolors=SURFACE, linewidths=1.0, zorder=3)
    ax.set_xlim(-0.6, len(libs) - 0.4)
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    group_separators(ax, libs)

def plot_prf(runs, out):
    tools, libs, by, style = axes_of(runs)
    metrics = [("precision", "Precision"), ("recall", "Recall"), ("f1", "F1"),
               ("specificity", "Specificity (decoys rejected)")]
    fig, axs = plt.subplots(2, 2, figsize=(0.5 * len(libs) * 2 + 2, 7), sharex=True)
    for ax, (m, title) in zip(axs.ravel(), metrics):
        dot_panel(ax, tools, libs, by, style, m)
        ax.set_title(title, loc="left"); ax.set_ylim(-0.03, 1.05)
    for ax in axs[1]: lib_ticks(ax, libs)
    tool_legend(fig, tools, style)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, out, "precision_recall_f1.png")

def plot_quant(runs, out):
    tools, libs, by, style = axes_of(runs)
    fig, ax = plt.subplots(figsize=(0.5 * len(libs) + 2, 3.8))
    dot_panel(ax, tools, libs, by, style, "quant_spearman")
    vals = [r.get("quant_spearman") for r in runs if r.get("quant_spearman") == r.get("quant_spearman")]
    ax.set_ylim(min([0.5] + [v - 0.05 for v in vals if v is not None]), 1.02)
    ax.set_ylabel("Spearman r")
    ax.set_title("Quantification: detected vs expected junction reads", loc="left")
    lib_ticks(ax, libs)
    missing = [TOOL_LABEL.get(t, t) for t in tools
               if all(by.get((t, l), {}).get("quant_spearman") != by.get((t, l), {}).get("quant_spearman")
                      or by.get((t, l), {}).get("quant_spearman") is None for l in libs)]
    if missing:
        ax.text(0.01, 0.04, "no read counts reported: " + ", ".join(missing),
                transform=ax.transAxes, fontsize=7, color=INK2)
    tool_legend(fig, tools, style)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save(fig, out, "quantification.png")

def plot_strand(runs, out):
    tools, libs, by, _ = axes_of(runs)
    fig, axs = plt.subplots(len(tools), 1, figsize=(0.5 * len(libs) + 2.5, 1.9 * len(tools) + 1.2),
                            sharex=True, squeeze=False)
    for ax, t in zip(axs[:, 0], tools):
        x = np.arange(len(libs))
        ok = np.array([by.get((t, l), {}).get("strand_correct", 0) for l in libs], float)
        bad = np.array([by.get((t, l), {}).get("strand_wrong", 0) for l in libs], float)
        unk = np.array([by.get((t, l), {}).get("strand_unknown", 0) for l in libs], float)
        # denominator is every true BSJ (P), so undetected BSJs show as the empty remainder
        P = np.array([by.get((t, l), {}).get("P", np.nan) for l in libs], float)
        miss = P - (ok + bad + unk)
        kw = dict(width=0.72, edgecolor=SURFACE, linewidth=1.0)
        ax.bar(x, ok / P, color=GOOD, label="correct strand", **kw)
        ax.bar(x, bad / P, bottom=ok / P, color=CRITICAL, label="wrong strand", **kw)
        ax.bar(x, unk / P, bottom=(ok + bad) / P, color=NONE, label="no strand reported", **kw)
        ax.bar(x, miss / P, bottom=(ok + bad + unk) / P, color=SURFACE, edgecolor=AXIS,
               linewidth=0.8, width=0.72 - 0.04, label="not detected")
        for j in range(len(libs)):
            if not np.isnan(P[j]):
                ax.text(j, 1.03, f"{ok[j]:.0f}/{P[j]:.0f}", ha="center", va="bottom", fontsize=6.5, color=INK2)
        ax.set_ylim(0, 1.18); ax.set_yticks([0, 0.5, 1])
        ax.set_title(TOOL_LABEL.get(t, t), loc="left")
        ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
        group_separators(ax, libs)
    lib_ticks(axs[-1, 0], libs)
    fig.supylabel("fraction of true BSJs", color=INK2, fontsize=8)
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.suptitle("Declared strand vs truth, over all true BSJs (labels: correct strand / true BSJs)",
                 x=0.02, y=1.03, ha="left", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, out, "strand_accuracy.png")

def load_truth(dataset):
    """[(bsj_id, category, strand)] in truth_bsj.tsv order, or [] if not found."""
    path = os.path.join(dataset, "truth_bsj.tsv")
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(hdr)}
        return [(f[ix["bsj_id"]], f[ix["category"]], f[ix["strand"]])
                for f in (ln.rstrip("\n").split("\t") for ln in fh) if len(f) > 1]

def plot_strand_heatmap(runs, truth, out):
    """Outcome of every true BSJ in every run. The run JSON lists missed / strand_wrong /
    strand_unknown ids; every other true BSJ was detected on the correct strand."""
    tools, libs, by, _ = axes_of(runs)
    OUTCOMES = [("correct strand", GOOD), ("wrong strand", CRITICAL),
                ("no strand reported", INK2), ("not detected", "#ecebe6")]
    cmap = ListedColormap([c for _, c in OUTCOMES])
    ids = [b for b, _, _ in truth]
    fig, axs = plt.subplots(1, len(tools), figsize=(0.2 * len(libs) * len(tools) + 3, 0.17 * len(ids) + 2.2),
                            sharey=True, squeeze=False)
    for ax, t in zip(axs[0], tools):
        M = np.full((len(ids), len(libs)), np.nan)
        for j, l in enumerate(libs):
            r = by.get((t, l))
            if r is None: continue
            miss, bad, unk = set(r.get("missed", [])), set(r.get("strand_wrong_ids", [])), set(r.get("strand_unknown_ids", []))
            for i, b in enumerate(ids):
                M[i, j] = 3 if b in miss else 1 if b in bad else 2 if b in unk else 0
        ax.imshow(M, cmap=cmap, vmin=-0.5, vmax=3.5, aspect="auto", interpolation="nearest")
        # cell gaps: surface-coloured gridlines between every row and column
        ax.set_xticks(np.arange(-0.5, len(libs)), minor=True); ax.set_yticks(np.arange(-0.5, len(ids)), minor=True)
        ax.grid(which="minor", color=SURFACE, lw=1.0); ax.tick_params(which="minor", length=0)
        cats = [c for _, c, _ in truth]
        for i in range(1, len(ids)):                   # heavier line between test-case categories
            if cats[i][0] != cats[i - 1][0]:
                ax.axhline(i - 0.5, color=INK2, lw=0.6)
        group_separators(ax, libs)
        for line in ax.lines: line.set_zorder(3)
        ax.set_title(TOOL_LABEL.get(t, t), loc="left")
        for sp in ax.spines.values(): sp.set_visible(False)
        ax.tick_params(length=0)
        lib_ticks(ax, libs)
    axs[0][0].set_yticks(range(len(ids)))
    axs[0][0].set_yticklabels([f"{b.replace('BSJ_', '')} ({s})" for b, _, s in truth], fontsize=6.5)
    fig.legend(handles=[Patch(facecolor=c, edgecolor=AXIS if k == 3 else c, label=n)
                        for k, (n, c) in enumerate(OUTCOMES)],
               loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.suptitle("Outcome per true BSJ (rows: truth id and strand; lines separate categories A–M)",
                 x=0.02, y=1.03, ha="left", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, out, "strand_heatmap.png")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default=None)
    ap.add_argument("--dataset", default="../dataset", help="for truth_bsj.tsv (strand_heatmap)")
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
    truth = load_truth(a.dataset)
    if truth:
        plot_strand_heatmap(runs, truth, out)
    else:
        print(f"no truth_bsj.tsv under {a.dataset}/ — skipping strand_heatmap.png")
    print(f"Wrote {5 + bool(truth)} plots for {len(runs)} run(s) -> {out}/")

if __name__ == "__main__":
    main()
