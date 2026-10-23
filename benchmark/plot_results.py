#!/usr/bin/env python3
"""
Plot benchmark results produced by evaluate.py (one JSON per tool x library x variant).

Every figure is small multiples with a shared library axis, so runs stay readable as
tools and libraries are added. Scaled libraries (scale > 1 in libraries.tsv) appear only in
runtime.png; summary.tsv still lists them. Strandedness variants (*_rev, *_unstr) sit next to their
base library. Tools keep a fixed color in every figure.

Outputs (each as PNG and SVG, into --out):
  recall_vs_coverage.png        recall per junction-read coverage bin: library x bin heatmap, one panel per tool
  per_category_recall.png       recall per test-case category (A-M): category x library heatmap, one panel per tool
  precision_recall_f1.png       precision / recall / F1 / specificity: one radar per metric, spokes = libraries
  quantification.png            Spearman of detected vs expected junction reads: radar, spokes = libraries
  precision_recall_f1_scatter.png
  quantification_scatter.png    the same two as dot plots: libraries on x, one dot (colour + marker) per tool
  runtime.png                   wall time, CPU time and peak memory per tool, SCALED libraries only (scale > 1
                                in libraries.tsv), from the GNU time files in <results>/timing/ (skipped if
                                there are none); numbers also go to <results>/runtime_summary.tsv
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

# Fixed tool identity, the same in every figure: colour (categorical slots) + marker. A preset
# variant of a tool (VARIANTS) keeps that tool's colour and marker but is drawn hollow / dashed,
# so "same tool, other settings" reads at a glance and no extra hue is needed.
TOOL_ORDER = ["bsj_finder", "bsj_finder_fast", "ciri2", "circexplorer2_bed", "circexplorer2_known"]
TOOL_LABEL = {"bsj_finder": "bsj_finder (sensitive)", "bsj_finder_fast": "bsj_finder (fast)",
              "ciri2": "CIRI2",
              "circexplorer2_bed": "CIRCexplorer2 (de novo)",
              "circexplorer2_known": "CIRCexplorer2 (known)"}
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
MARKERS = ["o", "s", "D", "^", "v", "P", "X", "h"]   # secondary encoding in the scatter plots
TOOL_STYLE = {"bsj_finder": (SLOTS[0], "o"), "bsj_finder_fast": (SLOTS[0], "o"),
              "ciri2": (SLOTS[1], "s"), "circexplorer2_bed": (SLOTS[2], "D"),
              "circexplorer2_known": (SLOTS[3], "^")}
VARIANTS = {"bsj_finder_fast"}

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
    extra = [t for t in tools if t not in TOOL_STYLE]          # unknown tools: next free slots
    style = {t: TOOL_STYLE.get(t) or (SLOTS[(4 + extra.index(t)) % len(SLOTS)], MARKERS[(4 + extra.index(t)) % len(MARKERS)])
             for t in tools}
    return tools, libs, by, style

def dot_kw(t, style, size=34):
    """scatter kwargs: filled marker, or hollow for a preset variant"""
    c, m = style[t]
    if t in VARIANTS:
        return dict(s=size, marker=m, facecolors=SURFACE, edgecolors=c, linewidths=1.6)
    return dict(s=size, marker=m, color=c, edgecolors=SURFACE, linewidths=1.0)

def line_kw(t, style):
    return dict(color=style[t][0], lw=2, ls="--" if t in VARIANTS else "-")

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

def tool_legend(fig, tools, style, ncol=None, notes=None, fontsize=8, markers=False):
    """legend in fixed tool order: line swatches, or marker swatches for the scatter plots;
    notes = {tool: text} appended to its label."""
    notes = notes or {}
    label = lambda t: TOOL_LABEL.get(t, t) + (f" ({notes[t]})" if t in notes else "")
    if markers:
        handles = [plt.Line2D([], [], ls="", marker=style[t][1], ms=7, label=label(t),
                              mfc=SURFACE if t in VARIANTS else style[t][0],
                              mec=style[t][0] if t in VARIANTS else SURFACE, mew=1.6 if t in VARIANTS else 1.2)
                   for t in tools]
    else:
        handles = [plt.Line2D([], [], color=style[t][0], lw=2.5, ls="--" if t in VARIANTS else "-",
                              label=label(t)) for t in tools]
    fig.legend(handles=handles, loc="upper center", ncol=ncol or len(tools), frameon=False,
               bbox_to_anchor=(0.5, 1.0), fontsize=fontsize)

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
    """write <name>.png and a vector <name>.svg alongside it."""
    path = os.path.join(out, name)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    fig.savefig(os.path.splitext(path)[0] + ".svg", bbox_inches="tight")
    plt.close(fig)

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

def _num(v):
    return v is not None and not (isinstance(v, float) and np.isnan(v))

RADAR_MIN = 0.25   # radial axis floor (the centre); every current value is above it

def radar_panel(ax, tools, libs, by, style, metric, title, fs=7):
    """radar: one spoke per library (clockwise from 12 o'clock), one closed outline per tool.
    Lines only (fills of 4 overlapping tools turn muddy). Radius RADAR_MIN..1; values below
    the floor are clipped to the centre. The first tool in the fixed order (bsj_finder) is
    drawn on top. Tools with no value at all for this metric are skipped; a missing single
    value leaves a gap in that outline. fs = base font size."""
    ang = np.linspace(0, 2 * np.pi, len(libs), endpoint=False)
    ax.set_theta_offset(np.pi / 2); ax.set_theta_direction(-1)
    for k, t in enumerate(tools):
        v = np.array([by.get((t, l), {}).get(metric) if _num(by.get((t, l), {}).get(metric)) else np.nan
                      for l in libs], float)
        if np.isnan(v).all():
            continue
        v = np.maximum(v, RADAR_MIN)
        a, r = np.append(ang, ang[0]), np.append(v, v[0])
        ax.plot(a, r, **line_kw(t, style), zorder=3 + len(tools) - k)   # earlier tools on top
    ax.set_xticks(ang); ax.set_xticklabels(libs, fontsize=fs, color=INK2)
    ax.tick_params(axis="x", pad=6)
    for lab, th in zip(ax.get_xticklabels(), ang):     # anchor each label on the side away from the ring
        x = np.sin(th)                                  # clockwise from 12 o'clock: +x is the right half
        lab.set_horizontalalignment("left" if x > 0.01 else "right" if x < -0.01 else "center")
    ax.set_ylim(RADAR_MIN, 1); ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0.25", "0.5", "0.75", "1"], fontsize=fs - 1, color=MUTED)
    ax.set_rlabel_position(180 / len(libs))           # radial labels between the first two spokes
    ax.grid(color=GRID, lw=0.6); ax.spines["polar"].set_color(AXIS)
    ax.set_title(title, loc="left", pad=18, fontsize=fs + 2)

def plot_prf(runs, out):
    tools, libs, by, style = axes_of(runs)
    metrics = [("precision", "Precision"), ("recall", "Recall"), ("f1", "F1"),
               ("specificity", "Specificity (decoys rejected)")]
    fig, axs = plt.subplots(2, 2, figsize=(16, 14), subplot_kw=dict(projection="polar"))
    for ax, (m, title) in zip(axs.ravel(), metrics):
        radar_panel(ax, tools, libs, by, style, m, title, fs=11)
    tool_legend(fig, tools, style, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96), h_pad=3, w_pad=8)
    save(fig, out, "precision_recall_f1.png")

def plot_quant(runs, out):
    tools, libs, by, style = axes_of(runs)
    fig, ax = plt.subplots(figsize=(6.5, 6.8), subplot_kw=dict(projection="polar"))
    radar_panel(ax, tools, libs, by, style, "quant_spearman",
                "Quantification: Spearman r, detected vs expected junction reads")
    missing = {t: "no read counts reported" for t in tools
               if not any(_num(by.get((t, l), {}).get("quant_spearman")) for l in libs)}
    tool_legend(fig, tools, style, ncol=2, notes=missing)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save(fig, out, "quantification.png")

def dot_panel(ax, tools, libs, by, style, metric):
    """one dot per tool per library (tools side by side within a library); bsj_finder on top."""
    n = len(tools); w = 0.6
    for k, t in enumerate(tools):
        xs, ys = [], []
        for j, l in enumerate(libs):
            v = by.get((t, l), {}).get(metric)
            if _num(v):
                xs.append(j + (k - (n - 1) / 2) * w / max(n - 1, 1)); ys.append(v)
        ax.scatter(xs, ys, **dot_kw(t, style), zorder=3 + len(tools) - k)
    ax.set_xlim(-0.6, len(libs) - 0.4)
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    group_separators(ax, libs)

def plot_prf_scatter(runs, out):
    tools, libs, by, style = axes_of(runs)
    metrics = [("precision", "Precision"), ("recall", "Recall"), ("f1", "F1"),
               ("specificity", "Specificity (decoys rejected)")]
    fig, axs = plt.subplots(2, 2, figsize=(0.5 * len(libs) * 2 + 2, 7), sharex=True)
    for ax, (m, title) in zip(axs.ravel(), metrics):
        dot_panel(ax, tools, libs, by, style, m)
        ax.set_title(title, loc="left"); ax.set_ylim(-0.03, 1.05)
    for ax in axs[1]: lib_ticks(ax, libs)
    tool_legend(fig, tools, style, markers=True)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, out, "precision_recall_f1_scatter.png")

def plot_quant_scatter(runs, out):
    tools, libs, by, style = axes_of(runs)
    fig, ax = plt.subplots(figsize=(0.5 * len(libs) + 2, 3.8))
    dot_panel(ax, tools, libs, by, style, "quant_spearman")
    vals = [r.get("quant_spearman") for r in runs if _num(r.get("quant_spearman"))]
    ax.set_ylim(min([0.5] + [v - 0.05 for v in vals]), 1.02)
    ax.set_ylabel("Spearman r")
    ax.set_title("Quantification: detected vs expected junction reads", loc="left")
    lib_ticks(ax, libs)
    missing = {t: "no read counts reported" for t in tools
               if not any(_num(by.get((t, l), {}).get("quant_spearman")) for l in libs)}
    tool_legend(fig, tools, style, ncol=2 if missing else None, notes=missing, markers=True)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    save(fig, out, "quantification_scatter.png")

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

# GNU time files from the Snakefile: <results>/timing/<tool>/<lib>.<rule>.tsv, one-off index
# builds as shared.<rule>.tsv; one line per timed command: wall_s user_s sys_s max_rss_kb.
BENCH_TOOL = {"bsj_finder": ("bsj_finder", "bsj_finder (sensitive)"),
              "bsj_finder_fast": ("bsj_finder_fast", "bsj_finder (fast)"), "ciri2": ("ciri2", "CIRI2 (BWA-MEM + CIRI2)"),
              "circexplorer2": ("circexplorer2_bed", "CIRCexplorer2 (STAR + parse + annotate)")}

def load_benchmarks(results):
    """{tool: {lib|'shared': {rule: {wall, cpu, rss}}}}; wall/cpu in s summed over the rule's
    timed commands, rss in MB (largest)."""
    out = {}
    for f in sorted(glob.glob(os.path.join(results, "timing", "*", "*.tsv"))):
        tool = os.path.basename(os.path.dirname(f))
        lib, rule = os.path.basename(f)[:-4].split(".", 1)
        rows = []
        for ln in open(f):
            p = ln.rstrip("\n").split("\t")
            try:
                rows.append([float(x) for x in p])        # skips "Command exited with ..." lines
            except ValueError:
                continue
        rows = [r for r in rows if len(r) == 4]
        if not rows:
            continue
        out.setdefault(tool, {}).setdefault(lib, {})[rule] = dict(
            wall=sum(r[0] for r in rows), cpu=sum(r[1] + r[2] for r in rows), rss=max(r[3] for r in rows) / 1024)
    return out

def scaled_libraries(dataset):
    """libraries with scale > 1 in libraries.tsv ({name: scale}); None if there is no manifest."""
    path = os.path.join(dataset, "libraries.tsv")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(hdr)}
        if "scale" not in ix:
            return {}
        rows = [ln.rstrip("\n").split("\t") for ln in fh]
    return {f[ix["library"]]: int(f[ix["scale"]]) for f in rows if len(f) > ix["scale"] and int(f[ix["scale"]]) > 1}

def is_scaled_fn(dataset):
    """predicate: is this a scaled (runtime-only) library? From libraries.tsv's scale column,
    or the _x<N> name suffix when there is no manifest."""
    scaled = scaled_libraries(dataset)
    if scaled is None:
        return lambda l: bool(re.search(r"_x\d+$", l))
    return lambda l: l in scaled

def plot_runtime(runs, bench, results, out, dataset):
    """Runtime on the SCALED libraries only (100x reads): on the small libraries every step takes
    seconds and measures start-up, not work. Per tool x library: wall and CPU (user + sys) time
    summed over the tool's steps, peak RSS = largest single process. One-off index builds are
    kept out of the per-library numbers and listed separately."""
    _, _, _, run_style = axes_of(runs)
    tools = [t for t in BENCH_TOOL if t in bench] + sorted(t for t in bench if t not in BENCH_TOOL)
    scaled = scaled_libraries(dataset)
    is_scaled = is_scaled_fn(dataset)
    libs = sorted({l for t in tools for l in bench[t] if l != "shared" and is_scaled(l)}, key=lib_key)
    if not libs:                       # keep the Snakemake output; say why it is empty
        fig, ax = plt.subplots(figsize=(6, 1.5)); ax.axis("off")
        ax.text(0.5, 0.5, "No scaled libraries timed.\nGenerate them with: python3 generate_dataset.py --scaled",
                ha="center", va="center", color=INK2)
        save(fig, out, "runtime.png")
        return
    key = lambda t: BENCH_TOOL.get(t, (t,))[0]          # timing dir -> the tool key used for styles
    style = {t: run_style.get(key(t)) or TOOL_STYLE.get(key(t)) or (SLOTS[(4 + i) % len(SLOTS)], MARKERS[i % 8])
             for i, t in enumerate(tools)}
    label = {t: BENCH_TOOL.get(t, (t, t))[1] for t in tools}
    agg = {}   # (tool, lib) -> wall, cpu, rss
    for t in tools:
        for l, steps in bench[t].items():
            st = list(steps.values())
            agg[(t, l)] = dict(wall=sum(x["wall"] for x in st), cpu=sum(x["cpu"] for x in st),
                               rss=max(x["rss"] for x in st))
    panels = [("wall", "Wall time (s)", False),
              ("cpu", "CPU time (user + system, s)", False),
              ("rss", "Peak memory (largest process RSS, MB)", False)]
    fig, axs = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 4.6))
    n = len(tools); w = 0.5
    for ax, (key, title, logy) in zip(axs, panels):
        for k, t in enumerate(tools):
            xs = [j + (k - (n - 1) / 2) * w / max(n - 1, 1) for j, l in enumerate(libs) if (t, l) in agg]
            ys = [agg[(t, l)][key] for l in libs if (t, l) in agg]
            ax.scatter(xs, ys, **dot_kw(t, style, 60), zorder=3 + n - k)
        if logy: ax.set_yscale("log")
        else: ax.set_ylim(0, ax.get_ylim()[1] * 1.08)
        ax.set_title(title, loc="left")
        ax.set_xlim(-0.6, len(libs) - 0.4)
        ax.grid(axis="y", which="major", color=GRID, lw=0.6); ax.set_axisbelow(True)
        group_separators(ax, libs)
        ax.set_xticks(range(len(libs))); ax.set_xticklabels(libs, rotation=0)
    handles = [plt.Line2D([], [], ls="", marker=style[t][1], ms=7, label=label[t],
                          mfc=SURFACE if t in VARIANTS else style[t][0],
                          mec=style[t][0] if t in VARIANTS else SURFACE, mew=1.6 if t in VARIANTS else 1.2)
               for t in tools]
    fig.legend(handles=handles, loc="upper center", ncol=len(tools), frameon=False,
               bbox_to_anchor=(0.5, 1.0), fontsize=8)
    shared = [f"{label[t].split(' (')[0]} {r.replace('_', ' ')}: {v['wall']:.1f} s wall, {v['cpu']:.1f} s CPU, "
              f"{v['rss']:.0f} MB" for t in tools for r, v in bench[t].get("shared", {}).items()]
    if shared:
        fig.text(0.01, -0.01, "One-off builds, not included above: " + "; ".join(shared),
                 fontsize=7, color=INK2, ha="left", va="top", wrap=True)
    ks = sorted({scaled[l] for l in libs}) if scaled else []
    scale = f" ({ks[0]}x the reads of their base library)" if len(ks) == 1 else ""
    fig.suptitle(f"Runtime on the scaled libraries{scale}; 8 threads, one job at a time",
                 x=0.01, y=0.9, ha="left", fontsize=9, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    save(fig, out, "runtime.png")
    # summary table: one row per tool x library, plus the one-off builds as library "shared"
    with open(os.path.join(results, "runtime_summary.tsv"), "w") as fh:
        fh.write("tool\tlibrary\twall_s\tcpu_s\tpeak_rss_mb\n")
        for t in tools:
            for l in libs + ["shared"]:
                if (t, l) in agg:
                    x = agg[(t, l)]
                    fh.write(f"{t}\t{l}\t{x['wall']:.2f}\t{x['cpu']:.2f}\t{x['rss']:.1f}\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default=None)
    ap.add_argument("--dataset", default="../dataset", help="for truth_bsj.tsv (strand_heatmap)")
    a = ap.parse_args()
    out = a.out or os.path.join(a.results, "plots")
    os.makedirs(out, exist_ok=True)
    all_runs = load(a.results)
    if not all_runs:
        print(f"No JSON results under {a.results}/ — run the pipelines first."); return
    # scaled libraries are a runtime workload: only runtime.png shows them
    is_scaled = is_scaled_fn(a.dataset)
    runs = [r for r in all_runs if not is_scaled(r["library"])]
    if not runs:                       # a results tree of only scaled libraries (e.g. another genome)
        runs = all_runs
    plot_recall_vs_coverage(runs, out)
    plot_category_heatmap(runs, out)
    plot_prf(runs, out)
    plot_quant(runs, out)
    plot_prf_scatter(runs, out)
    plot_quant_scatter(runs, out)
    plot_strand(runs, out)
    truth = load_truth(a.dataset)
    if truth:
        plot_strand_heatmap(runs, truth, out)
    else:
        print(f"no truth_bsj.tsv under {a.dataset}/ — skipping strand_heatmap.png")
    bench = load_benchmarks(a.results)
    if bench:
        plot_runtime(all_runs, bench, a.results, out, a.dataset)
    print(f"Wrote {7 + bool(truth) + bool(bench)} plots for {len(runs)} run(s) "
          f"(+ {len(all_runs) - len(runs)} scaled, runtime only) -> {out}/")

if __name__ == "__main__":
    main()
