#!/usr/bin/env python3
"""
Runtime / memory scaling across library size and genome size, from the GNU time files of
several results trees (one per genome; see the Snakefile's `genome` / `results` config).

Outputs (PNG + SVG, into --out):
  runtime_by_library_size   wall and CPU time per tool vs reads per library (log-log), on the
                            --runtime-genome tree; paired-end and single-end panels
  memory_by_genome_size     peak memory per tool vs genome size: the per-library runs at the
                            --memory-scale library, and the one-off index builds
A --variant adds another run of one tool (e.g. bsj_finder's `fast` preset) as an extra line:
same colour as the tool, dashed, hollow markers. It names one results tree per genome size, so
it appears on both plots; --rename relabels the base tool's line to match.
Numbers also go to <out>/scaling_summary.tsv.

Usage:
  python3 plot_scaling.py --tree results:1.5e6 --tree results_genome_1gbp:1e9 \
      --runtime-genome results_genome_1gbp --dataset ../dataset --out results/plots \
      --rename "bsj_finder=bsj_finder (sensitive)" \
      --variant "bsj_finder (fast)|bsj_finder|results_k21:1.5e6,results_genome_1gbp_k21:1e9"
"""
import argparse, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_results import (load_benchmarks, scaled_libraries, save, BENCH_TOOL, SLOTS, MARKERS,
                          INK2, GRID, SURFACE)

BASE_READS = 3568          # reads (SE) / pairs (PE) in every base library
TOOLS = [t for t in BENCH_TOOL]                     # fixed order: bsj_finder, ciri2, circexplorer2
STYLE = {"bsj_finder": (SLOTS[0], MARKERS[0]), "ciri2": (SLOTS[1], MARKERS[1]),
         "circexplorer2": (SLOTS[2], MARKERS[2])}
LABEL = {t: BENCH_TOOL[t][1] for t in TOOLS}
INDEX = {"bsj_finder": "bsj_finder_index", "ciri2": "bwa_index", "circexplorer2": "star_index"}

def per_lib(bench):
    """{(tool, lib): wall, cpu, rss} summed / maxed over each tool's per-library steps"""
    out = {}
    for t, libs in bench.items():
        for l, steps in libs.items():
            if l == "shared": continue
            st = list(steps.values())
            out[(t, l)] = dict(wall=sum(x["wall"] for x in st), cpu=sum(x["cpu"] for x in st),
                               rss=max(x["rss"] for x in st))
    return out

def fmt_bp(n):
    return f"{n / 1e9:g} Gbp" if n >= 1e9 else f"{n / 1e6:g} Mbp"

def style_ax(ax, title):
    ax.set_title(title, loc="left")
    ax.grid(which="major", color=GRID, lw=0.6); ax.set_axisbelow(True)

def line(ax, t, xs, ys, variant=False):
    c, m = STYLE[t]
    ax.plot(xs, ys, color=c, marker=m, ms=7, lw=2, ls="--" if variant else "-",
            mfc=SURFACE if variant else c, mec=c if variant else SURFACE, mew=1.5)

def legend(fig, variants=()):
    h = [plt.Line2D([], [], color=STYLE[t][0], marker=STYLE[t][1], ms=7, mec=SURFACE, lw=2, label=LABEL[t])
         for t in TOOLS]
    h += [plt.Line2D([], [], color=STYLE[t][0], marker=STYLE[t][1], ms=7, mfc=SURFACE, mec=STYLE[t][0], mew=1.5,
                     lw=2, ls="--", label=lab) for lab, t, _ in variants]
    fig.legend(handles=h, loc="upper center", ncol=len(h), frameon=False, bbox_to_anchor=(0.5, 1.0))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", action="append", required=True,
                    help="results_dir:genome_size_bp, once per genome (e.g. results_genome_1gbp:1e9)")
    ap.add_argument("--runtime-genome", required=True, help="results dir used for runtime vs library size")
    ap.add_argument("--memory-scale", type=int, default=100, help="library scale compared across genomes")
    ap.add_argument("--dataset", default="../dataset")
    ap.add_argument("--out", required=True)
    ap.add_argument("--variant", action="append", default=[],
                    help='"label|tool|dir:genome_bp[,dir:genome_bp...]" -- another run of `tool`, '
                         "one results tree per genome size")
    ap.add_argument("--rename", action="append", default=[], help='"tool=label" for a base tool line')
    a = ap.parse_args()
    for r in a.rename:
        t, lab = r.split("=", 1); LABEL[t] = lab
    variants = []                                       # (label, tool, {dir: genome_bp})
    for v in a.variant:
        lab, t, trees_ = v.split("|", 2)
        variants.append((lab, t, {d: float(sz) for d, sz in (x.rsplit(":", 1) for x in trees_.split(","))}))
    vbench = {d: load_benchmarks(d) for _, _, tv in variants for d in tv}
    vdata = {d: per_lib(b) for d, b in vbench.items()}
    def vtree(tv, size):                                # the variant's tree for one genome size
        return next((d for d, sz in tv.items() if sz == size), None)
    os.makedirs(a.out, exist_ok=True)
    trees = [(d, float(sz)) for d, sz in (t.rsplit(":", 1) for t in a.tree)]
    scaled = scaled_libraries(a.dataset) or {}
    data = {d: (load_benchmarks(d), per_lib(load_benchmarks(d)), size) for d, size in trees}
    rows = []

    # ---- runtime vs library size (one genome)
    bench, agg, gsize = data[a.runtime_genome]
    fig, axs = plt.subplots(2, 2, figsize=(10, 7.5), sharex=True)
    for col, lay in enumerate(["pe", "se"]):
        libs = sorted((l for l in scaled if l.startswith(lay)), key=lambda l: scaled[l])
        for row, (key, name) in enumerate([("wall", "Wall time (s)"), ("cpu", "CPU time, user + system (s)")]):
            ax = axs[row][col]
            for t in TOOLS:
                pts = [(scaled[l] * BASE_READS, agg[(t, l)][key]) for l in libs if (t, l) in agg]
                if pts:
                    line(ax, t, *zip(*pts))
            for _, t, tv in variants:
                d = vtree(tv, gsize)
                pts = [(scaled[l] * BASE_READS, vdata[d][(t, l)][key]) for l in libs if d and (t, l) in vdata[d]]
                if pts:
                    line(ax, t, *zip(*pts), variant=True)
            ax.set_xscale("log"); ax.set_yscale("log")
            style_ax(ax, f"{name} — {'paired-end' if lay == 'pe' else 'single-end'}")
            if row == 1:
                ax.set_xlabel("reads per library" if lay == "se" else "read pairs per library")
        for l in libs:
            for t in TOOLS:
                if (t, l) in agg:
                    rows.append((os.path.basename(a.runtime_genome.rstrip("/")), gsize, LABEL[t], l, scaled[l] * BASE_READS,
                                 agg[(t, l)]["wall"], agg[(t, l)]["cpu"], agg[(t, l)]["rss"]))
            for lab, t, tv in variants:
                d = vtree(tv, gsize)
                if d and (t, l) in vdata[d]:
                    x = vdata[d][(t, l)]
                    rows.append((os.path.basename(d.rstrip("/")), gsize, lab, l, scaled[l] * BASE_READS,
                                 x["wall"], x["cpu"], x["rss"]))
    legend(fig, variants)
    fig.suptitle(f"Runtime vs library size on the {fmt_bp(gsize)} genome (8 threads, one job at a time; "
                 "one-off index builds excluded)", x=0.01, y=0.955, ha="left", fontsize=9, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, a.out, "runtime_by_library_size.png")

    # ---- peak memory vs genome size (per-library runs at one scale, and the index builds)
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.6))
    sizes = sorted(size for _, size in trees)
    panels = [(f"pe150_lowerr_x{a.memory_scale}", f"Per library: paired-end x{a.memory_scale}"),
              (f"se150_lowerr_x{a.memory_scale}", f"Per library: single-end x{a.memory_scale}"),
              ("shared", "One-off index build")]
    for ax, (lib, title) in zip(axs, panels):
        for t in TOOLS:
            pts = []
            for d, (bench, agg, size) in data.items():
                if lib == "shared":
                    v = bench.get(t, {}).get("shared", {}).get(INDEX[t])
                    if v: pts.append((size, v["rss"] / 1024)); rows.append((d, size, LABEL[t], f"index:{INDEX[t]}", 0,
                                                                            v["wall"], v["cpu"], v["rss"]))
                elif (t, lib) in agg:
                    pts.append((size, agg[(t, lib)]["rss"] / 1024))
                    if d != a.runtime_genome:
                        x = agg[(t, lib)]
                        rows.append((d, size, LABEL[t], lib, scaled.get(lib, 0) * BASE_READS, x["wall"], x["cpu"], x["rss"]))
            if pts:
                line(ax, t, *zip(*sorted(pts)))
        for lab, t, tv in variants:         # one point per genome size the variant was run on
            pts = []
            for d, size in tv.items():
                if lib == "shared":
                    v = vbench[d].get(t, {}).get("shared", {}).get(INDEX[t])
                    if v:
                        pts.append((size, v["rss"] / 1024))
                        rows.append((d, size, lab, f"index:{INDEX[t]}", 0, v["wall"], v["cpu"], v["rss"]))
                elif (t, lib) in vdata[d]:
                    pts.append((size, vdata[d][(t, lib)]["rss"] / 1024))
                    if size != gsize:
                        x = vdata[d][(t, lib)]
                        rows.append((d, size, lab, lib, scaled.get(lib, 0) * BASE_READS, x["wall"], x["cpu"], x["rss"]))
            if pts:
                line(ax, t, *zip(*sorted(pts)), variant=True)
        ax.set_xscale("log"); ax.set_ylim(0, None)
        ax.set_xticks(sizes); ax.set_xticklabels([fmt_bp(s) for s in sizes]); ax.minorticks_off()
        ax.set_xlabel("genome size"); style_ax(ax, title)
    axs[0].set_ylabel("peak memory (GB, largest process)")
    legend(fig, variants)
    fig.suptitle("Peak memory vs genome size", x=0.01, y=0.92, ha="left", fontsize=9, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save(fig, a.out, "memory_by_genome_size.png")

    with open(os.path.join(a.out, "scaling_summary.tsv"), "w") as fh:
        fh.write("results\tgenome_bp\ttool\tlibrary\treads\twall_s\tcpu_s\tpeak_rss_mb\n")
        for r in rows:
            fh.write("\t".join([r[0], f"{r[1]:.0f}", r[2], r[3], str(r[4]), f"{r[5]:.2f}", f"{r[6]:.2f}",
                                f"{r[7]:.1f}"]) + "\n")
    print(f"wrote runtime_by_library_size + memory_by_genome_size (.png/.svg) and scaling_summary.tsv -> {a.out}")

if __name__ == "__main__":
    main()
