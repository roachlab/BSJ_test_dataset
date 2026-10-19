#!/usr/bin/env python3
"""
Evaluate circRNA / BSJ detection output against the tiny test-dataset ground truth.

Tool-agnostic: a built-in parser normalises each tool's native output to
(chrom, left_bp, right_bp, junction_reads, strand), then matches against truth_bsj.tsv with a
per-junction tolerance window and scores the call set.

Matching is strand-blind (a BSJ is found if its breakpoints are found); the strand the tool
declared for each matched call is then scored separately (strand_* metrics), so a caller that
detects the junction but reports the wrong strand -- e.g. the read strand of a reverse-stranded
library instead of the transcript strand -- is caught without changing recall/precision.

Supported --tool parsers:
  circexplorer2_bed    CIRCexplorer2 'parse' back_spliced_junction.bed (de novo BSJs)
  circexplorer2_known  CIRCexplorer2 'annotate' circularRNA_known.txt (annotation-guided)
  ciri2                CIRI2 .ciri output
  bed                  generic BED6 (chrom,start,end,name,reads,strand) -- use for YOUR pipeline

Usage:
  evaluate.py --tool ciri2 --calls out.ciri --dataset ../dataset --library pe150_lowerr \
              --out-prefix results/ciri2_pe150_lowerr
"""
import argparse, json, os, math, sys

# ---------------------------------------------------------------- parsers
def _ints(*xs): return [int(x) for x in xs]

def _strand(f, i):
    """strand from column i, or '.' when the tool reports none / something else."""
    return f[i] if len(f) > i and f[i] in ("+", "-") else "."

def parse_bed(path):
    calls = []
    for ln in open(path):
        if not ln.strip() or ln.startswith(("#", "track")): continue
        f = ln.rstrip("\n").split("\t")
        chrom, start, end = f[0], int(f[1]), int(f[2])
        reads = 0
        if len(f) > 4:
            try: reads = int(float(f[4]))
            except ValueError: reads = 0
        strand = _strand(f, 5)
        calls.append((chrom, start, end - 1, reads, strand))   # left_bp 0based, right_bp 0based
    return calls

def parse_circexplorer2_known(path):
    calls = []
    for ln in open(path):
        if not ln.strip() or ln.startswith("#"): continue
        f = ln.rstrip("\n").split("\t")
        chrom, start, end = f[0], int(f[1]), int(f[2])
        reads = int(f[12]) if len(f) > 12 and f[12].isdigit() else 0
        calls.append((chrom, start, end - 1, reads, _strand(f, 5)))
    return calls

def parse_ciri2(path):
    calls = []
    for ln in open(path):
        if ln.startswith("circRNA_ID") or not ln.strip(): continue
        f = ln.rstrip("\n").split("\t")
        # circRNA_ID, chr, circRNA_start(1b), circRNA_end(1b), #junction_reads, ..., strand (col 11)
        chrom = f[1]; start1 = int(f[2]); end1 = int(f[3])
        reads = int(f[4]) if f[4].isdigit() else 0
        calls.append((chrom, start1 - 1, end1 - 1, reads, _strand(f, 10)))
    return calls

PARSERS = {
    "circexplorer2_bed": parse_bed,
    "circexplorer2_known": parse_circexplorer2_known,
    "ciri2": parse_ciri2,
    "bed": parse_bed,
}

# ---------------------------------------------------------------- truth
def load_genome(dataset):
    path = os.path.join(dataset, "genome.fa")
    g = {}; name = None
    try:
        with open(path) as fh:
            for ln in fh:
                if ln.startswith(">"):
                    name = ln[1:].split()[0]; g[name] = []
                else:
                    g[name].append(ln.strip())
    except FileNotFoundError:
        return {}
    return {k: "".join(v) for k, v in g.items()}

def left_align(chrom, left, right, genome):
    """Canonicalize a back-splice breakpoint pair by shifting it as far left as the
    microhomology allows. (left,right) and (left-1,right-1) describe the same junction
    when genome[left-1]==genome[right]; applied identically to truth and calls, this makes
    matching independent of which end of a microhomology a tool reports. left/right are
    0-based, inclusive breakpoint positions (left<right)."""
    s = genome.get(chrom)
    if not s:
        return left, right
    while left > 0 and right < len(s) and s[left - 1] == s[right]:
        left -= 1; right -= 1
    return left, right

def load_truth(dataset):
    bsj = []
    with open(os.path.join(dataset, "truth_bsj.tsv")) as fh:
        hdr = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            left = int(f[ix["circle_start_0based"]])
            right = int(f[ix["circle_end_0based"]]) - 1
            bsj.append(dict(
                id=f[ix["bsj_id"]], category=f[ix["category"]], chrom=f[ix["chrom"]],
                left=left, right=right, strand=f[ix["strand"]],
                expected=int(f[ix["expected_junction_reads"]]),
                tol=int(f[ix["breakpoint_tolerance_bp"]])))
    negs = []
    with open(os.path.join(dataset, "truth_negatives.tsv")) as fh:
        hdr = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            negs.append(dict(id=f[ix["neg_id"]], category=f[ix["category"]], chrom=f[ix["chrom"]],
                             start=int(f[ix["start_0based"]]), end=int(f[ix["end_0based"]])))
    return bsj, negs

def library_strandedness(dataset, library):
    """fwd / rev / unstranded from dataset/libraries.tsv, or 'NA' if not listed."""
    try:
        with open(os.path.join(dataset, "libraries.tsv")) as fh:
            hdr = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(hdr)}
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if f[ix["library"]] == library: return f[ix["strandedness"]]
    except FileNotFoundError:
        pass
    return "NA"

# ---------------------------------------------------------------- matching
def match(calls, bsj, negs, base_tol):
    truth_found = {b["id"]: None for b in bsj}
    truth_strand = {}  # truth_id -> strand declared by the matching call
    call_status = []   # per call: ('TP', truth_id) | ('FP_decoy', neg_id) | ('FP_offtarget', None)
    used_truth = set()
    for (chrom, l, r, reads, strand) in calls:
        best = None; bestd = None
        for b in bsj:
            if b["chrom"] != chrom or b["id"] in used_truth: continue
            tol = max(base_tol, b["tol"])
            if abs(l - b["left"]) <= tol and abs(r - b["right"]) <= tol:
                d = abs(l - b["left"]) + abs(r - b["right"])
                if bestd is None or d < bestd: best, bestd = b, d
        if best is not None:
            used_truth.add(best["id"]); truth_found[best["id"]] = reads
            truth_strand[best["id"]] = strand
            call_status.append(("TP", best["id"], reads)); continue
        hit = None
        for n in negs:
            if n["chrom"] == chrom and (n["start"] - base_tol) <= l <= (n["end"] + base_tol) \
               and (n["start"] - base_tol) <= r <= (n["end"] + base_tol):
                hit = n["id"]; break
        call_status.append((("FP_decoy" if hit else "FP_offtarget"), hit, reads))
    return truth_found, truth_strand, call_status

# ---------------------------------------------------------------- correlation
def _rank(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    rank = [0.0] * len(xs); i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[order[j+1]] == xs[order[i]]: j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1): rank[order[k]] = avg
        i = j + 1
    return rank

def _pearson(x, y):
    n = len(x)
    if n < 2: return float("nan")
    mx, my = sum(x)/n, sum(y)/n
    sx = math.sqrt(sum((a-mx)**2 for a in x)); sy = math.sqrt(sum((b-my)**2 for b in y))
    if sx == 0 or sy == 0: return float("nan")
    return sum((a-mx)*(b-my) for a, b in zip(x, y)) / (sx*sy)

def spearman(x, y):
    if len(x) < 2: return float("nan")
    return _pearson(_rank(x), _rank(y))

# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", required=True, choices=list(PARSERS))
    ap.add_argument("--tool-name", default=None,
                    help="name to report in output/JSON (default: the --tool parser name). "
                         "Use e.g. --tool bed --tool-name bsj_finder to label your pipeline.")
    ap.add_argument("--calls", required=True)
    ap.add_argument("--dataset", default="../dataset")
    ap.add_argument("--library", default="NA")
    ap.add_argument("--base-tol", type=int, default=2,
                    help="baseline breakpoint tolerance (added to per-junction tol)")
    ap.add_argument("--out-prefix", default=None)
    a = ap.parse_args()

    tool_name = a.tool_name or a.tool
    calls = PARSERS[a.tool](a.calls) if os.path.exists(a.calls) else []
    bsj, negs = load_truth(a.dataset)

    # Microhomology-aware normalization: left-align both truth and call breakpoints to a
    # canonical position, so a call reported at the other end of a microhomology run still
    # matches (see left_align). Convention-independent and fair across tools.
    genome = load_genome(a.dataset)
    if genome:
        for b in bsj:
            b["left"], b["right"] = left_align(b["chrom"], b["left"], b["right"], genome)
        calls = [(chrom, *left_align(chrom, l, r, genome), reads, strand)
                 for (chrom, l, r, reads, strand) in calls]

    truth_found, truth_strand, call_status = match(calls, bsj, negs, a.base_tol)
    strandedness = library_strandedness(a.dataset, a.library)

    P = len(bsj); Nn = len(negs)
    TP = sum(1 for v in truth_found.values() if v is not None)
    FN = P - TP
    FP_decoy = sum(1 for s in call_status if s[0] == "FP_decoy")
    FP_off = sum(1 for s in call_status if s[0] == "FP_offtarget")
    FP = FP_decoy + FP_off
    decoys_hit = set(s[1] for s in call_status if s[0] == "FP_decoy")
    TN = Nn - len(decoys_hit)

    precision = TP/(TP+FP) if (TP+FP) else float("nan")
    recall = TP/P if P else float("nan")
    f1 = 2*precision*recall/(precision+recall) if precision==precision and recall==recall and (precision+recall) else float("nan")
    specificity = TN/Nn if Nn else float("nan")
    accuracy = (TP+TN)/(P+Nn) if (P+Nn) else float("nan")
    fdr = FP/(TP+FP) if (TP+FP) else float("nan")

    # per-category recall
    cats = {}
    for b in bsj:
        c = cats.setdefault(b["category"], [0, 0])
        c[1] += 1
        if truth_found[b["id"]] is not None: c[0] += 1
    # coverage-stratified recall
    bins = [("1-2", 1, 2), ("3-9", 3, 9), ("10-49", 10, 49), ("50+", 50, 10**9)]
    covrec = {}
    for name, lo, hi in bins:
        sub = [b for b in bsj if lo <= b["expected"] <= hi]
        det = sum(1 for b in sub if truth_found[b["id"]] is not None)
        covrec[name] = (det, len(sub))
    # quantification (matched TP only, where call reads > 0)
    qx, qy = [], []
    for b in bsj:
        v = truth_found[b["id"]]
        if v is not None and v > 0:
            qx.append(b["expected"]); qy.append(v)
    pear, spear = _pearson(qx, qy), spearman(qx, qy)
    # strand of matched calls vs truth ('.' = tool declared no strand). Split by truth strand
    # because the truth set is mostly '+': a caller that always says '+' still scores 42/48.
    found = [b for b in bsj if truth_found[b["id"]] is not None]
    s_ok = [b for b in found if truth_strand[b["id"]] == b["strand"]]
    s_bad = [b for b in found if truth_strand[b["id"]] not in (b["strand"], ".")]
    s_unk = [b for b in found if truth_strand[b["id"]] == "."]
    strand_acc = len(s_ok)/TP if TP else float("nan")
    recall_strand = len(s_ok)/P if P else float("nan")
    strand_by_truth = {s: f"{sum(1 for b in s_ok if b['strand'] == s)}/"
                          f"{sum(1 for b in found if b['strand'] == s)}" for s in ("+", "-")}

    res = dict(tool=tool_name, library=a.library, strandedness=strandedness, calls=len(calls),
               P=P, N=Nn, TP=TP, FN=FN, FP=FP, FP_on_decoy=FP_decoy, FP_offtarget=FP_off, TN=TN,
               precision=precision, recall=recall, f1=f1, specificity=specificity,
               accuracy=accuracy, fdr=fdr,
               quant_pearson=pear, quant_spearman=spear,
               strand_correct=len(s_ok), strand_wrong=len(s_bad), strand_unknown=len(s_unk),
               strand_accuracy=strand_acc, recall_strand=recall_strand,
               strand_accuracy_by_truth_strand=strand_by_truth,
               per_category_recall={k: f"{v[0]}/{v[1]}" for k, v in sorted(cats.items())},
               coverage_recall={k: f"{v[0]}/{v[1]}" for k, v in covrec.items()},
               missed=[b["id"] for b in bsj if truth_found[b["id"]] is None],
               decoys_called=sorted(decoys_hit),
               strand_wrong_ids=[b["id"] for b in s_bad],
               strand_unknown_ids=[b["id"] for b in s_unk])

    def fmt(x): return "nan" if x != x else f"{x:.3f}"
    print(f"\n=== {tool_name}  [{a.library}, strandedness={strandedness}] ===")
    print(f"calls={len(calls)}  P(true)={P}  N(decoy)={Nn}")
    print(f"TP={TP}  FN={FN}  FP={FP} (decoy={FP_decoy}, off-target={FP_off})  TN={TN}")
    print(f"precision={fmt(precision)}  recall={fmt(recall)}  F1={fmt(f1)}")
    print(f"specificity={fmt(specificity)}  accuracy={fmt(accuracy)}  FDR={fmt(fdr)}")
    print(f"quant corr: Pearson={fmt(pear)}  Spearman={fmt(spear)}")
    print(f"strand (of TP): correct={len(s_ok)}  wrong={len(s_bad)}  none={len(s_unk)}  "
          f"accuracy={fmt(strand_acc)}  recall_strand={fmt(recall_strand)}  "
          f"by truth strand: +={strand_by_truth['+']} -={strand_by_truth['-']}")
    print("per-category recall: " + "  ".join(f"{k}={v[0]}/{v[1]}" for k, v in sorted(cats.items())))
    print("coverage recall:     " + "  ".join(f"{k}={covrec[k][0]}/{covrec[k][1]}" for k, _, _ in bins))

    if a.out_prefix:
        os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
        with open(a.out_prefix + ".json", "w") as fh: json.dump(res, fh, indent=2)
    return res

if __name__ == "__main__":
    main()
