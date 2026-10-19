#!/usr/bin/env python3
"""
Convert bsj_finder's per-read BSJ table into an aggregated BED6 for the benchmark.

bsj_finder (workflow/scripts/detect_bsj_cor.py) emits one row PER READ that spans a
back-splice junction, with genomic breakpoint coordinates:

    q_id  bsj_coord  chromosome  s_start  s_end  s_strand  overlap_trimmed_nt

`s_start` / `s_end` are the genomic min / max breakpoints of the circle (already
ordered low->high regardless of strand), on the same 0-based convention as the test
dataset's truth_bsj.tsv (`circle_start_0based` / `circle_end_0based`) and
truth_bsj.bed ([min(bp), max(bp)+1)). They therefore map directly to a BED6:

    chrom   s_start   s_end   name   junction_reads   strand

which evaluate.py --tool bed reads as (chrom, left=start, right=end-1, reads).

Reads supporting the same junction rarely share byte-identical breakpoints (fuzzy
breakpoints, microhomology, sequencing error), so calls whose BOTH breakpoints fall
within --merge-tol bp are clustered into a single junction, the representative
coordinate is the most frequent exact (s_start,s_end) in the cluster, and
junction_reads is the number of supporting reads (feeds the quantification metric).

Clustering ignores strand: in a reverse-stranded or unstranded library (and for R2 of a
pair) reads of one junction align to both strands, and splitting on strand would turn
one junction into two calls. The reported strand is the majority s_strand of the
cluster's reads ('.' on a tie), so it is the transcript strand only if bsj_finder's
s_strand is (i.e. after bsj_finder corrects for library strandedness); evaluate.py
scores it against the truth strand.

Dependency-free (stdlib only) so it runs anywhere the benchmark does.

Usage:
    bsj_finder_to_bed.py --in <sample>.bsj_cor.tsv.gz --out calls.bed [--merge-tol 5]
    # or read from stdin:
    zcat <sample>.bsj_cor.tsv.gz | bsj_finder_to_bed.py --in - --out calls.bed
"""
import argparse, gzip, io, sys
from collections import Counter


def _open(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, encoding="utf-8")


def read_rows(path):
    """Yield (chrom, start, end, strand) per read; tolerant of header/blank/short lines."""
    fh = _open(path)
    try:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            f = ln.split("\t")
            # skip header line
            if f[0] == "q_id":
                continue
            if len(f) < 6:
                continue
            chrom = f[2]
            try:
                s_start = int(f[3]); s_end = int(f[4])
            except ValueError:
                continue
            strand = f[5] if f[5] in ("+", "-") else "."
            if s_start > s_end:
                s_start, s_end = s_end, s_start
            if s_start < 0:
                s_start = 0
            yield (chrom, s_start, s_end, strand)
    finally:
        if fh is not sys.stdin:
            fh.close()


def cluster(rows, merge_tol):
    """Greedy single-linkage-ish clustering per chrom (strand-blind); merge when BOTH
    breakpoints are within merge_tol of a cluster's current representative."""
    by_group = {}
    for (chrom, s, e, strand) in rows:
        by_group.setdefault(chrom, []).append((s, e, strand))

    calls = []  # (chrom, start, end, reads, strand)
    for chrom, pts in by_group.items():
        pts.sort()
        clusters = []  # each: {"rep": (s,e), "members": Counter(), "strands": Counter()}
        for (s, e, strand) in pts:
            placed = False
            for c in clusters:
                rs, re = c["rep"]
                if abs(s - rs) <= merge_tol and abs(e - re) <= merge_tol:
                    c["members"][(s, e)] += 1
                    c["strands"][strand] += 1
                    # keep representative as the running mode (ties -> smallest)
                    c["rep"] = c["members"].most_common(1)[0][0]
                    placed = True
                    break
            if not placed:
                clusters.append({"rep": (s, e), "members": Counter({(s, e): 1}),
                                 "strands": Counter({strand: 1})})
        for c in clusters:
            reads = sum(c["members"].values())
            rs, re = c["members"].most_common(1)[0][0]
            n_plus, n_minus = c["strands"]["+"], c["strands"]["-"]
            strand = "+" if n_plus > n_minus else "-" if n_minus > n_plus else "."
            calls.append((chrom, rs, re, reads, strand))
    # stable, readable ordering
    calls.sort(key=lambda x: (x[0], x[1], x[2]))
    return calls


def main():
    ap = argparse.ArgumentParser(description="bsj_finder bsj_cor.tsv(.gz) -> aggregated BED6")
    ap.add_argument("--in", dest="inp", required=True, help="bsj_cor.tsv[.gz], or '-' for stdin")
    ap.add_argument("--out", default="-", help="output BED6 path, or '-' for stdout")
    ap.add_argument("--merge-tol", type=int, default=5,
                    help="merge calls whose both breakpoints are within this many bp (default 5)")
    ap.add_argument("--name-prefix", default="BSJF", help="BED feature name prefix")
    a = ap.parse_args()

    calls = cluster(read_rows(a.inp), a.merge_tol)

    out = sys.stdout if a.out == "-" else open(a.out, "w", encoding="utf-8")
    try:
        for i, (chrom, s, e, reads, strand) in enumerate(calls, 1):
            out.write(f"{chrom}\t{s}\t{e}\t{a.name_prefix}_{i}\t{reads}\t{strand}\n")
    finally:
        if out is not sys.stdout:
            out.close()
    sys.stderr.write(f"[bsj_finder_to_bed] {len(calls)} junctions written "
                     f"(merge_tol={a.merge_tol}) -> {a.out}\n")


if __name__ == "__main__":
    main()
