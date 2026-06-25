import gzip, os, collections
D = "dataset"

def fa(path):
    seqs, name = {}, None
    for ln in open(path):
        if ln.startswith(">"): name = ln[1:].strip().split()[0]; seqs[name] = []
        else: seqs[name].append(ln.strip())
    return {k: "".join(v) for k, v in seqs.items()}

COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]

errs, oks = [], []
def check(cond, msg):
    (oks if cond else errs).append(msg)

# 1. genome
g = fa(os.path.join(D, "genome.fa"))
check(len(g) == 3, f"3 contigs ({len(g)})")
check(len(g["chrT1"]) == 700000 and len(g["chrT2"]) == 600000 and len(g["chrP"]) == 200000,
      "contig sizes 700k/600k/200k")

# 2. truth_bsj: coords in-bounds + motif matches genome
rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(D, "truth_bsj.tsv"))]
hdr = rows[0]; ix = {c: i for i, c in enumerate(hdr)}
bsj = rows[1:]
motif_mismatch = 0; oob = 0
for r in bsj:
    c = r[ix["chrom"]]; acc = int(r[ix["acceptor_pos_0based"]]); don = int(r[ix["donor_pos_0based"]])
    strand = r[ix["strand"]]; motif = r[ix["splice_motif"]]
    if not (0 <= acc < len(g[c]) and 0 <= don < len(g[c])): oob += 1; continue
    if strand == "+":
        dn = g[c][don+1:don+3] + "/" + g[c][acc-2:acc]
    else:
        dn = rc(g[c][don-2:don]) + "/" + rc(g[c][acc+1:acc+3])
    if dn != motif: motif_mismatch += 1
check(oob == 0, f"all BSJ coords in-bounds (oob={oob})")
check(motif_mismatch == 0, f"all motifs match genome.fa (mismatch={motif_mismatch})")
canon = sum(1 for r in bsj if r[ix["splice_motif"]] == "GT/AG")
check(canon >= 30, f"canonical GT/AG circles = {canon}")

# 3. junction read counts per BSJ == expected (pe_lowerr library)
exp = {r[ix["bsj_id"]]: int(r[ix["expected_junction_reads"]]) for r in bsj}
cnt = collections.Counter()
for l in open(os.path.join(D, "read_origin.tsv")):
    p = l.rstrip("\n").split("\t")
    if p[0] == "read_id": continue
    rid, lib, origin, rclass = p[0], p[1], p[2], p[3]
    if lib == "pe_lowerr" and rclass == "circle_junction":
        cnt[origin] += 1
mism = [(k, exp[k], cnt.get(k, 0)) for k in exp if cnt.get(k, 0) != exp[k]]
check(not mism, f"junction read counts match expected (mismatches={len(mism)})")
if mism: print("  e.g.", mism[:5])

# 4. FASTQ integrity
def fq_stats(path):
    n = 0; lens = set()
    with gzip.open(path, "rt") as fh:
        for i, ln in enumerate(fh):
            if i % 4 == 1: lens.add(len(ln.strip())); n += 1
    return n, lens
for lib, rl in [("pe_lowerr", 150), ("pe_higherr", 150)]:
    n1, l1 = fq_stats(os.path.join(D, f"reads_{lib}_R1.fastq.gz"))
    n2, l2 = fq_stats(os.path.join(D, f"reads_{lib}_R2.fastq.gz"))
    check(n1 == n2, f"{lib}: R1==R2 read count ({n1}=={n2})")
    check(l1 == {rl} and l2 == {rl}, f"{lib}: all reads {rl}bp ({l1|l2})")
for lib, rl in [("se_lowerr", 100), ("se_higherr", 100)]:
    n, l = fq_stats(os.path.join(D, f"reads_{lib}.fastq.gz"))
    check(l == {rl}, f"{lib}: {n} reads all {rl}bp ({l})")

# 5. GTF omissions: category-C genes annotated with one exon dropped
gtf = [l.split("\t") for l in open(os.path.join(D, "annotation.gtf")) if not l.startswith("#")]
exon_by_gene = collections.Counter()
genes_present = set()
for f in gtf:
    if len(f) < 9: continue
    if f[2] == "exon":
        gid = f[8].split('gene_id "')[1].split('"')[0]
        exon_by_gene[gid] += 1
    if f[2] == "gene":
        gid = f[8].split('gene_id "')[1].split('"')[0]
        genes_present.add(gid)
c_ok = all(exon_by_gene[f"geneC{n}"] == 3 for n in range(1, 5))  # built with 4 exons, 1 hidden
check(c_ok, f"category-C genes show 3/4 exons (hidden exon omitted): "
            + ",".join(f'geneC{n}={exon_by_gene[f"geneC{n}"]}' for n in range(1,5)))
# intronic D circles must NOT be annotated as features
check(not any("BSJ_D" in gid for gid in genes_present), "no intronic-circle gene features in GTF")

# 6. negatives present
negs = [l for l in open(os.path.join(D, "truth_negatives.tsv"))][1:]
cats = set(l.split("\t")[1] for l in negs)
check(len([c for c in cats if c.startswith("DEC")]) >= 9, f"decoy categories present = {len(cats)}")

print("\n".join("OK   " + m for m in oks))
print("\n".join("FAIL " + m for m in errs))
print(f"\n{len(oks)} checks passed, {len(errs)} failed")
