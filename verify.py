import gzip, os, sys, collections
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

# 3. junction read counts per BSJ == expected, in EVERY library
libs = [l.rstrip("\n").split("\t") for l in open(os.path.join(D, "libraries.tsv"))]
lix = {c: i for i, c in enumerate(libs[0])}
ALL = {r[lix["library"]]: dict(layout=r[lix["layout"]], read_len=int(r[lix["read_len"]]),
                               strandedness=r[lix["strandedness"]], base=r[lix["derived_from"]],
                               scale=int(r[lix["scale"]]) if "scale" in lix else 1)
       for r in libs[1:]}
LIBS = {t: L for t, L in ALL.items() if L["scale"] == 1}        # libraries in read_origin.tsv
SCALED = {t: L for t, L in ALL.items() if L["scale"] != 1}      # dataset/scaled/, optional
PE_LIBS = [(t, L["read_len"]) for t, L in LIBS.items() if L["layout"] == "pe"]
SE_LIBS = [(t, L["read_len"]) for t, L in LIBS.items() if L["layout"] == "se"]
ALL_LIBS = list(LIBS)

exp = {r[ix["bsj_id"]]: int(r[ix["expected_junction_reads"]]) for r in bsj}
cnt = collections.defaultdict(collections.Counter)
lib_totals = collections.Counter()
for l in open(os.path.join(D, "read_origin.tsv")):
    p = l.rstrip("\n").split("\t")
    if p[0] == "read_id": continue
    rid, lib, origin, rclass = p[0], p[1], p[2], p[3]
    lib_totals[lib] += 1
    if rclass == "circle_junction":
        # isoforms sharing one junction (BSJ_J2_1/_2) are scored as one truth row (BSJ_J2)
        cnt[lib][origin if origin in exp else origin.rsplit("_", 1)[0]] += 1
bad = []
for lib in ALL_LIBS:
    bad += [(lib, k, exp[k], cnt[lib].get(k, 0)) for k in exp if cnt[lib].get(k, 0) != exp[k]]
check(not bad, f"junction read counts match expected in all {len(ALL_LIBS)} libraries (mismatches={len(bad)})")
if bad: print("  e.g.", bad[:5])

check(set(lib_totals) == set(ALL_LIBS),
      f"read_origin.tsv covers exactly the {len(ALL_LIBS)} libraries ({sorted(lib_totals)})")
check(len(set(lib_totals.values())) == 1,
      f"all libraries have equal read counts ({dict(lib_totals)})")

# 4. FASTQ integrity
def fq_stats(path):
    n = 0; lens = set()
    with gzip.open(path, "rt") as fh:
        for i, ln in enumerate(fh):
            if i % 4 == 1: lens.add(len(ln.strip())); n += 1
    return n, lens
for lib, rl in PE_LIBS:
    n1, l1 = fq_stats(os.path.join(D, f"reads_{lib}_R1.fastq.gz"))
    n2, l2 = fq_stats(os.path.join(D, f"reads_{lib}_R2.fastq.gz"))
    check(n1 == n2, f"{lib}: R1==R2 read count ({n1}=={n2})")
    check(l1 == {rl} and l2 == {rl}, f"{lib}: all reads {rl}bp ({l1|l2})")
for lib, rl in SE_LIBS:
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

# 6. strandedness: variants are the base library's fragments with only orientation changed
def fq_seqs(path):
    with gzip.open(path, "rt") as fh:
        return [ln.strip() for i, ln in enumerate(fh) if i % 4 == 1]
def lib_reads(lib):
    if LIBS[lib]["layout"] == "pe":
        return list(zip(fq_seqs(os.path.join(D, f"reads_{lib}_R1.fastq.gz")),
                        fq_seqs(os.path.join(D, f"reads_{lib}_R2.fastq.gz"))))
    return [(s, None) for s in fq_seqs(os.path.join(D, f"reads_{lib}.fastq.gz"))]
flip = lambda r: (r[1], r[0]) if r[1] is not None else (rc(r[0]), None)
for lib, L in LIBS.items():
    if L["base"] == ".": continue
    var, base = lib_reads(lib), lib_reads(L["base"])
    # (some DEC6 low-complexity SE reads are their own reverse complement, so test each
    #  read for either orientation rather than summing the two counts)
    n_rev = sum(v == flip(b) for v, b in zip(var, base))
    ok = len(var) == len(base) and all(v in (b, flip(b)) for v, b in zip(var, base))
    if L["strandedness"] == "rev":
        check(ok and n_rev == len(var), f"{lib}: every read is the strand-flipped {L['base']} read")
    else:
        check(ok and 0.4 < n_rev / len(var) < 0.6,
              f"{lib}: each fragment is {L['base']} as-is or flipped (flipped={n_rev}/{len(var)})")

# ...and the label matches the reads: in a fwd library R1 of a junction read is circle-sense,
# in a rev library it is antisense (checked on the first 25-mer against the truth circle).
circ = {k.split("|")[0]: v for k, v in fa(os.path.join(D, "truth_bsj_sequences.fa")).items()
        if k.endswith("|full_circle")}
origin_by_lib = collections.defaultdict(list)   # read_origin rows, in FASTQ order
for l in open(os.path.join(D, "read_origin.tsv")):
    p = l.rstrip("\n").split("\t")
    if p[0] != "read_id": origin_by_lib[p[1]].append(p)
for lib, L in LIBS.items():
    if L["strandedness"] == "unstranded": continue
    sense = anti = 0
    for r, (r1, _) in zip(origin_by_lib[lib], lib_reads(lib)):
        if r[3] != "circle_junction" or r[2] not in circ: continue
        c = circ[r[2]] * 3; k = r1[:25]
        sense += k in c; anti += rc(k) in c
    want = "sense" if L["strandedness"] == "fwd" else "antisense"
    got = sense if want == "sense" else anti
    check(got > 0.9 * (sense + anti),
          f"{lib}: labelled {L['strandedness']}, R1 junction reads circle-{want} ({got}/{sense + anti})")

# 6b. scaled runtime libraries: read count = scale x base, read length. Only with --scaled
#     (reading 36M-read files takes a few minutes).
for lib, L in (SCALED.items() if "--scaled" in sys.argv else []):
    files = ([os.path.join(D, "scaled", f"reads_{lib}_R{i}.fastq.gz") for i in (1, 2)]
             if L["layout"] == "pe" else [os.path.join(D, "scaled", f"reads_{lib}.fastq.gz")])
    if not all(os.path.exists(f) for f in files):
        print(f"SKIP {lib}: not generated (python3 generate_dataset.py --scaled)"); continue
    want = L["scale"] * lib_totals[L["base"]]
    stats = [fq_stats(f) for f in files]
    check(all(n == want for n, _ in stats) and all(l == {L["read_len"]} for _, l in stats),
          f"{lib}: {stats[0][0]} reads ({L['scale']}x {L['base']}), all {L['read_len']}bp")

# 6c. filler-padded genomes (only with --padded-genomes): real contigs unchanged, size as named
if "--padded-genomes" in sys.argv:
    for f in sorted(os.listdir(os.path.join(D, "genomes"))) if os.path.isdir(os.path.join(D, "genomes")) else []:
        if not f.endswith(".fa"): continue
        pg = fa(os.path.join(D, "genomes", f))
        same = all(pg.get(c) == g[c] for c in g)
        size = sum(len(v) for v in pg.values())
        check(same, f"genomes/{f}: real contigs identical to genome.fa")
        check(f != "genome_1gbp.fa" or size == 1_000_000_000, f"genomes/{f}: {size:,} bp in {len(pg)} contigs")

# 7. negatives present
negs = [l for l in open(os.path.join(D, "truth_negatives.tsv"))][1:]
cats = set(l.split("\t")[1] for l in negs)
check(len([c for c in cats if c.startswith("DEC")]) >= 9, f"decoy categories present = {len(cats)}")

print("\n".join("OK   " + m for m in oks))
print("\n".join("FAIL " + m for m in errs))
print(f"\n{len(oks)} checks passed, {len(errs)} failed")
