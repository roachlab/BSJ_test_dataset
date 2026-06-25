#!/usr/bin/env python3
"""
Tiny BSJ (back-splice junction / circRNA) test-dataset generator.

Produces a small, fully-labelled dataset for validating circRNA detection pipelines:
  - genome.fa            ~1.5 Mbp, 3 contigs
  - annotation.gtf       gene models (with deliberate omissions)
  - reads_*.fastq.gz     PE 2x150 + SE 100, low- and high-error
  - truth_bsj.tsv/.bed   true back-splice junctions
  - truth_bsj_sequences.fa
  - truth_negatives.tsv  decoy loci (expected: no circle call)
  - read_origin.tsv      per-read provenance for scoring
  - README.md

Everything is driven by a single seeded RNG -> fully reproducible.
"""

import os, gzip, random, textwrap

SEED = 20260624
rng = random.Random(SEED)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dataset")
os.makedirs(OUT, exist_ok=True)

BASES = "ACGT"
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")
def revcomp(s): return s.translate(COMP)[::-1]
def randseq(n): return "".join(rng.choice(BASES) for _ in range(n))

# --------------------------------------------------------------------------
# Genome scaffolds
# --------------------------------------------------------------------------
CONTIG_SIZES = {"chrT1": 700_000, "chrT2": 600_000, "chrP": 200_000}
contigs = {name: list(randseq(size)) for name, size in CONTIG_SIZES.items()}
cursors = {name: 12_000 for name in contigs}        # leave a clean flank at the start

def write_region(contig, start, seq):
    arr = contigs[contig]
    for i, b in enumerate(seq):
        arr[start + i] = b

def alloc(contig, length, gap=4_000):
    s = cursors[contig]
    cursors[contig] = s + length + gap
    return s

# --------------------------------------------------------------------------
# Gene model builder
# --------------------------------------------------------------------------
class Gene:
    def __init__(self, gid, contig, strand, exons, exon_seqs, sense_len):
        self.gid = gid
        self.contig = contig
        self.strand = strand
        self.exons = exons          # list of (gstart, gend) 0-based half-open, genome coords
        self.exon_seqs = exon_seqs  # sense (transcript 5'->3') exon sequences
        self.sense_len = sense_len
    def mature(self):               # spliced linear transcript, sense
        return "".join(self.exon_seqs)

def make_gene(gid, contig, strand, exon_lens, intron_lens,
              noncanon=False, microhom=0, retain_intron_after=None,
              dup_exon=False):
    """Build a gene in transcript (sense) orientation, then place into genome."""
    assert len(intron_lens) == len(exon_lens) - 1
    exon_seqs = [randseq(l) for l in exon_lens]

    if dup_exon:  # make exon[2] a near-identical copy of exon[1] (one SNP)
        if len(exon_seqs) > 2:
            s = list(exon_seqs[1]); s[len(s)//2] = rng.choice(BASES)
            exon_seqs[2] = "".join(s)[:exon_lens[2]] if len(s) >= exon_lens[2] else "".join(s)

    if microhom:  # share `microhom` bases between last exon's 3' end and first exon's 5' start
        mh = exon_seqs[-1][-microhom:]
        exon_seqs[0] = mh + exon_seqs[0][microhom:]

    don, acc = ("AT", "AC") if noncanon else ("GT", "AG")
    introns = []
    for k, il in enumerate(intron_lens):
        if retain_intron_after is not None and k == retain_intron_after:
            introns.append(None)            # marker: this intron is retained in the circle
            introns[-1] = (don + randseq(il - 4) + acc)  # still real sequence in genome
        else:
            introns.append(don + randseq(il - 4) + acc)

    # flanking introns so that terminal exons also have canonical (or matching) motifs
    flank5 = don + randseq(196) + acc        # upstream intron, ends in `acc`
    flank3 = don + randseq(196) + acc        # downstream intron, starts with `don`
    # assemble sense sequence and record exon offsets within it
    sense = [flank5]
    offsets = []
    pos = len(flank5)
    for k, es in enumerate(exon_seqs):
        offsets.append((pos, pos + len(es)))
        sense.append(es); pos += len(es)
        if k < len(introns):
            sense.append(introns[k]); pos += len(introns[k])
    sense.append(flank3)
    sense = "".join(sense)
    sense_len = len(sense)

    start = alloc(contig, sense_len)
    if strand == "+":
        write_region(contig, start, sense)
        exons = [(start + o0, start + o1) for (o0, o1) in offsets]
    else:
        write_region(contig, start, revcomp(sense))
        exons = [(start + (sense_len - o1), start + (sense_len - o0)) for (o0, o1) in offsets]
        exons = exons  # genomic intervals (already 5'->3' in transcript order maps high->low genome)

    return Gene(gid, contig, strand, exons, exon_seqs, sense_len)

# --------------------------------------------------------------------------
# Circle helper
# --------------------------------------------------------------------------
def circle_from_gene(gene, i, j, extra_seq=""):
    """Circle over exons i..j (inclusive). Returns mature sense seq + breakpoints."""
    mature = "".join(gene.exon_seqs[i:j+1]) + extra_seq
    if gene.strand == "+":
        acceptor = gene.exons[i][0]            # 5' start of circle (genomic)
        donor = gene.exons[j][1] - 1           # 3' end of circle (genomic)
    else:
        acceptor = gene.exons[i][1] - 1
        donor = gene.exons[j][0]
    gstart = min(gene.exons[k][0] for k in range(i, j+1))
    gend = max(gene.exons[k][1] for k in range(i, j+1))
    return mature, acceptor, donor, gstart, gend

# --------------------------------------------------------------------------
# Read jobs: collected first, then rendered into each library
# --------------------------------------------------------------------------
read_jobs = []   # dicts describing how to emit reads
def add_circle_job(origin_id, seq, n_junction, n_body, note):
    read_jobs.append(dict(kind="circle", origin=origin_id, seq=seq,
                          n_junction=n_junction, n_body=n_body,
                          spans=True, role="positive", note=note))
def add_linear_job(origin_id, seq, n, role="negative", note="", spans=False):
    read_jobs.append(dict(kind="linear", origin=origin_id, seq=seq,
                          n=n, spans=spans, role=role, note=note))
def add_chimera_job(origin_id, seq, n, note):
    # pre-built chimeric reads (decoys): seq is a pool string longer than read; we slice windows
    read_jobs.append(dict(kind="chimera", origin=origin_id, seq=seq,
                          n=n, spans=True, role="negative", note=note))

# --------------------------------------------------------------------------
# Registries for truth output
# --------------------------------------------------------------------------
annotated_genes = []   # genes that go into the GTF
bsj_truth = []         # true back-splice junctions
neg_truth = []         # decoy loci

def reg_bsj(bid, category, gene, acceptor, donor, gstart, gend, motif, ctype,
            cov, tol=0, note=""):
    bsj_truth.append(dict(id=bid, category=category, chrom=gene.contig,
                          acceptor=acceptor, donor=donor, gstart=gstart, gend=gend,
                          strand=gene.strand, motif=motif, ctype=ctype,
                          expected_reads=cov, tol=tol, note=note))

def reg_neg(nid, category, contig, gstart, gend, strand, trap):
    neg_truth.append(dict(id=nid, category=category, chrom=contig,
                          gstart=gstart, gend=gend, strand=strand, trap=trap))

# ==========================================================================
#  POSITIVE TEST CASES
# ==========================================================================

# ---- A: canonical multi-exon circles (5) ----
A_cov = [200, 100, 40, 12, 3]
for n, cov in enumerate(A_cov, 1):
    g = make_gene(f"geneA{n}", "chrT1", "+", [220, 180, 160, 200], [600, 700, 650])
    annotated_genes.append(g)
    seq, acc, don, gs, ge = circle_from_gene(g, 1, 2)   # circle over exons 2..3
    reg_bsj(f"BSJ_A{n}", "A_canonical_multiexon", g, acc, don, gs, ge, "GT/AG",
            "exonic", cov, note="canonical multi-exon circle")
    add_circle_job(f"BSJ_A{n}", seq, cov, max(2, cov//3), "A canonical")
    # host linear transcript reads (forward splicing) for realism
    add_linear_job(f"geneA{n}_linear", g.mature(), max(4, cov//4),
                   role="host_linear", note="linear host transcript", spans=True)

# ---- B: single-exon self-circles (4) ----
B_cov = [80, 30, 8, 2]
for n, cov in enumerate(B_cov, 1):
    g = make_gene(f"geneB{n}", "chrT1", "+", [400, 200], [600])
    annotated_genes.append(g)
    seq, acc, don, gs, ge = circle_from_gene(g, 0, 0)   # single exon circle
    reg_bsj(f"BSJ_B{n}", "B_single_exon", g, acc, don, gs, ge, "GT/AG",
            "single_exon", cov, note="single-exon self-circularization")
    add_circle_job(f"BSJ_B{n}", seq, cov, max(2, cov//3), "B single-exon")

# ---- C: BSJ on an UN-ANNOTATED exon (4) ----
C_cov = [60, 20, 6, 2]
for n, cov in enumerate(C_cov, 1):
    g = make_gene(f"geneC{n}", "chrT1", "+", [200, 150, 180, 170], [600, 650, 600])
    # annotate the gene but DROP exon index 2 (E3) from the GTF
    g.hidden_exons = {2}
    annotated_genes.append(g)
    seq, acc, don, gs, ge = circle_from_gene(g, 2, 2)   # circle uses the hidden exon
    reg_bsj(f"BSJ_C{n}", "C_unannotated_exon", g, acc, don, gs, ge, "GT/AG",
            "exonic_unannotated", cov, note="circle exon absent from GTF")
    add_circle_job(f"BSJ_C{n}", seq, cov, max(2, cov//3), "C unannotated exon")

# ---- D: intronic ciRNA (4) ----
D_cov = [25, 8, 3, 1]
for n, cov in enumerate(D_cov, 1):
    g = make_gene(f"geneD{n}", "chrT1", "+", [180, 160], [1200])
    annotated_genes.append(g)
    # circle lies entirely inside the intron (between the two exons)
    i_start = g.exons[0][1] + 30
    i_len = 220
    circ = "".join(contigs[g.contig][i_start:i_start+i_len])
    acc, don = i_start, i_start + i_len - 1
    reg_bsj(f"BSJ_D{n}", "D_intronic_ciRNA", g, acc, don, i_start, i_start+i_len,
            "GT/AG", "intronic", cov, note="circular intronic RNA, no exon overlap")
    add_circle_job(f"BSJ_D{n}", circ, cov, max(1, cov//4), "D intronic")

# ---- E: EIciRNA (retained intron in circle body) (3) ----
E_cov = [40, 10, 3]
for n, cov in enumerate(E_cov, 1):
    g = make_gene(f"geneE{n}", "chrT1", "+", [200, 180, 160], [300, 600],
                  retain_intron_after=0)
    annotated_genes.append(g)
    # circle = exon1 + retained intron1 + exon2  (genomic contiguous span E1..E2)
    gs = g.exons[0][0]; ge = g.exons[1][1]
    circ = "".join(contigs[g.contig][gs:ge])
    acc = gs; don = ge - 1
    reg_bsj(f"BSJ_E{n}", "E_EIciRNA", g, acc, don, gs, ge, "GT/AG",
            "exon_intron", cov, note="exon-intron circle (retained intron)")
    add_circle_job(f"BSJ_E{n}", circ, cov, max(2, cov//3), "E EIciRNA")

# ---- F: minus-strand circles (4) ----
F_cov = [120, 40, 10, 2]
for n, cov in enumerate(F_cov, 1):
    g = make_gene(f"geneF{n}", "chrT2", "-", [220, 180, 200], [600, 650])
    annotated_genes.append(g)
    seq, acc, don, gs, ge = circle_from_gene(g, 0, 1)
    reg_bsj(f"BSJ_F{n}", "F_minus_strand", g, acc, don, gs, ge, "GT/AG",
            "exonic", cov, note="minus-strand circle")
    add_circle_job(f"BSJ_F{n}", seq, cov, max(2, cov//3), "F minus-strand")

# ---- G: non-canonical motif (3) ----
G_cov = [30, 8, 2]
for n, cov in enumerate(G_cov, 1):
    g = make_gene(f"geneG{n}", "chrT2", "+", [200, 180, 190], [600, 650], noncanon=True)
    annotated_genes.append(g)
    seq, acc, don, gs, ge = circle_from_gene(g, 0, 1)
    reg_bsj(f"BSJ_G{n}", "G_noncanonical", g, acc, don, gs, ge, "AT/AC",
            "exonic", cov, note="non-canonical back-splice motif")
    add_circle_job(f"BSJ_G{n}", seq, cov, max(2, cov//3), "G non-canonical")

# ---- H: fuzzy breakpoint / microhomology (4) ----
H_cov = [60, 20, 6, 2]
H_mh = [6, 5, 4, 3]
for n, (cov, mh) in enumerate(zip(H_cov, H_mh), 1):
    g = make_gene(f"geneH{n}", "chrT2", "+", [210, 190, 180], [600, 650], microhom=mh)
    annotated_genes.append(g)
    seq, acc, don, gs, ge = circle_from_gene(g, 0, 1)
    reg_bsj(f"BSJ_H{n}", "H_fuzzy_breakpoint", g, acc, don, gs, ge, "GT/AG",
            "exonic", cov, tol=mh, note=f"{mh}bp microhomology -> breakpoint ambiguous +/-{mh}")
    add_circle_job(f"BSJ_H{n}", seq, cov, max(2, cov//3), "H fuzzy")
    # DEC-10 companion: SNP near the junction of this circle
    if n <= 3:
        snp_seq = list(seq); snp_seq[3] = rng.choice([b for b in BASES if b != seq[3]])
        add_linear_job(f"DEC10_snp_H{n}", "".join(snp_seq), 6,
                       role="negative", note="SNP 3bp from breakpoint", spans=True)
        reg_neg(f"DEC10_{n}", "DEC10_snp_near_junction", g.contig, gs, ge, g.strand,
                "variant near breakpoint must not create a spurious isoform/shift")

# ---- I: sub-read-length / rolling-circle (4) ----
I_cov = [50, 15, 5, 2]
I_len = [90, 70, 120, 60]   # all shorter than/near a 150 or 100bp read
for n, (cov, L) in enumerate(zip(I_cov, I_len), 1):
    g = make_gene(f"geneI{n}", "chrT2", "+", [L, 150], [600])
    annotated_genes.append(g)
    seq, acc, don, gs, ge = circle_from_gene(g, 0, 0)   # tiny single-exon circle
    reg_bsj(f"BSJ_I{n}", "I_subreadlen_rolling", g, acc, don, gs, ge, "GT/AG",
            "single_exon_short", cov, note=f"circle {L}bp < read len; reads cross BSJ multiple times")
    add_circle_job(f"BSJ_I{n}", seq, cov, 0, "I rolling-circle")

# ---- J: alternative back-splicing (5 circles across 2 genes) ----
# gene J1: two circles share the acceptor (exon0 start) but differ in donor (exon1 vs exon2)
gj1 = make_gene("geneJ1", "chrT2", "+", [200, 170, 160, 180], [600, 650, 600])
annotated_genes.append(gj1)
for n, (j_end, cov) in enumerate([(1, 80), (2, 25)], 1):
    seq, acc, don, gs, ge = circle_from_gene(gj1, 0, j_end)
    reg_bsj(f"BSJ_J1_{n}", "J_alt_backsplice", gj1, acc, don, gs, ge, "GT/AG",
            "exonic_altBS", cov, note="alternative 3' back-splice site (shared acceptor)")
    add_circle_job(f"BSJ_J1_{n}", seq, cov, max(2, cov//3), "J alt-backsplice shared acceptor")
# gene J2: same BSJ, internal alternative splicing (full vs exon-skipped body)
gj2 = make_gene("geneJ2", "chrT2", "+", [200, 150, 150, 180], [600, 650, 600])
annotated_genes.append(gj2)
seqA, acc, don, gs, ge = circle_from_gene(gj2, 0, 3)               # full body
seqB = gj2.exon_seqs[0] + gj2.exon_seqs[3]                          # exon-skipped body, same BSJ
for n, (sq, cov) in enumerate([(seqA, 40), (seqB, 10)], 1):
    reg_bsj(f"BSJ_J2_{n}", "J_alt_backsplice", gj2, acc, don, gs, ge, "GT/AG",
            "exonic_altIntern", cov, note="same BSJ, internal alternative splicing")
    add_circle_job(f"BSJ_J2_{n}", sq, cov, max(2, cov//3), "J alt internal splicing")

# ---- K: nested / overlapping circles (3 from one gene) ----
gk = make_gene("geneK1", "chrT2", "+", [200, 170, 160, 150, 180], [600, 650, 600, 650])
annotated_genes.append(gk)
for n, (i, j, cov) in enumerate([(0, 4, 50), (1, 3, 18), (2, 2, 5)], 1):
    seq, acc, don, gs, ge = circle_from_gene(gk, i, j)
    reg_bsj(f"BSJ_K{n}", "K_nested", gk, acc, don, gs, ge, "GT/AG",
            "exonic_nested", cov, note="nested/overlapping circle from shared gene")
    add_circle_job(f"BSJ_K{n}", seq, cov, max(2, cov//3), "K nested")

# ---- L: indel / non-templated bases at the junction (3) ----
L_cov = [30, 9, 3]
for n, cov in enumerate(L_cov, 1):
    g = make_gene(f"geneL{n}", "chrT2", "+", [200, 180, 190], [600, 650])
    annotated_genes.append(g)
    seq, acc, don, gs, ge = circle_from_gene(g, 0, 1)
    ins = randseq(rng.choice([2, 3, 4]))            # non-templated bases at the BSJ
    seq_ins = seq + ins                              # inserted right at the wrap point
    reg_bsj(f"BSJ_L{n}", "L_indel_junction", g, acc, don, gs, ge, "GT/AG",
            "exonic_indel", cov, tol=len(ins),
            note=f"{len(ins)} non-templated bases inserted at junction")
    add_circle_job(f"BSJ_L{n}", seq_ins, cov, max(2, cov//3), "L indel junction")

# ---- M: intergenic & antisense circles (4) ----
# two fully intergenic (no gene model at all)
M_cov = [25, 10]
for n, cov in enumerate(M_cov, 1):
    contig = "chrT1"
    s = alloc(contig, 800)
    e1s, e1e = s + 50, s + 280
    # build two pseudo-exons separated by a canonical intron, but DO NOT annotate
    write_region(contig, e1e, "GT" + "".join(contigs[contig][e1e+2:e1e+498]) + "AG")
    e2s, e2e = e1e + 500, e1e + 740
    circ = "".join(contigs[contig][e1s:e1e]) + "".join(contigs[contig][e2s:e2e])
    acc, don = e1s, e2e - 1
    class _Tmp:  # lightweight gene-like holder
        strand = "+"; contig = "chrT1"
    reg_bsj(f"BSJ_M_ig{n}", "M_intergenic", _Tmp, acc, don, e1s, e2e, "GT/AG",
            "intergenic", cov, note="circle from fully un-annotated intergenic locus")
    bsj_truth[-1]["chrom"] = contig
    add_circle_job(f"BSJ_M_ig{n}", circ, cov, max(2, cov//3), "M intergenic")
# two antisense circles overlapping an annotated sense gene
M2_cov = [4, 2]
for n, cov in enumerate(M2_cov, 1):
    host = make_gene(f"geneMsense{n}", "chrT1", "+", [300, 250, 280], [600, 650])
    annotated_genes.append(host)
    # antisense circle: take the antisense strand of host's exon span
    gs = host.exons[0][0]; ge = host.exons[2][1]
    circ = revcomp("".join(contigs[host.contig][gs:ge]))
    acc, don = ge - 1, gs   # on the - strand
    class _TmpM:
        strand = "-"; contig = "chrT1"
    _TmpM.contig = host.contig
    reg_bsj(f"BSJ_M_as{n}", "M_antisense", _TmpM, acc, don, gs, ge, "GT/AG",
            "antisense", cov, note="antisense circle overlapping a sense gene")
    bsj_truth[-1]["chrom"] = host.contig
    add_circle_job(f"BSJ_M_as{n}", circ, cov, max(1, cov//3), "M antisense")

# ==========================================================================
#  DECOYS / NEGATIVES
# ==========================================================================

def junction_chimera(a, b):
    """A read pool that, aligned to the genome, looks back-spliced: b then a (out of order)."""
    return b + a

# ---- DEC-1: tandem duplication (4) ----
for n in range(1, 5):
    contig = "chrT2"
    s = alloc(contig, 1200)
    seg = "".join(contigs[contig][s:s+500])
    write_region(contig, s + 500, seg)          # exact tandem copy
    # rotation read across the repeat boundary -> apparent out-of-order alignment
    pool = seg[250:] + seg[:250]
    add_chimera_job(f"DEC1_tandem{n}", pool, 25, "tandem duplication boundary mimics BSJ")
    reg_neg(f"DEC1_{n}", "DEC1_tandem_duplication", contig, s, s+1000, "+",
            "linear tandem duplication; rotated reads look back-spliced (GT/AG flanked)")

# ---- DEC-2: linear-only gene (5) ----
for n in range(1, 6):
    g = make_gene(f"geneDEC2_{n}", "chrT2", "+", [220, 180, 200], [600, 650])
    annotated_genes.append(g)
    add_linear_job(f"DEC2_linear{n}", g.mature(), 30, role="negative",
                   note="forward splicing only; no circle", spans=True)
    reg_neg(f"DEC2_{n}", "DEC2_linear_only", g.contig, g.exons[0][0], g.exons[-1][1], "+",
            "normal forward splice junctions must not be called as BSJ")

# ---- DEC-3: paralog pairs (4) ----  copy on chrP and chrT2 (cross-contig multimapping)
for n in range(1, 5):
    g = make_gene(f"geneDEC3a_{n}", "chrT2", "+", [240, 200], [600])
    annotated_genes.append(g)
    # near-identical paralog on chrP
    ps = alloc("chrP", g.sense_len)
    para = list("".join(contigs[g.contig][g.exons[0][0]:g.exons[-1][1]]))
    para[len(para)//2] = rng.choice(BASES)
    write_region("chrP", ps, "".join(para))
    gp = Gene(f"geneDEC3b_{n}", "chrP", "+",
              [(ps, ps+240), (ps+ (g.exons[1][0]-g.exons[0][0]), ps+(g.exons[1][1]-g.exons[0][0]))],
              g.exon_seqs, g.sense_len)
    annotated_genes.append(gp)
    add_linear_job(f"DEC3_paralog{n}", g.mature(), 25, role="negative",
                   note="multi-mapping between paralogous copies", spans=True)
    reg_neg(f"DEC3_{n}", "DEC3_paralog", g.contig, g.exons[0][0], g.exons[-1][1], "+",
            "multi-mapping reads between paralogs can fake chimeric BSJ")

# ---- DEC-4: inverted repeat / palindrome (3) ----
for n in range(1, 4):
    contig = "chrT2"; s = alloc(contig, 700)
    seg = "".join(contigs[contig][s:s+250])
    write_region(contig, s+300, revcomp(seg))
    pool = seg + revcomp(seg)
    add_chimera_job(f"DEC4_invrep{n}", pool, 20, "inverted repeat: read aligns to both strands")
    reg_neg(f"DEC4_{n}", "DEC4_inverted_repeat", contig, s, s+550, "+",
            "inverted repeat read mimics chimeric/back-splice junction")

# ---- DEC-5: readthrough / trans-splice (3) ----
for n in range(1, 4):
    g1 = make_gene(f"geneDEC5a_{n}", "chrT2", "+", [200, 180], [600])
    g2 = make_gene(f"geneDEC5b_{n}", "chrT2", "+", [200, 180], [600])
    annotated_genes += [g1, g2]
    pool = g1.mature()[-120:] + g2.mature()[:120]   # linear readthrough across two genes
    add_chimera_job(f"DEC5_readthrough{n}", pool, 18, "adjacent-gene readthrough looks chimeric")
    reg_neg(f"DEC5_{n}", "DEC5_readthrough", g1.contig, g1.exons[0][0], g2.exons[-1][1], "+",
            "linear readthrough/trans-splice vs true back-splice")

# ---- DEC-6: low-complexity (3) ----
for n in range(1, 4):
    contig = "chrT1"; s = alloc(contig, 400)
    unit = rng.choice(["A", "AT", "CAG", "GGC"])
    lc = (unit * 400)[:300]
    write_region(contig, s, lc)
    add_chimera_job(f"DEC6_lowcplx{n}", lc + lc[:80], 15, "low-complexity ambiguous soft-clips")
    reg_neg(f"DEC6_{n}", "DEC6_low_complexity", contig, s, s+300, "+",
            "homopolymer/STR produces ambiguous anchors")

# ---- DEC-7: repeat within a spliced, expressed gene (4) ----  hard negative
for n in range(1, 5):
    g = make_gene(f"geneDEC7_{n}", "chrT1", "+", [220, 300, 200], [600, 650])
    annotated_genes.append(g)
    # embed a tandem repeat inside exon 2 (middle exon)
    es, ee = g.exons[1]
    rep = "".join(contigs[g.contig][es+40:es+140])
    write_region(g.contig, es+140, rep)         # tandem copy inside the exon
    add_linear_job(f"DEC7_fwdSJ{n}", g.mature(), 25, role="host_linear",
                   note="true forward SJ reads (must be retained)", spans=True)
    add_chimera_job(f"DEC7_pseudoBSJ{n}", rep[50:] + rep[:50], 15,
                    "intra-exon repeat boundary mimics BSJ within an expressed gene")
    reg_neg(f"DEC7_{n}", "DEC7_repeat_in_gene", g.contig, g.exons[0][0], g.exons[-1][1], "+",
            "JOINT: retain forward SJs AND reject repeat-induced pseudo-BSJ")

# ---- DEC-8: duplicated exon within a gene (4) ----  hard negative
for n in range(1, 5):
    g = make_gene(f"geneDEC8_{n}", "chrT1", "+", [200, 180, 180, 200], [600, 650, 600],
                  dup_exon=True)
    annotated_genes.append(g)
    add_linear_job(f"DEC8_fwdSJ{n}", g.mature(), 25, role="host_linear",
                   note="true forward SJ reads (must be retained)", spans=True)
    # read joining downstream exon copy back to the upstream near-identical exon
    pool = g.exon_seqs[2] + g.exon_seqs[1]
    add_chimera_job(f"DEC8_pseudoBSJ{n}", pool, 15,
                    "near-identical duplicated exons fake a within-gene BSJ")
    reg_neg(f"DEC8_{n}", "DEC8_duplicated_exon", g.contig, g.exons[0][0], g.exons[-1][1], "+",
            "JOINT: retain forward SJs AND reject duplicated-exon pseudo-BSJ")

# ---- DEC-9: RT template-switching chimera (4) ----
for n in range(1, 5):
    contig = "chrT1"; s = alloc(contig, 900)
    a = "".join(contigs[contig][s:s+200])
    b = "".join(contigs[contig][s+400:s+600])
    mh = a[-5:]                                  # short microhomology drives the switch
    bb = mh + b[5:]
    pool = a + bb                                # RT jumps from a into b at the microhomology
    add_chimera_job(f"DEC9_templateswitch{n}", pool, 18,
                    "RT template switch at microhomology mimics BSJ (library artifact)")
    reg_neg(f"DEC9_{n}", "DEC9_template_switching", contig, s, s+600, "+",
            "reverse-transcriptase template switch, not a real circle")

# ---- DEC-11: high-expression linear gene (2) ----
for n in range(1, 3):
    g = make_gene(f"geneDEC11_{n}", "chrT1", "+", [400, 350, 300], [600, 650])
    annotated_genes.append(g)
    add_linear_job(f"DEC11_highexpr{n}", g.mature(), 400, role="host_linear",
                   note="very high linear coverage (dynamic-range swamp)", spans=True)
    reg_neg(f"DEC11_{n}", "DEC11_high_expression", g.contig, g.exons[0][0], g.exons[-1][1], "+",
            "high linear background must not drown low-coverage circles")

# ==========================================================================
#  READ RENDERING
# ==========================================================================
def circ_read(seq, start, length):
    L = len(seq)
    return "".join(seq[(start + i) % L] for i in range(length))

def crosses_junction(seq, start, length):
    L = len(seq)
    return ((start + length) // L) - (start // L) >= 1 if L else False

def apply_errors(seq, sub_rate, indel_rate, want_len):
    out = []
    i = 0
    while len(out) < want_len and i < len(seq):
        b = seq[i]
        r = rng.random()
        if indel_rate and r < indel_rate / 2:           # insertion
            out.append(rng.choice(BASES)); continue     # don't advance i
        elif indel_rate and r < indel_rate:             # deletion
            i += 1; continue
        if rng.random() < sub_rate:
            b = rng.choice([x for x in BASES if x != b])
        out.append(b); i += 1
    while len(out) < want_len:                            # pad if we ran short
        out.append(rng.choice(BASES))
    return "".join(out[:want_len])

LIBS = [
    ("pe_lowerr",  True,  0.002, 0.000),
    ("pe_higherr", True,  0.020, 0.003),
    ("se_lowerr",  False, 0.002, 0.000),
    ("se_higherr", False, 0.020, 0.003),
]
RLEN = {"pe": 150, "se": 100}
FRAG = 330

origin_rows = []   # (read_id, library, origin, kind, spans, role, note)

def emit_fastq():
    for libtag, is_pe, sub, indel in LIBS:
        rlen = RLEN["pe" if is_pe else "se"]
        if is_pe:
            f1 = gzip.open(os.path.join(OUT, f"reads_{libtag}_R1.fastq.gz"), "wt")
            f2 = gzip.open(os.path.join(OUT, f"reads_{libtag}_R2.fastq.gz"), "wt")
        else:
            f1 = gzip.open(os.path.join(OUT, f"reads_{libtag}.fastq.gz"), "wt")
            f2 = None
        counter = 0
        def write_pair(rid, s1, s2):
            q1 = "I" * len(s1)
            f1.write(f"@{rid}/1\n{s1}\n+\n{q1}\n" if is_pe else f"@{rid}\n{s1}\n+\n{q1}\n")
            if is_pe:
                q2 = "I" * len(s2)
                f2.write(f"@{rid}/2\n{s2}\n+\n{q2}\n")

        for job in read_jobs:
            seq = job["seq"]
            if len(seq) < 20:
                continue
            if job["kind"] == "circle":
                njunc, nbody = job["n_junction"], job["n_body"]
                short_circle = len(seq) < rlen
                for k in range(njunc):
                    counter += 1
                    rid = f"{libtag}:{job['origin']}:J{k}"
                    if short_circle:
                        start = rng.randrange(len(seq))
                    else:
                        anchor = rng.choice([5, 6, 7, 8] + list(range(15, min(rlen-15, len(seq)-15))))
                        start = (len(seq) - anchor) % len(seq)
                    raw1 = circ_read(seq, start, rlen + 8)
                    s1 = apply_errors(raw1, sub, indel, rlen)
                    s2 = ""
                    if is_pe:
                        raw2 = circ_read(seq, start + FRAG - rlen, rlen + 8)
                        s2 = revcomp(apply_errors(raw2, sub, indel, rlen))
                    write_pair(rid, s1, s2)
                    origin_rows.append((rid, libtag, job["origin"], "circle_junction", "TRUE", job["role"], job["note"]))
                for k in range(nbody):
                    if short_circle:
                        break
                    counter += 1
                    rid = f"{libtag}:{job['origin']}:B{k}"
                    start = rng.randrange(0, max(1, len(seq) - rlen))
                    s1 = apply_errors(seq[start:start+rlen+8], sub, indel, rlen)
                    s2 = ""
                    if is_pe:
                        b2 = max(0, start + FRAG - rlen)
                        s2 = revcomp(apply_errors(seq[b2:b2+rlen+8], sub, indel, rlen))
                    write_pair(rid, s1, s2)
                    origin_rows.append((rid, libtag, job["origin"], "circle_body", "FALSE", job["role"], job["note"]))
            elif job["kind"] == "linear":
                n = job["n"]
                for k in range(n):
                    counter += 1
                    rid = f"{libtag}:{job['origin']}:L{k}"
                    start = rng.randrange(0, max(1, len(seq) - rlen))
                    s1 = apply_errors(seq[start:start+rlen+8], sub, indel, rlen)
                    s2 = ""
                    if is_pe:
                        b2 = max(0, start + FRAG - rlen)
                        s2 = revcomp(apply_errors(seq[b2:b2+rlen+8], sub, indel, rlen))
                    write_pair(rid, s1, s2)
                    origin_rows.append((rid, libtag, job["origin"], "linear", str(job["spans"]).upper(), job["role"], job["note"]))
            elif job["kind"] == "chimera":
                n = job["n"]
                for k in range(n):
                    counter += 1
                    rid = f"{libtag}:{job['origin']}:C{k}"
                    start = rng.randrange(0, max(1, len(seq) - rlen))
                    s1 = apply_errors(seq[start:start+rlen+8], sub, indel, rlen)
                    s2 = ""
                    if is_pe:
                        b2 = max(0, start + FRAG - rlen)
                        s2 = revcomp(apply_errors(seq[b2:b2+rlen+8], sub, indel, rlen))
                    write_pair(rid, s1, s2)
                    origin_rows.append((rid, libtag, job["origin"], "chimera", "TRUE", job["role"], job["note"]))
        f1.close()
        if f2: f2.close()

# ==========================================================================
#  OUTPUT: genome, GTF, truth files
# ==========================================================================
def wrap60(s): return "\n".join(s[i:i+60] for i in range(0, len(s), 60))

def write_genome():
    with open(os.path.join(OUT, "genome.fa"), "w") as fh:
        for name in CONTIG_SIZES:
            fh.write(f">{name}\n{wrap60(''.join(contigs[name]))}\n")

def write_gtf():
    with open(os.path.join(OUT, "annotation.gtf"), "w") as fh:
        fh.write("# Tiny BSJ test dataset annotation. NOTE: category-C circle exons, all\n")
        fh.write("# intronic/EIciRNA features, intergenic and antisense circles are\n")
        fh.write("# DELIBERATELY OMITTED (see README).\n")
        for g in annotated_genes:
            hidden = getattr(g, "hidden_exons", set())
            gstart = min(e[0] for e in g.exons) + 1
            gend = max(e[1] for e in g.exons)
            attr_g = f'gene_id "{g.gid}";'
            fh.write(f"{g.contig}\tsim\tgene\t{gstart}\t{gend}\t.\t{g.strand}\t.\t{attr_g}\n")
            attr_t = f'gene_id "{g.gid}"; transcript_id "{g.gid}.t1";'
            fh.write(f"{g.contig}\tsim\ttranscript\t{gstart}\t{gend}\t.\t{g.strand}\t.\t{attr_t}\n")
            for ei, (s, e) in enumerate(g.exons):
                if ei in hidden:
                    continue
                attr_e = f'gene_id "{g.gid}"; transcript_id "{g.gid}.t1"; exon_number "{ei+1}";'
                fh.write(f"{g.contig}\tsim\texon\t{s+1}\t{e}\t.\t{g.strand}\t.\t{attr_e}\n")

def write_truth():
    # BSJ TSV
    with open(os.path.join(OUT, "truth_bsj.tsv"), "w") as fh:
        cols = ["bsj_id","category","chrom","acceptor_pos_0based","donor_pos_0based",
                "circle_start_0based","circle_end_0based","strand","splice_motif",
                "circle_type","expected_junction_reads","breakpoint_tolerance_bp","notes"]
        fh.write("\t".join(cols) + "\n")
        for b in bsj_truth:
            fh.write("\t".join(str(x) for x in [
                b["id"], b["category"], b["chrom"], b["acceptor"], b["donor"],
                b["gstart"], b["gend"], b["strand"], b["motif"], b["ctype"],
                b["expected_reads"], b["tol"], b["note"]]) + "\n")
    # BSJ BED (two breakpoints as a block)
    with open(os.path.join(OUT, "truth_bsj.bed"), "w") as fh:
        for b in bsj_truth:
            lo, hi = min(b["acceptor"], b["donor"]), max(b["acceptor"], b["donor"]) + 1
            fh.write("\t".join(str(x) for x in [
                b["chrom"], lo, hi, b["id"], b["expected_reads"], b["strand"]]) + "\n")
    # negatives
    with open(os.path.join(OUT, "truth_negatives.tsv"), "w") as fh:
        fh.write("neg_id\tcategory\tchrom\tstart_0based\tend_0based\tstrand\texpected\ttrap\n")
        for d in neg_truth:
            fh.write("\t".join(str(x) for x in [
                d["id"], d["category"], d["chrom"], d["gstart"], d["gend"],
                d["strand"], "no_circle_call", d["trap"]]) + "\n")

def write_bsj_sequences():
    with open(os.path.join(OUT, "truth_bsj_sequences.fa"), "w") as fh:
        for job in read_jobs:
            if job["kind"] != "circle":
                continue
            seq = job["seq"]
            # junction window: 50bp upstream of BSJ (3' end) + 50bp downstream (5' start)
            w = min(50, len(seq))
            jwin = seq[-w:] + seq[:w]
            fh.write(f">{job['origin']}|junction_window\n{wrap60(jwin)}\n")
            fh.write(f">{job['origin']}|full_circle\n{wrap60(seq)}\n")

def write_read_origin():
    with open(os.path.join(OUT, "read_origin.tsv"), "w") as fh:
        fh.write("read_id\tlibrary\torigin\tread_class\tspans_junction\texpected_role\tnote\n")
        for r in origin_rows:
            fh.write("\t".join(r) + "\n")

# ==========================================================================
#  Derive each BSJ's actual splice motif from the genome (truth == sequence)
# ==========================================================================
def genome_sub(contig, a, b):
    if a < 0 or b > len(contigs[contig]):
        return "NN"
    return "".join(contigs[contig][a:b])

for b in bsj_truth:
    c, acc, don, strand = b["chrom"], b["acceptor"], b["donor"], b["strand"]
    if strand == "+":
        don_dn = genome_sub(c, don + 1, don + 3)
        acc_dn = genome_sub(c, acc - 2, acc)
    else:
        don_dn = revcomp(genome_sub(c, don - 2, don))
        acc_dn = revcomp(genome_sub(c, acc + 1, acc + 3))
    b["motif"] = f"{don_dn}/{acc_dn}"

# ==========================================================================
#  RUN
# ==========================================================================
write_genome()
write_gtf()
emit_fastq()
write_truth()
write_bsj_sequences()
write_read_origin()

print("contigs:", {k: len(v) for k, v in contigs.items()})
print("annotated genes:", len(annotated_genes))
print("true BSJs:", len(bsj_truth))
print("decoy loci:", len(neg_truth))
print("read jobs:", len(read_jobs))
print("total read records (origin rows):", len(origin_rows))
