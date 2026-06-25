# Tiny BSJ Test Dataset — Design Plan

A small but statistically meaningful dataset for validating a back-splice-junction
(BSJ / circRNA) detection pipeline. Every read's origin is known, so the pipeline can be
scored for **sensitivity** (recall of true BSJs), **specificity** (rejection of look-alikes),
and **quantification** (junction-read counts). With **3–5 instances per test case**, each
category yields a real denominator for per-category sensitivity/specificity, not just a single
pass/fail. All outputs in tool-agnostic standard formats (FASTA / GTF / FASTQ + truth tables),
usable with CIRCexplorer2, CIRI2, find_circ, DCC, etc.

---

## 1. Design principles

- **Small enough for GitHub, large enough for statistics.** ~1.5 Mbp genome across a few
  contigs (a 1.5 Mbp FASTA is ≈1.5 MB; gzipped reads stay well within repo limits).
- **Replicated cases.** Each positive and decoy category is instantiated 3–5 times at
  independent loci, so per-category recall/precision is computed over a meaningful N.
- **Ground truth at the read level.** Every read records its origin (which circle / which
  linear transcript / genomic background / which decoy) and whether it spans a junction.
- **Engineered splice sites.** Real GT…AG canonical introns by default, with deliberate
  exceptions (non-canonical motif, minus strand) to probe motif filters.
- **Seeded, reproducible generation.** A single seeded RNG so genome, transcripts, and reads
  regenerate identically.

---

## 2. Reference genome layout

Target **~1.5 Mbp** total across three contigs, GC≈45% background, genes spaced generously so
flanking introns / repeats don't overlap between loci:

- `chrT1` (~700 kb) — most positive-case genes.
- `chrT2` (~600 kb) — remaining positives + most decoys.
- `chrP` (~200 kb) — paralog / multi-mapping decoys (near-identical copies placed here and on
  `chrT2` to force cross-contig multi-mapping).

Background is random sequence; exon/intron boundaries, repeats, and junction motifs are placed
deliberately. Where biogenesis realism matters, flanking introns carry long **inverted
complementary repeats** (the classic circRNA-promoting element). Loci are laid out
programmatically from a locus table, so adding/spacing instances is automatic.

---

## 3. Positive test cases (true BSJs — should be detected)

Each category gets **3–5 independent loci**, with coverage spread across orders of magnitude
*within* each category so the sensitivity floor can be mapped per case type.

| Category | Instances | What it tests | Structure | Coverage spread (junction reads) |
|----------|-----------|---------------|-----------|----------------------------------|
| A — Canonical multi-exon | 5 | Positive control / typical circle | Back-splice donor(downstream exon 3′)→acceptor(upstream exon 5′), exonic body, GT/AG | ~200, 100, 40, 12, 3 |
| B — Single-exon circle | 4 | Self-circularization of one large exon | Exon 3′ end → its own 5′ start, GT/AG | ~80, 30, 8, 2 |
| C — Un-annotated exon | 4 | BSJ where the exon is **absent from the GTF** | Real exon exists in genome + reads but omitted from annotation; circle uses it | ~60, 20, 6, 2 |
| D — Intronic (ciRNA) | 4 | Circle inside an intron, no exon overlap | Intron-contained circle (lariat-style) | ~25, 8, 3, 1 |
| E — EIciRNA (exon-intron) | 3 | Retained intron stays in the circle body | Back-splice on exons, intron retained between them | ~40, 10, 3 |
| F — Minus-strand | 4 | Strand handling / motif on − strand | As A but − strand | ~120, 40, 10, 2 |
| G — Non-canonical motif | 3 | Motif filters (keep vs drop) | Back-splice flanked by non-GT/AG (e.g. AT/AC) | ~30, 8, 2 |
| H — Fuzzy breakpoint (microhomology) | 4 | Breakpoint ambiguity | Donor/acceptor share 2–6 identical bases, so the exact breakpoint is shiftable; truth records a tolerance window | ~60, 20, 6, 2 |
| I — Sub-read-length circle (rolling) | 4 | Reads crossing the BSJ ≥2× | Circle shorter than the read; reads tile the loop and contain tandem copies | ~50, 15, 5, 2 |
| J — Alternative back-splicing | 5 (2 genes) | Isoform resolution / read attribution | Same gene → multiple circles sharing one breakpoint with a different partner, plus internal alternative splicing under one BSJ | ~80/40, 25/10, 4 |
| K — Nested / overlapping circles | 3 (1 gene) | Breakpoint disambiguation | Two circles from one gene whose coordinates nest inside each other | ~50, 18, 5 |
| L — Indel / non-templated bases at junction | 3 | Exact-vs-fuzzy calling | A few inserted/deleted bases exactly at the back-splice | ~30, 9, 3 |
| M — Intergenic & antisense circles | 4 | No annotation / strand assignment | 2 circles from fully un-annotated intergenic loci; 2 antisense circles overlapping a sense gene | ~25, 10, 4, 2 |

**≈42 true BSJs total.** The deliberately low-coverage tail (down to 1–2 reads) probes
single-read detection limits; the high-coverage members anchor quantification.

---

## 4. Decoys / negatives (look like BSJs — should NOT be called)

Each decoy category gets **3–5 loci**, giving a real false-positive denominator for
specificity. **Critically, the false breakpoints are engineered onto canonical GT/AG motifs**
so the trap cannot be escaped by motif filtering alone — a naive caller that trusts any
out-of-order read flanked by GT/AG will fire on these.

**Hard negatives — genuine splicing *and* a pseudo-BSJ in the same gene.** The most important
group: actively transcribed, correctly spliced genes that *also* emit reads mimicking a
back-splice. Here the pipeline must accept the real forward splice junctions while rejecting
the false BSJ from the same locus — it can't just blocklist the gene.

| Category | Instances | The trap |
|----------|-----------|----------|
| DEC-1 — Tandem duplication | 4 | Two identical genomic copies in a row; reads spanning the copy boundary map "out of order" and mimic a back-splice (canonical motif), but it is linear DNA. |
| DEC-2 — Linear-only gene | 5 | Normally spliced gene with strong **forward** SJ reads and no circle; forward splicing must not be miscalled as BSJ. |
| DEC-3 — Paralog pair | 4 | Near-identical copies across `chrP`/`chrT2` → multi-mapping reads that can produce spurious chimeric calls. |
| DEC-4 — Inverted repeat / palindrome | 3 | Short inverted repeat where a read aligns to both strands, mimicking a chimeric junction. |
| DEC-5 — Readthrough / trans-splice | 3 | Adjacent genes whose junction looks chimeric but is linear — distinguishing BSJ from fusion. |
| DEC-6 — Low-complexity region | 3 | Homopolymer / simple repeat producing ambiguous soft-clips and anchor artifacts. |
| **DEC-7 — Repeat within a spliced gene** | 4 | An expressed, canonically spliced gene that contains a **tandem repeat inside an exon/intron**. Reads carry **valid forward SJs** *and* reads spanning the repeat boundary look like a **valid BSJ** (out-of-order, GT/AG-flanked). Pipeline must keep the real splicing, reject the repeat-induced back-splice. |
| **DEC-8 — Duplicated exon within a gene** | 4 | A spliced gene with two **near-identical exons**; a read from the downstream copy can appear to back-splice to the upstream copy, producing a within-gene pseudo-BSJ that even sits on canonical exon boundaries — while the gene's true forward junctions are all correct. |
| **DEC-9 — RT template-switching chimera** | 4 | Reverse-transcriptase jumps between two nearby sequences sharing short microhomology, producing a chimeric read that mimics a BSJ — a genuine library artifact distinct from genomic repeats; the realistic false-positive source. |
| DEC-10 — SNP/variant near junction | 3 | A mismatch a few bp from a real breakpoint that must not be misread as a junction shift or a spurious second isoform (paired with category-H/J loci). |
| DEC-11 — High-expression linear gene | 2 | Very high linear coverage to swamp the background, so low-coverage circles must still be found (dynamic-range stress). |

**≈39 decoy loci total.** DEC-7 and DEC-8 each generate *both* truth-positive forward-SJ reads
(annotated, expected) and truth-negative pseudo-BSJ reads (expected: no circle call) from the
same locus, recorded separately in `read_origin.tsv`.

---

## 5. Read simulation strategy

A small custom Python simulator (not a heavyweight tool) — needed for controlled wrap-around
reads and exact per-read ground truth.

- **Circular reads:** build each mature circle (exons concatenated in circle order), sample
  start positions uniformly around the loop, and let reads/mate-pairs **wrap past the BSJ**.
  Track anchor length on each side; deliberately include some **very short anchors (5–8 bp)**
  to stress anchor thresholds.
- **Linear reads:** host-transcript reads for each gene so the circular-to-linear ratio is
  realistic and specificity is exercised.
- **Background + decoy reads** from genomic and decoy loci.

**Variables swept (same underlying truth across all):**

- **Library type:** a **PE** set (2×150) *and* an **SE** set (100 bp). PE adds the inverted
  mate-pair orientation signal; SE relies on junction-spanning reads only.
- **Error model:** a **low-error** replicate (~0.2% substitution) and a **high-error**
  replicate (~2%, optional indels) to stress junction mismatch tolerance.
- **Coverage:** highly variable per BSJ per §3, down to 1–2 reads.

Output sets, e.g.: `reads_pe_lowerr.fastq.gz`, `reads_pe_higherr.fastq.gz`,
`reads_se_lowerr.fastq.gz`, `reads_se_higherr.fastq.gz` (each as R1/R2 for PE).

---

## 6. Ground-truth / output files

| File | Contents |
|------|----------|
| `genome.fa` | `chrT1` + `chrT2` + `chrP` (~1.5 Mbp). |
| `annotation.gtf` | Gene models — **deliberately omits** the category-C exons and all intronic-circle features (documented). |
| `reads_*.fastq.gz` | PE/SE × low/high-error sets. |
| `truth_bsj.tsv` + `truth_bsj.bed` | One row per true BSJ: id, category, chrom, donor_pos, acceptor_pos, strand, splice_motif, circle_type, expected_junction_reads, notes. |
| `truth_bsj_sequences.fa` | Junction sequence (±N bp around the back-splice) and full circle sequence per circle. |
| `truth_negatives.tsv` | Each decoy: id, category, coordinates, the trap, expected result = no call. |
| `read_origin.tsv` | read_id → origin (circle id / linear gene / genomic / decoy) + spans-junction flag, for per-read scoring. |
| `README.md` | Coordinate conventions, file manifest, scoring recipe, expected per-category metrics. |

**Coordinate convention (documented in README):** report each BSJ as its two genomic
breakpoints — acceptor (upstream 5′ of circle) and donor (downstream 3′ of circle); the
back-splice joins donor→acceptor. Provide both 1-based inclusive and 0-based BED, and state
strand handling explicitly.

---

## 7. Scoring recipe

With replicated cases, metrics are reported **per category and overall**:

- **Sensitivity** = detected / true BSJs, breakpoints within a tolerance window — per category
  (A–M) and as a coverage curve. The tolerance window is **required** for category H/L, whose
  truth rows carry an explicit allowed breakpoint range.
- **Specificity / precision** = decoy loci correctly *not* called, and matched true calls vs
  off-target calls — per decoy category (DEC-1…8). For DEC-7/DEC-8, scored as a **joint**
  test: the forward SJs must be retained *and* the pseudo-BSJ rejected at the same locus.
- **Quantification** = correlation of detected junction reads vs `expected_junction_reads`.
- **Sensitivity floor** = lowest coverage / shortest anchor still recovered, swept across the
  PE/SE and low/high-error replicates.

---

## 8. Build order (when you give the go-ahead)

1. Generate `chrT1`/`chrT2`/`chrP` (~1.5 Mbp) from a seeded locus table, placing exons,
   introns, repeats, and junction motifs.
2. Emit `annotation.gtf` with the intentional omissions.
3. Build circle and linear transcript sequences for every locus.
4. Simulate reads: PE/SE × low/high error × the per-locus coverage profile.
5. Emit all truth files + README.
6. **Verify:** re-derive every junction sequence from genome coordinates, assert read counts
   match targets, and sanity-map a handful of reads back to confirm the back-splice signature.
