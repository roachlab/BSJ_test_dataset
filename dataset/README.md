# Tiny BSJ / circRNA test dataset

A small, fully-labelled dataset for validating back-splice-junction (BSJ / circRNA)
detection pipelines. Every read's origin is recorded, so you can score **sensitivity**,
**specificity**, and **quantification**. Reproducible from a single seed
(`generate_dataset.py`, SEED=20260624).

## Files

| File | Description |
|------|-------------|
| `genome.fa` | Reference, 3 contigs: `chrT1` (700 kb), `chrT2` (600 kb), `chrP` (200 kb). |
| `annotation.gtf` | Gene models. **Deliberately omits**: category-C circle exons, all intronic (D) / EIciRNA (E) features, intergenic & antisense (M) circles. |
| `reads_{pe,se}{150,100}_{lowerr,higherr}*.fastq.gz` | 8 forward-stranded read libraries — see below. |
| `reads_{pe,se}150_lowerr_{rev,unstr}*.fastq.gz` | 4 strandedness variants (reverse-stranded, unstranded) — see below. |
| `libraries.tsv` | Per library: layout, read length, error rates, strandedness, base library. |
| `truth_bsj.tsv` | One row per true BSJ (see columns below). |
| `truth_bsj.bed` | BED6 of the two breakpoints per circle. |
| `truth_bsj_sequences.fa` | Per circle: 100 bp junction window (50 bp each side of the BSJ) + full circle sequence. |
| `truth_negatives.tsv` | Decoy loci; expected result = no circle call. |
| `read_origin.tsv` | Per-read provenance for scoring. |

## Read libraries

Full factorial: **{PE, SE} x {150, 100} bp x {low, high error} = 8 libraries**. Layout and
read length vary independently, so a PE-vs-SE difference can be attributed to the mate-pair
signal rather than to read length (and vice versa). Every library carries the **same 3568
reads from the same underlying truth** — only layout, length and error differ.

| Library | Layout | Read length | Error | Files |
|---------|--------|-------------|-------|-------|
| `pe150_lowerr` | paired | 2x150 | ~0.2% sub | `reads_pe150_lowerr_R1/R2.fastq.gz` |
| `pe150_higherr` | paired | 2x150 | ~2% sub + ~0.3% indel | `reads_pe150_higherr_R1/R2.fastq.gz` |
| `pe100_lowerr` | paired | 2x100 | ~0.2% sub | `reads_pe100_lowerr_R1/R2.fastq.gz` |
| `pe100_higherr` | paired | 2x100 | ~2% sub + ~0.3% indel | `reads_pe100_higherr_R1/R2.fastq.gz` |
| `se150_lowerr` | single | 150 | ~0.2% sub | `reads_se150_lowerr.fastq.gz` |
| `se150_higherr` | single | 150 | ~2% sub + ~0.3% indel | `reads_se150_higherr.fastq.gz` |
| `se100_lowerr` | single | 100 | ~0.2% sub | `reads_se100_lowerr.fastq.gz` |
| `se100_higherr` | single | 100 | ~2% sub + ~0.3% indel | `reads_se100_higherr.fastq.gz` |

Paired libraries use a 330 bp fragment. To isolate one axis, hold the others fixed —
e.g. `pe150_lowerr` vs `se150_lowerr` isolates paired-ness; `pe150_lowerr` vs `pe100_lowerr`
isolates read length.

### Strandedness

All eight libraries above are **forward-stranded**: R1 (or the single read) is the
transcript-sense sequence (ligation-type protocols; featureCounts `-s 1`, HISAT2 `FR`,
Salmon `ISF`/`SF`). Most total-RNA-seq (dUTP, e.g. TruSeq Stranded Total RNA) is
**reverse-stranded**: R1 is antisense (featureCounts `-s 2`, HISAT2 `RF`, Salmon `ISR`/`SR`).
Four extra libraries re-emit the **exact same fragments** as `pe150_lowerr` / `se150_lowerr`
with only the strand protocol changed, so any difference against the base library is due to
strandedness alone:

| Library | Derived from | Strandedness | How |
|---------|--------------|--------------|-----|
| `pe150_lowerr_rev` | `pe150_lowerr` | reverse | R1/R2 swapped |
| `pe150_lowerr_unstr` | `pe150_lowerr` | unstranded | each fragment fwd or swapped (50/50) |
| `se150_lowerr_rev` | `se150_lowerr` | reverse | read reverse-complemented |
| `se150_lowerr_unstr` | `se150_lowerr` | unstranded | each read as-is or reverse-complemented (50/50) |

`libraries.tsv` lists every library's layout, read length, error model and strandedness.
The evaluator scores the **strand each caller declares** for every detected BSJ against
`truth_bsj.tsv`. That catches a caller that finds the junction but reports the read-alignment
strand instead of the transcript strand, which is invisible in a forward-stranded library and
flips every call in a reverse-stranded one.

**Read-length-dependent case:** category-I circles are shorter than a read, so reads wrap the
loop and cross the BSJ more than once. `BSJ_I3` is 120 bp — between the two read lengths — so
it rolls in the 150 bp libraries but behaves as an ordinary circle in the 100 bp ones. That is
deliberate, and comparing I3 across the two lengths measures read-length dependence directly.
`BSJ_I1`/`I2`/`I4` (90/70/60 bp) roll in every library.

## Coordinate conventions

- All positions in `truth_bsj.tsv` are **0-based**.
- `acceptor_pos` = genomic coordinate of the circle's 5' start; `donor_pos` = genomic
  coordinate of the circle's 3' end. The back-splice joins **donor -> acceptor**.
- For `+` strand, `acceptor_pos < donor_pos`; for `-` strand, `acceptor_pos > donor_pos`.
- `truth_bsj.bed` reports `[min(bp), max(bp)+1)` as the feature span (0-based, BED standard).
- `splice_motif` is derived directly from `genome.fa` at the two breakpoints (strand-aware),
  so it always matches the sequence.

## truth_bsj.tsv columns

`bsj_id, category, chrom, acceptor_pos_0based, donor_pos_0based, circle_start_0based,
circle_end_0based, strand, splice_motif, circle_type, expected_junction_reads,
breakpoint_tolerance_bp, notes`

- `expected_junction_reads` — junction-spanning reads emitted for this circle **per library**.
- `breakpoint_tolerance_bp` — allowed +/- shift when matching a call to truth. Nonzero for
  fuzzy-breakpoint (H, microhomology) and indel-junction (L) cases; use it as the match window.

## Test-case categories

**48 true back-splice junctions** (`truth_bsj.tsv`, column `category`) and **39 decoy loci**
(`truth_negatives.tsv`). Instances per category:

Positives — A canonical multi-exon (5), B single-exon (4), C un-annotated exon (4),
D intronic ciRNA (4), E EIciRNA (3), F minus-strand (4), G non-canonical motif (3),
H fuzzy breakpoint/microhomology (4), I sub-read-length rolling circle (4),
J alternative back-splicing (4), K nested/overlapping (3), L indel at junction (3),
M intergenic (2) & antisense (2).

Decoys — must NOT be called: DEC1 tandem duplication (4), DEC2 linear-only gene (5),
DEC3 paralog (4), DEC4 inverted repeat (3), DEC5 readthrough/trans-splice (3),
DEC6 low-complexity (3), DEC7 repeat within a spliced gene* (4),
DEC8 duplicated exon within a gene* (4), DEC9 RT template-switching (4),
DEC10 SNP near junction (3), DEC11 high-expression linear gene (2).

\* DEC7/DEC8 are **joint** tests: the locus emits BOTH true forward-splice-junction reads
(`expected_role=host_linear`, must be retained) and pseudo-BSJ reads (`expected_role=negative`,
must be rejected).

## read_origin.tsv columns

`read_id, library, origin, read_class, spans_junction, expected_role, note`

- `read_class` — `circle_junction`, `circle_body`, `linear`, or `chimera`.
- `expected_role` — `positive` (true circle), `host_linear` (real linear transcript),
  or `negative` (decoy; should not yield a circle call).
- `read_id` is unique and encodes `library:origin:tag`.

## Scoring recipe

1. **Sensitivity** = detected / total true BSJs, breakpoints within `breakpoint_tolerance_bp`
   (>=0). Report per category and as a coverage curve using `expected_junction_reads`.
   Note: the former `BSJ_J2_1` / `BSJ_J2_2` (two internal-splicing isoforms of one
   back-splice) share identical breakpoints and are collapsed into a single junction
   `BSJ_J2` — this benchmark scores back-splice-junction discovery, not isoform resolution.
2. **Specificity / precision** = decoy loci correctly not called. For DEC7/DEC8, score jointly:
   forward SJ retained AND pseudo-BSJ rejected.
3. **Quantification** = correlation of detected junction reads vs `expected_junction_reads`.
4. Run all eight libraries (PE/SE x 150/100 bp x low/high error) to expose error-, length-
   and anchor-sensitivity, and to separate the read-length effect from the mate-pair effect.
5. **Strand** = for each detected BSJ, does the reported strand equal `truth_bsj.tsv` `strand`?
   Compare `*_rev` / `*_unstr` against their base library: detection should be unchanged and
   the declared strand should still be the transcript strand.

## Regenerating

```
python3 generate_dataset.py   # writes everything into ./dataset/
```
