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
| `reads_pe_lowerr_R1/R2.fastq.gz` | Paired-end 2x150, ~0.2% substitution error. |
| `reads_pe_higherr_R1/R2.fastq.gz` | Paired-end 2x150, ~2% substitution + ~0.3% indel error. |
| `reads_se_lowerr.fastq.gz` | Single-end 100 bp, low error. |
| `reads_se_higherr.fastq.gz` | Single-end 100 bp, high error. |
| `truth_bsj.tsv` | One row per true BSJ (see columns below). |
| `truth_bsj.bed` | BED6 of the two breakpoints per circle. |
| `truth_bsj_sequences.fa` | Per circle: 100 bp junction window (50 bp each side of the BSJ) + full circle sequence. |
| `truth_negatives.tsv` | Decoy loci; expected result = no circle call. |
| `read_origin.tsv` | Per-read provenance for scoring. |

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

Positives (`truth_bsj.tsv`, column `category`):
A canonical multi-exon, B single-exon, C un-annotated exon, D intronic ciRNA, E EIciRNA,
F minus-strand, G non-canonical motif, H fuzzy breakpoint/microhomology, I sub-read-length
rolling circle, J alternative back-splicing, K nested/overlapping, L indel at junction,
M intergenic & antisense.

Decoys (`truth_negatives.tsv`, column `category`) — must NOT be called:
DEC1 tandem duplication, DEC2 linear-only gene, DEC3 paralog, DEC4 inverted repeat,
DEC5 readthrough/trans-splice, DEC6 low-complexity, DEC7 repeat within a spliced gene*,
DEC8 duplicated exon within a gene*, DEC9 RT template-switching, DEC10 SNP near junction,
DEC11 high-expression linear gene.

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
2. **Specificity / precision** = decoy loci correctly not called. For DEC7/DEC8, score jointly:
   forward SJ retained AND pseudo-BSJ rejected.
3. **Quantification** = correlation of detected junction reads vs `expected_junction_reads`.
4. Run all four libraries (PE/SE x low/high error) to expose error- and anchor-sensitivity.

## Regenerating

```
python3 generate_dataset.py   # writes everything into ./dataset/
```
