# Flagellar gene order: a signed, overlap-aware synteny distance

## Method

### Why the breakpoint distance was replaced

The earlier analysis used the normalised breakpoint distance. Two genomes are
compared over the flagellar genes they share. Every pair of neighbouring genes
in one genome that is not also a neighbouring pair in the other is a
breakpoint (Watterson et al. 1982; Sankoff & Blanchette 1998):

$$
d_{\text{BP}}(A,B) = 1 - \frac{\left|\,\{\{g,h\} : g,h \text{ adjacent in } A\} \cap \{\{g,h\} : g,h \text{ adjacent in } B\}\,\right|}{n-1}
$$

This distance ignores two features that are biologically meaningful in
prokaryotic operons:

- **Strand.** Inverting a single gene changes its transcriptional direction,
  so it can leave the operon. But the gene keeps the same neighbours, so the
  breakpoint distance does not count it.
- **Gene overlap.** Overlapping stop and start codons are common in
  prokaryotes. About a third of microbial genes overlap a neighbour. Overlaps
  are conserved and associated with translational coupling (Johnson &
  Chisholm 2004). A neighbouring pair that overlaps in one genome but not in
  another is the same pair to the breakpoint distance.

### Definition

The distance is implemented in `scripts/common/synteny_distance.py`.

**Shared genes.** For genomes *A* and *B*, let *S* be the flagellar genes
present once on the main flagellar contig of both, with *n* = |*S*|. Genes
missing from either genome are removed first, so gene loss is not counted as
rearrangement. Pairs sharing fewer than 4 genes are not scored.

**Gene ends.** Each gene *g* has two ends, a tail *g*ᵗ and a head *g*ʰ. Read
left to right along the contig:

- a gene on the + strand runs *g*ᵗ → *g*ʰ;
- a gene on the − strand runs *g*ʰ → *g*ᵗ.

**Adjacencies.** An adjacency is the unordered pair of facing ends of two
consecutive genes: {right end of gene *i*, left end of gene *i*+1}. 𝒜_A is
the set of adjacencies in *A*. On a closed circular chromosome, the last and
first genes also form an adjacency. This is the signed adjacency, or
extremity, representation of genome rearrangement theory (Bergeron, Mixtacki
& Stoye 2006; Feijão & Meidanis 2011).

**Overlap.** For each adjacency *x*, *o*_A(*x*) is true when the annotated
coordinates of its two genes intersect in *A*.

**Distance.**

$$
d(A,B) \;=\; 1 \;-\; \frac{\displaystyle\sum_{x \,\in\, \mathcal{A}_A \cap \mathcal{A}_B} w(x)}{P},
\qquad
w(x)=\begin{cases}1 & o_A(x)=o_B(x)\\[2pt] \tfrac12 & o_A(x)\neq o_B(x)\end{cases},
\qquad
P=\begin{cases}n-1 & \text{linear order}\\ n & \text{closed circular chromosome}\end{cases}
$$

*d* = 0 when order, orientation and overlap status are all identical. *d* = 1
when no adjacency is shared.

**Properties:**

- Inverting a single gene breaks both of its adjacencies.
- Inverting a block of genes breaks only the two adjacencies at its ends,
  because the block's internal adjacencies are unchanged.
- Flipping a whole contig (the arbitrary orientation of an assembly) has no
  effect.
- An adjacency conserved with a change in overlap status counts as half
  conserved.
- Because a conserved signed adjacency implies a conserved neighbour pair,
  and *w* ≤ 1, *d* ≥ *d*_BP for any pair of genomes.

**Novel elements.** Signed adjacencies come from the established extremity
framework. Weighting conserved adjacencies by agreement in overlap status is
a design choice made for this work, motivated by the conservation of
overlapping genes. It is not a published metric.

### Worked example

| genome | gene order (strand) | overlap |
|---|---|---|
| *A* | *fliA*+ *fliM*+ *fliN*+ | *fliM*/*fliN* overlap |
| *B* | *fliA*+ *fliM*− *fliN*+ | none |
| *C* | *fliA*+ *fliM*+ *fliN*+ | none |

- 𝒜_A = {fliAʰ·fliMᵗ, fliMʰ·fliNᵗ}.
- 𝒜_B = {fliAʰ·fliMʰ, fliMᵗ·fliNᵗ}, because *fliM* is inverted.
- *A* vs *B*: no adjacency is shared, so *d* = 1 − 0/2 = **1.0**. The
  breakpoint distance is **0**, because the neighbours are unchanged.
- *A* vs *C*: both adjacencies are shared, but one differs in overlap status,
  so *d* = 1 − (1 + ½)/2 = **0.25**. The breakpoint distance is **0**.

### Comparison with phylogeny

**Data.**

- Gene order came from the combined OrthoFinder + eggNOG flagellar gene
  calls.
- The *reliable* set is restricted to genomes whose flagellar region lies on
  a closed chromosome: a native complete assembly or a same-species RagTag
  scaffold (see M3). For this set the circular form of the distance was used.
- The *all genomes* and *scaffolded* sets use the linear form on each
  genome's main flagellar contig.

**Test.**

- Pairwise synteny distances were compared with patristic distances on the
  GTDB reference tree.
- The test was a Mantel test: Pearson *r* over all scorable pairs, 999
  permutations of genome labels, seed 0 (Mantel 1967).
- The minimum attainable *p* is 0.001.
- The old breakpoint distance was recomputed on exactly the same genomes and
  pairs (`--breakpoint`; circular for the reliable set).

## Results

### The new distance tracks phylogeny more closely in every genome set

| genome set | genomes | pairs | Mantel *r*, breakpoint | Mantel *r*, signed + overlap | *p* |
|---|---|---|---|---|---|
| all genomes | 164 | 10,749 | 0.367 | **0.471** | ≤ 0.001 |
| all genomes, RagTag-scaffolded | 164 | 11,266 | 0.383 | **0.487** | ≤ 0.001 |
| reliable (closed chromosome, circular) | 103 | 5,248 | 0.615 | **0.665** | ≤ 0.001 |

**Reliable set.** It started from 164 genomes. 52 were excluded for
incomplete flagellar coverage on a closed chromosome, and 9 because they
shared fewer than 4 genes.

**What the table shows:**

- Gene order is strongly correlated with phylogenetic distance in all three
  sets, under both distances.
- Adding strand and overlap information raised the correlation in every set:
  by 0.10 for all genomes and 0.05 for the reliable set.
- So orientation and overlap changes are not noise. They accumulate with
  phylogenetic distance in the same way as changes in gene order.
- Restricting to closed chromosomes gave the largest gain (0.47 → 0.67).
  Much of the scatter in the full set comes from fragmented assemblies, not
  from biology.

### What the new distance detects that the breakpoint distance does not

These figures are for the reliable set, 5,248 pairs.

| | breakpoint | signed + overlap |
|---|---|---|
| median distance | 0.649 | 0.744 |
| mean distance | 0.623 | 0.712 |
| pairs with *d* = 0 | 113 (2.2 %) | 30 (0.6 %) |
| pairs with *d* = 1 | 13 | 47 |

- **83 pairs** had identical neighbour order (*d*_BP = 0) but differed in
  gene orientation or overlap status. These differences are invisible to the
  breakpoint distance.
- The new distance was higher for 5,010 pairs and equal for 238, and never
  lower. This is expected from its definition.
- The higher correlation is therefore not a rescaling. The extra differences
  the new distance detects are distributed in line with phylogeny.

### Caveats

- **Start-codon calls.** Overlap status is read from annotated coordinates,
  so it depends on start-codon calls. A mis-called start can create or remove
  an apparent overlap. The half weight limits the influence of any one call.
- **Mantel test.** Mantel tests on distance matrices derived from a tree have
  inflated type I error and are criticised as tests of phylogenetic signal
  (Harmon & Glor 2010; Guillot & Rousset 2013). *r* is used here to compare
  the two distances on identical pairs, not as standalone evidence. The
  tree-based analyses (RF, gCF) carry that argument.
- **Pairwise, not ancestral.** The distance is pairwise and does not
  reconstruct the rearrangements themselves, so it cannot say where on the
  tree a rearrangement happened.

## Reproduce

```
# signed + overlap (default)
python3 scripts/07_combined_synteny/orthology_synteny_preview.py gene_order/gene_order_combined.json gene_order/synteny_analysis_combined.json
python3 scripts/07_combined_synteny/orthology_synteny_preview.py gene_order/gene_order_combined_scaffolded.json gene_order/synteny_analysis_combined_scaffolded.json
python3 scripts/04_synteny/synteny_analysis_full_contigs.py gene_order/gene_order_combined_scaffolded.json gene_order/synteny_analysis_combined_reliable.json --closed
# breakpoint comparison: add --breakpoint and write to a different output
```

## References

- Bergeron A, Mixtacki J, Stoye J (2006). A unifying view of genome
  rearrangements. *WABI 2006*, LNBI 4175:163–173. doi:10.1007/11851561_16
- Feijão P, Meidanis J (2011). SCJ: a breakpoint-like distance that
  simplifies several rearrangement problems. *IEEE/ACM TCBB* 8:1318–1329.
  doi:10.1109/TCBB.2011.34
- Guillot G, Rousset F (2013). Dismantling the Mantel tests. *Methods in
  Ecology and Evolution* 4:336–344. doi:10.1111/2041-210x.12018
- Harmon LJ, Glor RE (2010). Poor statistical performance of the Mantel test
  in phylogenetic comparative analyses. *Evolution* 64:2173–2178.
  doi:10.1111/j.1558-5646.2010.00973.x
- Johnson ZI, Chisholm SW (2004). Properties of overlapping genes are
  conserved across microbial genomes. *Genome Research* 14:2268–2272.
  doi:10.1101/gr.2433104
- Mantel N (1967). The detection of disease clustering and a generalized
  regression approach. *Cancer Research* 27:209–220.
- Sankoff D, Blanchette M (1998). Multiple genome rearrangement and
  breakpoint phylogeny. *Journal of Computational Biology* 5:555–570.
  doi:10.1089/cmb.1998.5.555
- Watterson GA, Ewens WJ, Hall TE, Morgan A (1982). The chromosome inversion
  problem. *Journal of Theoretical Biology* 99:1–7.
  doi:10.1016/0022-5193(82)90384-8
