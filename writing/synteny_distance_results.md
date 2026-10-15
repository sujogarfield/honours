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
prokaryotic operons, whose gene order and orientation are conserved because
co-transcription and co-translation depend on them (Rogozin et al. 2002;
Rocha 2008):

- **Strand.** Inverting a single gene changes its transcriptional direction,
  so it can leave the operon. But the gene keeps the same neighbours, so the
  breakpoint distance does not count it.
- **Gene overlap.** Overlapping stop and start codons are common in
  prokaryotes. About a third of microbial genes overlap a neighbour. Overlaps
  are conserved and associated with translational coupling (Johnson &
  Chisholm 2004). A neighbouring pair that overlaps in one genome but not in
  another is the same pair to the breakpoint distance.

### How the equation was built, step by step

The goal (11 Oct 2026) was a distance that, unlike the breakpoint distance,
counts **inversions** and **changes in gene overlap**. Rather than inventing
a new kind of distance, it builds on the breakpoint distance and changes one
thing at a time. Each step was checked on hand-built gene orders before being
used on real data.

**Step 1. Start from the breakpoint distance.** A gene-order distance should
be 0 for identical orders and grow with rearrangement. The proportion of
neighbour pairs that are not shared does that, and it is the standard
starting point (Watterson et al. 1982; Sankoff & Blanchette 1998). Two
details from the old method were kept:

- *Shared genes only.* Gene gain and loss are a separate process from
  rearrangement, so genes missing from either genome are removed before
  comparing.
- *Normalise by the number of possible pairs.* This makes genomes with
  different numbers of flagellar genes comparable.

**Step 2. Make adjacencies signed, so inversions count.**

- **First idea, rejected:** attach the strand to each gene, giving pairs like
  (fliA+, fliM+). This breaks the rule that the orientation of an assembly is
  arbitrary. The same chromosome submitted as its reverse complement reads
  (fliM−, fliA−), which no longer matches, so every adjacency would look
  broken.
- **Solution: record which ends face each other.** Give each gene a tail and a
  head, and record an adjacency as the pair of gene ends that face each other.
  This is the extremity representation of genome rearrangement theory
  (Bergeron, Mixtacki & Stoye 2006).
  - Reading the chromosome backwards still gives the same facing ends, so whole
    orientation has no effect.
  - Inverting one gene changes which of its ends face its neighbours, so both
    its adjacencies break.
  - Inverting a block of genes leaves the ends facing each other *inside* the
    block unchanged, so only the two adjacencies at the block's edges break.
    That matches the rearrangement event: one inversion, two breakpoints. The
    same extremity model underlies the double-cut-and-join distance
    (Yancopoulos, Attie & Friedberg 2005), in which an inversion is one
    operation.

**Step 3. Treat closed chromosomes as circular.** Bacterial chromosomes are
circular, so the last flagellar gene is followed by the first. On a closed
chromosome the wrap-around pair is a real adjacency, and the denominator
becomes *n* instead of *n* − 1. It also makes the distance independent of
where the assembly's coordinate origin falls. Draft contigs are not known to
be complete, so they stay linear. Large inversions in bacteria are frequently
symmetric about the replication origin (Eisen et al. 2000), so they can
change a gene's orientation relative to distant flagellar genes.

**Step 4. Add overlap as half credit.** A conserved adjacency in which the two
genes overlap in one genome but not the other has not been rearranged. The
same ends still face each other, so scoring it 0 (fully broken) would put a
start-codon shift on the same footing as an inversion. But the relationship
has changed (for example, translational coupling through overlapping stop and
start codons, first shown in the *trp* operon by Oppenheim & Yanofsky 1980,
can be gained or lost), so scoring
it 1 would ignore the change. Half credit, *w* = ½, is the midpoint between
"conserved" and "broken". It was fixed **in advance** and not tuned to the
data (see "Sensitivity to the overlap weight" in Results).

**Step 5. Check on hand-built cases** before running on real genomes. The
reference order is *a*+ *b*+ *c*+ *d*+ with no overlaps, giving 4 adjacencies
when circular and 3 when linear.

| test case | expected behaviour | circular *d* | linear *d* |
|---|---|---|---|
| *b* flipped: *a*+ *b*− *c*+ *d*+ | both of *b*'s adjacencies break | 1 − 2/4 = 0.50 | 1 − 1/3 = 0.67 |
| block *b c* inverted: *a*+ *c*− *b*− *d*+ | only the 2 edge adjacencies break; *b*·*c* survives | 0.50 | 0.67 |
| whole order reversed: *d*− *c*− *b*− *a*− | 0, assembly orientation is arbitrary | 0 | 0 |
| rotated: *c d a b* | 0 when circular, the coordinate origin is arbitrary | 0 | 0.33 |
| *a* and *b* now overlap | half of one adjacency lost | 1 − 3.5/4 = 0.125 | 1 − 2.5/3 = 0.17 |

Under the breakpoint distance, the flip, the block inversion and the overlap
change all score 0. The block row shows the key property: inverting two genes
together costs the same as inverting one, because both are a single inversion
event. The `--breakpoint` option was also checked to reproduce the old output
files byte for byte.

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

### Relationship to published distances

**Without the overlap term (*w* = 1), the distance is a published one,
rescaled.** Let 𝒜_A and 𝒜_B be the sets of extremity adjacencies over the
shared genes.

- **SCJ distance (Feijão & Meidanis 2011).** It counts the adjacencies
  present in one genome but not the other:

  $$d_{\text{SCJ}} = |\mathcal{A}_A \setminus \mathcal{A}_B| + |\mathcal{A}_B \setminus \mathcal{A}_A|$$

  With the same genes on one chromosome, both genomes have *P* adjacencies
  (*P* = *n* circular, *n* − 1 linear), so *d*_SCJ = 2(*P* − |𝒜_A ∩ 𝒜_B|).
  Therefore

  $$d_{\text{strand-only}} = \frac{d_{\text{SCJ}}}{2P}$$

- **Breakpoint distance for signed genomes (Tannier, Zheng & Sankoff 2009).**
  It is *n* − (*A* + *T*/2), where *A* is the number of shared adjacencies and
  *T* the number of shared telomeres. A circular chromosome has no
  telomeres, so for the reliable set this is *n* − *A*, and the strand-only
  distance is exactly this divided by *n*. For linear contigs the two differ
  only in how chromosome ends are treated: here the ends are ignored and the
  denominator is *n* − 1.
- **The old breakpoint distance** used in Thesis B is the classic *unsigned*
  breakpoint distance (Watterson et al. 1982; Sankoff & Blanchette 1998).

**Not the same as the DCJ distance.** The DCJ distance (Yancopoulos, Attie &
Friedberg 2005; Bergeron, Mixtacki & Stoye 2006) uses the same extremity
adjacencies. But it counts the *minimum number of operations*, *n* − (*C* +
*I*/2), computed from cycles and paths in the adjacency graph. The
breakpoint/SCJ family counts adjacencies lost. That is simpler and does not
assume a most-parsimonious scenario, which makes it more robust to the
incomplete gene sets compared here.

**The overlap term has no published equivalent.** The closest idea is
rearrangement models that track intergenic sizes as well as gene order, such
as the wDCJ distance (Fertin, Jean & Tannier 2017). An overlap is an
intergenic size below zero. Those models track the exact size of each
intergenic region and need an NP-hard optimisation. This distance keeps only
a binary property of the intergenic region (overlapping or not) and adds it
as a weight on a conserved adjacency.

**How to state it in the methods:** "the normalised single-cut-or-join
(signed breakpoint) distance (Feijão & Meidanis 2011; Tannier et al. 2009),
modified so that conserved adjacencies whose overlap status differs count as
half conserved."

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

### Real examples from the data

Both pairs below are from the reliable set (closed chromosomes, circular).
Their breakpoint distance is 0, because the order of shared flagellar genes
is identical, but the new distance is not.

**1. A change in overlap.** *Campylobacter insulaenigrae* (GCF_000816185.1)
vs *C. ornithocola* (GCF_013201605.1), with 45 shared genes, so 45 possible
adjacencies.

| genome | *fliA* | *fliM* | relationship |
|---|---|---|---|
| *C. insulaenigrae* | + 202,453–203,139 | + 203,141–204,232 | 1 bp gap, no overlap |
| *C. ornithocola* | + 216,617–217,303 | + 217,303–218,394 | share 1 bp, overlap |

- The fliAʰ·fliMᵗ adjacency is conserved, but its overlap status differs, so
  it scores ½. All other 44 adjacencies match.
- *d* = 1 − 44.5/45 = **0.011**. The breakpoint distance is 0.
- This is a real case of the start-codon caveat. A difference of a few bases
  in where *fliM*'s start is called is enough to switch overlap on or off.
  The half weight keeps this from counting as a rearrangement.

**2. An inversion.** *Aliarcobacter lanthieri* (GCF_013201625.1) vs
*A. vitoriensis* (GCF_003316695.1), with 29 shared genes. The shared genes
are in the same circular order, but *cheB* is on the opposite strand
relative to its neighbours.

| genome | ... *motA* → *cheB* → *fliP* ... | facing ends |
|---|---|---|
| *A. lanthieri* | *motA*− *cheB*+ *fliP*+ | motAᵗ·cheBᵗ, cheBʰ·fliPᵗ |
| *A. vitoriensis* | *motA*− *cheB*− *fliP*+ (*motA* is the last gene, *cheB* the first: wrap-around) | motAᵗ·cheBʰ, cheBᵗ·fliPᵗ |

- Both of *cheB*'s adjacencies are broken. The other 27 are conserved with
  the same overlap status.
- *d* = 1 − 27/29 = **0.069**. The breakpoint distance is 0.
- "Adjacent" here means consecutive among the shared flagellar genes, not
  physically next to each other. In *A. lanthieri*, *cheB* sits about 0.8 Mb
  from *motA*. So this records that *cheB* changed orientation relative to
  the rest of the flagellar genes, through an inversion of the region
  containing it. It does not record a change within an operon.

**How common each case is.** Among all 83 reliable-set pairs that have
identical shared-gene order but non-zero new distance:

| difference | pairs |
|---|---|
| overlap only | 41 |
| orientation only | 18 |
| both | 24 |

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

### Contribution of each step

These are on the reliable set: 103 genomes, the same 5,248 pairs throughout,
circular.

| distance | Mantel *r* |
|---|---|
| breakpoint (step 1 + circular) | 0.615 |
| + signed adjacencies (strand only, step 2) | 0.633 |
| + overlap half credit (full distance, step 4) | **0.665** |

Both additions raise the correlation with phylogeny. Strand alone gives
+0.018; overlap adds a further +0.032. Strand-only does not depend on
start-codon calls, so the improvement over the breakpoint distance holds
without the overlap term.

### Sensitivity to the overlap weight

| weight *w* for an overlap mismatch | Mantel *r* |
|---|---|
| 1 (overlap ignored; = strand only) | 0.633 |
| 0.75 | 0.650 |
| **0.5 (used)** | **0.665** |
| 0.25 | 0.677 |
| 0 (overlap change = broken adjacency) | 0.688 |

- *r* rises steadily as overlap changes are penalised more. Any weight
  below 1 improves on strand-only.
- *w* was **not** chosen to maximise *r*. Picking the weight by its fit to the
  GTDB tree would be circular, because the same correlation is then used to
  argue that gene order tracks phylogeny.
- ½ was fixed in advance, on the biological argument in step 4: an overlap
  change is not a rearrangement.
- The trend does show that overlap status is strongly conserved along the
  phylogeny, which agrees with Johnson & Chisholm (2004).
- An alternative explanation is consistent annotation between close
  relatives: start codons predicted the same way give the same overlap
  calls. This cannot be fully separated from biology here.

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
  an apparent overlap. Start-codon (translation initiation site) prediction is
  the least accurate part of prokaryotic gene calling (Hyatt et al. 2010).
  All genomes here are RefSeq annotations from the same pipeline (PGAP;
  Tatusova et al. 2016), which keeps calls consistent but means shared
  errors between close relatives cannot be ruled out. The half weight limits
  the influence of any one call.
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

† = citation details checked online (15–16 Oct 2026). ‡ = details checked,
but confirm the DOI in Zotero before submission. The rest are standard
citations; import them by DOI to check.

**Gene-order distances and the extremity (signed adjacency) model**

- Watterson GA, Ewens WJ, Hall TE, Morgan A (1982). The chromosome inversion
  problem. *Journal of Theoretical Biology* 99:1–7.
  doi:10.1016/0022-5193(82)90384-8 (breakpoint idea)
- Sankoff D, Blanchette M (1998). Multiple genome rearrangement and
  breakpoint phylogeny. *Journal of Computational Biology* 5:555–570.
  doi:10.1089/cmb.1998.5.555 † (breakpoint distance)
- Yancopoulos S, Attie O, Friedberg R (2005). Efficient sorting of genomic
  permutations by translocation, inversion and block interchange.
  *Bioinformatics* 21:3340–3346. doi:10.1093/bioinformatics/bti535 ‡
  (double-cut-and-join; gene extremities; an inversion = one operation)
- Bergeron A, Mixtacki J, Stoye J (2006). A unifying view of genome
  rearrangements. In: *Algorithms in Bioinformatics (WABI 2006)*, LNBI
  4175:163–173. doi:10.1007/11851561_16 † (adjacencies of gene extremities;
  linear and circular chromosomes)
- Feijão P, Meidanis J (2011). SCJ: a breakpoint-like distance that
  simplifies several rearrangement problems. *IEEE/ACM Transactions on
  Computational Biology and Bioinformatics* 8:1318–1329.
  doi:10.1109/TCBB.2011.34 (adjacency-set breakpoint distance; closest
  published relative of this distance)

- Tannier E, Zheng C, Sankoff D (2009). Multichromosomal median and
  halving problems under different genomic distances. *BMC Bioinformatics*
  10:120. doi:10.1186/1471-2105-10-120 † (signed breakpoint distance with
  adjacencies and telomeres)
- Fertin G, Jean G, Tannier E (2017). Algorithms for computing the double cut
  and join distance on both gene order and intergenic sizes. *Algorithms for
  Molecular Biology* 12:16. doi:10.1186/s13015-017-0107-y † (wDCJ; closest
  relative of the overlap term)

**Why gene order and orientation carry signal in bacteria**

- Dandekar T, Snel B, Huynen M, Bork P (1998). Conservation of gene order: a
  fingerprint of proteins that physically interact. *Trends in Biochemical
  Sciences* 23:324–328. doi:10.1016/S0968-0004(98)01274-2 †
- Tamames J (2001). Evolution of gene order conservation in prokaryotes.
  *Genome Biology* 2:research0020. doi:10.1186/gb-2001-2-6-research0020 †
- Rogozin IB, Makarova KS, Murvai J, Czabarka E, Wolf YI, Tatusov RL,
  Székely LA, Koonin EV (2002). Connected gene neighborhoods in prokaryotic
  genomes. *Nucleic Acids Research* 30:2212–2223. doi:10.1093/nar/30.10.2212 †
- Rocha EPC (2008). The organization of the bacterial genome. *Annual Review
  of Genetics* 42:211–233. doi:10.1146/annurev.genet.42.110807.091653 †
  (operons, strand and replication-related organisation)
- Eisen JA, Heidelberg JF, White O, Salzberg SL (2000). Evidence for
  symmetric chromosomal inversions around the replication origin in
  bacteria. *Genome Biology* 1:research0011.
  doi:10.1186/gb-2000-1-6-research0011 † (inversions as a common bacterial
  rearrangement)

**Gene overlap (the half-credit term)**

- Oppenheim DS, Yanofsky C (1980). Translational coupling during expression
  of the tryptophan operon of *Escherichia coli*. *Genetics* 95:785–795.
  doi:10.1093/genetics/95.4.785 † (overlapping stop/start codons couple
  translation)
- Johnson ZI, Chisholm SW (2004). Properties of overlapping genes are
  conserved across microbial genomes. *Genome Research* 14:2268–2272.
  doi:10.1101/gr.2433104 † (~⅓ of microbial genes overlap; overlaps are
  conserved)

**Annotation (start-codon caveat)**

- Hyatt D, Chen GL, LoCascio PF, Land ML, Larimer FW, Hauser LJ (2010).
  Prodigal: prokaryotic gene recognition and translation initiation site
  identification. *BMC Bioinformatics* 11:119. doi:10.1186/1471-2105-11-119 †
- Tatusova T, DiCuccio M, Badretdin A, Chetvernin V, Nawrocki EP, Zaslavsky
  L, Lomsadze A, Pruitt KD, Borodovsky M, Ostell J (2016). NCBI prokaryotic
  genome annotation pipeline. *Nucleic Acids Research* 44:6614–6624.
  doi:10.1093/nar/gkw569 †

**Statistics**

- Mantel N (1967). The detection of disease clustering and a generalized
  regression approach. *Cancer Research* 27:209–220.
- Harmon LJ, Glor RE (2010). Poor statistical performance of the Mantel test
  in phylogenetic comparative analyses. *Evolution* 64:2173–2178.
  doi:10.1111/j.1558-5646.2010.00973.x ‡
- Guillot G, Rousset F (2013). Dismantling the Mantel tests. *Methods in
  Ecology and Evolution* 4:336–344. doi:10.1111/2041-210x.12018 ‡
