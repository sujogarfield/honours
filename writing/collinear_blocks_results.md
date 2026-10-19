# Collinear flagellar gene clusters

Script: `scripts/07_combined_synteny/collinear_blocks.py`.
Outputs: `gene_order/collinear_blocks/`.

## Method

**Genomes.** These are the 103 reliable genomes used for the synteny
distance: every flagellar gene lies on one closed chromosome, either a native
complete assembly or a same-species RagTag scaffold. Genes present in more
than one copy are excluded, because their order is ambiguous.

**Physical adjacency.**

- Two consecutive flagellar genes are counted as adjacent when the
  intergenic gap is ≤ 300 bp.
- The gap distribution is clearly split. 35 % of gaps between neighbouring
  flagellar genes are ≤ 50 bp, typical of genes within an operon (Salgado et
  al. 2000). Only 5 % fall between 50 and 300 bp, after which the
  distribution flattens out.
- A gap of ≤ 300 bp is too short to hide an unannotated gene, so adjacent
  genes are truly next to each other.
- Adjacencies are signed: they record which gene ends face each other, as
  in the synteny distance. Two genes therefore count as collinear only if
  they are in the same order *and* the same relative orientation.

**Lineages.** These are the GTDB groups used throughout:

- *Campylobacter* (30 genomes)
- Arcobacteraceae (29)
- *Helicobacter* + *Wolinella* (16)
- *Sulfurimonas* + *Sulfuricurvum* (12)
- Nitratiruptoraceae (7)
- *Sulfurospirillum* (5)
- Nautiliales (3)
- Desulfurellia outgroup (1 genome, *Hippea*)

A lineage **holds** an adjacency when at least 50 % of its genomes have it.
Only lineages with at least 3 genomes vote.

**Core collinear blocks.**

- These are built from adjacencies held by **at least two lineages**.
- Adjacencies are added greedily, from the most to the least frequent
  overall. Each gene end is used once, and no cycles are allowed, so every
  block is a single linear chain.
- Adjacencies held by exactly one lineage are reported separately as
  **lineage-specific**. Requiring two lineages stops the largest lineages
  (*Campylobacter* 30, Arcobacteraceae 29) from defining the consensus on
  their own.

**Lineage consensus loci.** The same chaining is run separately for each
lineage, using every adjacency that lineage holds. The result is the set of
collinear clusters carried by a typical member of the lineage.

**Status of each block in each genome:**

| status | meaning |
|---|---|
| **intact** | all block genes present, every block adjacency present |
| **partial** | at least one block adjacency present, but not all |
| **dispersed** | at least 2 block genes present, none adjacent as in the block |
| **absent** | fewer than 2 block genes present |

## Results

### Ten core collinear blocks

| block | genes (order and strand) | intact | partial | dispersed | absent |
|---|---|---|---|---|---|
| C1 | *fliF*+ *fliG*+ *fliH*+ (MS ring, rotor, export) | **96** | 1 | 5 | 1 |
| C2 | *flgB*+ *flgC*+ *fliE*+ (proximal rod) | 72 | 1 | 28 | 2 |
| C3 | *flaG*+ *fliD*+ *fliS*+ (filament cap, chaperone) | 56 | 28 | 15 | 4 |
| C4 | *fliA*+ *fliM*+ *fliY*+ (σ28, switch) | 67 | 2 | 4 | 30 |
| C5 | *cheV*+ *cheA*+ *cheW*+ (chemotaxis) | 50 | 17 | 1 | 35 |
| C6 | *flgI*+ *flgJ*+ *flgM*+ (P ring, rod cap, anti-σ28) | 34 | 33 | 30 | 6 |
| C7 | *motA*+ *motB*+ (stator) | 56 | 0 | 1 | 46 |
| C8 | *fliK*+ *flgD*+ (hook length, hook cap) | 39 | 0 | 19 | 45 |
| C9 | *flgQ*+ *flgR*+ | 24 | 0 | 11 | 68 |
| C10 | *flaB*+ *maf_2356*+ (flagellin, motility-accessory factor) | 18 | 0 | 7 | 78 |

Block IDs match the figures and TSV outputs. Blocks are numbered by length,
then by adjacency support.

Many "absent" cells reflect gene calls, not loss. For example, *motA*/*motB*
and *cheV* were not called in most Arcobacteraceae, and *flgQ*/*flgR* are
named in few genomes outside *Campylobacter*.

### Genomes with each block intact, by lineage

| block | *Campylobacter* | *Helicobacter* | *Sulfurospirillum* | Arcobacteraceae | *Sulfurimonas* | Nautiliales | Nitratiruptoraceae | *Hippea* |
|---|---|---|---|---|---|---|---|---|
| C1 *fliFGH* | 30/30 | 16/16 | 4/5 | 28/29 | 11/12 | 3/3 | 3/7 | 1/1 |
| C2 *flgBC fliE* | 30/30 | 15/16 | 5/5 | **0/29** | 11/12 | 3/3 | 7/7 | 1/1 |
| C3 *flaG fliD fliS* | 30/30 | 15/16 | 4/5 | **0/29** | 0/12 | 3/3 | 4/7 | 0/1 |
| C4 *fliA fliM fliY* | 30/30 | 15/16 | 5/5 | **0/29** | 11/12 | 3/3 | 3/7 | 0/1 |
| C5 *cheV cheA cheW* | 30/30 | 12/16 | 5/5 | 0/29 | 0/12 | 0/3 | 3/7 | 0/1 |
| C6 *flgI flgJ flgM* | 16/30 | 0/16 | 4/5 | 0/29 | 11/12 | 0/3 | 3/7 | 0/1 |
| C7 *motAB* | 25/30 | 11/16 | 5/5 | 14/29 | 1/12 | 0/3 | 0/7 | 0/1 |
| C8 *fliK flgD* | 12/30 | 8/16 | 4/5 | 0/29 | 9/12 | 3/3 | 3/7 | 0/1 |
| C9 *flgQR* | 19/30 | 0/16 | 5/5 | 0/29 | 0/12 | 0/3 | 0/7 | 0/1 |
| C10 *flaB maf* | 0/30 | 14/16 | 4/5 | 0/29 | 0/12 | 0/3 | 0/7 | 0/1 |

**Main patterns.**

1. **The flagellar system is split into small, conserved collinear blocks,
   not one operon.** Blocks have 2–3 genes. The most conserved, *fliF fliG
   fliH*, is intact in 96 of 103 genomes, across every lineage except
   *Nitrosophilus*/*Nitratiruptor*. The proximal rod (*flgB flgC fliE*) and
   the switch (*fliA fliM fliY*) are intact in every lineage except
   Arcobacteraceae.
2. **Gene order follows the phylogeny.** In Figure 1, collinearity changes
   between lineages, not within them. Each lineage carries its own
   combination of blocks and lineage-specific joins, and that combination is
   shared by almost all of its members. This agrees with the vertical signal
   in the gene trees and in the synteny-distance Mantel test.
3. **Arcobacteraceae have a rearranged, more consolidated organisation, shared
   by the whole family.** In 23–29 of 29 genomes:
   - *flgB* is joined to *fliF* instead of *flgC*, and the MS-ring block is
     extended to give *flgB fliF fliG fliH fliN flgD*.
   - *flgI* is joined to *fliM* instead of *flgJ*, with *fliA* not adjacent.
   - *fliE* now sits *before* *flgC* (*fliE flgC*).
   - New clusters appear: *flgH flgK fliD fliS*, *fliI flgL flhA* and *fliR
     flhB*.

   Every one of these arrangements is shared across *Arcobacter*,
   *Aliarcobacter*, *Malaciobacter*, *Halarcobacter* and *Poseidonibacter*.
   That fits one reorganisation in the family's ancestor, inherited
   vertically afterwards, rather than repeated rearrangement or transfer.
4. **Smaller lineage-specific joins**, each held by one lineage:
   - *Campylobacter*: *flgP flgQ flgR* (30/30) and *flgM flgN flgK* (16–29/30).
   - *Helicobacter*: *flgI flgJ flgS* (13/16).
   - *Helicobacter* and *Sulfurospirillum* share *flaB*–*maf*. This is a
     two-lineage block (C10), and it is absent from *Campylobacter*.
5. **The order in *Nitrosophilus*/*Nitratiruptor* resembles the outgroup's,
   not the phylum's.** These 4 genomes:
   - lack the near-universal *fliF fliG fliH* adjacency;
   - instead carry *flgB flgC fliE fliF*, *flgG flgA flgH flgI*, *fliP fliQ
     fliR*, *flhA flhF* and *flgK flgL*.

   Within the dataset, 6 of their 8 lineage-specific adjacencies are found
   **only** in *Hippea maritima*, the Desulfurellia outgroup. Two
   explanations fit, and gene order alone can't separate them:
   - **Ancestral arrangement kept.** These resemble the operon arrangements
     common across bacteria (e.g. *fliE*–*fliF*, *fliPQR*), so the lineage may
     have kept an ancestral order that the rest of the phylum lost.
   - **Acquired from outside.** The flagellar system may have come from
     outside Campylobacterota.

   Either way, it supports the gene-tree result (RESULTS.md §4b) that the
   flagellar genes of this lineage are unlike the rest of the phylum. The
   donor search proposed there would distinguish the two.
   *Hydrogenimonas*, the other Nitratiruptoraceae genus, follows the
   phylum-wide pattern (C1, C2, C4, C5, C6 and C8 intact). That is why Nitratiruptoraceae
   shows 3/7 in several rows.

### Figures

- **`fig_collinear_blocks.png`:**
  - left: the GTDB tree with lineage colours;
  - middle: the status of each core block in each genome;
  - right: lineage-specific adjacencies, coloured by the lineage that holds
    each one.
- **`fig_lineage_loci.png`:** the consensus collinear loci of each lineage,
  drawn as arrows. Arrows are coloured by core block, so rearranged blocks
  can be followed between lineages (for example, the red *flgB flgC fliE*
  block split in Arcobacteraceae). Grey marks genes not in a core block.

Tables (TSV):

- `core_blocks.tsv`
- `group_summary.tsv`
- `genome_blocks.tsv` (per genome)
- `lineage_specific.tsv`
- `lineage_loci.tsv`

## Caveats

- **Gene-call completeness.** "Absent" mixes true absence with genes that
  weren't called. *Arcobacter acticola* has only 16 flagellar genes called
  (12 single-copy), against about 40 in its relatives, so most of its blocks
  are absent. Check whether it is non-motile or its calls are incomplete
  before interpreting that row.
- **Unnamed genes.** Some genes (*cheY*, *fliO*, *fliT*, *flgF*) are usually
  unnamed in NCBI and were not called. Blocks containing them may be longer
  than shown.
- **Flagellar genes only.** Adjacency is measured between flagellar genes.
  The 300 bp threshold ensures no other *gene* lies between two adjacent
  genes. But blocks separated by non-flagellar genes are not joined, even if
  they form one transcriptional unit.
- **Greedy chaining.** Chaining gives one linear consensus per lineage. A
  genome with an alternative arrangement held by <50 % of its lineage is
  shown only as "partial" or "dispersed".
- **Small lineages.** Lineages with 3–7 genomes (Nautiliales,
  Nitratiruptoraceae, *Sulfurospirillum*) give coarse percentages. The
  outgroup is a single genome, so it can't vote.

## Reproduce

```
python3 scripts/07_combined_synteny/collinear_blocks.py
```

## References

† = details checked online; ‡ = confirm the DOI before submission.

- Rogozin IB, Makarova KS, Murvai J, et al. (2002). Connected gene
  neighborhoods in prokaryotic genomes. *Nucleic Acids Research*
  30:2212–2223. doi:10.1093/nar/30.10.2212 † (conserved gene neighbourhoods
  and rearranged arrays)
- Tamames J (2001). Evolution of gene order conservation in prokaryotes.
  *Genome Biology* 2:research0020. doi:10.1186/gb-2001-2-6-research0020 †
- Rocha EPC (2008). The organization of the bacterial genome. *Annual Review
  of Genetics* 42:211–233. doi:10.1146/annurev.genet.42.110807.091653 †
- Salgado H, Moreno-Hagelsieb G, Smith TF, Collado-Vides J (2000). Operons in
  *Escherichia coli*: genomic analyses and predictions. *PNAS*
  97:6652–6657. doi:10.1073/pnas.120163297 ‡ (short intergenic distances
  within operons)
- Bergeron A, Mixtacki J, Stoye J (2006). A unifying view of genome
  rearrangements. *WABI 2006*, LNBI 4175:163–173. doi:10.1007/11851561_16 †
  (signed adjacencies)
