# Flagellar trees vs GTDB: results (Katana run 2, 27 Sep 2026)

37 flagellar/chemotaxis genes × 149 Campylobacterota genomes. Per-gene IQ-TREE
trees (MAFFT L-INS-i, trimAl -automated1, ModelFinder, 1000 UFBoot), a
partitioned concatenated tree (9,927 aa sites, 10.8 % missing data), gene and
site concordance factors (gCF, sCFL) on both the flagellar and GTDB trees, and
targeted AU tests. Figures in `figures/` (`plot_results.py`).

Species tree: GTDB **R232** bac120 tree (15 Apr 2026, 189,801 genomes; file
pinned in `02_species_tree/gtdb_tree.py`), pruned to our genomes by exact
accession with branch lengths preserved — the relationships GTDB inferred among
these genomes within the full tree, not a tree re-inferred from them alone.

Run 1 (34 genes / 153 genomes, 26 Sep) is superseded; its headline numbers are
within 1–2 % of run 2's (table in §1).

## 0. Gene set and verification

Selection (`prepare_gene_sets.py`, `gene_selection.json`): gene in ≥ 100 of 164
genomes, ≤ 10 % of those with ambiguous extra copies; genomes kept if they carry
≥ 60 % of the selected genes (≥ 22 of 37).

Copies are resolved with OrthoFinder orthogroups: of several same-named copies,
the one in the gene's main orthogroup is used (192 calls recovered); a genome's
only copy is dropped if it sits in a paralog orthogroup that elsewhere occurs
alongside the main one (47 calls, in flgE, fliK, fliW) and kept if its
orthogroup is lineage-specific (8 Hippea/Desulfurella/Nitrosophilus/
Nitratiruptor FlgG: own orthogroup, but 47–62 % identical to FlgG across the
phylum vs ≤ 32 % to the other rod proteins in the same genome).

All 4,964 tree-input proteins were checked against five reference proteomes
spanning the phylum (`verify/`, DIAMOND very-sensitive, ≥ 50 % query cover):
94.3 % agree by majority; the 6 flagged were reviewed by sequence and are
correct calls (`verify/REVIEW.md`). flgJ/flgR, which the references don't name,
are covered by the original OrthoFinder + eggNOG confirmation.

Not included (22 genes, all < 100 genomes; flaA/flaB/maf_2356 also ambiguous
multi-copy). Some low counts are detection gaps, not absence: flgF and flhG
calls were dropped by the confirmed-only eggNOG rule, cheY/fliO/fliT are
usually unnamed in NCBI, and **flagellin detection is incomplete** — 58 genomes
with an NCBI "flagellin" protein have no flaA/flaB/flaC call (their flagellins,
mostly Arcobacteraceae and Sulfurimonas, fall in orthogroups other than the
C. jejuni-seeded OG0000010). Presence/absence claims for these genes should
read "not detected".

## 1. Flagellar genes carry a strong vertical (species-tree) signal

| Measure | Run 2 | (run 1) |
|---|---|---|
| GTDB branches recovered by the flagellar tree | **108 / 146 (74.0 %)** | 112/150 (74.7 %) |
| …of which flagellar UFBoot ≥ 95 | 101 | 105 |
| Random-shuffle null (1000 tip relabellings) | mean 0.19, max 2, **p = 0.001** | 0.22 |
| Robinson–Foulds | 76 (normalised 0.26) | 0.25 |
| Mantel, patristic distances | **r = 0.57, p = 0.0001** | 0.59 |
| gCF on flagellar tree | mean 62 %, none below 10 % | 63 % |
| gCF on GTDB tree | mean 55 %, 81/146 ≥ 50 %, 15 below 10 %, 6 at 0 % | 56 % |

The flagellar system was mostly inherited vertically, as one unit: the genes
agree with each other and recover three quarters of species-tree branches,
against ~0 expected by chance (the synteny-UPGMA tree from Thesis B recovered
23.8 %). The shuffle null is a weak baseline (random trees share almost
nothing), so RF and concordance carry the argument; the Mantel test is
supporting only (Mantel tests on tree-derived distances are criticised: Harmon
& Glor 2010; Guillot & Rousset 2013).

## 2. Most per-gene disagreement is weak signal, not transfer

Normalised RF to GTDB falls with trimmed alignment length: **Spearman
ρ = −0.84** across 37 genes (`fig2_signal_vs_length.png`). Most discordant are
the shortest — flaG (51 sites, 0.64), fliE (67), flgJ (78), fliS (110), fliN
(72), fliQ (86); least discordant are long — pflB (589 sites, 0.29), fliF, flgR,
flgS, fliR, flhB. A single gene disagreeing with GTDB is expected at this
length and is not evidence of HGT.

## 3. GTDB branches no flagellar gene recovers are unresolved, not contradicted

`fig3_concordance_gtdb.png`. For the 9 branches with gCF ≤ 5 %, the gene-tree
disagreement is mostly **gDFP** (44–100 %: gene trees messy around the clade)
rather than one consistent alternative (largest one-sided gDF 33 %), and the
**sCF is 39–86 %, above both alternatives for all 9** — at site level the
alignment still leans towards GTDB. They sit on the deep backbone
(Arcobacteraceae / Nautiliales / Desulfurellia), the Nitratiruptoraceae, and
internal Sulfurimonas branches.

The deep thermophile grouping (*Hippea*/*Desulfurella*, the GTDB root, attaching
next to *Nitrosophilus*/*Nitratiruptor* rather than Nautiliales; 6 gene trees)
looks like a composition/long-branch artefact: these are exactly the genomes
failing IQ-TREE's amino-acid composition test most often (*Nitrosophilus* ×3,
*Nitratiruptor* 7 genes each; *Hippea*, *Desulfurella* 6; *Lebetimonas* 5; 171
of 4,964 sequence-gene tests overall). To test: site-heterogeneous model
(LG+C20+F+G / PMSF) or Dayhoff-6 recoding.

## 4. HGT candidates

### Screen
`screen_hgt_candidates.py`: 501 strongly supported (UFBoot ≥ 95) conflicts
with GTDB in 178 distinct groups (`hgt_candidates.tsv`); most single-gene,
shallow, or in short genes.

### Shortlist
Explicit rule, applied to run 2 (it formalises the criteria used informally on
run 1, so it is a screening rule, not an independent test):

1. in ≥ 10 gene trees **and** the concatenated tree;
2. one-sided: at the contradicted GTDB branch, ≥ 50 % of gene trees favour this
   alternative and ≥ 3× as many as the other alternative;
3. the GTDB branch is not unusually short (≥ 25th percentile, 0.0086 subst./site).

**Primary** (pass all three): *C. mucosalis* + *C. suis*; *H. saguini* +
*H. didelphidarum*; *P. lekithochrous* + *A. roscoffensis*;
*C. pinnipediorum* / *gastrosuis* / *majalis*.
**Exploratory** (fail one criterion, tested for a stated reason):
*Hydrogenimonas* (GTDB branch 0.0079, just short; most recurrent, 20 genes),
*Sulfurimonas* subgroup (not one-sided; deepest within-genus conflict),
*Nitrosophilus* (7 genes; tests the thermophile question), *S. tamanense*
(5 genes; deep and one-sided, GTDB splits the genus).
The rule selects nothing that was not already on the list.

### AU tests
`prepare_au_tests.py` → `katana/03_au_tests.pbs` → `summarise_au.py`. For each
candidate, the constraint forces only the GTDB branch(es) the grouping
contradicts; best of 3 constrained searches (different seeds) vs the
unconstrained tree; AU with 10,000 RELL replicates, on the concatenated
alignment and every informative gene. Constrained and free concatenated trees
differ only at the forced branch(es) (1 split; 2 for the Campylobacter pair;
4 for Sulfurimonas), and the three constrained replicates agree to < 0.01 logL.

| Candidate | | Concat p (GTDB) | Holm (8) | ΔlogL | Genes rejecting GTDB, p < 0.05 / BH |
|---|---|---|---|---|---|
| *P. lekithochrous* + *A. roscoffensis* | primary | **< 10⁻⁴** | 0.0008 | 262 | 6 / 0 of 25 |
| *H. saguini* + *H. didelphidarum* | primary | **< 10⁻⁴** | 0.0008 | 181 | 6 / 0 of 36 |
| *C. mucosalis* + *C. suis* | primary | **< 10⁻⁴** | 0.0008 | 700 | 16 / 11 of 37 |
| *C. pinnipediorum* / *gastrosuis* / *majalis* | primary | **< 10⁻⁴** | 0.0008 | 700 | 12 / 9 of 36 |
| *Hydrogenimonas* | exploratory | **< 10⁻⁴** | 0.0008 | 286 | 7 / 3 of 35 |
| *Sulfurimonas* subgroup | exploratory | **< 10⁻⁴** | 0.0008 | 908 | 24 / 23 of 36 |
| *S. tamanense* | exploratory | **< 10⁻⁴** | 0.0008 | 148 | 3 / 0 of 35 |
| *Nitrosophilus* | exploratory | 0.0044 | 0.0044 | 41 | 2 / 0 of 23 |

"< 10⁻⁴": the GTDB-constrained tree won none of the 10,000 RELL replicates
(KH and SH p also 0). IQ-TREE's AU value is then an unstable extrapolation
(identical tree pairs gave 0.037 and 8 × 10⁻⁷), so 10⁻⁴ is reported and used
for Holm. No gene significantly preferred GTDB over its free tree.

**Reading.** At all eight spots the flagellar genes, taken together,
significantly reject the GTDB arrangement. Individually few genes do (except
*Sulfurimonas*: 23/36 after BH, and the *Campylobacter* pair: 9–11): the signal
is spread across many genes, i.e. the flagellar gene set shares one history
there that differs from GTDB — consistent with the flagellar system moving
(or recombining) as a unit between close relatives, or with GTDB being wrong at
those nodes. The two *Campylobacter* candidates are one conflict (their
constrained trees are identical).

**Caveats.** (1) The candidates were chosen from these same trees, so the tests
confirm the conflicts are real features of the flagellar data (not search
noise), but the p-values are post-selection and are not independent
confirmation. (2) Rejecting GTDB is not by itself HGT (see next steps).

## 5. Next steps

1. **Is GTDB wrong there?** For each candidate, check whether GTDB's own
   bac120 marker genes (or a whole-genome tree) support the GTDB branch or the
   flagellar one. Markers agreeing with the flagellar genes → species-tree
   uncertainty, not HGT. Markers agreeing with GTDB → the flagellar genes have a
   genuinely different history → HGT/recombination candidate.
2. **Synteny link**: flagellar gene order and nearby mobile elements for the
   candidate genomes (existing synteny outputs).
3. **Thermophile artefact check**: site-heterogeneous model / recoding (§3).
4. Donor search (DIAMOND vs RefSeq) only for candidates surviving 1.
5. **Planned sensitivity analysis — motility-essential genes below the
   100-genome cut-off:** motA (85 genomes) and motB (88) at a lower threshold;
   flaA/flaB after re-detecting flagellin properly (§0) and a copy rule for the
   tandem duplicates. Rerun trees + concordance and check the main results.

## Caveats

- The pruned GTDB tree's node labels are GTDB's support for the corresponding
  splits in the full 189,801-genome tree, not support for a 149-genome tree;
  they are not used in any analysis here.
- UFBoot ≥ 95 as "strong" (UFBoot is slightly optimistic under model
  misspecification).
- trimAl -automated1 trims some genes heavily (e.g. flaG → 51 sites); a less
  aggressive trimmer (BMGE/ClipKIT) would be a useful sensitivity check.
- 15 genomes excluded (< 22 of 37 genes), including *C. hominis* and
  *C. gracilis* (no selected genes) and 4 Arcobacter at 20–21 genes.
- Per-gene constrained searches: best of 3 used; 72 of 271 tests had
  replicates > 2 logL apart, and in 8 gene tests the constrained tree beat the
  original free gene tree by 0.6–7 logL (gene-tree searches are not perfectly
  optimal). Concatenated tests were stable.
- Per-gene trees did not save UFBoot trees (`--wbtl`), so reconciliation with
  ALE would need the gene-tree step rerun.
