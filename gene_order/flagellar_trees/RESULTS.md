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

**Summary.** The 37 flagellar genes recover 74 % of GTDB species-tree branches
(chance ≈ 0 %); per-gene disagreement tracks gene length (ρ = −0.84) and deep
conflicts are unresolved rather than contradicted — the flagellar system is
inherited largely vertically, as a unit. Of 8 local spots where the flagellar
genes significantly reject the GTDB tree (AU p < 10⁻⁴), GTDB's own marker
genes show 4 to be species-tree problems, 3 remain unresolved, and **1 is a
genuine discordance: *Poseidonibacter lekithochrous* + *Arcobacter
roscoffensis***, where the genome backbone rejects the grouping the flagellar
genes support (p = 0.006; Holm 0.048) — non-vertical inheritance of the
flagellar genes. A single-gene screen (73 cases) adds 5 clean cases of
individual genes with a history different from the genome (fliG, flgD, and the
FlgS/FlgR regulator pair), none recent.

Run 1 (34 genes / 153 genomes, 26 Sep) is superseded; its headline numbers are
within 1–2 % of run 2's (table in §1).

## Methods

Scripts are in `scripts/` (numbered by step; `scripts/README.md`); run from the
project root. Paths below are relative to it.

### M1. Genomes and reference species tree
- **Genome set.** GTDB species-representative genomes of phylum
  Campylobacterota with a RefSeq (GCF_) accession, from a GTDB advanced search
  (`campy_fetched/gtdb-adv-search.tsv`, `01_fetch/extract_tsv.py`): 166 genomes.
  Genome sequences, RefSeq GFF annotations and proteomes were downloaded from
  NCBI (`01_fetch/fetch.py`, `fetch_proteins.py`).
- **Reference tree.** GTDB R232 bac120 tree (`bac120_r232.tree`; inferred by
  GTDB with FastTree 2.1.10 under WAG from the concatenated, masked alignment of
  120 bacterial marker genes; 189,801 genomes). Pruned to the study genomes by
  exact accession with ete3 `prune(preserve_branch_length=True)`
  (`02_species_tree/gtdb_tree.py`): all 166 found; checked against an
  independent re-prune (RF = 0, pairwise distances preserved to < 10⁻⁶).
  Rooting is GTDB's (Desulfurellia as outgroup); used only for figures — all
  comparisons are unrooted.

### M2. Flagellar and chemotaxis gene calling
1. **Annotation names** (`03_gene_order/gene_order.py`): RefSeq GFF `gene`
   features whose `gene=` tag exactly matches one of 73 flagellar, motor,
   chemotaxis, regulatory and glycosylation gene names (fl*, mot*, che*, maf*,
   pfl*). No pseudogenes were matched.
2. **Orthology extension** (`06_orthology/orthology_crossref*.py`). OrthoFinder
   3.1.5 on all 166 proteomes (DIAMOND search, FAMSA alignments, otherwise
   default settings). For each gene name, the orthogroups containing its
   name-matched proteins were taken as seeds, and every other member of those
   orthogroups became a candidate call for that gene. Flagellin (flaA/flaB),
   almost never named in RefSeq, was seeded externally: DIAMOND search of
   C. jejuni NCTC 11168 FlaA (Cj1339c) and FlaB (Cj1338c) against all
   proteomes, top hits (≥ 92 % identity, 100 % coverage) all in OG0000010.
3. **eggNOG confirmation** (`06_orthology/orthology_eggnog_confirm*.py`).
   Candidates annotated with eggNOG-mapper (emapper 6.0.17 as recorded in the
   output headers; eggNOG DB 5.0.2; DIAMOND mode, one round in HMMER mode
   against Bacteria, taxid 2). A candidate was **confirmed** if eggNOG's
   Preferred_name matched the seed gene or the gene symbol appeared as a whole
   word in its Description; "reclassified", "unconfirmed" and "no hit" calls
   were excluded.
4. **Combined gene set** (`06_orthology/build_combined_gene_order.py` →
   `gene_order/gene_order_combined.json`): all name-matched calls plus
   eggNOG-confirmed orthology calls (pflA and yuxG excluded entirely),
   deduplicated by locus. 6,640 calls in 164 genomes.

### M3. Gene order and RagTag scaffolding (synteny analyses)
Draft assemblies that split the flagellar genes across contigs were scaffolded
with RagTag 2.1.0 (`ragtag.py scaffold`, minimap2 2.31) against a
**same-species** Complete Genome or Chromosome-level NCBI reference only
(`05_ragtag/check_ragtag_feasibility.py`, `run_ragtag_scaffold.py`), so that no
other species' gene order is imposed. Gene coordinates were lifted onto the
scaffolds from RagTag's AGP (`lift_ragtag_coordinates.py`); no re-annotation.
Adjacencies across RagTag joins are inferred from the reference, so scaffolded
results are reported alongside unscaffolded ones as a sensitivity analysis.
(Gene order is not used by the tree analyses below.)

### M4. Gene selection for phylogenetics (`11_flagellar_trees/prepare_gene_sets.py`)
- **Copy resolution with OrthoFinder.** Each gene's main orthogroup is the one
  most of its proteins fall in. A genome with several same-named copies keeps
  the single copy in the main orthogroup; 2+ copies in it = ambiguous. A
  genome's only copy in another orthogroup is kept if that orthogroup never
  occurs in the same genome as the main one (lineage-specific divergent copy),
  dropped if it does (paralog family). (Added after run 1, following a check
  of the gene calls — not in response to tree results; run 1 is superseded.)
- **Selection rule** (thresholds set before any tree was inspected; modelled on GTDB's
  marker criteria, Parks et al. 2017, loosened because flagellar genes are
  absent from non-motile lineages): gene present in ≥ 100 of 164 genomes and
  ambiguous in ≤ 10 % of them; genome kept if it carries ≥ 60 % of the selected
  genes. → 37 genes × 149 genomes.
- Protein sequences from the RefSeq proteomes (terminal `*` removed); one
  FASTA per gene, header = genome accession.

### M5. Verification of the tree inputs (`11_flagellar_trees/verify_gene_calls.py`)
Every tree-input protein (4,964) was searched with DIAMOND 2.2.8 blastp
(`--very-sensitive`, e ≤ 10⁻⁵, hit covering ≥ 50 % of the query) against the
complete proteomes of five reference genomes spanning the phylum
(C. jejuni, H. pylori, A. butzleri, S. denitrificans, N. profundicola). The
best hit's own RefSeq gene name (gene tag or gene symbol in the product) was
compared with the call; a protein was flagged if more named best hits disagreed
than agreed. Flagged proteins were reviewed by pairwise sequence comparison
(`verify/REVIEW.md`). Genes the references do not name (flgJ, flgR) rely on the
M2 eggNOG confirmation.

### M6. Alignment and tree inference (Katana; `katana/01_gene_trees.pbs`, `02_concat_concordance.pbs`)
- Per gene: MAFFT 7.526 L-INS-i (`--localpair --maxiterate 1000`), trimAl 1.5.0
  `-automated1`, IQ-TREE 2.4.0 with ModelFinder (`-m MFP`, BIC) and 1000
  ultrafast bootstraps (`-B 1000`).
- Concatenated tree: IQ-TREE on the directory of trimmed alignments
  (`-p trimmed/`, edge-proportional partition model, ModelFinder per partition,
  `-B 1000`): 149 taxa, 37 partitions, 9,927 sites, 10.8 % missing data.

### M7. Concordance factors
Gene concordance (gCF; `-t <tree> --gcf loci.treefile`) and likelihood-based
site concordance (sCFL; `-te <tree> -p concat.best_model.nex --scfl 100`,
Minh et al. 2020; Mo et al. 2022), computed separately (IQ-TREE 2.4 does not
combine them) on (a) the concatenated flagellar tree and (b) the pruned GTDB
tree, using the 37 per-gene ML trees.

### M8. Comparing the flagellar and GTDB trees (`11_flagellar_trees/compare_to_gtdb.py`)
On the shared taxa, unrooted: shared non-trivial bipartitions (also counting
only flagellar branches with UFBoot ≥ 95); a null distribution from 1,000
random relabellings of the flagellar tree's tips (topology kept); Robinson–
Foulds distance, raw and normalised by 2(n − 3); Mantel test (Pearson, 9,999
permutations) between patristic distance matrices; per-gene normalised RF to
GTDB and to the flagellar tree (each pruned to the gene's taxa), related to
parsimony-informative sites by Spearman correlation.

### M9. HGT screen and shortlist (`11_flagellar_trees/screen_hgt_candidates.py`)
A conflict is a branch in a gene tree or the concatenated tree with UFBoot ≥ 95
that is **incompatible** with at least one GTDB bipartition (not merely absent
from it). Conflicts were reduced to minimal groups (no smaller conflicting
group nested inside) and annotated with the smallest GTDB clade containing the
group, the gCF of the contradicted GTDB branches, alignment length and
recurrence across trees. Shortlist rule (§4): ≥ 10 gene trees + the
concatenated tree; favoured alternative gDF ≥ 50 % and ≥ 3× the other; GTDB
branch ≥ 25th percentile length. Four further candidates were added as
exploratory with stated reasons.

### M10. AU topology tests (`prepare_au_tests.py`, `katana/03_au_tests.pbs`, `au_one.sh`, `summarise_au.py`)
For each candidate, a constraint tree = the GTDB tree collapsed to only the
branch(es) incompatible with the candidate grouping (everything else free).
Constrained ML searches (`iqtree2 -g`, 3 seeds, best log-likelihood kept) on the
concatenated alignment (models fixed from `concat.best_model.nex`) and on every
gene alignment where the constraint remains informative (gene's best-fit
model). AU test (Shimodaira 2002) on {unconstrained tree, best constrained
tree} with 10,000 RELL replicates (`-z -n 0 -zb 10000 -au`). Where the
constrained tree wins no RELL replicate, p is reported as < 10⁻⁴ (IQ-TREE's AU
value is then an unstable extrapolation). Holm correction across candidates
(concatenated), Benjamini–Hochberg across genes within a candidate.

### M11. Species-marker test (`prepare_marker_test.py`, `katana/05–07_*.pbs`, `summarise_marker_test.py`)
GTDB R232's masked concatenated bac120 marker alignment for representative
genomes (`bac120_msa_reps_r232.faa.gz`), reduced to the 149 genomes (5,010
sites; `markers/bac120_msa_r232.faa`). Unconstrained IQ-TREE tree (ModelFinder → Q.yeast+F+R7,
1000 UFBoot); per candidate, best trees with (a) the flagellar grouping forced
and (b) the GTDB branch(es) forced (3 seeds each); AU on {free, flagellar-
forced, GTDB-forced}, 10,000 RELL; Holm across candidates. If the markers
significantly reject the flagellar grouping, the flagellar genes and the genome
backbone have different histories at that spot (non-vertical inheritance
candidate); if not, the species tree is uncertain there.

### M13. Single-gene screen (`prepare_single_gene_tests.py`, `katana/08–10_*.pbs`, `single_one.sh`, `summarise_single_gene.py`)
Looks for genes that moved on their own. The 37 gene trees were re-screened
(M9 logic) against the ML marker tree (M11), rooted on Desulfurellia like GTDB,
instead of GTDB, since GTDB was shown to be wrong or uncertain at several
spots. Cases kept: UFBoot ≥ 95; not in the concatenated tree and in ≤ 3 gene
trees (recurring conflicts were covered by M10–M11); gene with ≥ 150
parsimony-informative sites; enclosing marker-tree clade ≥ 3 genomes larger
than the group. → 73 cases in 25 genes, 58 distinct groups (filters are
screening choices, set before testing). Per case: (a) gene test — the gene
alignment, unconstrained gene tree vs best tree forcing the marker-tree
branch(es) the group contradicts; (b) marker test (per group) — the marker
alignment, unconstrained marker tree vs best tree forcing the gene's grouping.
3 seeds per constrained search, AU with 10,000 RELL, bp-RELL = 0 → p < 10⁻⁴;
Benjamini–Hochberg across cases (gene tests) and groups (marker tests). A case
is a candidate only if both q < 0.05; groups containing genomes that fail the
composition test in ≥ 5 genes are labelled with a composition caveat.
Each candidate was then checked for (i) search stability (3 replicates),
(ii) how far the constrained trees differ from the free ones, (iii) whether the
group's genomes had their gene copy resolved by orthogroup (paralogy risk), and
(iv) pairwise protein identity within the group for that gene vs the same
genome pairs across the other genes (trimmed alignments; a recent transfer
would show anomalously high identity).

### M12. Software and computing
| Tool | Version | Use |
|---|---|---|
| OrthoFinder | 3.1.5 | orthogroups (DIAMOND, FAMSA) |
| eggNOG-mapper / eggNOG DB | emapper 6.0.17 (as recorded) / 5.0.2 | gene-call confirmation |
| DIAMOND | 2.2.8 | flagellin seeding, verification |
| RagTag / minimap2 | 2.1.0 / 2.31 | same-species scaffolding |
| MAFFT | 7.526 | alignment |
| trimAl | 1.5.0 | trimming |
| IQ-TREE | 2.4.0 | ModelFinder, trees, UFBoot, gCF/sCFL, constrained searches, AU |
| ete3 / Python | 3.1.3 / 3.9 (local), 3.11 (Katana) | tree handling, analysis scripts |
| GTDB | R232 | reference tree, marker alignment |

Tree inference ran on UNSW Katana (PBS Pro; conda environment
`scripts/11_flagellar_trees/katana/env.yml`); gene calling, verification,
comparisons and figures ran locally. Full rerun: `katana/submit.sh`, then
`katana/submit_au.sh`, then `katana/submit_markers.sh`, then
`katana/submit_single.sh` (after `prepare_single_gene_tests.py`); locally
`plot_results.py`, `verify_gene_calls.py`.

# Results

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

### Species-marker test (is it the flagellar genes or GTDB?)
Same 8 spots, tested on GTDB R232's own bac120 marker alignment for the 149
genomes (M11). Tree 1 = unconstrained ML marker tree, 2 = flagellar grouping
forced, 3 = GTDB branch(es) forced.

| Candidate | Marker tree has | p, flagellar grouping (Holm, 8) | p, GTDB arrangement | Verdict |
|---|---|---|---|---|
| *P. lekithochrous* + *A. roscoffensis* | GTDB branch (UFBoot 100) | **0.006 (0.048)** | 0.52 | **markers reject the flagellar grouping → non-vertical inheritance of the flagellar genes** |
| *C. mucosalis* + *C. suis* | GTDB branches (UFBoot 96) | 0.048 (0.33) | 0.59 | suggestive only |
| *C. pinnipediorum* / *gastrosuis* / *majalis* | GTDB branch (UFBoot 97) | 0.050 (0.33) | 0.61 | suggestive only (same conflict) |
| *Nitrosophilus* | GTDB branch (UFBoot 100) | 0.10 (0.51) | 0.61 | not supported |
| *Hydrogenimonas* | flagellar grouping | 0.55 | 0.43 | markers can't separate → GTDB uncertain |
| *H. saguini* + *H. didelphidarum* | flagellar grouping | 0.52 | 0.43 | markers can't separate → GTDB uncertain |
| *Sulfurimonas* subgroup | flagellar grouping | 0.48 | **< 10⁻⁴** | markers reject GTDB → GTDB tree likely wrong here |
| *S. tamanense* | flagellar grouping | 0.55 | **0.011** | markers reject GTDB → GTDB tree likely wrong here |

Where the markers side with GTDB, the unconstrained marker tree equals the
GTDB-forced tree (ΔlogL ≈ 0) and forcing the flagellar grouping costs 74–87
logL; where they side with the flagellar genes, the reverse. With only the 4
primary candidates in the Holm correction, *Poseidonibacter* is p = 0.024.

**Reading.** Of the 8 spots where the flagellar genes reject the GTDB tree:
- **4 are species-tree problems, not flagellar ones.** An ML tree built from
  GTDB's own markers on these 149 genomes agrees with the flagellar genes
  (Hydrogenimonas, H. saguini/didelphidarum) or even rejects GTDB's
  arrangement (Sulfurimonas; S. tamanense, where the flagellar and marker genes
  keep Sulfurospirillum monophyletic). GTDB's tree is built with FastTree
  across ~190,000 genomes, so short, local branches can differ from a focused
  ML analysis. At these spots the flagellar genes follow the genome.
- **1 is a genuine discordance:** *P. lekithochrous* + *A. roscoffensis*. The
  genome backbone firmly supports the GTDB arrangement (UFBoot 100) and
  significantly rejects the grouping the flagellar genes support (13 gene
  trees + concatenated tree; AU p < 10⁻⁴ on the flagellar data). The flagellar
  gene set of one of these genomes has a history different from its genome —
  non-vertical inheritance (HGT or homologous recombination) of the flagellar
  system. Borderline after correcting across all 8 (Holm 0.048).
- **3 are unresolved:** the *Campylobacter mucosalis/suis/pinnipediorum* group
  (one conflict; markers lean against the flagellar grouping, p ≈ 0.05, not
  significant after correction) and *Nitrosophilus* (p = 0.10).

### Single-gene screen (genes that moved on their own)
M13. 73 single-gene conflicts with the ML marker tree tested; 58 not supported
(the gene cannot reject the marker arrangement), 1 half-supported, **14 with
both tests significant** (BH q < 0.05). All 14 had stable constrained searches
(replicate spread ≤ 1.0 logL) and marker-side constrained trees differing only
at the forced branches (2–4 splits). None shows the signature of a *recent*
transfer: the group's identity in the gene is not anomalously high relative to
the same genomes' other genes. After the checks:

| Status | Gene | Group | Note |
|---|---|---|---|
| **clean** | fliG | *H. turcicus* / *ibis* / *winghamensis* | same group also in flgE (below) |
| **clean** | flgR | *C. rectus* / *showae* / *curvus* / *massiliensis* | with flgS below: the FlgS–FlgR flagellar two-component regulator, both in the *C. concisus/rectus* group |
| **clean** | flgS | *C. concisus* / *rectus* / *massiliensis* | |
| **clean** | flgS | *C. subantarcticus* / *peloridis* | |
| **clean** | flgD | *Sulfurimonas hydrogeniphila* / *indica* | |
| paralogy caveat | flgE ×2 | *H. rodentium*/*valdiviensis*; *H. turcicus*/*ibis*/*winghamensis* | these genomes had two flgE copies, resolved by orthogroup |
| alignment caveat | fliD ×4 | *C. peloridis*/*ornithocola*; *S. denitrificans*/*hongkongensis*; *S. gotlandica*/*marisnigri*; *Sulfurospirillum cavolei*/*oryzae* | fliD (filament cap) is the fastest-evolving gene here (identities rank 28–36 of 36–37 for these pairs) with a hard-to-align variable middle; ΔlogL up to 714 with 11–14 rearranged branches fits alignment/model problems better than transfer |
| composition caveat | flgI, fliI, fliD | Hippea/Desulfurella/Nitrosophilus/Nitratiruptor; Nautiliales | thermophile artefact (§3) |

**Reading.** Single-gene non-vertical inheritance is also rare: 5 clean
cases in 5 genes out of 73 tested, none recent. The most coherent is the
flagellar regulatory pair FlgS/FlgR disagreeing with the genome in the same
*C. concisus/rectus* clade. Caveats as for the AU tests: cases were selected
from the same gene trees (post-selection p-values), and 150–600-site genes
carry limited signal.

## 5. Next steps

1. ~~Is GTDB wrong there?~~ Done (species-marker test, §4): 1 genuine
   discordance (*P. lekithochrous* + *A. roscoffensis*), 4 GTDB problems,
   3 unresolved.
2. ~~Single-gene screen~~ Done (M13, §4): 5 clean single-gene cases, 9
   caveated. Next for these: check FlgS/FlgR gene neighbourhoods in the
   *C. concisus/rectus* clade.
3. **Characterise the Poseidonibacter/Arcobacter case**: which genome's
   flagellar genes moved (compare branch lengths / placement in each gene
   tree), flagellar gene order vs. both relatives (existing synteny outputs),
   mobile elements and GC/codon usage around the flagellar clusters.
   Optionally a whole-genome core-gene tree as a second backbone.
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
- The marker test uses GTDB's marker genes, which can themselves carry a little
  HGT; the marker tree is re-inferred with IQ-TREE on 149 genomes, so it can
  differ from GTDB's FastTree tree built in the context of 189,801 genomes.
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
