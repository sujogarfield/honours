# Tools, methods and references

Every tool, database and method used in the project, with the version
actually run, what it was used for, and its citation.

- **How this list was built:** compiled 2026-10-15 from the code (`scripts/`,
  `katana/env.yml`, tool logs), so it reflects what was run, not what was
  planned.
- **Checking the citations:** import the DOIs into Zotero (or similar) to pull
  the full records. Entries marked † were checked online on 2026-10-15; the
  rest are standard citations from the tools' own documentation and
  well-known records.
- **Methods written for this project:** the Mantel test, Spearman, UPGMA and
  classical MDS were implemented by hand in `04_synteny/synteny_analysis.py`
  (no SciPy). So cite the **methods**, not a statistics package.

---

## 1. Genomes and reference taxonomy

| what | version / detail | used for | cite |
|---|---|---|---|
| **GTDB** (Genome Taxonomy Database) | Release 232 (R11-RS232), `bac120_r232.tree`, pinned (md5 `a8b589ff…`) | species list, reference species tree, outgroup (*Desulfurellia*: *Hippea*, *Desulfurella*) | Parks et al. 2025 †; Parks et al. 2018 (original GTDB) |
| **NCBI RefSeq** / Assembly | GCF accessions | genome FASTA + GFF annotation + proteomes | O'Leary et al. 2016 |
| **NCBI Entrez** via **Biopython** | Biopython 1.85 | accession → assembly lookup and download | Cock et al. 2009 |

- Parks DH, Chaumeil P-A, Mussig AJ, Rinke C, Chuvochina M, Hugenholtz P (2025). GTDB release 10: a complete and systematic taxonomy for 715 230 bacterial and 17 245 archaeal genomes. *Nucleic Acids Research*. doi:10.1093/nar/gkaf1040 †
- Parks DH, Chuvochina M, Waite DW, Rinke C, Skarshewski A, Chaumeil P-A, Hugenholtz P (2018). A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life. *Nature Biotechnology* 36:996–1004. doi:10.1038/nbt.4229
- O'Leary NA et al. (2016). Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. *Nucleic Acids Research* 44:D733–D745. doi:10.1093/nar/gkv1189
- Cock PJA et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25:1422–1423. doi:10.1093/bioinformatics/btp163

## 2. Finding the flagellar genes (orthology)

| what | version | used for | cite |
|---|---|---|---|
| **OrthoFinder** | 3.1.5 | orthogroups across all proteomes; cross-referencing name-matched genes to find missed/renamed orthologs | Emms et al. 2026 † (OrthoFinder 3); Emms & Kelly 2019 † |
| **DIAMOND** | 2.2.8 | protein searches inside OrthoFinder and eggNOG-mapper; reference-protein searches (regulon set, `verify_gene_calls.py`) | Buchfink et al. 2015; Buchfink et al. 2021 |
| **MCL** | (inside OrthoFinder) | orthogroup clustering | Van Dongen 2000 |
| **eggNOG-mapper** | 2.1.13 (logs say "emapper-6.0.17", an internal build label) | confirming candidate orthologs (preferred name, KEGG KO); only "confirmed" calls kept | Cantalapiedra et al. 2021 † |
| **eggNOG database** | 5.0.2 (the version this eggNOG-mapper expects) | orthology database behind eggNOG-mapper | Huerta-Cepas et al. 2019 |
| **KEGG** | flagellar assembly pathway ko02040; KO annotations from eggNOG | defining the flagellar gene set | Kanehisa & Goto 2000 |

**Citing OrthoFinder 3.1.5:** OrthoFinder's own `Citation.txt` asks for Emms
& Kelly 2019, and Emms & Kelly 2015 for the original method. Version 3 now has
its own paper, so cite **Emms et al. 2026** as well.

**Citing eggNOG:** cite eggNOG-mapper **v2** and the eggNOG **5.0** database,
not eggNOG 6 or 7. The installed package (2.1.13) uses eggNOG 5.0.2.

- Emms DM, Liu Y, Belcher L, Holmes J, Kelly S (2026). OrthoFinder: improved phylogenetic orthology inference with enhanced accuracy and scalability. *Nature Methods* 23. (bioRxiv 2025 preprint: "OrthoFinder: scalable phylogenetic orthology inference for comparative genomics") †
- Emms DM, Kelly S (2019). OrthoFinder: phylogenetic orthology inference for comparative genomics. *Genome Biology* 20:238. doi:10.1186/s13059-019-1832-y †
- Emms DM, Kelly S (2015). OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. *Genome Biology* 16:157. doi:10.1186/s13059-015-0721-2
- Buchfink B, Xie C, Huson DH (2015). Fast and sensitive protein alignment using DIAMOND. *Nature Methods* 12:59–60. doi:10.1038/nmeth.3176
- Buchfink B, Reuter K, Drost H-G (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nature Methods* 18:366–368. doi:10.1038/s41592-021-01101-x
- Van Dongen S (2000). *Graph clustering by flow simulation*. PhD thesis, University of Utrecht.
- Cantalapiedra CP, Hernández-Plaza A, Letunic I, Bork P, Huerta-Cepas J (2021). eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Molecular Biology and Evolution* 38:5825–5829. doi:10.1093/molbev/msab293 †
- Huerta-Cepas J et al. (2019). eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. *Nucleic Acids Research* 47:D309–D314. doi:10.1093/nar/gky1085
- Kanehisa M, Goto S (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research* 28:27–30. doi:10.1093/nar/28.1.27

## 3. Assembly scaffolding (gene order on draft genomes)

| what | version | used for | cite |
|---|---|---|---|
| **RagTag** (`scaffold`) | 2.1.0 | ordering draft contigs against a same-species complete genome so gene order is readable on one chromosome (13 genomes in the reliable set) | Alonge et al. 2022 † |
| **minimap2** | 2.31-r1302 | aligner used by RagTag | Li 2018 |

- Alonge M, Lebeigle L, Kirsche M, Jenike K, Ou S, Aganezov S, Wang X, Lippman ZB, Schatz MC, Soyk S (2022). Automated assembly scaffolding using RagTag elevates a new tomato system for high-throughput genome editing. *Genome Biology* 23:258. doi:10.1186/s13059-022-02823-7 †
- Li H (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34:3094–3100. doi:10.1093/bioinformatics/bty191

## 4. Gene order and synteny

### 4a. The synteny distance (new, Thesis C)

`scripts/common/synteny_distance.py`, used in `orthology_synteny_preview.py`
and `synteny_analysis_full_contigs.py --closed`.

**Definition:**

1. Each gene has two ends (tail, head). A gene on the + strand reads tail →
   head left to right; on the − strand, head → tail.
2. An **adjacency** is the pair of gene *ends* that face each other between
   two consecutive genes. This is a *signed* adjacency, the standard
   "extremity" representation of genome rearrangement theory.
   - Inverting one gene breaks both of its adjacencies.
   - Inverting a whole block breaks only the two adjacencies at its ends.
   - The arbitrary orientation of an assembly has no effect.
3. **Overlap.** An adjacency conserved in both genomes counts 1 if both
   genomes agree on whether the two genes' coordinates overlap, and 0.5 if
   only one has the overlap.
4. **Distance** = 1 − conserved / possible, over the genes the two genomes
   share:
   - possible = *n* − 1 adjacencies for a linear order;
   - possible = *n* for a closed circular chromosome (the last gene is adjacent
     to the first).
5. **Minimum** of 4 shared genes. The old unsigned breakpoint distance is kept
   behind `--breakpoint` for comparison.

**What to cite for each part:**

| part | basis | cite |
|---|---|---|
| breakpoint (adjacency-loss) distance | classic gene-order distance | Watterson et al. 1982; Sankoff & Blanchette 1998 |
| signed adjacencies of gene extremities; circular chromosomes | the adjacency/extremity framework | Bergeron, Mixtacki & Stoye 2006 †; Feijão & Meidanis 2011 (adjacency-based SCJ distance) |
| why gene order carries evolutionary/functional signal in prokaryotes | gene-order conservation | Dandekar et al. 1998 †; Tamames 2001 † |
| why overlaps are scored | overlapping genes are common (~⅓ of microbial genes), conserved, and linked to regulation (translational coupling) | Johnson & Chisholm 2004 † |

**Novel in this work:**
- combining signed adjacencies with **half credit for a change in overlap
  status**;
- normalising over shared genes, with circular handling for closed
  chromosomes.

This is a design choice, motivated by the papers above, not a published
metric. Say that explicitly in the methods.

- Watterson GA, Ewens WJ, Hall TE, Morgan A (1982). The chromosome inversion problem. *Journal of Theoretical Biology* 99:1–7. doi:10.1016/0022-5193(82)90384-8
- Sankoff D, Blanchette M (1998). Multiple genome rearrangement and breakpoint phylogeny. *Journal of Computational Biology* 5:555–570. doi:10.1089/cmb.1998.5.555 †
- Bergeron A, Mixtacki J, Stoye J (2006). A unifying view of genome rearrangements. In: *Algorithms in Bioinformatics (WABI 2006)*, LNBI 4175:163–173. Springer. doi:10.1007/11851561_16 †
- Feijão P, Meidanis J (2011). SCJ: a breakpoint-like distance that simplifies several rearrangement problems. *IEEE/ACM Transactions on Computational Biology and Bioinformatics* 8:1318–1329. doi:10.1109/TCBB.2011.34
- Dandekar T, Snel B, Huynen M, Bork P (1998). Conservation of gene order: a fingerprint of proteins that physically interact. *Trends in Biochemical Sciences* 23:324–328. doi:10.1016/S0968-0004(98)01274-2 †
- Tamames J (2001). Evolution of gene order conservation in prokaryotes. *Genome Biology* 2:research0020. doi:10.1186/gb-2001-2-6-research0020 †
- Johnson ZI, Chisholm SW (2004). Properties of overlapping genes are conserved across microbial genomes. *Genome Research* 14:2268–2272. doi:10.1101/gr.2433104 †

### 4b. Comparing gene order with phylogeny

**Statistics settings:**
- *Synteny (`04_synteny/synteny_analysis.py`):* Mantel test with Pearson r,
  999 permutations, seed 0.
- *GTDB vs flagellar patristic distances (`compare_to_gtdb.py`):* Mantel test
  with Pearson r, 9,999 permutations.
- *Agreement score:* Spearman rank correlation per genome.

| method | used for | cite |
|---|---|---|
| **Mantel test** | synteny distance vs GTDB distance; flagellar vs GTDB patristic distance | Mantel 1967 |
| **Spearman rank correlation** | per-genome agreement score; informative sites vs RF distance | Spearman 1904 |
| **UPGMA** | clustering genomes by synteny distance (dashboard) | Sokal & Michener 1958 |
| **Classical MDS / PCoA** | 2-D synteny map (dashboard) | Torgerson 1952; Gower 1966 |

- Mantel N (1967). The detection of disease clustering and a generalized regression approach. *Cancer Research* 27:209–220.
- Spearman C (1904). The proof and measurement of association between two things. *American Journal of Psychology* 15:72–101. doi:10.2307/1412159
- Sokal RR, Michener CD (1958). A statistical method for evaluating systematic relationships. *University of Kansas Science Bulletin* 38:1409–1438.
- Torgerson WS (1952). Multidimensional scaling: I. Theory and method. *Psychometrika* 17:401–419. doi:10.1007/BF02288916
- Gower JC (1966). Some distance properties of latent root and vector methods used in multivariate analysis. *Biometrika* 53:325–338. doi:10.1093/biomet/53.3-4.325

## 5. Phylogenetics (UNSW Katana)

| what | version | used for | cite |
|---|---|---|---|
| **MAFFT** L-INS-i | 7.526 | per-gene protein alignment | Katoh & Standley 2013 |
| **trimAl** `-automated1` | 1.5.0 | alignment trimming | Capella-Gutiérrez et al. 2009 |
| **IQ-TREE 2** | 2.4.0 | gene trees, concatenated partitioned tree, concordance, AU tests | Minh et al. 2020a |
| **ModelFinder** (`-m MFP`) | in IQ-TREE | substitution model selection per gene / partition | Kalyaanamoorthy et al. 2017 |
| **UFBoot2** (`-B 1000`) | in IQ-TREE | branch support | Hoang et al. 2018 |
| partition model (`-p`) | in IQ-TREE | concatenated analysis with per-gene models | Chernomor et al. 2016 |
| **gCF** (`--gcf`) | in IQ-TREE | gene concordance factors on GTDB and flagellar trees | Minh, Hahn & Lanfear 2020b †; Baum 2007 (concept) |
| **sCF, likelihood version** (`--scfl 100`) | in IQ-TREE ≥ 2.2.2 | site concordance factors (the only sCF used; robust to homoplasy) | Mo et al. 2023 †; Minh, Hahn & Lanfear 2020b (original sCF) |
| **AU test** (`-au -zb 10000`, RELL) | in IQ-TREE | GTDB vs gene-tree placement of HGT candidates; single-gene and marker tests | Shimodaira 2002; Kishino, Miyata & Hasegawa 1990 (RELL) |
| **ete3** | 3.1.3 | tree handling, pruning, RF, gCF re-implementation check, figures | Huerta-Cepas, Serra & Bork 2016 |
| **UNSW Katana** | HPC (PBS) | all tree inference | UNSW Research Technology Services, doi:10.26190/669X-A286 † |

**Interpreting concordance factors** (for the discussion of gCF on short
genes; see Lanfear & Hahn 2024):
- gene-tree error biases gCF *downward*;
- it is unstable at ≤ ~30 loci.

- Katoh K, Standley DM (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution* 30:772–780. doi:10.1093/molbev/mst010
- Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics* 25:1972–1973. doi:10.1093/bioinformatics/btp348
- Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R (2020a). IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Molecular Biology and Evolution* 37:1530–1534. doi:10.1093/molbev/msaa015
- Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS (2017). ModelFinder: fast model selection for accurate phylogenetic estimates. *Nature Methods* 14:587–589. doi:10.1038/nmeth.4285
- Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS (2018). UFBoot2: improving the ultrafast bootstrap approximation. *Molecular Biology and Evolution* 35:518–522. doi:10.1093/molbev/msx281
- Chernomor O, von Haeseler A, Minh BQ (2016). Terrace aware data structure for phylogenomic inference from supermatrices. *Systematic Biology* 65:997–1008. doi:10.1093/sysbio/syw037
- Minh BQ, Hahn MW, Lanfear R (2020b). New methods to calculate concordance factors for phylogenomic datasets. *Molecular Biology and Evolution* 37:2727–2733. doi:10.1093/molbev/msaa106 †
- Baum DA (2007). Concordance trees, concordance factors, and the exploration of reticulate genealogy. *Taxon* 56:417–426.
- Mo YK, Lanfear R, Hahn MW, Minh BQ (2023). Updated site concordance factors minimize effects of homoplasy and taxon sampling. *Bioinformatics* 39:btac741. doi:10.1093/bioinformatics/btac741 †
- Lanfear R, Hahn MW (2024). The meaning and measure of concordance factors in phylogenomics. *Molecular Biology and Evolution* 41:msae214. doi:10.1093/molbev/msae214 †
- Kishino H, Miyata T, Hasegawa M (1990). Maximum likelihood inference of protein phylogeny and the origin of chloroplasts. *Journal of Molecular Evolution* 31:151–160. doi:10.1007/BF02109483
- Shimodaira H (2002). An approximately unbiased test of phylogenetic tree selection. *Systematic Biology* 51:492–508. doi:10.1080/10635150290069913
- Huerta-Cepas J, Serra F, Bork P (2016). ETE 3: reconstruction, analysis, and visualization of phylogenomic data. *Molecular Biology and Evolution* 33:1635–1638. doi:10.1093/molbev/msw046

## 6. Comparing flagellar and species trees; HGT

| method | used for | cite |
|---|---|---|
| **Robinson–Foulds distance** (normalised) | flagellar vs GTDB tree; per-gene discordance | Robinson & Foulds 1981 |
| **phylogenetic incongruence as HGT evidence** | rationale for the HGT screen (strongly supported conflict, not just lack of resolution) | Ravenhall et al. 2015 † |
| **tanglegrams** | GTDB vs flagellar tree figures (child rotation to reduce crossings; topology unchanged) | Scornavacca, Zickmann & Huson 2011 |
| **bipartition compatibility** | HGT screen: gene-tree branch with UFBoot ≥ 95 incompatible with a GTDB branch | standard split-compatibility (e.g. Ravenhall et al. 2015 for the approach) |
| **planted-transfer positive control** | sensitivity of the HGT screen (`gene_order/hgt_spikein/`) | this work; detection rate with exact binomial CI (Clopper & Pearson 1934) |
| **within-genus gCF** (pruning to one genus) | separating root/attachment conflict from internal conflict | this work, using the gCF definition of Minh et al. 2020b (re-implementation matches IQ-TREE on all 145 branches) |

- Robinson DF, Foulds LR (1981). Comparison of phylogenetic trees. *Mathematical Biosciences* 53:131–147. doi:10.1016/0025-5564(81)90043-2
- Ravenhall M, Škunca N, Lassalle F, Dessimoz C (2015). Inferring horizontal gene transfer. *PLoS Computational Biology* 11:e1004095. doi:10.1371/journal.pcbi.1004095 †
- Scornavacca C, Zickmann F, Huson DH (2011). Tanglegrams for rooted phylogenetic trees and networks. *Bioinformatics* 27:i248–i256. doi:10.1093/bioinformatics/btr210
- Clopper CJ, Pearson ES (1934). The use of confidence or fiducial limits illustrated in the case of the binomial. *Biometrika* 26:404–413. doi:10.1093/biomet/26.4.404

## 7. Motility phenotypes

| what | used for | cite |
|---|---|---|
| **BacDive** (curated motility; genome-based prediction as fallback) | motility labels per species | Schober et al. 2025 †; Reimer et al. 2022 † |
| **Madin et al. trait synthesis** | second motility source (15/15 agreement with BacDive where both present) | Madin et al. 2020 † |

- Schober I, Koblitz J, Sardà Carbasse J, Ebeling C, Schmidt ML, Podstawka A, Gupta R, Ilangovan V, Chamanara J, Overmann J, Reimer LC (2025). BacDive in 2025: the core database for prokaryotic strain data. *Nucleic Acids Research* 53(D1). †
- Reimer LC, Sardà Carbasse J, Koblitz J, Ebeling C, Podstawka A, Overmann J (2022). BacDive in 2022: the knowledge base for standardized bacterial and archaeal data. *Nucleic Acids Research* 50:D741–D746. doi:10.1093/nar/gkab961 †
- Madin JS et al. (2020). A synthesis of bacterial and archaeal phenotypic trait data. *Scientific Data* 7:170. doi:10.1038/s41597-020-0497-4 †

## 8. Software libraries and visualisation

| what | version | used for | cite |
|---|---|---|---|
| Python | 3.9.6 (local), 3.11 (Katana) | all scripts | — (or Van Rossum & Drake 2009) |
| NumPy | 2.0.2 | distance matrices, permutations, MDS | Harris et al. 2020 |
| Matplotlib | 3.9.4 | static figures (tanglegrams, gCF, spike-in) | Hunter 2007 |
| D3.js | v7 | interactive tree and synteny dashboard (`frontend/`, `docs/`) | Bostock, Ogievetsky & Heer 2011 |

- Harris CR et al. (2020). Array programming with NumPy. *Nature* 585:357–362. doi:10.1038/s41586-020-2649-2
- Hunter JD (2007). Matplotlib: a 2D graphics environment. *Computing in Science & Engineering* 9:90–95. doi:10.1109/MCSE.2007.55
- Bostock M, Ogievetsky V, Heer J (2011). D³: data-driven documents. *IEEE Transactions on Visualization and Computer Graphics* 17:2301–2309. doi:10.1109/TVCG.2011.185

---

## Methods developed in this project (describe, don't cite)

These have no single source. Describe them in the methods and justify them
with the references above:

1. **Reliable gene order** (`common/reliability.py`): all flagellar genes on
   one closed chromosome, from either a native complete assembly or a RagTag
   scaffold against the same species, with ≥ 4 genes. That gives 103 genomes
   (90 native + 13 RagTag).
2. **Lifting orthology calls onto RagTag scaffolds** through the AGP
   transform, verified by mapping every gene back with an identical sequence.
3. **Circular display rotated to the most-shared gene** (fliM), and the
   overlap/adjacency bars (interactive tree).
4. **Signed-adjacency + overlap synteny distance** (section 4a).
5. **HGT screen** (minimal conflicting groups, GTDB span, recurrence) and
   **planted-transfer positive control**.
6. **Within-genus gCF** comparison.
