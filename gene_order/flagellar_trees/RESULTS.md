# Flagellar trees vs GTDB: results (Katana run 1, 26 Sep 2026)

> **Superseded — rerun pending.** Run 1 used 34 genes / 153 genomes. The gene
> inputs were then corrected (27 Sep): when a genome has several copies named
> as the same gene, the copy in the gene's main OrthoFinder orthogroup is used
> (recovering 192 calls that were previously excluded); a genome's only copy
> is dropped if it belongs to a paralog orthogroup that elsewhere sits
> alongside the main one (47 calls, in flgE, fliK, fliW), and kept if its
> orthogroup is lineage-specific (the 8 Hippea/Desulfurella/Nitrosophilus/
> Nitratiruptor FlgG: a divergent but genuine FlgG, 47–62 % identical to FlgG
> elsewhere). Result: 37 genes (+flgE, +fliK, +fliW) × 149 genomes. Numbers
> below are from run 1 and will be replaced after the rerun.

34 flagellar/chemotaxis genes × 153 Campylobacterota genomes (selection rule in
`gene_selection.json`). Per-gene IQ-TREE trees (ModelFinder + 1000 UFBoot), a
partitioned concatenated tree (9,427 aa sites, 10.8 % missing data), gene and site concordance
factors on both the flagellar and GTDB trees. Figures in `figures/`, regenerate
with `python scripts/11_flagellar_trees/plot_results.py`.

## 1. Flagellar genes carry a strong vertical (species-tree) signal

| Measure | Value |
|---|---|
| GTDB branches recovered by the flagellar tree | **112 / 150 (74.7 %)** |
| …of which flagellar UFBoot ≥ 95 | 105 |
| Random-shuffle null (1000 tip relabellings) | mean 0.22 shared, max 2, **p = 0.001** |
| Robinson–Foulds | 76 (normalised 0.25) |
| Mantel, patristic distances | **r = 0.59, p = 0.0001** (9999 perms, 11,628 pairs) |
| gCF on flagellar tree | mean 63 %, no branch below 10 % |
| gCF on GTDB tree | mean 56 %, 83/150 branches ≥ 50 %, 18 below 10 %, 7 at 0 % |

The flagellar system was mostly inherited vertically, as one unit: the genes
agree with each other (gCF 63 %) and three quarters of species-tree branches
are recovered — against ~0 expected by chance. Compare the synteny-UPGMA tree
from Thesis B, which recovered 23.8 %.

## 2. Most per-gene disagreement is weak signal, not transfer

Normalised RF to GTDB falls with alignment length: **Spearman ρ = −0.82**
across the 34 genes (`fig2_signal_vs_length.png`). The most discordant genes
are the shortest after trimming — fliE (67 sites), flgJ (74), flaG (63),
fliQ (86), fliN (87), fliS (110). The long genes (pflB, cheA, fliF, flhA,
flgK) are the closest to GTDB. A gene tree disagreeing with GTDB on its own is
therefore not evidence of HGT here; it is expected for 60–110-site genes.

## 3. The ten GTDB branches no flagellar gene recovers are mostly unresolved, not contradicted

`fig3_concordance_gtdb.png`. For every branch with gCF ≤ 5 %:

- the gene-tree disagreement is almost all **gDFP** (77–100 %) — gene trees
  that are simply messy around a large clade — rather than one consistent
  alternative (gDF1/gDF2 mostly < 20 %). The exception is the *Poseidonibacter*
  branch, where 75 % of gene trees agree on one alternative (shortlist #1, §4);
- the **sCF is 37–76 %, above the ⅓ no-signal line for all ten**, and above
  both alternatives (sDF1/sDF2) for all ten (narrowly for the 12-genome
  Sulfurimonas branch: 37 % vs 36 %). At site level the flagellar
  alignment still leans towards the GTDB branch.

Where they sit: the deep backbone (the 56-, 41- and 29-genome clades around
Arcobacteraceae/Sulfurimonas/Nautiliales), the attachment of the GTDB root
(Desulfurellia: *Hippea*, *Desulfurella*), and internal Sulfurimonas branches.

**The deep thermophile grouping looks like an artefact.** In the flagellar
trees *Hippea*/*Desulfurella* (the GTDB root) attach next to *Nitrosophilus*/
*Nitratiruptor* instead of Nautiliales, and *Hydrogenimonas* separates from
*Nitrosophilus* (`fig1_tanglegram.png`). These are all deep-sea-vent
thermophiles, and they are exactly the genomes failing IQ-TREE's amino-acid
composition test most often (*Nitrosophilus* ×3 and *Nitratiruptor* fail in
7 genes each; *Hippea*, *Lebetimonas*, *Desulfurella* 4–5; 167 failures of
4703 sequence-gene tests overall). Thermophile proteome composition bias plus a
long outgroup branch is a textbook source of spurious grouping. Treat as
artefact unless it survives recoding / site-heterogeneous models (below).

## 4. HGT candidates

The screen (`hgt_candidates.tsv`) found 458 strongly supported (UFBoot ≥ 95)
conflicts in 178 distinct groups; most are single-gene, shallow, or in short
genes. Two kinds are worth following up.

### A. Consistent alternative histories — shortlist

The same grouping in many gene trees *and* the concatenated tree, with the gene
trees favouring **one** specific alternative (gDF strongly one-sided; noise
would split roughly evenly). `fig4_recurrent_conflicts.png`.

| # | Group (flagellar genes put together) | Gene trees | gDF of the favoured alternative | GTDB branch length | Read |
|---|---|---|---|---|---|
| 1 | *Poseidonibacter lekithochrous* + *Arcobacter roscoffensis* | 13 + concat | 75 % (other alt. 0 %) | 0.0136 (median) | **Top candidate** — GTDB branch is not short, genes very one-sided |
| 2 | *Hydrogenimonas urashimensis* + *H. cancrithermarum* | 19 + concat | 73 % (3 %) | 0.0079 | Strong, but within-genus |
| 3 | *Helicobacter saguini* + *H. didelphidarum* | 15 + concat | 61 % (3 %) | 0.0086 | Strong |
| 4 | *Campylobacter pinnipediorum* + *C. gastrosuis* + *C. majalis* | 7 + concat | 55 % (0 %) | 0.0187 (long) | Worth testing — long GTDB branch |
| 5 | *Nitrosophilus* spp. (vs *Nitratiruptor*) | 7 + concat | 64 % (9 %) | 0.0131 | Thermophile caveat (§3) |
| 6 | *Sulfurimonas denitrificans* / *crateris* / *baltica* / *marisnigri* | 9 + concat | gCF 0 on 4 GTDB branches | 0.005 | Deepest recurrent within-genus conflict (GTDB span 14) |
| — | *Sulfurospirillum diekertiae* + *S. oryzae* | 13 + concat | 58 % | **0.0031** (short) | GTDB itself uncertain here — low priority |

These are shallow (close relatives), so the realistic interpretations are:
the flagellar operon moved as a unit between close relatives (HGT or
homologous recombination), or GTDB's placement of a short branch is wrong.
Distinguishing the two needs the AU test and a check of GTDB's own marker
genes at that node.

### B. Deeper single-/few-gene placements

- ***Sulfurospirillum tamanense***: GTDB places it as sister to *Campylobacter*,
  apart from the other eight *Sulfurospirillum* (so GTDB's *Sulfurospirillum*
  is not monophyletic). The flagellar genes put it back with its genus, which
  is sister to *Campylobacter* (fliA, fliL, fliR, flgS at UFBoot 95–100, and the
  concatenated tree). Here the flagellar genes follow the named genus and GTDB
  is the odd one out, so it may be the species tree rather than the flagellar
  genes that needs explaining. Moderate priority.
- **Enterohepatic *Helicobacter*** (11 spp.) move in flhB (UFBoot 98, span 31).
- ***Arcobacter nitrofigilis*** pairs with *Halarcobacter* in fliQ/fliI/fliN —
  but fliQ and fliN are two of the shortest genes. Low priority.

## 5. Next steps

1. **Targeted AU tests** (set up: `prepare_au_tests.py`, `katana/submit_au.sh`)
   on shortlist 1–6, *C. mucosalis* + *C. suis* and *S. tamanense*: for each, the flagellar
   alignment (concatenated and per-gene) with the group constrained to its GTDB
   position vs unconstrained (`iqtree2 -g constraint.nwk`, then
   `-z trees -n 0 -zb 10000 -au`).
2. **Artefact check for the thermophiles**: rerun the concatenated tree with a
   site-heterogeneous model (`LG+C20+F+G`, or PMSF) and/or Dayhoff-6 recoding;
   see whether *Hippea*/*Desulfurella* move back next to Nautiliales.
3. **Is GTDB wrong there?** For the shallow candidates, check whether GTDB's
   own bac120 marker gene trees support the GTDB branch or the flagellar one.
   If the markers split too, it's a species-tree uncertainty, not HGT.
4. **Synteny link**: check flagellar gene order and nearby mobile elements for
   the shortlisted genomes (existing synteny outputs).
5. Donor search (DIAMOND vs RefSeq) only if a candidate survives 1–3.
6. **Planned sensitivity analysis — motility-essential genes below the
   100-genome cut-off:** motA (85 genomes) and motB (88) as ordinary genes at a
   lower threshold; flaA/flaB (tandem flagellin duplicates, 38 % of genomes
   with 2+ copies in one orthogroup) need a copy rule first, e.g. one copy per
   genome when the copies are near-identical. Rerun trees + concordance and
   check whether the main results change.

## Caveats

- UFBoot ≥ 95 was the "strong" threshold; UFBoot is known to be slightly
  optimistic under model misspecification.
- Genomes with a multi-copy call for a gene are missing from that gene tree,
  which lowers gCF for clades containing them.
- 11 genomes were excluded (< 20 of 34 genes), including *C. hominis* and
  *C. gracilis* (no selected genes) — results say nothing about those.
- Per-gene trees did not save bootstrap trees (`--wbtl`), so reconciliation
  with ALE would need the gene-tree step rerun.
