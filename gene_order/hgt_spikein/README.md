# Planted-transfer positive control for the HGT screen

## Why

The flagellar gene trees mostly agree with the GTDB species tree, which supports
"predominantly vertical inheritance". That conclusion only means something if the
pipeline *would* have found a transfer had one happened. Here, transfers are
planted on purpose, where the right answer is known, and the unmodified
pipeline is run to see how often it finds them. It is also run to find **where
detection stops working**.

## Design (run 2: 72 tests)

Every test is built from copies of a real gene file in
`gene_order/flagellar_trees_ext/genes/`. Those source files are only read; the
prepare script checks their md5 sums before and after.

**Genes:** one short, one medium and one long. Signal scales with length; the
trimmed alignment lengths come from the 42-gene run.

| gene | length class | trimmed sites |
|---|---|---|
| *fliE* | short | 67 |
| *flhB* | medium (≈ median of the 42) | 297 |
| *flhA* | long | 556 |

**Recipients:** three, in different parts of the tree, because the recipient's
own position and branch length affect power.

| code | recipient | family |
|---|---|---|
| cjej | *Campylobacter jejuni* | Campylobacteraceae |
| hfel | *Helicobacter felis* | Helicobacteraceae (long-branched) |
| abut | *Aliarcobacter butzleri* | Arcobacteraceae |

*H. felis* was chosen over other *Helicobacter* because it is the one with a
genuinely close donor group (GTDB distance 0.055).

**Donor levels:** four, defined by position in the GTDB tree. The GTDB
distance (substitutions/site) is reported alongside, because branch lengths
differ a lot between families. *Helicobacter* is long-branched, so a "sister
family" donor is further from *H. felis* in substitutions than a "distant"
donor is from *C. jejuni*.

| level | rule | cjej donor (dist) | hfel donor (dist) | abut donor (dist) |
|---|---|---|---|---|
| very close | nearest GTDB group that leaves ≥ 2 of the donor's relatives once the donor is removed | *C. coli* (0.037) | *H. salomonis* (0.055) | *A. vitoriensis* (0.061) |
| within family | other half of the recipient's family | *C. fetus* (0.247) | *H. pullorum* (0.447) | *Malaciobacter marinus* (0.158) |
| sister family | GTDB clade sister to the recipient's family | *Sulfurospirillum deleyianum* (0.332) | *S. deleyianum* (0.590) | *Nitrosophilus alvini* (0.391) |
| distant | a family further away | *Sulfurimonas denitrificans* (0.543) | *S. denitrificans* (0.672) | *C. fetus* (0.569) |

**Donor relatives:** for each donor, its closest GTDB group with at least 2
genomes left once the donor is removed, and not containing the recipient.
These are listed in `design.tsv`. They are where a transferred copy should
land. Every donor has at least 2 in all three gene files, so the test is
always possible. (For *C. jejuni*, *S. tamanense* was rejected as the sister
family donor: it is the lone sister of *Campylobacter* and leaves no
relatives.)

**Conditions:** each (gene, recipient, donor) combination gets two files.

- **spike**: the recipient's protein sequence is replaced by the donor's
  sequence. The header is kept, so the recipient now "carries the donor's
  gene". The donor itself is removed from the file. Otherwise the planted copy
  would be identical to a sampled sequence and trivially sit next to it. With
  the donor gone, the copy has to group with the donor's relatives. That is
  what a gene acquired from that lineage some time ago looks like.
- **control**: the donor is removed and the recipient is left untouched. This
  separates the effect of the transfer from the effect of dropping one genome,
  and catches pre-existing conflicts around the recipient that would otherwise
  look like detections.

That gives 3 genes × 3 recipients × 4 donor levels × 2 conditions = **72 test
genes** (36 planted transfers). `design.tsv` lists every test.

## Pipeline (identical to the real run)

1. **Gene trees, on Katana.** This uses `katana/01_gene_trees.pbs` unchanged:
   MAFFT L-INS-i, then trimAl `-automated1`, then IQ-TREE ModelFinder with
   1000 UFBoot.
2. **HGT screen.** This uses `screen_hgt_candidates.py` unchanged. A candidate is
   a gene-tree branch with UFBoot ≥ 95 that is incompatible with a GTDB branch,
   keeping only the minimal conflicting groups.
3. **Scoring.** `summarise_hgt_spikein.py`, for each tree:
   - **placement**: where the recipient lands. Its sister group is the smaller
     subtree it joins in the unrooted gene tree.
     - `donor`: the sister lies on the donor's side of the recipient–donor
       split in GTDB and includes at least one donor relative.
     - `home`: the sister lies on the recipient's own side.
     - `elsewhere`: anything else.
     - The UFBoot of that branch is also reported.
   - **screen hit**: the screen reports a candidate group made of the recipient
     plus donor-side genomes only, including at least one donor relative.
   - **rank**: that candidate's position among all of the tree's candidates,
     ordered by depth (GTDB span minus group size), then UFBoot. Rank 1 means
     it is the deepest conflict in the tree. Each real gene tree carries 9–18
     background conflicts, so rank shows whether the transfer would stand out
     without knowing where to look.
   - **detected**: the spike tree has a screen hit and its control does not.
   - The overall rate is reported with an exact (Clopper–Pearson) 95%
     confidence bound.
4. **Figures.** `plot_hgt_spikein.py` draws one tanglegram per planted
   transfer, plus an overview per gene × recipient (controls on top, planted
   below, four donor levels across), with the recipient's line in red.

## Dry run (scoring check, not a result)

To check the scoring before running on Katana, the real 42-gene trees were
used as stand-ins. In each spike stand-in, the donor's leaf was relabelled as
the recipient; in each control stand-in, the donor was pruned. That gave 24/36
detected and 1/36 controls flagged:

| donor level | detected |
|---|---|
| very close | 3/9 |
| within family | 5/9 |
| sister family | 9/9 |
| distant | 7/9 |

A stand-in is only detected if the real gene tree already places the donor
with UFBoot ≥ 95 in a group that conflicts with GTDB once relabelled. So the
misses show where single flagellar genes do not resolve a donor's position
strongly enough. Mostly these are very close donors, and within-family donors
in the short *fliE*.

The flagged control was *flhB* for *H. felis*: in the real tree, *H. felis*
already groups with *H. labacensis* and *H. vulpis* differently from GTDB.
That is a pre-existing conflict, and it is why the control is needed. The real
run realigns and re-infers every tree, so its numbers will differ.

## Expected outcome and interpretation

- **Sister family and distant transfers** should be detected in most genes.
- **Very close transfers** are expected to be missed often. A *C. coli* gene in
  *C. jejuni* barely differs from *C. jejuni*'s own.
- **The detection rate per level × gene length** sets the resolution limit of
  the "predominantly vertical" conclusion. A known failure point (for example,
  "transfers between sister species are not detectable with single flagellar
  genes") makes the positive result above it more convincing.

## Limits

- **Whole-gene replacement.** Partial (recombined) transfers would be harder to
  detect.
- **Sampled donor lineage.** The donor's relatives are in the dataset. A
  transfer from an unsampled lineage would show up as a long, poorly placed
  branch instead.
- **One donor per level per recipient.** With 3 recipients and 3 genes this is
  9 tests per level: a spot check across the tree, not a full power curve.

## How to run

```
# local (done; genes/ and design.tsv are committed)
python scripts/11_flagellar_trees/prepare_hgt_spikein.py
# Katana, from the repo root, after git pull (72 array elements)
rm -rf gene_order/hgt_spikein/{aligned,trimmed,gene_trees}   # clear run-1 outputs
qsub -v FLAG_DIR=gene_order/hgt_spikein -J 1-72 scripts/11_flagellar_trees/katana/01_gene_trees.pbs
# local, after rsync of gene_order/hgt_spikein/gene_trees/ back
FLAG_DIR=gene_order/hgt_spikein python scripts/11_flagellar_trees/screen_hgt_candidates.py
python scripts/11_flagellar_trees/summarise_hgt_spikein.py   # writes results.tsv
python scripts/11_flagellar_trees/plot_hgt_spikein.py        # writes figures/
```

`aligned/`, `trimmed/`, `gene_trees/` and the screen outputs are gitignored.
`results.tsv` is kept.

## Run 2 results (2026-10-15): 36 planted transfers

**Detected 22/36 (61%; exact 95% CI lower bound 43%). Controls flagged 1/36.**
The planted copy was placed with the donor's relatives in 32/36 trees. A
detected transfer was the deepest conflict in its tree in 18/22, so it would
stand out without knowing where to look.

Cell = detected (rank among the tree's conflicts) or the reason it was missed.

| recipient | donor level (donor, GTDB dist) | *fliE* (67) | *flhB* (297) | *flhA* (556) |
|---|---|---|---|---|
| *C. jejuni* | very close (*C. coli*, 0.037) | missed, not with donor | missed, not with donor | missed, not with donor |
| | within family (*C. fetus*, 0.247) | DET 1/11 | DET 2/14 | missed, UFBoot 85 |
| | sister family (*S. deleyianum*, 0.332) | DET 1/10 | DET 2/17 | DET 1/17 |
| | distant (*S. denitrificans*, 0.543) | missed, UFBoot 34 | DET 1/14 | DET 1/17 |
| *H. felis* | very close (*H. salomonis*, 0.055) | missed, UFBoot 88 | control flagged | missed, no topological conflict |
| | within family (*H. pullorum*, 0.447) | missed, UFBoot 93 | DET 1/14 | DET 1/15 |
| | sister family (*S. deleyianum*, 0.590) | DET 1/9 | DET 1/15 | DET 1/17 |
| | distant (*S. denitrificans*, 0.672) | missed, masked (UFBoot 99) | DET 1/14 | DET 1/17 |
| *A. butzleri* | very close (*A. vitoriensis*, 0.061) | DET 2/10 | missed, UFBoot 47 | DET 1/16 |
| | within family (*M. marinus*, 0.158) | missed, UFBoot 81 | DET 3/14 | DET 1/18 |
| | sister family (*N. alvini*, 0.391) | DET 1/14 | DET 1/13 | DET 1/15 |
| | distant (*C. fetus*, 0.569) | missed, UFBoot 55 | DET 1/14 | missed, UFBoot 89 |

| | very close | within family | sister family | distant | all |
|---|---|---|---|---|---|
| detected | 2/9 | 6/9 | 9/9 | 5/9 | 22/36 |

| gene (sites) | detected | detected, excluding very close |
|---|---|---|
| *fliE* (67) | 5/12 | 4/9 |
| *flhB* (297) | 9/12 | 9/9 |
| *flhA* (556) | 8/12 | 7/9 |

**Why transfers were missed (14):**

- **Very close donors (7).** For *C. jejuni* ← *C. coli*, the planted copy
  stays inside *Campylobacter* next to other *C. jejuni*-group species and
  never forms a group with *C. coli*'s relatives. For *H. felis* ←
  *H. salomonis* in *flhA*, landing with the donor's relatives is not a
  conflict once *H. salomonis* is removed, so there is nothing to detect. The
  *flhB* *H. felis* control was already flagged: a pre-existing
  *H. felis*/*H. labacensis* conflict, the same one seen in the dry run.
- **Placed correctly but with UFBoot < 95 (6).** Four of these are in the
  67-site *fliE*. The gene tree puts the copy in the right place but without
  enough support to pass the screen.
- **Masked by a nested conflict (1).** In *fliE* *H. felis* ← *Sulfurimonas*,
  the copy groups with 13 *Sulfurimonas* at UFBoot 99. But the screen keeps
  only minimal conflicting groups, and a smaller unrelated conflict inside
  *Sulfurimonas* (UFBoot 96) is nested within it, so the outer group was
  dropped. This is a limitation of the minimal-group rule, not of the tree.

**Interpretation.**

- **Above the very-close level, length decides.** Medium and long genes
  detected 16/18 transfers from within-family or more distant donors (all 9
  for *flhB*). The 67-site *fliE* detected 4/9.
- **Distant is lower than sister family only because of short genes.** Three
  of the four distant misses are in *fliE*, and the fourth is *flhA*
  *A. butzleri* (UFBoot 89). It is not a trend with distance.
- **Very close transfers (between sister species) are mostly not detectable**
  by single-gene topology. That sets the lower resolution limit of the
  "predominantly vertical" conclusion.
- **False-positive rate:** 1/36 controls. That hit was a conflict already in
  the real tree, not one created by removing the donor.
- **Compared with the dry run (24/36):** the counts per level shifted,
  because real re-inference changes support, but the pattern is the same.

## Run 1 (2026-10-15): *C. jejuni* only, 9 planted transfers

These results are kept in `run1/` (its `design.tsv`, `results.tsv`, screen
output and figures). They used the earlier labels: "close genus" =
*H. pylori*, "distant genus" = *S. denitrificans*. In GTDB distance,
*H. pylori* is actually the further of the two (0.672 vs 0.543), which is why
run 2 uses tree position plus distance instead.

**Detected 9/9, 0/9 controls flagged.** The planted transfer was the deepest
conflict in its tree in 7/9.

| gene (sites) | *C. lari* (0.116) | *H. pylori* (0.672) | *S. denitrificans* (0.543) |
|---|---|---|---|
| *fliE* (67) | detected, rank 1/10 | detected, rank 1/11 | detected, rank 1/13 |
| *flhB* (297) | detected, rank 4/15 | detected, rank 1/15 | detected, rank 2/17 |
| *flhA* (556) | detected, rank 1/18 | detected, rank 1/17 | detected, rank 1/16 |
