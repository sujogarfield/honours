# Planted-transfer positive control for the HGT screen

## Why

The flagellar gene trees mostly agree with the GTDB species tree, which supports
"predominantly vertical inheritance". That conclusion only means something if the
pipeline *would* have found a transfer had one happened. Here, transfers are
planted on purpose, where the right answer is known, and the unmodified
pipeline is run to see how often it finds them.

## Design

Every test is built from copies of a real gene file in
`gene_order/flagellar_trees_ext/genes/`. Those source files are only read; the
prepare script checks their md5 sums before and after.

**Recipient (all tests):** *Campylobacter jejuni* (GCF_001457695.1)

**Genes:** one short, one medium and one long. Signal scales with length; the
trimmed alignment lengths come from the 42-gene run.

| gene | length class | trimmed sites |
|---|---|---|
| *fliE* | short | 67 |
| *flhB* | medium (≈ median of the 42) | 297 |
| *flhA* | long | 556 |

**Donors:** three transfer distances.

| transfer distance | donor | accession |
|---|---|---|
| within genus | *Campylobacter lari* (other half of *Campylobacter*) | GCF_013372185.1 |
| close genus | *Helicobacter pylori* | GCF_900478295.1 |
| distant genus | *Sulfurimonas denitrificans* | GCF_000012965.1 |

**Conditions:** each (gene, donor) pair gets two files.

- **spike**: the recipient's protein sequence is replaced by the donor's
  sequence. The header is kept, so *C. jejuni* now "carries the donor's gene".
  The donor itself is removed from the file. Otherwise the planted copy would be
  identical to a sampled sequence and trivially sit next to it. With the donor
  gone, the copy has to group with the donor's relatives. That is what a gene
  acquired from that lineage some time ago looks like.
- **control**: the donor is removed and the recipient is left untouched. This
  separates the effect of the transfer from the effect of dropping one genome,
  and catches pre-existing conflicts around *C. jejuni* that would otherwise
  look like detections.

That gives 3 genes × 3 distances × 2 conditions = **18 test genes**.
`design.tsv` lists every test.

## Pipeline (identical to the real run)

1. **Gene trees, on Katana.** This uses `katana/01_gene_trees.pbs` unchanged:
   MAFFT L-INS-i, then trimAl `-automated1`, then IQ-TREE ModelFinder with
   1000 UFBoot.
2. **HGT screen.** This uses `screen_hgt_candidates.py` unchanged. A candidate is
   a gene-tree branch with UFBoot ≥ 95 that is incompatible with a GTDB branch,
   keeping only the minimal conflicting groups.
3. **Scoring.** `summarise_hgt_spikein.py`, for each tree:
   - **placement**: where *C. jejuni* lands. Its sister group is the smaller
     subtree it joins in the unrooted gene tree. Placement is `donor` if every
     sister is in the donor lineage (the GTDB clade on the donor's side of the
     recipient–donor split), `home` if every sister is in its own lineage, and
     `elsewhere` otherwise. The UFBoot of that branch is also reported.
   - **screen hit**: the screen reports a candidate group made of *C. jejuni*
     plus only donor-lineage genomes.
   - **detected**: the spike tree has a screen hit and its control does not.

## Expected outcome and interpretation

- Distant and close-genus transfers in medium and long genes should be
  detected.
- Within-genus transfers in the short gene are the hardest case and may be
  missed.
- The detection rate per length × distance cell sets the resolution limit of
  the "predominantly vertical" conclusion. For example: "transfers between
  genera are detected in genes ≥ ~300 sites, so their absence in the real data
  is informative; transfers within *Campylobacter* in short genes could go
  unseen."
- Controls flagged by the screen would indicate false positives around the
  recipient. The answer is no in a dry run on the real 42-gene trees.

## Limits

- **One recipient.** The recipient's own branch length and position affect
  power, so this is a spot check of the method, not a full power curve. More
  recipients or AliSim-simulated replicates would extend it.
- **Whole-gene replacement.** The planted transfer replaces the whole gene.
  Partial (recombined) transfers would be harder to detect.
- **Sampled donor lineage.** The donor's relatives are in the dataset. A
  transfer from an unsampled lineage would show up as a long, poorly placed
  branch instead.

## How to run

```
# local (already done; the files in genes/ are committed)
python scripts/11_flagellar_trees/prepare_hgt_spikein.py
# Katana, from the repo root, after git pull (18 array elements)
qsub -v FLAG_DIR=gene_order/hgt_spikein -J 1-18 scripts/11_flagellar_trees/katana/01_gene_trees.pbs
# local, after rsync of gene_order/hgt_spikein/gene_trees/ back
FLAG_DIR=gene_order/hgt_spikein python scripts/11_flagellar_trees/screen_hgt_candidates.py
python scripts/11_flagellar_trees/summarise_hgt_spikein.py   # writes results.tsv
```

`aligned/`, `trimmed/`, `gene_trees/` and the screen outputs are gitignored.
`results.tsv` is kept.
