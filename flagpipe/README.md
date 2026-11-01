# flagpipe

The thesis flagellar-gene pipeline for **any GTDB taxon and gene list**. It
uses the same logic as the *Campylobacterota* analysis:

```
GTDB TSV ─► genomes ─► string match (gene= tags)
                         │
                         ├─► OrthoFinder orthogroups ◄── internal seeds (string-matched proteins)
                         │                          ◄── external seeds (reference genome, DIAMOND)
                         ├─► eggNOG-mapper confirmation ─► combined gene set
                         ├─► copy rules ─► per-gene trees + concatenated tree (IQ-TREE 2)
                         │                 ├─► gCF / sCF, comparison with GTDB
                         │                 └─► HGT screen ─► shortlist rule ─► AU tests
                         └─► synteny distance vs GTDB distance (Mantel)
Outputs: interactive tree, synteny scatter, tanglegram
```

The thesis scripts in `scripts/` are not changed. flagpipe reuses
`scripts/common/` (synteny distance, reliability, gene colours),
`scripts/10_visualisation/interactive_tree.py` and the tree-drawing helpers
from `scripts/11_flagellar_trees/plot_results.py`.

---

## On Katana: one-off setup (once per group)

```bash
cd /srv/scratch/$USER
git clone https://github.com/sujogarfield/honours.git
cd honours

# 1. conda environment (OrthoFinder, eggNOG-mapper 2.1.13, DIAMOND, MAFFT, trimAl, IQ-TREE 2.4.0)
qsub flagpipe/katana/setup_env.pbs              # -> /srv/scratch/baker/flagpipe_env

# 2. eggNOG database, ~50 GB, shared by everyone
qsub -v EGGNOG_DATA_DIR=/srv/scratch/baker/eggnog_data flagpipe/katana/download_eggnog.pbs
```

The defaults are `FLAGPIPE_ENV=/srv/scratch/baker/flagpipe_env` and
`CONDA_BASE=/srv/scratch/baker/miniforge3`, set in `katana/common.sh`.
Anyone in the group with read access to those folders can skip both steps.

If Katana compute nodes have no internet, the eggNOG download job will fail.
In that case run `download_eggnog_data.py -y --data_dir <dir>` on the login
node inside `screen`.

## A run

Every step takes a **run directory**. Put it in scratch, e.g.
`/srv/scratch/$USER/flagpipe_runs/listeriaceae`.

```bash
export PATH=/srv/scratch/baker/flagpipe_env/bin:$PATH
RUN=/srv/scratch/$USER/flagpipe_runs/listeriaceae

# Login node (light work, needs internet)
python flagpipe/flagpipe.py init    $RUN   # questions -> config.json
python flagpipe/flagpipe.py genomes $RUN   # filter + cap the TSV, pick the outgroup
python flagpipe/flagpipe.py fetch   $RUN   # NCBI FASTA + GFF + proteome
python flagpipe/flagpipe.py tree    $RUN   # GTDB tree pruned to the genomes
python flagpipe/flagpipe.py check   $RUN   # which gene names are found -> edit config.json "genes"

# Katana jobs
bash flagpipe/katana/submit.sh $RUN orthology  # OrthoFinder, eggNOG, gene sets, synteny, report
#   ...check $RUN/trees/gene_selection.json and $RUN/report/ ...
bash flagpipe/katana/submit.sh $RUN trees      # gene trees, concatenated tree, gCF/sCF, HGT screen, tanglegram
#   ...check $RUN/trees/hgt_shortlist.tsv ...
bash flagpipe/katana/submit.sh $RUN hgt        # AU tests on the shortlisted HGT candidates
```

**Getting the GTDB TSV:** on gtdb.ecogenomic.org, open Advanced Search, search
for the taxon (e.g. GTDB Taxonomy contains `f__Listeriaceae`) and download the
results as TSV. It needs the default columns: Accession, Organism Name, GTDB
Taxonomy, GTDB Representative of Species, GTDB Type Material.

**Choosing genes (`init` + `check`):**
- `init` offers the presets in `presets/`:
  - `flagellar_motility`: the thesis list;
  - `flagellar_structural`: core flagellar genes found in most flagellated bacteria.
- You can also give your own file or type names, then add or remove names.
- `check` runs after `fetch`. It shows how many genomes each name is found in.
  It also lists flagellar-looking gene names in the GFFs that are *not* on your
  list (on *Campylobacterota* it finds *fliI* in 148 genomes, which the thesis
  list missed). Edit `"genes"` in `config.json` and rerun `check` until the list
  looks right.

**The seed reference genome:**
- Pick a well-annotated genome inside the taxon, e.g. *B. subtilis* 168
  (`GCF_000009045.1`) for Bacillota, or *C. jejuni* NCTC 11168
  (`GCF_000009085.1`) for Campylobacterota.
- Its proteins whose product contains a keyword (flagell, chemotax, ...) seed
  extra orthogroups. That is how genes with no gene= tag, such as unnamed
  flagellins, are recovered.
- **Choose a reference inside the family if you can.** DIAMOND runs in
  `--ultra-sensitive` mode, but short, fast-evolving proteins (FliJ/FliK/FliS/FliT
  type, regulators) still won't hit across genera. In the *Listeria* test,
  *B. subtilis* 168 seeded 30 orthogroups, but 16 of its 59 flagellar/chemotaxis
  proteins had no hit. A reference from another phylum will find few orthogroups.
- **Expect to rely on the reference outside Campylobacterota.** RefSeq gives many
  flagellar genes no gene= tag in some taxa. In 4 *Listeria* genomes only 13 of
  38 core names were found by string matching (e.g. no *fliG*, *fliF*, *flhA*),
  even though the proteins are annotated "flagellar motor switch protein FliG" etc.
- **Name collisions:** *pflA*/*pflB* (pyruvate formate-lyase) and *maf* (septum
  formation) mean unrelated genes outside Campylobacterota. `check` warns if they
  are on your list. The `flagellar_motility` preset is Campylobacterota-specific;
  use `flagellar_structural` for other taxa.

### Outputs (`$RUN/report/`)

| File | After | What it shows |
|---|---|---|
| `interactive_tree.html` | stage 1 | pruned GTDB tree with each genome's gene track |
| `synteny_scatter.html` | stage 1 | synteny vs phylogenetic distance per genome pair, Mantel r; all genomes / reliable toggle |
| `tanglegram.png` | stage 2 | GTDB tree (branches shaded by gCF) vs the concatenated tree |
| `hgt_shortlist.tsv` | stage 3 | every conflicting group of genomes and whether it passed the shortlist rule (also in `trees/`) |
| `au_summary.tsv` | stage 3 | AU test per shortlisted candidate: does the flagellar alignment reject the GTDB arrangement? |

The HTML files are self-contained. Copy them to a laptop (`rsync`/FileZilla)
and open them in a browser.

Other results:
- `gene_order_combined.json`
- `orthology/candidates_confirmed.json` (eggNOG status of every candidate)
- `trees/gene_selection.json` (why each gene and genome was kept or dropped)
- `trees/comparison.json` (RF, shared branches, Mantel and gCF/sCF summaries)
- `trees/hgt_candidates.tsv` (every gene-tree branch with UFBoot ≥ 95 that contradicts GTDB)
- `trees/au/au_results.tsv` (AU result per candidate × alignment)
- `synteny/*.json`

### Config (`config.json`)

| Key | Default | Meaning |
|---|---|---|
| `taxon_filter` | — | keep TSV rows whose GTDB taxonomy contains this |
| `representatives_only` | true | one genome per GTDB species |
| `max_genomes`, `max_per_genus` | 300, 40 | caps (OrthoFinder is all-vs-all); type material is kept first |
| `outgroup` | auto | `auto` = 3 RefSeq genomes spread across the sister clade in the GTDB tree; or a list of accessions; or `none` |
| `genes` | preset | gene names to string-match |
| `exclude_genes` | [] | names never taken from orthology (string matches are kept, as in the thesis) |
| `seed_reference` | — | `{accession, keywords}` |
| `dominance` | 0.80 | eggNOG name agreement needed to label an external-seed orthogroup |
| `selection` | 0.60 / 0.10 / 0.60 | min fraction of genomes with the gene / max fraction ambiguous multi-copy / min fraction of selected genes per genome |
| `synteny` | 4 / 999 | min shared genes per pair / Mantel permutations |
| `hgt` | auto / 0.50 / 3 / 25 | shortlist rule: gene trees needed (auto = 25% of genes, 10 of 37 in the thesis) / one-sided fraction / one-sided ratio / GTDB branch-length percentile |

---

## HGT candidates (stages 2 and 3)

The thesis method (`scripts/11_flagellar_trees`, RESULTS.md section 4), automated:

1. **Screen** (`hgt-screen`, end of stage 2): `screen_hgt_candidates.py`. Every
   branch with UFBoot ≥ 95 in a gene tree or the concatenated tree that is
   incompatible with GTDB, reduced to minimal groups of genomes.
2. **Shortlist** (`hgt-shortlist`, end of stage 2). The thesis rule, which the
   thesis applied by hand, now in `fp/hgt.py`. A group passes if:
   1. it is a branch in ≥ `min_gene_trees` gene trees **and** the concatenated tree;
   2. it is one-sided: at the GTDB branch it contradicts, ≥ 50% of the gene
      trees that can resolve that branch support the group's alternative, and
      ≥ 3× as many as the other alternative;
   3. that GTDB branch is not unusually short (≥ 25th percentile).
   `trees/hgt_shortlist.tsv` lists every group with these numbers, so groups
   that just fail can be added by hand, as the thesis did for its four
   exploratory candidates.
3. **AU tests** (`submit.sh $RUN hgt`, stage 3): for each shortlisted group,
   constrained searches forcing only the GTDB branch(es) it contradicts
   (3 seeds, best kept), then AU (10,000 RELL) on the concatenated alignment
   and every informative gene (`prepare_au_tests.py`, `au_one.sh`). Summarised
   by `summarise_au.py`: Holm across candidates, BH across genes.

**Adding a candidate by hand:** add its name to `trees/au/candidates.txt` and
create `trees/au/<name>/` the same way `hgt-shortlist` does. The simplest way
is to copy the code for one candidate from `fp/hgt.py:shortlist`.

**Not included:** the species-marker test and the core-gene backbone test
(thesis M11/M12). They need the GTDB bac120 marker alignment for the genomes
and are run separately.

**Reading the result:** rejecting GTDB means the flagellar genes share a
history at that spot that differs from the species tree. That is consistent
with HGT, but also with GTDB being wrong there. The candidates are picked from
the same trees they are tested on, so the p-values are post-selection.

## Differences from the thesis run

1. **Coverage threshold as a fraction.** 0.60 of genomes. That is the thesis's
   100 of 166.
2. **Each candidate is annotated by eggNOG once**, in DIAMOND mode. The thesis
   ran eggNOG in five batches, some in HMMER mode, and the regulon batch
   annotated only proteins not seen before.
3. **External-seed naming** combines the two thesis rules:
   - the dominant eggNOG name if ≥ 80% of an orthogroup's named members agree
     (the regulon sweep);
   - otherwise, if every named member has a name on the gene list (a paralog
     family of target genes, like flaA/flaB), each protein keeps its own name
     (the flagellin script).
4. **No lineage-specific special cases.** The thesis *maf* (`maf_2356`) rule
   is not included.
5. **No RagTag.** The reliable synteny set counts only natively closed
   chromosomes.
6. **Fixed IQ-TREE seed** (`--seed 1`), so trees can be reproduced.

## Test on the thesis data (2026-10-22)

Run locally on the 166 *Campylobacterota* genomes, using the thesis OrthoFinder
results and eggNOG annotations (no Katana).

| Step | flagpipe | thesis |
|---|---|---|
| genomes after filters | 166 | 166 |
| string-matched genes | 3,895 | 3,895 |
| external seed: reference proteins → orthogroups | 50 → 47 | 50 → 47 |
| internal-seed genes kept | 1,826 | 1,826 |
| external-seed genes kept | 993 | 919 |
| genes selected for trees | 38 | 37 |
| genomes in trees | 149 (identical set) | 149 |
| gene FASTA identical to thesis | 31 of 37 | — |
| synteny, all genomes, on the thesis gene set | r = 0.471, 164 genomes, 10,749 pairs | r = 0.471, 164 genomes, 10,749 pairs |

**Why the external-seed counts differ:**
- **+145 *flgG2*:** flagpipe confirms OG0000613 as *flgG2*. The thesis never did,
  because those proteins were already annotated in an earlier batch.
- **+34 *flgM*, +33 *fliY* and a few others:** the same batch effect.
- **−166 *maf_2356*:** the thesis special case (difference 4).
- **No conflicts:** every locus found by both is given the same gene name.

Those extra copies change the selection:
- *flgG2* and *flgM* now pass;
- *fliK* now fails the multi-copy rule.

So 6 gene FASTA files differ: *cheV*, *flgD*, *flgE*, *flgL*, *fliK*, *fliY*.

**Bug found and fixed by this test:** organism names containing `=` were
altered by ete3 when written to newick, which broke their link to the tree.
Labels are now cleaned, and the `tree` step checks they survive a round trip.

## Test of the HGT stage on the thesis trees (2026-11-02)

The thesis gene trees, concatenated tree and gCF output were copied into a
scratch run (the thesis folder was not touched).

| Check | Result |
|---|---|
| `hgt-screen` output | byte-identical to `gene_order/flagellar_trees/hgt_candidates.tsv` (178 groups) |
| 25th-percentile GTDB branch length | 0.0086, as in RESULTS.md |
| groups passing the rule | exactly the 4 thesis primary candidates |
| 4 exploratory candidates | fail for the reasons in RESULTS.md: *Hydrogenimonas* branch 0.0079; *Sulfurimonas* not one-sided (14%); *Nitrosophilus* 7 genes; *Sulfurospirillum* 5 genes |
| AU constraint trees (concatenated and per gene) | byte-identical to the thesis files |
| `hgt-summary` on the thesis AU outputs | same p, ΔlogL and gene counts as the thesis `au_summary.tsv` |
| `au_one.sh` run locally on one gene (*fliK*) | completes; values differ slightly from the thesis because the laptop has IQ-TREE 3.1.3 (Katana env: 2.4.0, as the thesis) |

Not yet run: the AU array job on Katana.

## Test of the download steps (2 *Listeria* species + auto outgroup)

- `genomes` → `fetch` → `tree` → `check` → `string-match` work against NCBI.
- OrthoFinder (run locally, 54 s) and `crossref` also work.
- The auto outgroup picked 2 RefSeq genomes from the sister clade, all 4 genomes
  were in the GTDB tree, and the *B. subtilis* reference downloaded.
- This test found the DIAMOND sensitivity problem (fixed: 19 → 30 seeded
  orthogroups; the Campylobacterota result was unchanged at 47) and the name
  collisions above.
- eggNOG-mapper wasn't run, because there is no local eggNOG database.

**Not yet tested on Katana:**
- `setup_env.pbs` and `download_eggnog.pbs`;
- the PBS job chain (including the new HGT stage);
- eggNOG-mapper on a new taxon.

Run a small family first (e.g. `f__Listeriaceae`) to check these.
