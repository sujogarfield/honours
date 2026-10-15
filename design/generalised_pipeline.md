# Design: a generalised gene-order + phylogeny pipeline

**Status:** planned, not started. Begin **after the *Campylobacter* analyses are
finished and frozen.**
**Written:** 2026-10-15

## Goal

Run the whole *Campylobacter* workflow on **any GTDB taxon and any gene set**
from a single config file. The workflow covers:

- genomes;
- gene search with eggNOG confirmation;
- gene order and synteny;
- the interactive tree and synteny dashboard;
- gene trees, concordance, the HGT screen, the planted-transfer control and the
  marker test.

Example target: flagellar genes (KEGG ko02040) across *Bacillota*, or a family
within it.

## Principles

1. **Refactor first, then generalise.**
   - Phase 1 changes *where values come from*, not *what the code does*. It
     must reproduce the frozen *Campylobacter* results.
   - Only then does phase 2 add new capability.
   - Otherwise a different answer on a new taxon can't be separated from a bug
     in the new code.
2. **One config per run.** Nothing taxon- or gene-specific is left in a script.
3. **The *Campylobacter* config is the regression test**, kept for good: any
   later change must still reproduce it.
4. **Reuse what already generalises.**
   - `scripts/11_flagellar_trees/` already runs from `FLAG_DIR`.
   - `common/synteny_distance.py` and `common/reliability.py` are already
     mostly taxon-neutral.

---

## Phase 0: freeze *Campylobacter*

- **Tag the final commit:** `git tag campy-final`. All phase 1 work branches
  from this tag.
- **Record the reference outputs** in `design/campy_reference_manifest.tsv`:
  path, md5 and the key numbers for every file listed under "Regression test"
  below.
- **Set IQ-TREE seeds.** The current Katana jobs don't pass `--seed`, so their
  trees can't be reproduced byte-for-byte. Rerun the reference tree steps with
  `--seed 1` before freezing. Or accept the tolerance-based checks below for
  those steps.

---

## Phase 1: move hardcoded values into config (no behaviour change)

Work on branch `pipeline-config`, off `campy-final`. Merge only when the
regression test passes.

### Config schema (draft)

```yaml
# config/campylobacterota_flagella.yaml   -- reproduces the thesis run
name: campylobacterota_flagella
data_root: data/campylobacterota_flagella      # replaces campy_fetched/, campy_orthologs/
genomes:
  accessions: campy_fetched/campy_refseq_accessions.txt   # phase 2: or `taxon:`
gtdb:
  release: 232
  tree: scripts/02_species_tree/bac120.tree
  outgroup_genera: [Hippea, Desulfurella]                 # phase 2: `auto`
genes:
  name_patterns: config/genes/campy_flagellar_names.txt   # current string-match list
  extra_sets:                                             # current crossref sets
    flagellin: config/genes/campy_flagellin.tsv
    maf: config/genes/campy_maf.tsv
    regulon: config/genes/campy_regulon.tsv
  excluded: [pfla, yuxg]
  anchor_genes: [flim, flhf]                              # circular rotation anchor
gene_order:
  reliable: complete_or_ragtag                            # complete_only | complete_or_ragtag
  min_shared_genes: 4
  ragtag: true
trees:
  run_dir: gene_order/flagellar_trees_ext
  ufboot_min: 95
  gcf_low: 5
spikein:
  genes: auto            # short / median / long from trimmed lengths
  recipients: [Campylobacter jejuni, Helicobacter felis, Aliarcobacter butzleri]
display:
  groups: config/display/campy_groups.tsv                 # replaces helpers.py and plot_results GROUPS
```

`scripts/common/config.py` loads it, with `load_config(path)` returning a
dict. Every script takes `--config`, defaulting to the *Campylobacter* config
so that existing commands keep working.

### Inventory of hardcoded values

As of 2026-10-15, 41 scripts contain hardcoded *Campylobacter* values.

| kind | where | moves to |
|---|---|---|
| **Data paths** (`campy_fetched/…`, `campy_orthologs/…`, `gene_order/…`) | nearly every script; heaviest in `06_orthology/*` (3–10 each), `09_hgt/build_concat_tree.py`, `11_flagellar_trees/prepare_*` | `data_root` + fixed sub-layout |
| **Genome list** | `01_fetch/fetch.py`, `fetch_proteins.py`, `02_species_tree/gtdb_tree.py` (read `campy_refseq_accessions.txt`) | `genomes.accessions` |
| **Gene name list** (string match) | `03_gene_order/gene_order.py` (~50 names), `gene_order_overlay.py`, `04_synteny/*` | `genes.name_patterns` |
| **Extra gene sets** (flagellin, *maf*, regulon) | the three `06_orthology/orthology_crossref_*.py` + `orthology_eggnog_confirm_*.py` pairs, which are near-duplicates | `genes.extra_sets`; the pairs collapse to one crossref + one confirm script |
| **Reference accessions** (GCF IDs used as query genomes) | `orthology_crossref_flagellin/maf/regulon.py`, `verify_gene_calls.py`, `prepare_gene_sets.py`, `build_concat_tree.py`, `export_tanglegram_newick.py` | per-gene-set reference in the gene-set TSV |
| **Taxon names → colours/groups** | `common/helpers.py` (genus `if/elif` chain), `plot_results.py` (`GROUPS`), `gene_order_dashboard.py`, `prepare_au_tests.py` | `display.groups`; phase 2: generated from GTDB taxonomy |
| **Outgroup** | `plot_results.py`, `plot_hgt_spikein.py` (`Hippea`, `Desulfurella`) | `gtdb.outgroup_genera` |
| **Anchor genes, thresholds** | `interactive_tree.py` (fliM/flhF), `reliability.py`, `synteny_distance.py` (`MIN_SHARED_GENES`), `screen_hgt_candidates.py` (`UFBOOT_MIN`, `GCF_LOW`) | `gene_order.*`, `trees.*` |
| **Spike-in donors/recipients** | `prepare_hgt_spikein.py` (hand-picked) | `spikein.recipients`; donors already follow a GTDB-position rule, so they can be chosen automatically |
| **Side branch** | `side_single_gene_trees/*` | leave as is (not part of the pipeline) |

### Regression test

Run the full pipeline with the *Campylobacter* config and compare against the
`campy-final` manifest.

| output | comparison |
|---|---|
| `gene_order/gene_order_combined*.json` | identical (md5) |
| `synteny_analysis_combined*.json` | identical; Mantel r 0.665 (reliable, new distance) and 0.471 (all) |
| interactive tree / dashboard HTML in `frontend/`, `docs/` | identical, or a diff limited to build timestamps |
| reliable set | the same 103 genomes (90 native + 13 RagTag) |
| gene trees, concat tree | same topology (RF = 0) with seeds set; otherwise RF ≤ the run-to-run noise measured on 2 reruns |
| `concordance/*.cf.stat` | identical with seeds; otherwise per-branch gCF within ±2 |
| `hgt_candidates.tsv`, spike-in `results.tsv` | identical sets of candidate groups and detections |

Script: `scripts/common/regression_check.py --manifest design/campy_reference_manifest.tsv`.

---

## Phase 2: generalise

### Stages

| # | stage | inputs → outputs | tools | reuses |
|---|---|---|---|---|
| 1 | **Select genomes** | `taxon:` (any GTDB rank) → accession list | GTDB metadata TSV (`bac120_metadata_r232`): species representatives, then CheckM completeness ≥ 90 / contamination ≤ 5 | `01_fetch/extract_tsv.py` |
| 2 | **Download** | accessions → genome FASTA, GFF, proteins | NCBI `datasets` CLI (batch; replaces per-genome Entrez calls) | `01_fetch/fetch*.py` |
| 3 | **Species tree** | GTDB tree → pruned reference + outgroup | ete3; outgroup = sister clade of the taxon in GTDB (`outgroup: auto`) | `02_species_tree/gtdb_tree.py` |
| 4 | **Find genes** | gene set + proteomes → candidate hits | DIAMOND (`--more-sensitive`) or HMMER against reference proteins / KO profiles | replaces OrthoFinder + crossref |
| 5 | **Confirm with eggNOG** | hits → confirmed orthologs (KO / preferred name) | eggNOG-mapper on **hits only**, not whole proteomes | `06_orthology/orthology_eggnog_confirm.py` (single generic version) |
| 6 | **Gene order** | confirmed genes + GFF → `gene_order_combined.json` | — | `03_gene_order`, `06_orthology/build_combined_gene_order*.py` |
| 7 | **Reliability / RagTag** | assemblies → reliable set | RagTag (optional; needs a same-species complete reference) | `common/reliability.py`, `05_ragtag/*` |
| 8 | **Synteny** | gene order + tree → distances, Mantel, clustering | — | `common/synteny_distance.py`, `07_combined_synteny/*` |
| 9 | **Frontend** | → interactive tree + dashboard HTML | — | `10_visualisation/*` |
| 10 | **Trees and tests** | confirmed genes → gene trees, concat, gCF/sCF, GTDB comparison, HGT screen, marker test, spike-in | MAFFT, trimAl, IQ-TREE 2 | `11_flagellar_trees/*` |
| 11 | **Motility traits** (optional) | taxon names → BacDive / Madin traits | — | `08_motility/*` |

### Gene set input

Three accepted forms, in order of preference:

1. **KEGG KO list or pathway** (`kegg:ko02040`). eggNOG assigns KOs, so this
   works for any phylum without reference sequences. This is the default.
2. **FASTA of reference proteins.** For genes without a KO, or for
   lineage-specific genes (as with the *maf* and regulon sets here).
3. **Name patterns.** The current string match; kept only for the regression
   test.

### Workflow runner

**Snakemake**, with a PBS profile for Katana (`profiles/katana/`). This
replaces the hand-chained `katana/submit*.sh` scripts. Snakemake already
provides:

- job arrays;
- dependencies;
- rerunning only missing outputs;
- per-rule resources.

```
snakemake --profile profiles/katana --configfile config/bacillaceae_flagella.yaml
```

Layout: `workflow/Snakefile`, `workflow/rules/{fetch,genes,synteny,trees,tests,frontend}.smk`,
`config/`, `profiles/`.

### Scaling

*Campylobacterota* was about 166 genomes. *Bacillota* has thousands of GTDB
species.

- **Sampling:**
  - one representative per species;
  - cap per genus (`max_per_genus`);
  - total cap (`max_genomes`, about 300), stratified so every genus keeps at
    least one genome.
- **OrthoFinder:** dropped. All-vs-all doesn't scale, and a fixed gene set
  doesn't need it.
- **eggNOG:** run on hits only (step 5), not whole proteomes.
- **IQ-TREE:** the per-gene array is fine. The concatenated tree at about 300
  taxa needs more walltime (≈ 6–12 h); the concordance step is cheap.
- **Start small:** run one family first (e.g. `f__Bacillaceae` or
  `f__Listeriaceae`) before a phylum.

### Taxon-specific choices that must become rules

| choice made by hand for *Campylobacter* | rule for any taxon |
|---|---|
| Colour groups (`helpers.py`, `GROUPS`) | colour by GTDB family (or genus if ≤ 8 families) |
| Outgroup *Hippea*/*Desulfurella* | GTDB sister clade of the target taxon, 2–4 genomes |
| Anchor gene fliM for circular display | most-shared gene in the set (already how `rank_anchor_genes` works) |
| Spike-in recipients | one recipient per largest family (up to 3), with a close donor group available |
| Spike-in genes | shortest, median and longest trimmed alignment |
| Marker test | GTDB bac120 markers (already generic) |

---

## Open questions (decide at the start of phase 2)

1. **Scope:** a whole phylum, or one or two families? This decides sampling and
   compute.
2. **Gene set for *Bacillota*:** flagellar genes only (ko02040), or also
   chemotaxis (ko02030)? *Bacillota* flagellar genes are in fewer, larger
   operons than in *Campylobacterota*, which suits the synteny analysis.
3. **RagTag:** keep it as an option, or restrict to complete genomes only? With
   a larger pool of complete genomes, "complete only" may be enough.
4. **Where it runs:** all on Katana (Snakemake + PBS), or fetch/frontend
   locally and trees on Katana as now.
5. **Thesis use:** generalisation as "future work" only, or a *Bacillota*
   comparison in the thesis? This decides whether phase 2 happens before or
   after submission.

## Checklist

- [ ] Phase 0: finish *Campylobacter*; set IQ-TREE seeds; tag `campy-final`; write manifest
- [ ] Phase 1: `common/config.py` + *Campylobacter* config
- [ ] Phase 1: replace hardcoded values, script group by script group (01 → 11)
- [ ] Phase 1: collapse the three crossref/confirm script pairs into one each
- [ ] Phase 1: `regression_check.py` passes → merge `pipeline-config`
- [ ] Phase 2: Snakemake skeleton wrapping the existing steps (*Campylobacter* still passes)
- [ ] Phase 2: stages 1–2 (GTDB selection + `datasets` download)
- [ ] Phase 2: stages 4–5 (DIAMOND/HMMER search + eggNOG confirm on hits)
- [ ] Phase 2: automatic rules (colours, outgroup, spike-in choices)
- [ ] Phase 2: first new taxon (one family), then scale up
