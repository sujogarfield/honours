# scripts/

Run every script **from the project root**, e.g.

```
.venv/bin/python scripts/04_synteny/synteny_analysis.py
```

All data paths inside the scripts (`gene_order/...`, `campy_fetched/...`,
`images_and_newick/...`, `frontend/...`) are relative to the project root.

| Folder | Step | Scripts |
|---|---|---|
| `common/` | shared imports | `helpers.py`, `gene_colors.py` |
| `01_fetch/` | GTDB accession list → NCBI genomes, GFFs, proteins | `extract_tsv.py`, `fetch.py`, `fetch_proteins.py`, `list_genomes.py` |
| `02_species_tree/` | prune GTDB bac120 tree to the 166 genomes | `gtdb_tree.py` |
| `03_gene_order/` | string-match flagellar genes from GFFs, static tree overlay | `gene_order.py`, `gene_order_overlay.py` |
| `04_synteny/` | gene order vs. phylogeny (Mantel, agreement), full-contig sensitivity check | `synteny_analysis.py`, `synteny_analysis_full_contigs.py`, `build_merged_synteny.py` |
| `05_ragtag/` | RagTag scaffolding of split assemblies, lift coordinates | `check_ragtag_feasibility.py`, `run_ragtag_scaffold.py`, `lift_ragtag_coordinates.py` |
| `06_orthology/` | OrthoFinder crossref + eggNOG confirmation → combined gene set | `orthology_crossref*.py`, `orthology_eggnog_confirm*.py`, `build_combined_gene_order*.py`, `orthology_overlay_preview.py` |
| `07_combined_synteny/` | synteny analysis on the combined gene set | `orthology_synteny_preview.py`, `build_merged_synteny_combined.py` |
| `08_motility/` | BacDive + Madin motility traits, gene content vs. motility | `fetch_motility.py`, `fetch_madin_traits.py`, `motility_clustering.py` |
| `09_hgt/` | gene-tree vs. species-tree tanglegrams, concatenated tree | `build_gene_tanglegram_data.py`, `export_tanglegram_newick.py`, `build_concat_tree.py` |
| `10_visualisation/` | interactive HTML outputs in `frontend/` | `interactive_tree.py`, `gene_order_dashboard.py` |
| `11_flagellar_trees/` | per-gene + concatenated IQ-TREE flagellar trees on Katana, gCF/sCF, comparison with GTDB, HGT candidate screen, AU tests, independent check of tree-input gene calls | `prepare_gene_sets.py`, `katana/submit.sh`, `compare_to_gtdb.py`, `screen_hgt_candidates.py`, `plot_results.py`, `prepare_au_tests.py`, `katana/submit_au.sh`, `summarise_au.py`, `verify_gene_calls.py` |
| `side_single_gene_trees/` | early single-gene phylogenies (side branch) | `gene_search.py`, `gene_extract.py`, `gene_aligner.py`, `etetest.py` |

Scripts that import from another folder (`common/`, or `06_orthology/` for
`build_concat_tree.py`) add it to `sys.path` at the top of the file.
