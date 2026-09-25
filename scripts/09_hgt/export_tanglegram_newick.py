#!/usr/bin/env python3
"""
export_tanglegram_newick.py
Plain-Newick export of the gene-tree/species-tree pair behind each HGT
tanglegram, so it's reproducible outside this project's own D3 dashboard --
e.g. in R with ape::cophyloplot() or phytools::cophylo(), which is what the
dashboard's tanglegram tab visually reimplements. Exports the same pairing
logic as build_gene_tanglegram_data.py, just as three flat files per
orthogroup instead of one nested JSON for D3.

Both trees' tip labels are set to the SAME identifiers (each gene-tree
leaf's own name, e.g. "GCF_000009085_1_WP_..." ; species-tree leaves with
paralogs are duplicated once per copy, one duplicate per gene-tree leaf,
under a zero-branch-length spacer node) specifically so R's default
label-matching (cophylo(tr1, tr2), no assoc table) works out of the box.
assoc.csv is also written as an explicit, robust fallback (two columns:
species-tree tip, gene-tree tip) in case any tip label needs re-quoting for
a given R Newick parser.

Usage: python3 export_tanglegram_newick.py <orthogroup_id> [<orthogroup_id> ...]
Output (per orthogroup, under gene_order/tanglegram_r_export/<og_id>/):
  gene_tree.nwk      -- OrthoFinder's resolved gene tree, original branch lengths
  species_tree.nwk   -- GTDB species tree, pruned to this orthogroup's genomes
                        and duplicated per paralog copy (zero-length spacer
                        nodes), tip labels matching gene_tree.nwk
  assoc.csv           -- species_tip,gene_tip pairing table (fallback for R)
  leaf_labels.csv      -- gene_tip,organism,protein_id (human-readable key)
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

from ete3 import Tree

from build_gene_tanglegram_data import load_gcf_prefix_map, load_orthogroup_tree

OUTPUT_DIR = Path("gene_order/tanglegram_r_export")
NAMED_NWK = "images_and_newick/thesis_b/gtdb_ref_tree.nwk"


def build_one(og_id: str, prefix_map: dict, species_tree: Tree):
    newick = load_orthogroup_tree(og_id)
    gtree = Tree(newick, format=1)

    leaf_organism_of = {}
    protein_of = {}
    unmatched = []
    sorted_prefixes = sorted(prefix_map, key=len, reverse=True)
    for leaf in gtree.get_leaves():
        matched = None
        for prefix in sorted_prefixes:
            if leaf.name.startswith(prefix + "_"):
                matched = prefix
                break
        if matched is None:
            unmatched.append(leaf.name)
            continue
        gcf_id, organism = prefix_map[matched]
        leaf_organism_of[leaf.name] = organism
        protein_of[leaf.name] = leaf.name[len(matched) + 1:]

    if unmatched:
        print(f"  {og_id}: {len(unmatched)} unmatched gene-tree leaves dropped (not one of our 166 genomes)")
        gtree.prune([leaf.name for leaf in gtree.get_leaves() if leaf.name not in unmatched],
                    preserve_branch_length=True)

    gene_leaves_by_organism = defaultdict(list)
    for gene_leaf_name, organism in leaf_organism_of.items():
        gene_leaves_by_organism[organism].append(gene_leaf_name)
    organisms_present = set(gene_leaves_by_organism)

    pruned = species_tree.copy()
    pruned.prune([leaf.name for leaf in pruned.get_leaves() if leaf.name in organisms_present],
                 preserve_branch_length=True)

    assoc_rows = []
    for leaf in pruned.get_leaves():
        copies = gene_leaves_by_organism[leaf.name]
        if len(copies) == 1:
            assoc_rows.append((leaf.name, copies[0]))
            leaf.name = copies[0]  # matching tip label for R's default cophylo() matching
        else:
            # duplicate under a zero-length spacer, one child per paralog copy
            leaf.dist = 0.0
            for i, gene_leaf_name in enumerate(copies):
                child = leaf.add_child(dist=0.0)
                child.name = gene_leaf_name
                assoc_rows.append((gene_leaf_name, gene_leaf_name))
            leaf.name = ""  # now an internal spacer, not a tip

    og_dir = OUTPUT_DIR / og_id
    og_dir.mkdir(parents=True, exist_ok=True)
    gtree.write(outfile=str(og_dir / "gene_tree.nwk"), format=1)
    pruned.write(outfile=str(og_dir / "species_tree.nwk"), format=1)

    with open(og_dir / "assoc.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["species_tip", "gene_tip"])
        w.writerows(assoc_rows)

    with open(og_dir / "leaf_labels.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene_tip", "organism", "protein_id"])
        for leaf_name, organism in leaf_organism_of.items():
            w.writerow([leaf_name, organism, protein_of[leaf_name]])

    print(f"  {og_id}: wrote {len(gtree.get_leaves())} gene-tree tips, "
          f"{len(pruned.get_leaves())} species-tree tips -> {og_dir}/")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    og_ids = sys.argv[1:]

    print("Loading GCF->organism map ...")
    prefix_map = load_gcf_prefix_map()
    print(f"Loading species tree from {NAMED_NWK} ...")
    species_tree = Tree(NAMED_NWK, format=1, quoted_node_names=True)

    for og_id in og_ids:
        print(f"Exporting {og_id} ...")
        build_one(og_id, prefix_map, species_tree)

    print(f"\nDone. See {OUTPUT_DIR}/<orthogroup>/ for gene_tree.nwk, species_tree.nwk, assoc.csv, leaf_labels.csv")
    print("R usage (phytools):")
    print('  library(phytools)')
    print('  gt <- read.tree("gene_tree.nwk"); st <- read.tree("species_tree.nwk")')
    print('  obj <- cophylo(st, gt)   # matches by identical tip labels automatically')
    print('  plot(obj)')
    print("R usage (ape, with explicit assoc table):")
    print('  library(ape)')
    print('  assoc <- as.matrix(read.csv("assoc.csv"))')
    print('  cophyloplot(read.tree("species_tree.nwk"), read.tree("gene_tree.nwk"), assoc=assoc, space=60)')


if __name__ == "__main__":
    main()
