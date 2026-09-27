#!/usr/bin/env python3 

import os
from Bio import SeqIO
from ete3 import Tree, TextFace, RectFace, faces
from ete3.treeview import TreeStyle, NodeStyle
import subprocess
from pathlib import Path
from collections import defaultdict
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from helpers import get_organism_from_fasta, validate_pruned_tree_integrity, colour_tree_by_species, validate_bootstrap_support

results = []

with open("campy_fetched/campy_refseq_accessions.txt") as f:
    gcf_accessions = [line.strip() for line in f if line.strip()]

# Pinned to GTDB R232 (released 15 Apr 2026). The committed bac120.tree is
# byte-identical to this file (md5 a8b589ffabbf41ebd7e8ab45795b8e6a). Don't use
# releases/latest: it would silently change the reference tree on a rerun.
GTDB_TREE_URL = "https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_r232.tree"
if not os.path.exists("scripts/02_species_tree/bac120.tree"):
    subprocess.run(["wget", "-O", "scripts/02_species_tree/bac120.tree", GTDB_TREE_URL], check=True)

tree = Tree("scripts/02_species_tree/bac120.tree", format=1, quoted_node_names=True)

target_leaves = []
found = []
not_found = []
for leaf in tree.iter_leaves():
    for acc in gcf_accessions:
        if acc in leaf.name:
            target_leaves.append(leaf.name)
            found.append(acc)
            break
        else:
            not_found.append(acc)

print(f"{len(found)} found out of {len(gcf_accessions)}")

output_dir = Path("images_and_newick/thesis_b")
if len(target_leaves) > 1:
    pruned_tree = tree.copy()
    pruned_tree.prune(target_leaves, preserve_branch_length=True)
    og_pruned = pruned_tree.copy()
    fasta_dir = Path("campy_fetched/campy_fas")

    for leaf in pruned_tree.get_leaves():
        for acc in gcf_accessions:
            if acc in leaf.name:
                matching_files = list(fasta_dir.glob(f"{acc}*.fna"))
                if matching_files:
                    fasta_file = matching_files[0]

                    if fasta_file.exists():
                        organism_name = get_organism_from_fasta(fasta_file)
                        leaf.name = organism_name
                        print(f"accession: {acc}, organism: {organism_name}")
                    else:
                        leaf.name = acc
                        print(f"accession: {acc}, organism: couldn't retrieve")
                    break

pruned_tree.write(outfile=str(output_dir / "gtdb_ref_tree.nwk"))
og_pruned.write(outfile=str(output_dir / "gtdb_ref_tree_og.nwk"))

tree = Tree(str(output_dir / "gtdb_ref_tree.nwk"), format=1, quoted_node_names=True)
og_tree = Tree(str(output_dir / "gtdb_ref_tree_og.nwk"), format=1, quoted_node_names=True)

ts = TreeStyle()

# ts.show_leaf_name = True

ts.show_leaf_name = False

for leaf in tree.iter_leaves():
    face = TextFace(leaf.name, fstyle="italic")
    leaf.add_face(face, column=0, position="branch-right")

for node in tree.traverse():
    node.dist = 1

nstyle = NodeStyle()
nstyle["size"] = 0

for node in tree.traverse():
    node.set_style(nstyle)

colour_tree_by_species(tree)

tree.render(str(output_dir / "gtdb_ref_tree.png"), tree_style=ts)

for leaf in og_tree.iter_leaves():
    face = TextFace(leaf.name, fstyle="italic")
    leaf.add_face(face, column=0, position="branch-right")

for node in og_tree.traverse():
    node.dist = 1

nstyle = NodeStyle()
nstyle["size"] = 0

for node in og_tree.traverse():
    node.set_style(nstyle)

og_tree.render(str(output_dir / "gtdb_ref_tree_og.png"), tree_style=ts)

validate_pruned_tree_integrity()
validate_bootstrap_support(
    small_nwk="images_and_newick/thesis_b/gtdb_ref_tree_og.nwk",
    output_png="images_and_newick/thesis_b/gtdb_ref_tree_support.png"
)
