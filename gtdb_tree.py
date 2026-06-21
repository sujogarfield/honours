#!/usr/bin/env python3 

import os
from Bio import SeqIO
from ete3 import Tree, TextFace, RectFace, faces
from ete3.treeview import TreeStyle, NodeStyle
import subprocess
from pathlib import Path
from collections import defaultdict
from helpers import get_organism_from_fasta, validate_pruned_tree_integrity, colour_tree_by_species, validate_bootstrap_support

results = []

with open("campy_refseq_accessions.txt") as f:
    gcf_accessions = [line.strip() for line in f if line.strip()]

subprocess.run([
    "wget",
    "-nc",
    "https://data.gtdb.ecogenomic.org/releases/latest/bac120.tree"
])

tree = Tree("bac120.tree", format=1, quoted_node_names=True)

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

output_dir = Path("phylo_trees")
if len(target_leaves) > 1:
    pruned_tree = tree.copy()
    pruned_tree.prune(target_leaves, preserve_branch_length=True)
    og_pruned = pruned_tree.copy()
    fasta_dir = Path("campy_fas")

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
    small_nwk="phylo_trees/gtdb_ref_tree_og.nwk",
    output_png="phylo_trees/gtdb_ref_tree_support.png"
)
