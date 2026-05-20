#!/usr/bin/env python3 

from ete3 import Tree, TextFace
from ete3.treeview import TreeStyle, NodeStyle
import subprocess
from pathlib import Path
from helpers import get_organism_from_fasta

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
            target_leaves.append(leaf)
            found.append(acc)
            break
        else:
            not_found.append(acc)

print(f"{len(found)} found out of {len(gcf_accessions)}")

output_dir = Path("phylo_trees")
if len(target_leaves) > 1:
    ancestor = tree.get_common_ancestor(target_leaves)
    leaves_to_keep = set(target_leaves)
    for leaf in ancestor.get_leaves():
        if leaf not in leaves_to_keep:
            leaf.delete()
    
    fasta_dir = Path("campy_fas")
    for leaf in ancestor.get_leaves():
        for acc in gcf_accessions:
            if acc in leaf.name:
                matching_files = list(fasta_dir.glob(f"{acc}*.fna"))
                if matching_files:
                    fasta_file = matching_files[0]
                
                if fasta_file.exists():
                    organism_name = get_organism_from_fasta(fasta_file)
                    leaf.name = organism_name
                    print(f"Renamed {acc} -> {organism_name}")
                else:
                    leaf.name = acc
                    print("keep")
                break

    ancestor.write(outfile=str(output_dir / "gtdb_ref_tree.nwk"))

    tree = Tree("phylo_trees/gtdb_ref_tree.nwk", format=1, quoted_node_names=True)

    ts = TreeStyle()

    # ts.show_leaf_name = True

    ts.show_leaf_name = False

    for leaf in tree.iter_leaves():
        face = TextFace(
            leaf.name,
            fstyle="italic"
        )
        leaf.add_face(face, column=0, position="branch-right")
    
    for node in tree.traverse():
        node.dist = 1

    nstyle = NodeStyle()
    nstyle["size"] = 0

    for node in tree.traverse():
        node.set_style(nstyle)

    tree.render(str(output_dir / "gtdb_ref_tree.png"), tree_style=ts)
