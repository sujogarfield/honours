#!/usr/bin/env python3

from Bio import AlignIO
import matplotlib
import matplotlib.pyplot as plt
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from pathlib import Path
import sys

matplotlib.use("Agg") 

if len(sys.argv) == 1:
    print("please enter one existing gene file in gene_alignments directory to generate a phylogenetic tree")
    exit(1)

input_fasta = sys.argv[1]
file_path = Path("gene_alignments") / f"{input_fasta}"

if not file_path.exists():
    print(f"{file_path} not found")
    sys.exit(1)

parts = input_fasta.strip().split("_")
gene = parts[0].strip()

tree_file = f"{gene}_tree.nwk"
image_file = f"{gene}_tree.png"

tree_dir = Path("gene_phylo_trees")
tree_dir.mkdir(exist_ok=True)
tree_path = tree_dir / tree_file
image_path = tree_dir / image_file

alignment = AlignIO.read(file_path, "fasta")
print(f"Alignment length: {alignment.get_alignment_length()}")

calculator = DistanceCalculator("identity")
constructor = DistanceTreeConstructor(calculator, "nj")
tree = constructor.build_tree(alignment)

Phylo.write(tree, tree_path, "newick")

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(1, 1, 1)

for clade in tree.get_nonterminals():
    clade.name = None

Phylo.draw(tree, axes=ax)
print("Number of clades:", len(tree.get_terminals()))
plt.savefig(image_path, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Tree saved as {tree_path}")
print(f"Image saved as {image_path}")
