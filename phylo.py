#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import sys, subprocess, os, re
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from pathlib import Path

matplotlib.use("Agg") 

if len(sys.argv) == 1:
    print("please enter one existing gene file to generate phylo tree")
    exit(1)

input_fasta = sys.argv[1]

if not os.path.exists(input_fasta):
    print(f"{input_fasta} not found")

pattern = r"^[A-Za-z0-9]+_extracted\.fna$"

if not re.match(pattern, input_fasta):
    print(f"usage: python3 phylo.py (gene)_extracted.fna")
    sys.exit(1)

parts = input_fasta.strip().split("_")
gene = parts[0].strip()

aligned_fasta = "aligned.fna"
tree_file = f"{gene}_tree.nwk"
image_file = f"{gene}_tree.png"

tree_dir = Path("gene_phylo_trees")
tree_dir.mkdir(exist_ok=True)
tree_path = tree_dir / tree_file
image_path = tree_dir / image_file

subprocess.run(
    ["muscle", "-align", input_fasta, "-output", aligned_fasta],
    check=True
)

alignment = AlignIO.read(aligned_fasta, "fasta")
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

if os.path.exists(aligned_fasta):
    os.remove(aligned_fasta)

print(f"Tree saved as {tree_path}")
print(f"Image saved as {image_path}")
