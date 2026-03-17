#!/usr/bin/env python3

import sys, subprocess, os, re, sys
from pathlib import Path

if len(sys.argv) == 1:
    print("please enter one existing gene file in gene_extractions directory to generate alignment fasta file")
    sys.exit(1)

input_fasta = sys.argv[1]
file_path = Path("gene_extractions") / f"{input_fasta}"

if not file_path.exists():
    print(f"{file_path} not found")
    sys.exit(1)

pattern = r"^[A-Za-z0-9]+_extracted\.fna$"

if not re.match(pattern, input_fasta):
    print(f"usage: python3 phylo.py (gene)_extracted.fna")
    sys.exit(1)

parts = input_fasta.strip().split("_")
gene = parts[0].strip()

align_dir = Path("gene_alignments")
align_dir.mkdir(exist_ok=True)

aligned_fasta = align_dir / f"{gene}_aligned.fna"
tree_file = f"{gene}_tree.nwk"
image_file = f"{gene}_tree.png"

tree_dir = Path("gene_phylo_trees")
tree_dir.mkdir(exist_ok=True)
tree_path = tree_dir / tree_file
image_path = tree_dir / image_file

subprocess.run(
    ["muscle", "-align", str(file_path), "-output", str(aligned_fasta)],
    check=True
)
