#!/usr/bin/env python3

from pathlib import Path
import sys

align_dir = "gene_alignments"

if len(sys.argv) == 1:
    print("usage: python2 align_cleaner.py [name of aligned file]")
    sys.exit(1)

file = sys.argv[1]
path = Path(align_dir) / file

if not path.exists():
    print(f"path {path} not found")
    sys.exit(1)

parts = file.split('.')
cleaned = parts[0].strip() + "_cleaned.fna"
output = Path(align_dir) / cleaned


with open(path) as infile, open(output, "w") as outfile:
    for line in infile:
        if line.startswith(">"):
            line = line.rstrip("\n")
            line = line.split(",")[0] + "\n"
        outfile.write(line)
