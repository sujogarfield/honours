#!/usr/bin/env python3

# Script that searches for gene presence in each of the
# Campylobacter species listed in campy_ann. If gene is present
# it returns start and end location and strand of gene.

import os, gzip

genes = ["fliI", "fliF", "fliC"]

for file in os.listdir("campy_ann"):
    if not file.endswith(".gff.gz"):
        continue

    print("-" * 70)
    print(f"\033[95mFile: {file}\033[0m")

    filepath = os.path.join("campy_ann", file)

    counter = 0
    with gzip.open(filepath, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")

            if len(parts) > 8 and parts[2] == "gene":
                attributes = parts[8]

                for gene in genes:
                    if gene.lower() in attributes.lower():
                        counter = 1
                        contig = parts[0]
                        start = parts[3]
                        end = parts[4]
                        strand = parts[6]

                        print(f"\nGene: {gene}")
                        print(f"Contig: {contig}")
                        print(f"Start: {start}")
                        print(f"End: {end}")
                        print(f"Strand: {strand}")
    if counter == 0:
        print("No matches!")
    print(f"\nSource: https://www.ncbi.nlm.nih.gov/nuccore/{contig}")
    print("-" * 70)
