#!/usr/bin/env python3

# Script that searches for gene presence in each of the
# Campylobacter species listed in campy_ann. If gene is present
# it returns start and end location and strand of gene.

import os, sys
from Bio import SeqIO
from helpers import get_organism_from_fasta

gff_dir = "campy_ann"
fasta_dir = "campy_fas"

if len(sys.argv) == 1:
    print("please enter at least one gene to search")
    exit(1)

genes = sys.argv[1:]

print(genes)

print(f"\033[95mSEARCH RESULTS:\033[0m")

for file in os.listdir(gff_dir):
    if not file.endswith(".gff"):
        continue

    print("-" * 70)
    print(f"\033[95mFile: {file}\033[0m")

    gff_path = os.path.join(gff_dir, file)
    fasta_name = file.replace(".gff", ".fna")
    fasta_path = os.path.join(fasta_dir, fasta_name)

    if not os.path.exists(fasta_path):
        print(f"No matching FASTA found for {file}")
        continue

    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    organism_name = get_organism_from_fasta(fasta_path)

    genes_found = {gene: False for gene in genes}

    contig = ""

    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue

            attributes = parts[8]
            contig = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]

            for gene in genes:
                if gene.lower() in attributes.lower():
                    seq = genome[contig].seq[start-1:end]
                    if strand == "-":
                        seq = seq.reverse_complement()

                    print(f"\nGene: {gene}")
                    print(f"Contig: {contig}")
                    print(f"Start: {start}")
                    print(f"End: {end}")
                    print(f"Strand: {strand}")
                    print(f"From: {organism_name}")
                    print(f"\nSource: https://www.ncbi.nlm.nih.gov/nuccore/{contig}")
                    
                    genes_found[gene] = True

    for gene, found in genes_found.items():
        if not found:
            print(f"\nNo matches found for {gene}.")
            print(f"From: {organism_name}")
            print(f"\nSource: https://www.ncbi.nlm.nih.gov/nuccore/{contig}")

    print("-" * 70)
    