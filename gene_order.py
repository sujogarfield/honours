#!/usr/bin/env python3

import os, sys
from Bio import SeqIO
from helpers import get_organism_from_fasta
from collections import defaultdict
from helpers import plot_species

gff_dir = "campy_ann"
fasta_dir = "campy_fas"

genes = ["flie", "flif", "flig", "flih", "flij", "flik",
         "flil", "flim", "flin", "flio", "flip", "fliq",
         "flir", "flhe", "flha", "flhb", "flga", "flgb",
         "flgc", "flgd", "flge", "flgf", "flgg", "flgh", 
         "flgi", "flgj", "flgk", "flgl", "flgo", "flgp",
         "flgq", "flgt", "motx", "moty", "flic", "flid",
         "flis", "flit", "mota", "motb", "motc", "motd",
         "flgm", "flgn", "fliy", "fliz"]

for file in os.listdir(gff_dir):
    gff_path = os.path.join(gff_dir, file)
    fasta_name = file.replace(".gff", ".fna")
    fasta_path = os.path.join(fasta_dir, fasta_name)
    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    organism_name = get_organism_from_fasta(fasta_path)

    species_genes = []

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
                    species_genes.append({
                        "gene": gene,
                        "contig": contig,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "organism": organism_name
                    })

    for thing in species_genes:
        print(thing)
    
    if species_genes:
        plot_species(organism_name, species_genes)

print("see gene_order directory for output")
