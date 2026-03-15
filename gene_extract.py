#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys, os
from pathlib import Path
from helpers import get_organism_from_fasta

gff_dir = "campy_ann"
fasta_dir = "campy_fas"
extract_dir = Path("gene_extractions")
extract_dir.mkdir(exist_ok=True)

if len(sys.argv) == 1:
    print("please enter one gene to extract")
    exit(1)

gene = sys.argv[1]

records = []

print(f"\033[36mEXTRACTION:\033[0m")

for file in os.listdir(gff_dir):
    if not file.endswith(".gff"):
        continue

    print("-" * 70)
    print(f"\033[36mFile: {file}\033[0m")

    gff_path = os.path.join(gff_dir, file)
    fasta_name = file.replace(".gff", ".fna")
    fasta_path = os.path.join(fasta_dir, fasta_name)

    if not os.path.exists(fasta_path):
        print(f"No matching FASTA found for {file}")
        continue

    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    organism_name = get_organism_from_fasta(fasta_path)

    found = False

    contig = ""

    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            attributes = ""
            if parts[2] == "gene":
                attributes = parts[8]

            if gene.lower() in attributes.lower():
                contig = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]

                print(f"\nGene: {gene}")
                print(f"Contig: {contig}")
                print(f"Start: {start}")
                print(f"End: {end}")
                print(f"Strand: {strand}")

                seq = genome[contig].seq[start-1:end]

                if strand == "-":
                    seq = seq.reverse_complement()

                record_id = organism_name

                print(f"From: {organism_name}")
                print(f"\nSource: https://www.ncbi.nlm.nih.gov/nuccore/{contig}")

                records.append(
                    SeqRecord(seq, id=record_id, description="")
                )

                found = True
    
    if not found:
        print("No matches found.")
        print(f"From: {organism_name}")
        print(f"\nSource: https://www.ncbi.nlm.nih.gov/nuccore/{contig}")

    print("-" * 70)

output_file = extract_dir / f"{gene}_extracted.fna"
SeqIO.write(records, output_file, "fasta")

print(f"\n\033[36mSaved extracted sequences to {extract_dir}/{gene}_extracted.fna\033[0m")
