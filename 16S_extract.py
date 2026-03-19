#!/usr/bin/env python3
import os
import subprocess
from Bio import SeqIO
from helpers import get_organism_from_fasta
from pathlib import Path

fasta_dir = "campy_fas"
output_dir = "campy_16S"
merged_fasta = Path("gene_extractions") / "16S_extracted.fna"

os.makedirs(output_dir, exist_ok=True)

all_records = []

for fasta_file in os.listdir(fasta_dir):
    if not fasta_file.endswith(".fna") and not fasta_file.endswith(".fa"):
        continue

    fasta_path = os.path.join(fasta_dir, fasta_file)
    basename = os.path.splitext(fasta_file)[0]
    print(f"\nProcessing {fasta_file} ...")

    outseq = os.path.join(output_dir, f"{basename}.16S.fna")
    subprocess.run([
        "barrnap",
        "--kingdom", "bac",
        "--threads", "1",
        "--outseq", outseq,
        fasta_path
    ], check=True)

    fai_file = fasta_path + ".fai"
    if os.path.exists(fai_file):
        os.remove(fai_file)

    records = list(SeqIO.parse(outseq, "fasta"))
    first_16S = records[0]
    species_name = get_organism_from_fasta(fasta_path)
    first_16S.id = species_name
    first_16S.description = ""
    all_records.append(first_16S)

SeqIO.write(all_records, merged_fasta, "fasta")
print(f"\nSaved {len(all_records)} 16S sequences to {merged_fasta}")
