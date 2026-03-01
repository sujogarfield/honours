#!/usr/bin/env python3

from Bio import SeqIO

def get_organism_from_fasta(fasta_path):
    record = next(SeqIO.parse(fasta_path, "fasta"))
    description = record.description

    parts = description.split(" ", 1)

    if len(parts) > 1:
        name_part = parts[1]

        name_part = name_part.split(" chromosome")[0]
        name_part = name_part.split(" complete genome")[0]

        return name_part.replace(" ", "_")

    return "Unknown"