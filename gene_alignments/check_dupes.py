#!/usr/bin/env python3

from collections import Counter
from Bio import SeqIO

names = [rec.id for rec in SeqIO.parse("fliI_aligned.fna", "fasta")]
dups = [name for name, count in Counter(names).items() if count > 1]

print("Duplicate IDs:", dups)