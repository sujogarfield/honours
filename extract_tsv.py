#!/usr/bin/env python3

import pandas as pd

tsv_file = "gtdb-adv-search.tsv"
meta = pd.read_csv(tsv_file, sep="\t")

campy = meta[
    (meta["GTDB Taxonomy"].str.contains("p__Campylobacterota")) &
    (meta["GTDB Representative of Species"] == True)
]

refseq_accessions = campy["Accession"].tolist()

print(f"{len(refseq_accessions)} Campylobacterota representative genomes")

with open("campy_refseq_accessions.txt", "w") as f:
    for acc in refseq_accessions:
        f.write(acc + "\n")
