#!/usr/bin/env python3

import pandas as pd

tsv_file = "campy_fetched/gtdb-adv-search.tsv"
meta = pd.read_csv(tsv_file, sep="\t")

campy = meta[
    (meta["GTDB Taxonomy"].str.contains("p__Campylobacterota")) &
    (meta["GTDB Representative of Species"] == True) &
    # RefSeq (GCF_) only -- GenBank-only (GCA_) assemblies have no RefSeq
    # annotation for fetch.py to download
    (meta["Accession"].str.startswith("GCF_"))
]

refseq_accessions = campy["Accession"].tolist()

print(f"{len(refseq_accessions)} Campylobacterota representative genomes")

with open("campy_fetched/campy_refseq_accessions.txt", "w") as f:
    for acc in refseq_accessions:
        f.write(acc + "\n")
