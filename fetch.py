#!/usr/bin/env python3

# Script that downloads fasta and annotation (gff) files
# from NCBI. Fasta is saved to campy_fas and gff is saved to
# campy_ann. Takes unique species from 200 results.

from Bio import Entrez
import os, urllib.request, gzip, shutil

Entrez.email = "suhanijones2@gmail.com"

search = Entrez.esearch(
    db="assembly",
    term="Campylobacter[Organism]",
    retmax=200
)

record = Entrez.read(search)
assembly_ids = record["IdList"]

summaries = Entrez.esummary(db="assembly", id=",".join(assembly_ids))
records = Entrez.read(summaries)

seen = set()
count = 0

fas_dir = "campy_fas"
ann_dir = "campy_ann"
os.makedirs(fas_dir, exist_ok=True)
os.makedirs(ann_dir, exist_ok=True)

for rec in records["DocumentSummarySet"]["DocumentSummary"]:
    species = rec["Organism"]

    if species not in seen:
        seen.add(species)
        ftp = rec["FtpPath_RefSeq"]

        if ftp:
            base = ftp.split("/")[-1]

            files_to_download = [
                base + "_genomic.fna.gz",
                base + "_genomic.gff.gz"
            ]

            for filename in files_to_download:
                url = ftp + "/" + filename
                
                if filename.endswith("_genomic.fna.gz"):
                    filepath = os.path.join(fas_dir, filename)
                elif filename.endswith("_genomic.gff.gz"):
                    filepath = os.path.join(ann_dir, filename)
                else:
                    continue

                unzippedpath = filepath[:-3]

                if os.path.exists(filepath) or os.path.exists(unzippedpath):
                    print(f"Skipping: {filename}")
                else:
                    print(f"Downloading: {filename}")
                    urllib.request.urlretrieve(url, filepath)
                    # unzip here
                    with gzip.open(filepath, "rb") as f_in:
                        with open(unzippedpath, "wb") as f_out:
                            shutil.copyfileobj(f_in, f_out)
                
                os.remove(filepath)

            count += 1
            if count == 10:
                break
