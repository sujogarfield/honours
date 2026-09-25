#!/usr/bin/env python3
"""
fetch_proteins.py
Downloads RefSeq protein FASTA (_protein.faa) for every genome already
present in campy_fetched/campy_ann + campy_fetched/campy_fas, needed as OrthoFinder/eggNOG-mapper input.

FTP path is derived directly from the GCF_xxx_ASMname filename base already
used for campy_fetched/campy_ann + campy_fetched/campy_fas (no Entrez lookup needed).
"""

import os
import urllib.request
import gzip
import shutil

ann_dir = "campy_fetched/campy_ann"
prot_dir = "campy_fetched/campy_prot"
os.makedirs(prot_dir, exist_ok=True)

bases = sorted(f.replace(".gff", "") for f in os.listdir(ann_dir) if f.endswith(".gff"))
print(f"{len(bases)} genomes to fetch protein FASTA for")

failed = []
for i, base in enumerate(bases, 1):
    faa_local = os.path.join(prot_dir, f"{base}.faa")
    if os.path.exists(faa_local):
        print(f"[{i}/{len(bases)}] {base} already exists: skipping")
        continue

    gcf, ver_and_rest = base.split("_", 1)[1].split(".", 1)
    version, _, asm = ver_and_rest.partition("_")
    digits = gcf  # 9-digit accession number, no version
    triplets = "/".join([digits[0:3], digits[3:6], digits[6:9]])
    ftp_path = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{triplets}/{base}"
    url = f"{ftp_path}/{base}_protein.faa.gz"

    try:
        gz_local = faa_local + ".gz"
        print(f"[{i}/{len(bases)}] downloading {url}")
        urllib.request.urlretrieve(url, gz_local)
        with gzip.open(gz_local, "rb") as f_in, open(faa_local, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(gz_local)
    except Exception as e:
        print(f"[{i}/{len(bases)}] FAILED {base}: {e}")
        failed.append(base)

print(f"\nDone. {len(bases) - len(failed)} succeeded, {len(failed)} failed.")
if failed:
    print("Failed:", failed)
