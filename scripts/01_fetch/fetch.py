#!/usr/bin/env python3

# Script that downloads fasta, annotation (gff) and proteome (faa) files
# from NCBI. Fasta is saved to campy_fetched/campy_fas, gff is saved to
# campy_fetched/campy_ann and protein fasta (OrthoFinder/eggNOG-mapper input)
# is saved to campy_fetched/campy_prot.

from Bio import Entrez
import os, urllib.request, gzip, shutil

fas_dir = "campy_fetched/campy_fas"
ann_dir = "campy_fetched/campy_ann"
prot_dir = "campy_fetched/campy_prot"
os.makedirs(fas_dir, exist_ok=True)
os.makedirs(ann_dir, exist_ok=True)
os.makedirs(prot_dir, exist_ok=True)

Entrez.email = "suhanijones2@gmail.com"

with open("campy_fetched/campy_refseq_accessions.txt") as f:
    gcf_accessions = [line.strip() for line in f if line.strip()]

for gcf in gcf_accessions:
    handle = Entrez.esearch(db="assembly", term=gcf)
    record = Entrez.read(handle)
    handle.close()

    if len(record["IdList"]) == 0:
        print(f"{gcf} not found")
        continue
    
    assembly_id = record["IdList"][0]

    summary_handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
    summary = Entrez.read(summary_handle)
    summary_handle.close()
    
    ftp_path = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
    if not ftp_path:
        print(f"no FTP path for {gcf}")
        continue
    
    fname_base = ftp_path.split("/")[-1]
    fna_local = os.path.join(fas_dir, f"{fname_base}.fna")
    gff_local = os.path.join(ann_dir, f"{fname_base}.gff")
    faa_local = os.path.join(prot_dir, f"{fname_base}.faa")
    fna_url = f"{ftp_path}/{fname_base}_genomic.fna.gz"
    gff_url = f"{ftp_path}/{fname_base}_genomic.gff.gz"
    faa_url = f"{ftp_path}/{fname_base}_protein.faa.gz"

    if not os.path.exists(fna_local):
        print(f"downloading fasta {fna_url} ...")
        fna_local_gz = fna_local + ".gz"
        urllib.request.urlretrieve(fna_url, fna_local_gz)
        with gzip.open(fna_local_gz, 'rb') as f_in, open(fna_local, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(fna_local_gz)
    else:
        print(f"{fna_local} already exists: skipping")

    if not os.path.exists(gff_local):
        print(f"downloading annotation {gff_url} ...")
        gff_local_gz = gff_local + ".gz"
        urllib.request.urlretrieve(gff_url, gff_local_gz)
        with gzip.open(gff_local_gz, 'rb') as f_in, open(gff_local, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(gff_local_gz)
    else:
        print(f"{gff_local} already exists: skipping")

    if not os.path.exists(faa_local):
        print(f"downloading proteome {faa_url} ...")
        faa_local_gz = faa_local + ".gz"
        urllib.request.urlretrieve(faa_url, faa_local_gz)
        with gzip.open(faa_local_gz, 'rb') as f_in, open(faa_local, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(faa_local_gz)
    else:
        print(f"{faa_local} already exists: skipping")

    print(f"finished {gcf}\n")
