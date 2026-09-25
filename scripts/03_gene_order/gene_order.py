#!/usr/bin/env python3

import os
import json
from Bio import SeqIO
from pathlib import Path
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from helpers import get_organism_from_fasta

gff_dir = "campy_fetched/campy_ann"
fasta_dir = "campy_fetched/campy_fas"

genes = ["flie", "flif", "flig", "flih", "flij", "flik",
         "flil", "flim", "flin", "flio", "flip", "fliq",
         "flir", "flhe", "flha", "flhb", "flga", "flgb",
         "flgc", "flgd", "flge", "flgf", "flgg", "flgh",
         "flgi", "flgj", "flgk", "flgl", "flgo", "flgp",
         "flgq", "flgt", "motx", "moty", "flic", "flid",
         "flis", "flit", "mota", "motb", "motc", "motd",
         "flgm", "flgn", "fliy", "fliz",
         # flg-prefixed two-component regulators (control flagellar gene expression,
         # not structural) + epsilonproteobacteria-specific motor disc proteins +
         # flagellin glycosylation island genes -- all missing from the structural
         # flg/fli/mot list above
         "flgr", "flgs", "pfla", "pflb",
         "maf1", "maf2", "maf3", "maf4", "maf5", "maf6", "maf7", "maf9",
         # flhF/flhG: polar flagellum placement/number (FtsY-like GTPase /
         # MinD-like ATPase) -- well-populated (144/15 gene= hits) but
         # missing from the structural list above.
         "flhf", "flhg",
         # fliW: flagellin translation feedback control, well-populated
         # (191 hits) but not flg/fli-structural, so missed until now.
         "fliw",
         # flaA/flaB/flaC: the ACTUAL flagellin gene names used in
         # Campylobacter/Helicobacter annotations (fliC, already in the
         # list above, has zero hits in this dataset -- the enteric-style
         # name simply isn't used here). flaA/flaB have ~zero gene= hits
         # too (structural flagellin is frequently left unnamed --
         # handled via external-reference orthology seeding instead, see
         # orthology_crossref_flagellin.py); flaC listed here since it has
         # some direct hits.
         "flaa", "flab", "flac",
         # che*: chemotaxis genes that gate flagellar motor switching --
         # real, non-trivial counts in this dataset's annotations
         # (cheB/cheP/cheQ/cheW/cheZ; cheA/cheR/cheV/cheY singletons).
         "chea", "cheb", "chep", "cheq", "cher", "chev", "chew", "chey", "chez"]

Path("gene_order").mkdir(exist_ok=True)
all_species_genes = {}

for file in os.listdir(gff_dir):
    gff_path = os.path.join(gff_dir, file)
    fasta_name = file.replace(".gff", ".fna")
    fasta_path = os.path.join(fasta_dir, fasta_name)
    gcf_id = file.replace(".gff", "")
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

            # Extract the gene= tag value exactly (e.g. gene=fliK)
            gene_tag = ""
            for field in attributes.split(";"):
                if field.lower().startswith("gene="):
                    gene_tag = field.split("=", 1)[1].strip().lower()
                    break

            if gene_tag in genes:
                species_genes.append({
                    "gene": gene_tag,
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "organism": organism_name
                })

    for thing in species_genes:
        print(thing)
    
    if species_genes:
        all_species_genes[gcf_id] = {
            "organism": organism_name,
            "genes": species_genes,
        }
        # plot_species(organism_name, species_genes)

with open("campy_orthologs/archive_json/gene_order.json", "w") as f:
    json.dump(all_species_genes, f, indent=2)

print("Gene order saved to gene_order/gene_order.json")
print("See gene_order directory for output")
