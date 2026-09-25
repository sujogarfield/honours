#!/usr/bin/env python3
"""
orthology_crossref_maf.py
Seeds orthogroups for the flagellin-glycosylation `maf` gene family from an
EXTERNAL reference (C. jejuni NCTC 11168, GCF_000009085.1) rather than from
this dataset's own string-matched genes, since string matching found zero
maf hits across all 166 genomes (3 of the 5 reference maf genes are
annotated "hypothetical protein" with no gene= tag even in NCTC 11168 --
there was never a name to match against).

Reference proteins (from GCF_000009085.1.gff, verified locus_tag lookup):
  maf2 (pseD) = YP_002344722.1   maf5 (pseE) = YP_002344725.1
  maf3         = YP_002344723.1   maf6         = YP_002344729.1
  maf7         = YP_002344730.1

Diamond search of these against all 166 proteomes (see /tmp/maf_ref_hits.tsv)
found unambiguous top hits, which land in two OrthoFinder orthogroups:
  OG0000011 <- maf2, maf3, maf5, maf6 (NOT separable -- see caveat below)
  OG0001020 <- maf7

CAVEAT: maf2/maf3/maf5/maf6 all cluster into the SAME orthogroup -- these
paralogs are too similar for OrthoFinder's clustering to separate at this
resolution. Members of OG0000011 are therefore output as gene="maf_2356"
(an unresolved paralog family), NOT assigned to an individual maf number.
Only OG0001020 (maf7) gets a specific, individual gene label.

Run AFTER orthofinder has produced Orthogroups.tsv.
"""

import csv
import json
from collections import defaultdict

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from orthology_crossref import build_gff_index, GFF_DIR
import glob
import os

ORTHOGROUPS_TSV = "campy_orthologs/orthofinder_results/Results_Sep23/Orthogroups/Orthogroups.tsv"
OUTPUT_JSON = "campy_orthologs/archive_json/gene_order_maf_candidates.json"

SEED_OGS = {
    "OG0000011": "maf_2356",  # maf2/maf3/maf5/maf6 -- unresolved paralog family
    "OG0001020": "maf7",
}


def main():
    print("Indexing GFF CDS <-> protein_id for all genomes ...")
    gff_bases = sorted(
        os.path.basename(p)[:-4] for p in glob.glob(os.path.join(GFF_DIR, "*.gff"))
    )
    protein_idx = {}
    for base in gff_bases:
        _, p2l = build_gff_index(os.path.join(GFF_DIR, f"{base}.gff"))
        protein_idx[base] = p2l
    print(f"Indexed {len(gff_bases)} genomes")

    print("Loading Orthogroups.tsv ...")
    with open(ORTHOGROUPS_TSV) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        species_cols = header[1:]
        og_members = defaultdict(lambda: defaultdict(list))
        for row in reader:
            og_id = row[0]
            if og_id not in SEED_OGS:
                continue
            for species, cell in zip(species_cols, row[1:]):
                if not cell.strip():
                    continue
                prot_ids = [p.strip() for p in cell.split(",") if p.strip()]
                og_members[og_id][species].extend(prot_ids)

    from helpers import get_organism_from_fasta
    out = defaultdict(lambda: {"organism": None, "genes": []})
    n = 0
    for og_id, gene_label in SEED_OGS.items():
        for gcf_id, prot_ids in og_members[og_id].items():
            if gcf_id not in protein_idx:
                continue
            for protein_id in prot_ids:
                locus = protein_idx[gcf_id].get(protein_id)
                if locus is None:
                    continue
                contig, start, end, strand = locus
                out[gcf_id]["genes"].append({
                    "gene": gene_label,
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "source": "orthology_external_seed",
                    "orthogroup": og_id,
                    "protein_id": protein_id,
                    "seed_note": "unresolved paralog family (maf2/3/5/6)" if og_id == "OG0000011"
                                 else "distinct orthogroup",
                })
                n += 1

    for gcf_id, data in out.items():
        fasta_path = os.path.join("campy_fetched/campy_fas", f"{gcf_id}.fna")
        data["organism"] = get_organism_from_fasta(fasta_path) if os.path.exists(fasta_path) else gcf_id
        for r in data["genes"]:
            r["organism"] = data["organism"]

    print(f"\n{n} maf-family candidate genes found across {len(out)} genomes")
    os.makedirs("gene_order", exist_ok=True)
    with open(OUTPUT_JSON, "w") as f:
        json.dump(dict(out), f, indent=2)
    print(f"Saved to {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
