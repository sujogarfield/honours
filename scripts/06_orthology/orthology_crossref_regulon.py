#!/usr/bin/env python3
"""
orthology_crossref_regulon.py
Comprehensive sweep: instead of adding flagellar/chemotaxis genes one at a
time as gaps get discovered (flgR/flgS/pflA/pflB, maf, flaA/flaB, flhF/
flhG/fliW/flaC/che*, ...), seed from EVERY gene in C. jejuni NCTC 11168
(GCF_000009085.1) whose product description matches flagell*/chemotax*/
chemorecept*/methyl-accepting/motility -- 50 reference proteins, found by
keyword-searching the PRODUCT field (not just gene= tags), which is what
catches genes NCTC11168 itself never gave a gene symbol to (6 of the 50:
multiple methyl-accepting chemotaxis paralogs, a second flagellin copy, a
generic "motility protein"). Also surfaced two genuine gaps in our own
gene name list that had nothing to do with missing gene= tags: fliI
(flagellum-specific ATP synthase) and flaG (flagellar length control) were
simply never added to gene_order.py's search list at all.

Each of the 50 reference proteins is diamond-searched (as query) against
all 166 genomes' proteomes (best hit only), then that hit's OrthoFinder
orthogroup is looked up -- 47 unique orthogroups result (some reference
proteins share an orthogroup, e.g. flaA/flaB both -> OG0000010). Every
member of every one of those 47 orthogroups, across all 166 genomes,
becomes a candidate -- comprehensively, in one pass, rather than
incrementally. Redundancy with already-known genes (e.g. mota/motb/flig,
which string matching alone already covered completely) is harmless: it
just re-derives already-known candidates, filtered out naturally at the
"already annotated" check before the next eggNOG batch.

Run AFTER building /tmp/regulon_ref_ogs.json (diamond search of the 50
NCTC11168 reference proteins against /tmp/campy_diamond_db, best hit's
orthogroup looked up in Orthogroups.tsv -- see conversation for the
one-off commands that produced it).
"""

import csv
import json
import os
import glob
from collections import defaultdict

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from orthology_crossref import build_gff_index, GFF_DIR

ORTHOGROUPS_TSV = "campy_orthologs/orthofinder_results/Results_Sep23/Orthogroups/Orthogroups.tsv"
REF_OGS_JSON = "campy_orthologs/regulon_ref_ogs.json"
OUTPUT_JSON = "campy_orthologs/archive_json/gene_order_regulon_candidates.json"


def main():
    ref_ogs = json.load(open(REF_OGS_JSON))
    seed_ogs = sorted({og for og in ref_ogs.values() if og != "NOT_IN_ANY_OG"})
    print(f"{len(seed_ogs)} unique orthogroups seeded from the NCTC 11168 regulon sweep")

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
    seed_set = set(seed_ogs)
    with open(ORTHOGROUPS_TSV) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        species_cols = header[1:]
        og_members = defaultdict(lambda: defaultdict(list))
        for row in reader:
            if row[0] not in seed_set:
                continue
            for species, cell in zip(species_cols, row[1:]):
                if not cell.strip():
                    continue
                og_members[row[0]][species].extend(p.strip() for p in cell.split(",") if p.strip())

    from helpers import get_organism_from_fasta
    out = defaultdict(lambda: {"organism": None, "genes": []})
    n = 0
    for og_id in seed_ogs:
        for gcf_id, prot_ids in og_members[og_id].items():
            if gcf_id not in protein_idx:
                continue
            for protein_id in prot_ids:
                locus = protein_idx[gcf_id].get(protein_id)
                if locus is None:
                    continue
                contig, start, end, strand = locus
                out[gcf_id]["genes"].append({
                    "gene": f"regulon_{og_id}",  # placeholder -- real identity comes from eggNOG
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "source": "orthology_external_seed",
                    "orthogroup": og_id,
                    "protein_id": protein_id,
                    "seed_note": "NCTC 11168 comprehensive regulon sweep",
                })
                n += 1

    for gcf_id, data in out.items():
        fasta_path = os.path.join("campy_fetched/campy_fas", f"{gcf_id}.fna")
        data["organism"] = get_organism_from_fasta(fasta_path) if os.path.exists(fasta_path) else gcf_id
        for r in data["genes"]:
            r["organism"] = data["organism"]

    print(f"\n{n} candidate genes found across {len(out)} genomes (before dedup against existing data)")
    os.makedirs("gene_order", exist_ok=True)
    with open(OUTPUT_JSON, "w") as f:
        json.dump(dict(out), f, indent=2)
    print(f"Saved to {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
