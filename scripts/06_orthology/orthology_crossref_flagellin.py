#!/usr/bin/env python3
"""
orthology_crossref_flagellin.py
Seeds the flagellin structural gene (flaA/flaB) from an EXTERNAL reference
(C. jejuni NCTC 11168, GCF_000009085.1) rather than this dataset's own
string-matched genes: fliC (the enteric-bacteria flagellin gene name,
already in gene_order.py's list) has ZERO hits anywhere in this dataset,
and flaA/flaB (the actual Campylobacter/Helicobacter naming convention)
have essentially zero direct gene= hits too -- flagellin, the core
structural filament protein, is present as an unnamed/hypothetical protein
in the large majority of these annotations.

Reference proteins (from GCF_000009085.1.gff, verified locus_tag lookup):
  flaA (Cj1339c) = YP_002344727.1   flaB (Cj1338c) = YP_002344726.1

Diamond search of these against all 166 proteomes found unambiguous top
hits (92%+ identity, 100% coverage -- see /tmp/flagellin_ref_hits.tsv),
landing in a single orthogroup:
  OG0000010 (209 proteins / 92 genomes) <- flaA AND flaB both land here --
  NOT separable, same as the maf2/3/5/6 situation: these are recent
  tandem-duplicate paralogs too similar for OrthoFinder's clustering to
  tell apart. Members output as gene="flaab" (unresolved paralog family),
  NOT assigned to an individual flaA/flaB identity.

Run AFTER orthofinder has produced Orthogroups.tsv.
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
OUTPUT_JSON = "campy_orthologs/archive_json/gene_order_flagellin_candidates.json"

SEED_OG = "OG0000010"
GENE_LABEL = "flaab"  # unresolved flaA/flaB paralog family


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
        og_members = defaultdict(list)
        for row in reader:
            if row[0] != SEED_OG:
                continue
            for species, cell in zip(species_cols, row[1:]):
                if not cell.strip():
                    continue
                og_members[species].extend(p.strip() for p in cell.split(",") if p.strip())

    from helpers import get_organism_from_fasta
    out = defaultdict(lambda: {"organism": None, "genes": []})
    n = 0
    for gcf_id, prot_ids in og_members.items():
        if gcf_id not in protein_idx:
            continue
        for protein_id in prot_ids:
            locus = protein_idx[gcf_id].get(protein_id)
            if locus is None:
                continue
            contig, start, end, strand = locus
            out[gcf_id]["genes"].append({
                "gene": GENE_LABEL,
                "contig": contig,
                "start": start,
                "end": end,
                "strand": strand,
                "source": "orthology_external_seed",
                "orthogroup": SEED_OG,
                "protein_id": protein_id,
                "seed_note": "unresolved paralog family (flaA/flaB)",
            })
            n += 1

    for gcf_id, data in out.items():
        fasta_path = os.path.join("campy_fetched/campy_fas", f"{gcf_id}.fna")
        data["organism"] = get_organism_from_fasta(fasta_path) if os.path.exists(fasta_path) else gcf_id
        for r in data["genes"]:
            r["organism"] = data["organism"]

    print(f"\n{n} flagellin (flaA/flaB) candidate genes found across {len(out)} genomes")
    os.makedirs("gene_order", exist_ok=True)
    with open(OUTPUT_JSON, "w") as f:
        json.dump(dict(out), f, indent=2)
    print(f"Saved to {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
