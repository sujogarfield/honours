#!/usr/bin/env python3
"""
build_combined_gene_order_scaffolded.py
Same strict inclusion rule as build_combined_gene_order.py (see that file's
docstring), based on gene_order_scaffolded.json instead of the raw
gene_order.json, so the 17 RagTag-scaffolded genomes are included with
their lifted (scaffold-coordinate) gene positions.

Orthology candidates for the 17 scaffolded genomes are on the original
draft contigs (campy_fetched/campy_ann GFF coordinates), so they are lifted
onto the scaffolds with the same AGP transform lift_ragtag_coordinates.py
applies to the string-match genes (lift_gene). Previously they were skipped,
which left these genomes with a string-match-only gene set -- and so made
some look fully scaffolded when an orthology-only gene sat on an unplaced
contig.

Run AFTER lift_ragtag_coordinates.py and orthology_eggnog_confirm.py.
"""

import json
import os
import sys
from collections import defaultdict, Counter
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "05_ragtag"))
from lift_ragtag_coordinates import SCAFFOLD_OUTPUT_DIR, parse_agp, lift_gene

GENE_ORDER_SCAFFOLDED_JSON = "campy_orthologs/archive_json/gene_order_scaffolded.json"
CANDIDATES_JSON = "campy_orthologs/archive_json/gene_order_orthology_candidates.json"
MAF_CANDIDATES_JSON = "campy_orthologs/archive_json/gene_order_maf_candidates.json"
FLAGELLIN_CANDIDATES_JSON = "campy_orthologs/archive_json/gene_order_flagellin_candidates.json"
REGULON_CANDIDATES_JSON = "campy_orthologs/archive_json/gene_order_regulon_candidates.json"
OUTPUT_JSON = "gene_order/gene_order_combined_scaffolded.json"

FULL_GENE_EXCLUSIONS = {"pfla", "yuxg"}


def main():
    with open(GENE_ORDER_SCAFFOLDED_JSON) as f:
        base = json.load(f)

    scaffolded_gcfs = {gcf for gcf, d in base.items() if d.get("scaffolded")}
    agp_maps = {gcf: parse_agp(SCAFFOLD_OUTPUT_DIR / gcf / "ragtag_out" / "ragtag.scaffold.agp")
                for gcf in scaffolded_gcfs}
    print(f"{len(scaffolded_gcfs)} genomes are RagTag-scaffolded -- orthology "
          f"candidates lifted onto their scaffolds via the RagTag AGP")

    def placed(gcf_id, rec):
        """rec as it should be stored: lifted to scaffold coordinates for
        RagTag-scaffolded genomes, unchanged otherwise."""
        if gcf_id not in scaffolded_gcfs:
            return dict(rec)
        contig, start, end, strand = lift_gene(rec, agp_maps[gcf_id])
        return {**rec, "contig": contig, "start": start, "end": end, "strand": strand}

    combined = defaultdict(lambda: {"organism": None, "genes": []})
    for gcf_id, data in base.items():
        combined[gcf_id]["organism"] = data["organism"]
        combined[gcf_id]["scaffolded"] = data.get("scaffolded", False)
        for rec in data["genes"]:
            r = dict(rec)
            r.setdefault("source", "string_match")
            combined[gcf_id]["genes"].append(r)

    with open(CANDIDATES_JSON) as f:
        cand = json.load(f)

    n_confirmed = n_excluded_status = n_excluded_full_gene = 0
    for gcf_id, data in cand.items():
        combined[gcf_id]["organism"] = data["organism"]
        for rec in data["genes"]:
            gene = rec["gene"]
            status = rec.get("eggnog_status")

            if gene in FULL_GENE_EXCLUSIONS:
                n_excluded_full_gene += 1
                continue

            if status == "confirmed":
                combined[gcf_id]["genes"].append(placed(gcf_id, rec))
                n_confirmed += 1
            else:
                n_excluded_status += 1

    print(f"Main candidates: {n_confirmed} confirmed, "
          f"{n_excluded_status} reclassified/unconfirmed/no_hit (dropped, confirmed-only rule), "
          f"{n_excluded_full_gene} pflA/yuxG (dropped)")

    n_maf_confirmed = n_maf_excluded = 0
    with open(MAF_CANDIDATES_JSON) as f:
        maf = json.load(f)
    for gcf_id, data in maf.items():
        combined[gcf_id]["organism"] = data["organism"]
        for rec in data["genes"]:
            if rec.get("eggnog_status") == "confirmed":
                combined[gcf_id]["genes"].append(placed(gcf_id, rec))
                n_maf_confirmed += 1
            else:
                n_maf_excluded += 1

    print(f"Maf candidates: {n_maf_confirmed} confirmed; {n_maf_excluded} excluded")

    n_flagellin_confirmed = n_flagellin_excluded = 0
    with open(FLAGELLIN_CANDIDATES_JSON) as f:
        flagellin = json.load(f)
    for gcf_id, data in flagellin.items():
        combined[gcf_id]["organism"] = data["organism"]
        for rec in data["genes"]:
            if rec.get("eggnog_status") == "confirmed":
                combined[gcf_id]["genes"].append(placed(gcf_id, rec))
                n_flagellin_confirmed += 1
            else:
                n_flagellin_excluded += 1

    print(f"Flagellin candidates: {n_flagellin_confirmed} confirmed; {n_flagellin_excluded} excluded")

    existing_loci = defaultdict(set)
    for gcf_id, data in combined.items():
        for rec in data["genes"]:
            existing_loci[gcf_id].add((rec["contig"], rec["start"], rec["end"], rec["strand"]))

    n_regulon_new = n_regulon_dup = n_regulon_excluded = 0
    with open(REGULON_CANDIDATES_JSON) as f:
        regulon = json.load(f)
    for gcf_id, data in regulon.items():
        combined[gcf_id]["organism"] = data["organism"]
        for rec in data["genes"]:
            if rec.get("eggnog_status") != "confirmed" or rec["gene"] in FULL_GENE_EXCLUSIONS:
                n_regulon_excluded += 1
                continue
            rec = placed(gcf_id, rec)
            locus = (rec["contig"], rec["start"], rec["end"], rec["strand"])
            if locus in existing_loci[gcf_id]:
                n_regulon_dup += 1
                continue
            combined[gcf_id]["genes"].append(rec)
            existing_loci[gcf_id].add(locus)
            n_regulon_new += 1

    print(f"Regulon sweep candidates: {n_regulon_new} confirmed and newly added; "
          f"{n_regulon_dup} confirmed but already present (deduped); "
          f"{n_regulon_excluded} unconfirmed/no_hit (dropped)")

    total_genes = sum(len(d["genes"]) for d in combined.values())
    print(f"\nCombined+scaffolded dataset: {total_genes} gene records across {len(combined)} genomes")

    with open(OUTPUT_JSON, "w") as f:
        json.dump(dict(combined), f, indent=2)
    print(f"Saved to {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
