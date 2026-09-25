#!/usr/bin/env python3
"""
orthology_eggnog_confirm.py
Cross-checks each orthology-recovered candidate gene (found by
orthology_crossref.py, seeded from OrthoFinder orthogroups) against its
independent eggNOG-mapper annotation.

A candidate is "confirmed" if eggNOG's Preferred_name matches the seed gene
name, or the seed gene's 4-letter code appears as a whole word in eggNOG's
Description -- i.e. two unrelated orthology methods (OrthoFinder's
genome-set-specific clustering vs eggNOG's precomputed reference OGs) agree
on the same protein's identity. Unconfirmed candidates are kept (not
dropped) but flagged, since eggNOG can simply lack a precise call for a
divergent protein -- that's not the same as OrthoFinder being wrong.

Run AFTER orthology_crossref.py and after eggnog-mapper has produced
eggnog_out/candidates.emapper.annotations.
"""

import json
import re
from collections import Counter

CANDIDATES_JSON = "campy_orthologs/archive_json/gene_order_orthology_candidates.json"
ANNOTATIONS_TSV = "campy_orthologs/annotations_merged_all.tsv"  # rounds 1-3 combined (see hgt_detection_notes.txt / findings summary)
OUTPUT_JSON = "campy_orthologs/archive_json/gene_order_orthology_candidates.json"  # updated in place


def load_annotations(path: str) -> dict:
    """protein_id -> {preferred_name, description}"""
    ann = {}
    with open(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            fields = line.rstrip("\n").split("\t")
            row = dict(zip(header, fields))
            ann[row["query"]] = {
                "preferred_name": row.get("Preferred_name", "-"),
                "description": row.get("Description", "-"),
            }
    return ann


def classify(gene: str, pref_name: str, description: str) -> str:
    """confirmed: eggNOG's own symbol or description names this exact gene.
    reclassified: eggNOG gives a *different*, specific symbol -- e.g. seed
      was flgF but eggNOG calls it flgG2 (a real rod-protein paralog pair
      OrthoFinder's clustering may not separate). Surfaced, not discarded.
    unconfirmed: an eggNOG hit exists but gives no symbol/description
      support either way (may still be correct -- eggNOG's own annotation
      can simply be incomplete for this protein family)."""
    has_pref = pref_name and pref_name not in ("-", "")
    if has_pref and pref_name.lower() == gene.lower():
        return "confirmed"
    if description and description != "-":
        if re.search(rf"\b{re.escape(gene)}\b", description, re.IGNORECASE):
            return "confirmed"
    if has_pref:
        return "reclassified"
    return "unconfirmed"


def main():
    with open(CANDIDATES_JSON) as f:
        candidates = json.load(f)

    annotations = load_annotations(ANNOTATIONS_TSV)
    print(f"Loaded {len(annotations)} eggNOG annotations")

    status_counts = Counter()
    per_gene = {}

    for gcf_id, data in candidates.items():
        for rec in data["genes"]:
            protein_id = rec["protein_id"]
            gene = rec["gene"]
            ann = annotations.get(protein_id)
            if ann is None:
                rec["eggnog_status"] = "no_hit"
                rec["eggnog_preferred_name"] = None
                rec["eggnog_description"] = None
                status_counts["no_hit"] += 1
            else:
                rec["eggnog_status"] = classify(gene, ann["preferred_name"], ann["description"])
                rec["eggnog_preferred_name"] = ann["preferred_name"]
                rec["eggnog_description"] = ann["description"]
                status_counts[rec["eggnog_status"]] += 1

            g = per_gene.setdefault(gene, Counter())
            g[rec["eggnog_status"]] += 1

    print("\nOverall:")
    for status, n in status_counts.most_common():
        print(f"  {status:12s} {n}")

    print("\nPer gene (confirmed / reclassified / unconfirmed / no_hit):")
    for gene in sorted(per_gene):
        c = per_gene[gene]
        print(f"  {gene:6s} {c.get('confirmed', 0):4d} / {c.get('reclassified', 0):4d} / "
              f"{c.get('unconfirmed', 0):4d} / {c.get('no_hit', 0):4d}")

    reclass_names = Counter()
    for gcf_id, data in candidates.items():
        for rec in data["genes"]:
            if rec.get("eggnog_status") == "reclassified":
                reclass_names[(rec["gene"], rec["eggnog_preferred_name"])] += 1
    if reclass_names:
        print("\nReclassified as (seed_gene -> eggNOG symbol : count):")
        for (gene, alt), n in reclass_names.most_common():
            print(f"  {gene} -> {alt} : {n}")

    with open(OUTPUT_JSON, "w") as f:
        json.dump(candidates, f, indent=2)
    print(f"\nUpdated {OUTPUT_JSON} with eggNOG confirmation status")


if __name__ == "__main__":
    main()
