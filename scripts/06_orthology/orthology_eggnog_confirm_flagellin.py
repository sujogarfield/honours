#!/usr/bin/env python3
"""
orthology_eggnog_confirm_flagellin.py
eggNOG confirmation for the flagellin (flaA/flaB) candidates, seeded
externally (see orthology_crossref_flagellin.py) since OrthoFinder's own
clustering couldn't separate the two paralogs (both landed in OG0000010,
labelled "flaab").

Unlike every other gene in this pipeline, eggNOG's broader reference
database CAN tell flaA and flaB apart even though our own 166-genome
OrthoFinder run couldn't -- so this script relabels each candidate to its
actual eggNOG-resolved identity (flaa/flab) rather than keeping the
lumped "flaab" placeholder, whenever eggNOG gives a specific answer.

Result (2026-09-25): 154 -> flaB, 52 -> flaA, 3 no_hit (dropped).

Run AFTER eggnog-mapper has annotated the flagellin candidate proteins
(merged into /tmp/annotations_merged_all.tsv alongside everything else).
"""

import json

CANDIDATES_JSON = "campy_orthologs/archive_json/gene_order_flagellin_candidates.json"
ANNOTATIONS_TSV = "campy_orthologs/annotations_merged_all.tsv"
OUTPUT_JSON = "campy_orthologs/archive_json/gene_order_flagellin_candidates.json"  # updated in place

RESOLVABLE_NAMES = {"flaa": "flaa", "flab": "flab"}


def load_annotations(path: str) -> dict:
    ann = {}
    header = None
    with open(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            ann[row["query"]] = {
                "preferred_name": row.get("Preferred_name", "-"),
                "description": row.get("Description", "-"),
            }
    return ann


def main():
    with open(CANDIDATES_JSON) as f:
        candidates = json.load(f)

    annotations = load_annotations(ANNOTATIONS_TSV)

    counts = {"flaa": 0, "flab": 0, "no_hit": 0, "unresolved": 0}
    for gcf_id, data in candidates.items():
        for rec in data["genes"]:
            ann = annotations.get(rec["protein_id"])
            if ann is None:
                rec["eggnog_status"] = "no_hit"
                counts["no_hit"] += 1
                continue
            pref = (ann["preferred_name"] or "").lower()
            rec["eggnog_preferred_name"] = ann["preferred_name"]
            rec["eggnog_description"] = ann["description"]
            if pref in RESOLVABLE_NAMES:
                rec["original_seed_gene"] = rec["gene"]  # "flaab"
                rec["gene"] = RESOLVABLE_NAMES[pref]
                rec["eggnog_status"] = "confirmed"
                counts[pref] += 1
            else:
                rec["eggnog_status"] = "unconfirmed"
                counts["unresolved"] += 1

    print("Flagellin resolution:")
    for k, v in counts.items():
        print(f"  {k}: {v}")

    with open(OUTPUT_JSON, "w") as f:
        json.dump(candidates, f, indent=2)
    print(f"\nUpdated {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
