#!/usr/bin/env python3
"""
orthology_eggnog_confirm_regulon.py
Confirms the comprehensive NCTC 11168 regulon-sweep candidates
(gene_order/gene_order_regulon_candidates.json, seeded by
orthology_crossref_regulon.py) against their eggNOG-mapper annotations
(eggnog_out5/candidates_regulon.emapper.annotations).

Unlike the earlier gene-by-gene confirmation scripts, these candidates were
never seeded with one specific expected gene name -- each carries only a
placeholder "regulon_<orthogroup>", since the sweep's whole point was to
catch genes the incremental approach never named. So there is nothing to
compare eggNOG's answer against; instead, identity comes entirely from
eggNOG, using the same DOMINANCE_THRESHOLD-style rule already used for
seed-name reclassification elsewhere in this pipeline: an orthogroup is
"confirmed" under eggNOG's Preferred_name only if that name is what >= 80%
of its own (has-a-name) members independently got from eggNOG. A mixed or
no-name-majority orthogroup is dropped rather than guessed at -- consistent
with the project's standing rule that both tools must agree, and that
ambiguous cases are excluded rather than force-labelled.

Members with no eggNOG Preferred_name at all (e.g. generic "methyl-accepting
chemotaxis protein" hits, common for the many individually-unnamed MCP
paralogs) are dropped under this strict standard -- real chemotaxis genes,
but without a specific two-tool-agreed identity to keep them under.
"""

import json
from collections import Counter, defaultdict

CANDIDATES_JSON = "campy_orthologs/archive_json/gene_order_regulon_candidates.json"
ANNOTATIONS_TSV = "campy_orthologs/eggnog_out5/candidates_regulon.emapper.annotations"
OUTPUT_JSON = "campy_orthologs/archive_json/gene_order_regulon_candidates.json"
DOMINANCE_THRESHOLD = 0.80


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
            fields = line.rstrip("\n").split("\t")
            row = dict(zip(header, fields))
            ann[row["query"]] = {
                "preferred_name": row.get("Preferred_name", "-"),
                "description": row.get("Description", "-"),
            }
    return ann


def main():
    candidates = json.load(open(CANDIDATES_JSON))
    annotations = load_annotations(ANNOTATIONS_TSV)
    print(f"Loaded {len(annotations)} eggNOG annotations for "
          f"{sum(len(d['genes']) for d in candidates.values())} regulon candidates")

    # Pass 1: attach raw eggNOG results, tally per-orthogroup name votes
    og_votes = defaultdict(Counter)
    for gcf_id, data in candidates.items():
        for rec in data["genes"]:
            ann = annotations.get(rec["protein_id"])
            if ann is None:
                rec["eggnog_preferred_name"] = None
                rec["eggnog_description"] = None
                continue
            rec["eggnog_preferred_name"] = ann["preferred_name"]
            rec["eggnog_description"] = ann["description"]
            name = ann["preferred_name"]
            if name and name not in ("-", ""):
                og_votes[rec["orthogroup"]][name.lower()] += 1

    # Pass 2: per-orthogroup dominance decision
    og_decision = {}  # og_id -> resolved gene name, or None (drop)
    for og_id, votes in og_votes.items():
        total = sum(votes.values())
        name, n = votes.most_common(1)[0]
        og_decision[og_id] = name if n / total >= DOMINANCE_THRESHOLD else None

    # Pass 3: apply, tag status, relabel
    status_counts = Counter()
    for gcf_id, data in candidates.items():
        for rec in data["genes"]:
            og_id = rec["orthogroup"]
            resolved = og_decision.get(og_id)
            if rec.get("eggnog_preferred_name") is None:
                rec["eggnog_status"] = "no_hit"
            elif resolved is None:
                rec["eggnog_status"] = "unconfirmed"
            else:
                rec["eggnog_status"] = "confirmed"
                rec["gene"] = resolved
            status_counts[rec["eggnog_status"]] += 1

    print("\nOverall:")
    for status, n in status_counts.most_common():
        print(f"  {status:12s} {n}")

    print(f"\n{len([o for o, d in og_decision.items() if d])} / {len(og_votes)} "
          f"orthogroups resolved to a dominant (>= {DOMINANCE_THRESHOLD:.0%}) eggNOG name")

    resolved_genes = Counter()
    for og_id, name in og_decision.items():
        if name:
            resolved_genes[name] += og_votes[og_id][name]
    print("\nConfirmed gene -> record count (top 40):")
    for name, n in resolved_genes.most_common(40):
        print(f"  {name:12s} {n}")

    with open(OUTPUT_JSON, "w") as f:
        json.dump(candidates, f, indent=2)
    print(f"\nUpdated {OUTPUT_JSON} with eggNOG confirmation status")


if __name__ == "__main__":
    main()
