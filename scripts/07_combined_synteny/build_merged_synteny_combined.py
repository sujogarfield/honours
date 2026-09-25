#!/usr/bin/env python3
"""
build_merged_synteny_combined.py
Same as build_merged_synteny.py, but for the combined (string match +
orthology) dataset: swaps in the reliable-population entry (coverage==1.0,
108 genomes including RagTag-scaffolded) for organisms that reach it, flags
every organism "reliable"/"scaffolded" -- this is what interactive_tree.py's
"reliable only" toggle needs.

Run AFTER orthology_synteny_preview.py and synteny_analysis_full_contigs.py
have both been run against gene_order/gene_order_combined_scaffolded.json.
"""
import json

ALL_JSON = "gene_order/synteny_analysis_combined_scaffolded.json"
RELIABLE_JSON = "gene_order/synteny_analysis_combined_reliable.json"
SCAFFOLDED_GENE_ORDER = "gene_order/gene_order_combined_scaffolded.json"
OUTPUT_JSON = "gene_order/synteny_analysis_combined_merged.json"


def main():
    all_data = json.load(open(ALL_JSON))
    reliable_data = json.load(open(RELIABLE_JSON))
    scaffolded_gene_order = json.load(open(SCAFFOLDED_GENE_ORDER))

    scaffolded_organisms = {
        data["organism"] for data in scaffolded_gene_order.values()
        if data.get("scaffolded")
    }

    merged = {}
    for org, entry in all_data["organisms"].items():
        merged[org] = dict(entry)
        merged[org]["reliable"] = False
        merged[org]["scaffolded"] = org in scaffolded_organisms

    for org, entry in reliable_data["organisms"].items():
        merged[org] = dict(entry)
        merged[org]["reliable"] = True
        merged[org]["scaffolded"] = org in scaffolded_organisms

    n_reliable = sum(1 for v in merged.values() if v["reliable"])
    n_scaffolded = sum(1 for v in merged.values() if v["scaffolded"])
    print(f"{len(merged)} organisms total, {n_reliable} reliable, {n_scaffolded} scaffolded (RagTag)")

    output = {"organisms": merged}
    with open(OUTPUT_JSON, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Saved: {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
