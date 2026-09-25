#!/usr/bin/env python3
"""
build_merged_synteny.py
Builds the per-organism dataset the interactive tree's "reliable only"
toggle needs: every organism from synteny_analysis.json, but for the 108
that reach full flagellar-gene coverage (93 originally + 15 recovered by
RagTag scaffolding), swaps in the entry from synteny_analysis_reliable.json
instead — computed among that reliable population, not diluted by the
fragmented ones. Every organism gets an explicit "reliable" flag (drives
the toggle) and "scaffolded" flag (drives the UI marker for gene order
that's inferred via RagTag rather than directly assembled).

Run AFTER synteny_analysis.py and synteny_analysis_full_contigs.py
(against gene_order_scaffolded.json -> synteny_analysis_reliable.json).
"""
import json

ALL_JSON = "campy_orthologs/archive_json/synteny_analysis.json"
RELIABLE_JSON = "campy_orthologs/archive_json/synteny_analysis_reliable.json"
SCAFFOLDED_GENE_ORDER = "campy_orthologs/archive_json/gene_order_scaffolded.json"
OUTPUT_JSON = "campy_orthologs/archive_json/synteny_analysis_merged.json"


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
