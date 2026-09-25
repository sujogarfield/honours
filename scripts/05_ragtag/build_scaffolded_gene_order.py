#!/usr/bin/env python3
"""
build_scaffolded_gene_order.py
Folds the RagTag-scaffolded gene coordinates (lifted from the original
NCBI GFF positions onto the new scaffolds via each genome's AGP file --
see the RagTag pilot) into a copy of gene_order.json.

For the 17 genomes RagTag was run on, replaces their gene records with the
lifted (scaffold-coordinate) versions. Every other genome is left exactly
as-is. Does NOT touch gene_order/gene_order.json itself -- writes a
separate file so the "as originally sequenced" dataset stays available.

Input: gene_order/gene_order.json, the RagTag pilot's liftover_results.json
Output: gene_order/gene_order_scaffolded.json
"""
import json
from pathlib import Path

GENE_ORDER_JSON = "campy_orthologs/archive_json/gene_order.json"
LIFTOVER_JSON = ("/private/tmp/claude-501/-Users-suhani-Desktop-honours-honours/"
                  "586a0baf-7a38-46f1-8eb9-b113de76f63a/scratchpad/ragtag_pilot/liftover_results.json")
OUTPUT_JSON = "campy_orthologs/archive_json/gene_order_scaffolded.json"


def main():
    gene_order = json.load(open(GENE_ORDER_JSON))
    liftover = json.load(open(LIFTOVER_JSON))

    patched = 0
    for entry in liftover:
        acc = entry["accession"]
        if acc not in gene_order:
            print(f"WARNING: {acc} not found in {GENE_ORDER_JSON}, skipping")
            continue
        gene_order[acc]["genes"] = entry["lifted_gene_records"]
        gene_order[acc]["scaffolded"] = True
        gene_order[acc]["scaffold_coverage"] = entry["after_coverage"]
        patched += 1

    Path("gene_order").mkdir(exist_ok=True)
    with open(OUTPUT_JSON, "w") as f:
        json.dump(gene_order, f, indent=2)
    print(f"Patched {patched} genomes with RagTag-scaffolded coordinates")
    print(f"Saved: {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
