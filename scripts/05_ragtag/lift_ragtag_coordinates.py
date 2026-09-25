#!/usr/bin/env python3
"""
lift_ragtag_coordinates.py
Lifts the original NCBI-GFF-derived flagellar gene coordinates (from
gene_order.json, on the ORIGINAL draft contigs) onto each RagTag-scaffolded
assembly, using the AGP file RagTag itself produced -- no re-annotation, no
new tool in the loop for gene positions, just a coordinate transform on
data already trusted (NCBI's own GFF) and RagTag's own placement record.

Writes gene_order/gene_order_scaffolded.json: a copy of gene_order.json
with the scaffolded organisms' gene records replaced by their lifted
(scaffold-coordinate) versions, each tagged "scaffolded": true and
"scaffold_coverage": <post-lift coverage>. Everything else is untouched.

Run AFTER run_ragtag_scaffold.py.
"""
import json
from collections import Counter
from pathlib import Path

FEASIBILITY_JSON = "gene_order/ragtag_feasibility.json"
GENE_ORDER_JSON = "campy_orthologs/archive_json/gene_order.json"
SCAFFOLD_OUTPUT_DIR = Path("campy_fetched/ragtag_scaffolding/output")
OUTPUT_JSON = "campy_orthologs/archive_json/gene_order_scaffolded.json"


def parse_agp(agp_path) -> dict:
    """Returns {orig_contig: (new_object, object_beg, comp_beg, comp_end, orientation)}."""
    mapping = {}
    for line in open(agp_path):
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9 or parts[4] != "W":
            continue
        obj, obj_beg, obj_end, part_no, ctype, comp_id, comp_beg, comp_end, orient = parts[:9]
        mapping[comp_id] = (obj, int(obj_beg), int(comp_beg), int(comp_end), orient)
    return mapping


def lift_gene(rec: dict, agp_map: dict) -> tuple:
    """(new_contig, new_start, new_end, new_strand). Genes on a contig RagTag
    left unplaced keep their original contig/coordinates unchanged."""
    contig = rec["contig"]
    if contig not in agp_map:
        return contig, rec["start"], rec["end"], rec["strand"]
    obj, obj_beg, comp_beg, comp_end, orient = agp_map[contig]
    if orient == "+":
        new_start = obj_beg + (rec["start"] - comp_beg)
        new_end = obj_beg + (rec["end"] - comp_beg)
        new_strand = rec["strand"]
    else:
        new_start = obj_beg + (comp_end - rec["end"])
        new_end = obj_beg + (comp_end - rec["start"])
        new_strand = "-" if rec["strand"] == "+" else "+"
    return obj, new_start, new_end, new_strand


def coverage_from_records(records: list) -> tuple:
    by_contig = {}
    for r in records:
        by_contig.setdefault(r["contig"], []).append(r)
    if not by_contig:
        return 0, 0
    primary = max(by_contig.values(), key=len)
    return len(primary), len(records)


def main():
    feasibility = json.load(open(FEASIBILITY_JSON))
    gene_order = json.load(open(GENE_ORDER_JSON))

    targets = {org: v for org, v in feasibility.items() if v["ragtag_possible"]}
    patched = 0
    print(f"{'organism':32s} {'before':>8s} {'after':>8s}")

    for org, info in sorted(targets.items()):
        acc = info["accession"]
        agp_path = SCAFFOLD_OUTPUT_DIR / acc / "ragtag_out" / "ragtag.scaffold.agp"
        if not agp_path.exists():
            print(f"{org:32s}  no scaffold output yet ({agp_path}), skipping")
            continue
        if acc not in gene_order:
            print(f"{org:32s}  {acc} not in {GENE_ORDER_JSON}, skipping")
            continue

        agp_map = parse_agp(agp_path)
        orig_records = gene_order[acc]["genes"]
        before_primary, before_total = coverage_from_records(orig_records)

        lifted = []
        for r in orig_records:
            new_contig, new_start, new_end, new_strand = lift_gene(r, agp_map)
            lifted.append({**r, "contig": new_contig, "start": new_start, "end": new_end, "strand": new_strand})

        after_primary, after_total = coverage_from_records(lifted)
        before_cov = before_primary / before_total if before_total else 0
        after_cov = after_primary / after_total if after_total else 0
        print(f"{org:32s} {before_cov:>7.1%} -> {after_cov:>6.1%}")

        gene_order[acc]["genes"] = lifted
        gene_order[acc]["scaffolded"] = True
        gene_order[acc]["scaffold_coverage"] = round(after_cov, 4)
        patched += 1

    Path("gene_order").mkdir(exist_ok=True)
    with open(OUTPUT_JSON, "w") as f:
        json.dump(gene_order, f, indent=2)
    print(f"\nPatched {patched} genomes with RagTag-scaffolded coordinates")
    print(f"Saved: {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
