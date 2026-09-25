#!/usr/bin/env python3
"""
build_combined_gene_order.py
Merges string-matched genes (gene_order/gene_order.json) with orthology
candidates, restricted to what BOTH OrthoFinder and eggNOG agree on.

Inclusion rule (strict -- confirmed-only, this feeds a quantitative Mantel
test):
  - String-matched genes: always included.
  - eggnog_status == "confirmed": included under the seed gene name. This
    is the ONLY orthology-derived status now included.
  - eggnog_status == "reclassified": always excluded, as of the
    confirmed-only tightening (see below) -- even a dominant (>= 80%)
    alternate identity is eggNOG correcting OrthoFinder's seed name, one
    inferential step removed from direct two-tool agreement. flgG2 and
    mind (flgF's and flhG's former reclassification targets) are dropped
    under this rule; they are real, defensible identifications
    (documented in orthology_findings_summary.txt Round 4) but excluded
    here for the stricter standard.
  - eggnog_status in ("unconfirmed", "no_hit"): always excluded.
  - EXPLICIT FULL EXCLUSION: pflA (orthogroup OG0000627's reclassified/
    unconfirmed members carry unrelated descriptions -- beta-lactamase,
    peptidoglycan synthesis, chaperone folding -- so even its "confirmed"
    calls aren't trustworthy without independent manual verification);
    yuxG (OG0001663, confirmed status but an unreviewed, unusual
    Bacillus-style gene symbol -- see Round 4 notes).

Run AFTER orthology_eggnog_confirm.py.
"""

import json
from collections import defaultdict, Counter

GENE_ORDER_JSON = "campy_orthologs/archive_json/gene_order.json"
CANDIDATES_JSON = "campy_orthologs/archive_json/gene_order_orthology_candidates.json"
MAF_CANDIDATES_JSON = "campy_orthologs/archive_json/gene_order_maf_candidates.json"
FLAGELLIN_CANDIDATES_JSON = "campy_orthologs/archive_json/gene_order_flagellin_candidates.json"
REGULON_CANDIDATES_JSON = "campy_orthologs/archive_json/gene_order_regulon_candidates.json"
OUTPUT_JSON = "gene_order/gene_order_combined.json"

FULL_GENE_EXCLUSIONS = {"pfla", "yuxg"}


def main():
    with open(GENE_ORDER_JSON) as f:
        base = json.load(f)

    combined = defaultdict(lambda: {"organism": None, "genes": []})
    for gcf_id, data in base.items():
        combined[gcf_id]["organism"] = data["organism"]
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
                combined[gcf_id]["genes"].append(dict(rec))
                n_confirmed += 1
            else:
                n_excluded_status += 1

    print(f"\nMain candidates: {n_confirmed} confirmed, "
          f"{n_excluded_status} reclassified/unconfirmed/no_hit (dropped, confirmed-only rule), "
          f"{n_excluded_full_gene} pflA/yuxG (dropped, full exclusion)")

    n_maf_confirmed = n_maf_excluded = 0
    with open(MAF_CANDIDATES_JSON) as f:
        maf = json.load(f)
    for gcf_id, data in maf.items():
        combined[gcf_id]["organism"] = data["organism"]
        for rec in data["genes"]:
            if rec.get("eggnog_status") == "confirmed":
                combined[gcf_id]["genes"].append(dict(rec))
                n_maf_confirmed += 1
            else:
                n_maf_excluded += 1

    print(f"Maf candidates: {n_maf_confirmed} confirmed (maf_2356); "
          f"{n_maf_excluded} excluded (maf7's whole orthogroup has 0% eggNOG support)")

    n_flagellin_confirmed = n_flagellin_excluded = 0
    with open(FLAGELLIN_CANDIDATES_JSON) as f:
        flagellin = json.load(f)
    for gcf_id, data in flagellin.items():
        combined[gcf_id]["organism"] = data["organism"]
        for rec in data["genes"]:
            if rec.get("eggnog_status") == "confirmed":
                combined[gcf_id]["genes"].append(dict(rec))
                n_flagellin_confirmed += 1
            else:
                n_flagellin_excluded += 1

    print(f"Flagellin candidates: {n_flagellin_confirmed} confirmed (flaA/flaB, eggNOG-resolved); "
          f"{n_flagellin_excluded} excluded (no eggNOG hit)")

    # Regulon sweep candidates need an explicit locus-dedup check: unlike
    # maf/flagellin (which had ~zero prior string-match coverage), the
    # sweep's 41 orthogroups re-derive many genes already covered above --
    # by design, since the point was to catch what the incremental approach
    # missed, not to replace it. Only genuinely new loci get added.
    existing_loci = defaultdict(set)
    for gcf_id, data in combined.items():
        for rec in data["genes"]:
            existing_loci[gcf_id].add((rec["contig"], rec["start"], rec["end"], rec["strand"]))

    n_regulon_confirmed_new = n_regulon_confirmed_dup = n_regulon_excluded = 0
    with open(REGULON_CANDIDATES_JSON) as f:
        regulon = json.load(f)
    for gcf_id, data in regulon.items():
        combined[gcf_id]["organism"] = data["organism"]
        for rec in data["genes"]:
            if rec.get("eggnog_status") != "confirmed" or rec["gene"] in FULL_GENE_EXCLUSIONS:
                n_regulon_excluded += 1
                continue
            locus = (rec["contig"], rec["start"], rec["end"], rec["strand"])
            if locus in existing_loci[gcf_id]:
                n_regulon_confirmed_dup += 1
                continue
            combined[gcf_id]["genes"].append(dict(rec))
            existing_loci[gcf_id].add(locus)
            n_regulon_confirmed_new += 1

    print(f"Regulon sweep candidates: {n_regulon_confirmed_new} confirmed and newly added, "
          f"{n_regulon_confirmed_dup} confirmed but already present (deduped by locus), "
          f"{n_regulon_excluded} unconfirmed/no_hit (dropped)")

    total_genes = sum(len(d["genes"]) for d in combined.values())
    print(f"\nCombined dataset: {total_genes} gene records across {len(combined)} genomes")
    print(f"(string-matched baseline had "
          f"{sum(len(d['genes']) for d in base.values())} records across {len(base)} genomes)")

    with open(OUTPUT_JSON, "w") as f:
        json.dump(dict(combined), f, indent=2)
    print(f"\nSaved to {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
