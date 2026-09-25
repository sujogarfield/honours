#!/usr/bin/env python3
"""
check_ragtag_feasibility.py
For every organism whose flagellar-gene coverage isn't 100% (coverage < 1.0
in gene_order/synteny_analysis.json — i.e. its genes are split across more
than one contig), checks NCBI for a real same-species reference assembly at
Complete Genome or Chromosome level (excluding our own accession).

This does NOT run RagTag or download anything — it only flags which
genomes have a legitimate reference to scaffold against, and records which
accession that is, so the actual scaffolding step (run_ragtag_scaffold.py)
can be pointed at exactly those without re-deriving this each time.

A same-species reference is required deliberately: scaffolding against a
different species would impose that species' gene order rather than
recovering this genome's own, which would be circular for a study that's
testing whether gene order agrees with phylogeny. See the project
discussion this was decided in.

Output: gene_order/ragtag_feasibility.json
    {organism: {accession, ragtag_possible, reference_accession,
                 reference_status, reference_n50}}

Run AFTER synteny_analysis.py. Safe to re-run — it doesn't mutate anything
else, and NCBI availability can change over time so re-checking periodically
is reasonable.
"""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

SYNTENY_JSON = "campy_orthologs/archive_json/synteny_analysis.json"
GENE_ORDER_JSON = "campy_orthologs/archive_json/gene_order.json"
OUTPUT_JSON = "gene_order/ragtag_feasibility.json"

REQUEST_DELAY = 0.5  # seconds between NCBI e-utils calls (no API key -> be conservative)


def fetch_json(url: str, params: dict, retries: int = 6) -> dict:
    q = url + "?" + urllib.parse.urlencode(params)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(q, timeout=20) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            if e.code == 429 and attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            raise


def best_reference(species: str, exclude_acc_prefix: str) -> dict | None:
    """Best same-species reference at Complete Genome / Chromosome level,
    excluding our own accession. None if nothing qualifies."""
    term = f'"{species}"[Organism]'
    res = fetch_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                      {"db": "assembly", "term": term, "retmode": "json", "retmax": 50})
    ids = res.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None
    time.sleep(REQUEST_DELAY)
    summ = fetch_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                       {"db": "assembly", "id": ",".join(ids), "retmode": "json"})
    result = summ.get("result", {})
    candidates = []
    for uid in result.get("uids", []):
        doc = result[uid]
        acc = doc.get("assemblyaccession", "")
        if acc.startswith(exclude_acc_prefix):
            continue
        status = doc.get("assemblystatus", "")
        if status not in ("Complete Genome", "Chromosome"):
            continue
        n50 = int(doc.get("contign50") or 0)
        ftp = doc.get("ftppath_refseq") or doc.get("ftppath_genbank")
        if not ftp:
            continue
        rank = 0 if status == "Complete Genome" else 1
        candidates.append((rank, -n50, acc, status, n50, ftp))
    if not candidates:
        return None
    candidates.sort()
    rank, negn50, acc, status, n50, ftp = candidates[0]
    return {"accession": acc, "status": status, "n50": n50, "ftp_path": ftp}


def main():
    synteny = json.load(open(SYNTENY_JSON))
    gene_order = json.load(open(GENE_ORDER_JSON))

    org_to_acc = {data["organism"]: acc for acc, data in gene_order.items()}

    not_full_coverage = sorted(
        org for org, v in synteny["organisms"].items() if v["coverage"] < 1.0
    )
    print(f"{len(not_full_coverage)} organisms below 100% coverage to check")

    results = {}
    if Path(OUTPUT_JSON).exists():
        results = json.load(open(OUTPUT_JSON))
        print(f"Resuming: {len(results)} already checked in {OUTPUT_JSON}")

    for i, org in enumerate(not_full_coverage, 1):
        if org in results:
            continue
        acc = org_to_acc.get(org)
        if not acc:
            print(f"[{i}/{len(not_full_coverage)}] {org}: no accession found, skipping")
            continue
        prefix = acc.split("_")[0] + "_" + acc.split("_")[1]
        ref = best_reference(org, prefix)
        time.sleep(REQUEST_DELAY)

        if ref:
            print(f"[{i}/{len(not_full_coverage)}] {org}: RagTag possible -> {ref['accession']} ({ref['status']})")
            results[org] = {
                "accession": acc,
                "ragtag_possible": True,
                "reference_accession": ref["accession"],
                "reference_status": ref["status"],
                "reference_n50": ref["n50"],
                "reference_ftp_path": ref["ftp_path"],
            }
        else:
            print(f"[{i}/{len(not_full_coverage)}] {org}: no usable reference in NCBI")
            results[org] = {
                "accession": acc,
                "ragtag_possible": False,
                "reference_accession": None,
                "reference_status": None,
                "reference_n50": None,
                "reference_ftp_path": None,
            }

        # Save incrementally so a rate-limit failure partway through doesn't
        # lose progress -- re-running resumes from here.
        with open(OUTPUT_JSON, "w") as f:
            json.dump(results, f, indent=2)

    n_possible = sum(1 for v in results.values() if v["ragtag_possible"])
    print(f"\n{n_possible} of {len(results)} flagged ragtag_possible")
    print(f"Saved: {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
