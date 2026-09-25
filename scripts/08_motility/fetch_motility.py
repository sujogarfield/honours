#!/usr/bin/env python3
"""
fetch_motility.py
Fetches motility phenotype data from BacDive (DSMZ) for every organism in
the pruned GTDB tree, to show alongside the flagellar gene-order analysis —
motility is the phenotype flagella actually produce.

BacDive's API needs no signup/API key (DSMZ dropped that requirement in
2026). Species are looked up by genus+species name directly; results are
aggregated across every matched strain record, since coverage/detail varies
a lot strain to strain (many returned records have empty morphology data).

Run AFTER gtdb_tree.py (reads images_and_newick/thesis_b/gtdb_ref_tree.nwk for the leaf list).
"""

import json
import re
import time
from pathlib import Path

import bacdive
import requests
from ete3 import Tree

NAMED_NWK = "images_and_newick/thesis_b/gtdb_ref_tree.nwk"
OUTPUT_JSON = "gene_order/motility.json"
REQUEST_DELAY = 0.3  # seconds between BacDive queries — be polite to a live API
PREDICTION_CONFIDENCE_THRESHOLD = 60  # only surfaced in the UI above this

# BacDive's "Genome-based predictions" section (an ML model over genome
# content, separate from curated/observed data) isn't exposed by the public
# API — only rendered server-side into the strain page's HTML — so this is
# scraped directly. Only used as a fallback for organisms with no empirical
# motility record; a different, independent model on the same pages
# sometimes disagrees with this one (confirmed by hand for one species), so
# this is a best-effort secondary signal, not treated as equal to BacDive's
# curated data.
PREDICTION_PATTERN = re.compile(
    r'<td class="border ">flagellated</td>\s*'
    r'<td class="border ">([^<]*)(?:<span[^>]*>.*?</span>)?</td>\s*'
    r'<td class="border ">([^<]*)</td>\s*'
    r'<td class="border ">([\d.]+)</td>',
    re.DOTALL
)


def fetch_flagellated_prediction(bacdive_id):
    """Scrape the 'flagellated' genome-based ML prediction for one strain.
    Returns {"answer": "yes"|"no", "confidence": float} or None if that
    strain's page has no such row."""
    try:
        resp = requests.get(f"https://bacdive.dsmz.de/strain/{bacdive_id}", timeout=20)
        m = PREDICTION_PATTERN.search(resp.text)
        if not m:
            return None
        return {"answer": m.group(2).strip(), "confidence": float(m.group(3))}
    except Exception:
        return None


def as_list(value):
    """BacDive sometimes returns a single dict, sometimes a list of dicts,
    for the same field (e.g. when multiple references disagree). Normalise."""
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def fetch_motility_for(client, organism: str) -> dict:
    """Query BacDive for one organism, aggregating motility across every
    matched strain record. Returns motility=None if nothing matched or no
    strain had morphology data — that's an expected outcome, not an error."""
    genus_species = " ".join(organism.split(" ")[:2])

    try:
        client.setSearchType("exact")
        n_matched = client.search(taxonomy=genus_species)
        strains = list(client.retrieve()) if n_matched else []
    except Exception as e:
        print(f"  {organism}: BacDive query failed ({e})")
        return {"motility": None, "flagellum_arrangement": None,
                "n_strains_matched": 0, "n_strains_with_motility_info": 0}

    motility_values = []
    flagellum_arrangement = None
    for strain in strains:
        for cm in as_list(strain.get("Morphology", {}).get("cell morphology")):
            if not isinstance(cm, dict):
                continue
            if cm.get("motility"):
                motility_values.append(cm["motility"])
            if flagellum_arrangement is None and cm.get("flagellum arrangement"):
                flagellum_arrangement = cm["flagellum arrangement"]

    if any(v == "yes" for v in motility_values):
        motility = "motile"
    elif any(v == "no" for v in motility_values):
        motility = "non-motile"
    else:
        motility = None

    predicted_motility = None
    if motility is None and strains:
        bacdive_id = strains[0].get("General", {}).get("BacDive-ID")
        if bacdive_id:
            predicted_motility = fetch_flagellated_prediction(bacdive_id)
            time.sleep(REQUEST_DELAY)

    return {
        "motility": motility,
        "flagellum_arrangement": flagellum_arrangement,
        "n_strains_matched": len(strains),
        "n_strains_with_motility_info": len(motility_values),
        "predicted_motility": predicted_motility,
    }


def main():
    Path("gene_order").mkdir(exist_ok=True)

    print(f"Loading {NAMED_NWK} ...")
    tree = Tree(NAMED_NWK, format=1, quoted_node_names=True)
    organisms = sorted(leaf.name for leaf in tree.get_leaves())
    print(f"{len(organisms)} organisms to look up")

    client = bacdive.BacdiveClient()

    results = {}
    for i, organism in enumerate(organisms, 1):
        info = fetch_motility_for(client, organism)
        results[organism] = info
        tag = info["motility"] or "unknown"
        pred = info["predicted_motility"]
        pred_note = f", predicted {pred['answer']} ({pred['confidence']}%)" if pred else ""
        print(f"[{i}/{len(organisms)}] {organism}: {tag} "
              f"({info['n_strains_matched']} strains, "
              f"{info['n_strains_with_motility_info']} with motility data{pred_note})")
        time.sleep(REQUEST_DELAY)

    motile = sum(1 for v in results.values() if v["motility"] == "motile")
    non_motile = sum(1 for v in results.values() if v["motility"] == "non-motile")
    unknown = sum(1 for v in results.values() if v["motility"] is None)
    predicted_motile_confident = sum(
        1 for v in results.values()
        if v["predicted_motility"] and v["predicted_motility"]["answer"] == "yes"
        and v["predicted_motility"]["confidence"] > PREDICTION_CONFIDENCE_THRESHOLD
    )
    predicted_non_motile_confident = sum(
        1 for v in results.values()
        if v["predicted_motility"] and v["predicted_motility"]["answer"] == "no"
        and v["predicted_motility"]["confidence"] > PREDICTION_CONFIDENCE_THRESHOLD
    )
    print(f"\nSummary: {motile} motile, {non_motile} non-motile, {unknown} unknown/no BacDive record")
    print(f"Of the unknowns, {predicted_motile_confident} have a genome-based "
          f"'motile' prediction above {PREDICTION_CONFIDENCE_THRESHOLD}% confidence, "
          f"and {predicted_non_motile_confident} have a genome-based "
          f"'non-motile' prediction above {PREDICTION_CONFIDENCE_THRESHOLD}% confidence")

    with open(OUTPUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved: {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
