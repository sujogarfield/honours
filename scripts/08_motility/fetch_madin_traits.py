#!/usr/bin/env python3
"""
fetch_madin_traits.py
Adds a second curated motility source to gene_order/motility.json: Madin et
al. 2020 ("A synthesis of bacterial and archaeal phenotypic trait data",
Scientific Data) -- a species-level aggregation of 26 phenotype databases
(itself including BacDive among its sources), downloaded from the project's
own GitHub repo (output/condensed_species_GTDB.csv, GTDB-taxonomy-matched,
so species names should line up with this project's own GTDB-based tree
without translation).

Coverage check against our 166 genomes: only 27/166 have any Madin motility
entry at all (this dataset is a 2019/2020 snapshot; most of our genomes are
from more recently described Campylobacterota species it predates), but of
those 27, 15 overlap with existing curated BacDive calls -- 15/15 agree,
zero disagreements, which is a real (if small) validation signal -- and 4
resolve genomes that had NO prior curated data at all:
  Aliarcobacter lanthieri, Campylobacter iguaniorum, Hippea maritima,
  Sulfurospirillum deleyianum (all "motile", high strain-agreement).

Adds a new "motility_madin" field per organism (raw value + supporting
strain count/proportion) to motility.json -- does NOT touch the existing
"motility" (BacDive curated) or "predicted_motility" fields. Downstream,
motility_clustering.py's classify_motility() should prefer, in order:
curated BacDive -> Madin et al. -> BacDive's own genome-based ML prediction
(the least trustworthy of the three, per fetch_motility.py's own caveat).

Run AFTER fetch_motility.py.
"""

import csv
import json
import urllib.request

MOTILITY_JSON = "gene_order/motility.json"
MADIN_CSV_URL = ("https://raw.githubusercontent.com/bacteria-archaea-traits/"
                  "bacteria-archaea-traits/master/output/condensed_species_GTDB.csv")
MADIN_CSV_CACHE = "gene_order/madin_condensed_species_GTDB.csv"


def normalize(value: str) -> str:
    """Madin's raw motility values -> this project's motile/non-motile
    convention. 'flagella' is folded into motile (that's exactly our
    structural-gene question); 'gliding'/'axial filament' are non-flagellar
    motility mechanisms -- kept as their own literal value rather than
    silently folded into "motile", since a gliding organism's motility
    doesn't depend on the flagellar gene set this whole analysis is about."""
    if value in ("yes", "flagella"):
        return "motile"
    if value == "no":
        return "non-motile"
    return value  # "gliding" / "axial filament" -- caller decides what to do


def main():
    try:
        with open(MADIN_CSV_CACHE) as f:
            pass
        print(f"Using cached {MADIN_CSV_CACHE}")
    except FileNotFoundError:
        print(f"Downloading {MADIN_CSV_URL} ...")
        urllib.request.urlretrieve(MADIN_CSV_URL, MADIN_CSV_CACHE)
        print(f"Cached to {MADIN_CSV_CACHE}")

    madin = {}
    with open(MADIN_CSV_CACHE) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["species"].strip()  # NOTE: Madin's "species" column is the full binomial, not just the epithet
            if row["motility"] != "NA":
                madin[name] = {
                    "motility": normalize(row["motility"]),
                    "raw_value": row["motility"],
                    "n_strains": int(row["motility.count"]) if row["motility.count"] else None,
                    "agreement_pct": float(row["motility.prop"]) if row["motility.prop"] else None,
                }

    with open(MOTILITY_JSON) as f:
        motility = json.load(f)

    n_matched = n_new = n_agree = n_disagree = 0
    for org, entry in motility.items():
        if org not in madin:
            continue
        n_matched += 1
        entry["motility_madin"] = madin[org]
        if entry.get("motility") is None:
            n_new += 1
        elif entry["motility"] == madin[org]["motility"]:
            n_agree += 1
        elif madin[org]["motility"] in ("motile", "non-motile"):
            n_disagree += 1
            print(f"DISAGREEMENT: {org} -- BacDive curated={entry['motility']!r}, "
                  f"Madin={madin[org]['motility']!r}")

    print(f"\n{n_matched}/{len(motility)} organisms matched in Madin et al.")
    print(f"  {n_new} newly resolved (had no curated BacDive record)")
    print(f"  {n_agree} agree with existing curated BacDive call")
    print(f"  {n_disagree} disagree (see above)")

    with open(MOTILITY_JSON, "w") as f:
        json.dump(motility, f, indent=2)
    print(f"\nUpdated {MOTILITY_JSON} with motility_madin field")


if __name__ == "__main__":
    main()
