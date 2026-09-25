#!/usr/bin/env python3
"""
motility_clustering.py
Prepares the data behind the "Gene content vs. motility" tab in the
synteny dashboard: for each gene in the combined (string-match +
OrthoFinder+eggNOG-verified) dataset, how often is it present in genomes
known to be motile vs. known to be non-motile -- and which genes cluster
together by having similar presence patterns across genomes.

Motility classification (gene_order/motility.json, from fetch_motility.py):
  - "curated": BacDive's curated motility field (motility: "motile" /
    "non-motile"). Used whenever available -- 65 genomes (57 motile, 8
    non-motile). The non-motile group is small; treat differences built on
    it as suggestive, not conclusive.
  - "predicted": BacDive's genome-based ML prediction, used ONLY as a
    fallback for genomes with no curated record, and only kept if the
    dashboard's toggle enables it (off by default -- fetch_motility.py's
    own docstring flags this as a secondary, sometimes-disagreeing signal).

Gene clustering: genes are ordered by average-linkage hierarchical
clustering on Jaccard distance between their presence patterns across
genomes with known motility (1 - |genomes_with_both| / |genomes_with_either|)
-- genes that tend to co-occur in the same genomes end up adjacent, so a
gene specifically tied to motile genomes clusters away from one that's
either universal or tied to non-motile genomes. Hand-rolled (no scipy in
the venv), same average-linkage approach as synteny_analysis.py's
upgma_tree, but here only the leaf ORDER is needed, not a rendered
dendrogram.

Run AFTER build_combined_gene_order_scaffolded.py and fetch_motility.py.
"""

import json
import itertools

GENE_ORDER_JSON = "gene_order/gene_order_combined_scaffolded.json"
MOTILITY_JSON = "gene_order/motility.json"
OUTPUT_JSON = "gene_order/motility_clustering.json"

PREDICTION_CONFIDENCE_THRESHOLD = 60


def classify_motility(entry: dict) -> tuple[str, str]:
    """Returns (status, source) where status in ("motile","non-motile","unknown")
    and source in ("curated","madin","predicted","none"). Priority: BacDive
    curated (highest trust) -> Madin et al. 2020 synthesis (also curated-
    style, but a broader multi-database aggregation -- see
    fetch_madin_traits.py; 15/15 agreement with BacDive where both have
    data) -> BacDive's own genome-based ML prediction (lowest trust, a
    secondary/sometimes-disagreeing signal per fetch_motility.py)."""
    curated = entry.get("motility")
    if curated in ("motile", "non-motile"):
        return curated, "curated"

    madin = entry.get("motility_madin")
    if madin and madin.get("motility") in ("motile", "non-motile"):
        return madin["motility"], "madin"

    pred = entry.get("predicted_motility")
    if pred and pred.get("confidence", 0) >= PREDICTION_CONFIDENCE_THRESHOLD:
        status = "motile" if pred["answer"] == "yes" else "non-motile" if pred["answer"] == "no" else None
        if status:
            return status, "predicted"
    return "unknown", "none"


def jaccard_distance(a: set, b: set) -> float:
    union = a | b
    if not union:
        return 0.0
    return 1 - len(a & b) / len(union)


def upgma_leaf_order(items: list[str], dist: dict) -> list[str]:
    """items: list of ids. dist: {(i,j): distance} for i<j. Average-linkage
    clustering; returns items reordered so adjacent items are most similar
    (in-order traversal of the resulting merge tree)."""
    n = len(items)
    if n <= 1:
        return items

    clusters = {i: {"members": [i], "leaf_order": [i]} for i in range(n)}
    active = list(range(n))

    def cdist(a, b):
        pairs = list(itertools.product(clusters[a]["members"], clusters[b]["members"]))
        vals = [dist[(min(i, j), max(i, j))] for i, j in pairs if i != j]
        return sum(vals) / len(vals) if vals else 0.0

    next_id = n
    while len(active) > 1:
        best = None
        for x, y in itertools.combinations(active, 2):
            d = cdist(x, y)
            if best is None or d < best[0]:
                best = (d, x, y)
        _, a, b = best
        clusters[next_id] = {
            "members": clusters[a]["members"] + clusters[b]["members"],
            "leaf_order": clusters[a]["leaf_order"] + clusters[b]["leaf_order"],
        }
        active = [c for c in active if c not in (a, b)] + [next_id]
        next_id += 1

    return [items[i] for i in clusters[active[0]]["leaf_order"]]


def main():
    with open(GENE_ORDER_JSON) as f:
        gene_order = json.load(f)
    with open(MOTILITY_JSON) as f:
        motility_raw = json.load(f)

    # organism -> set of gene names present (deduped; paralogs count once)
    presence: dict[str, set[str]] = {}
    for gcf_id, data in gene_order.items():
        org = data["organism"]
        genes = {g["gene"].lower() for g in data["genes"]}
        presence.setdefault(org, set()).update(genes)

    organisms = sorted(presence)
    all_genes = sorted({g for genes in presence.values() for g in genes})

    org_status = {}
    for org in organisms:
        entry = motility_raw.get(org, {})
        status, source = classify_motility(entry)
        org_status[org] = {"status": status, "source": source}

    # Madin et al. is itself curated-style data (validated 15/15 agreement
    # with BacDive where both have entries -- see fetch_madin_traits.py), so
    # it's grouped with "curated" for the trusted/default population. Only
    # BacDive's own genome-based ML prediction sits behind the opt-in toggle.
    TRUSTED_SOURCES = {"curated", "madin"}

    n_trusted_motile = sum(1 for v in org_status.values() if v["status"] == "motile" and v["source"] in TRUSTED_SOURCES)
    n_trusted_nonmotile = sum(1 for v in org_status.values() if v["status"] == "non-motile" and v["source"] in TRUSTED_SOURCES)
    n_madin_only = sum(1 for v in org_status.values() if v["source"] == "madin")
    n_predicted_motile = sum(1 for v in org_status.values() if v["status"] == "motile" and v["source"] == "predicted")
    n_predicted_nonmotile = sum(1 for v in org_status.values() if v["status"] == "non-motile" and v["source"] == "predicted")
    print(f"Motility classification: {n_trusted_motile} trusted-motile (curated+madin), "
          f"{n_trusted_nonmotile} trusted-non-motile (curated+madin), of which {n_madin_only} newly from madin; "
          f"{n_predicted_motile} predicted-motile, {n_predicted_nonmotile} predicted-non-motile, "
          f"{len(organisms) - n_trusted_motile - n_trusted_nonmotile - n_predicted_motile - n_predicted_nonmotile} unknown")

    # Per-gene presence stats, computed separately for the trusted-only and
    # trusted+predicted groupings (the dashboard toggles between them).
    def gene_stats(gene: str, include_predicted: bool):
        motile_orgs = [o for o in organisms if org_status[o]["status"] == "motile"
                       and (org_status[o]["source"] in TRUSTED_SOURCES or include_predicted)]
        nonmotile_orgs = [o for o in organisms if org_status[o]["status"] == "non-motile"
                          and (org_status[o]["source"] in TRUSTED_SOURCES or include_predicted)]
        n_motile = len(motile_orgs)
        n_nonmotile = len(nonmotile_orgs)
        pct_motile = 100 * sum(1 for o in motile_orgs if gene in presence[o]) / n_motile if n_motile else None
        pct_nonmotile = 100 * sum(1 for o in nonmotile_orgs if gene in presence[o]) / n_nonmotile if n_nonmotile else None
        return pct_motile, pct_nonmotile, n_motile, n_nonmotile

    gene_records = {}
    for gene in all_genes:
        pm_c, pn_c, nm_c, nn_c = gene_stats(gene, include_predicted=False)
        pm_p, pn_p, nm_p, nn_p = gene_stats(gene, include_predicted=True)
        gene_records[gene] = {
            "curated": {"pct_motile": pm_c, "pct_nonmotile": pn_c, "n_motile": nm_c, "n_nonmotile": nn_c},
            "extended": {"pct_motile": pm_p, "pct_nonmotile": pn_p, "n_motile": nm_p, "n_nonmotile": nn_p},
        }

    # Cluster gene rows by Jaccard distance of presence pattern across
    # genomes with KNOWN motility (curated only, to keep the clustering
    # itself independent of the lower-confidence predictions).
    known_orgs = [o for o in organisms if org_status[o]["source"] in TRUSTED_SOURCES]
    gene_to_orgset = {g: {o for o in known_orgs if g in presence[o]} for g in all_genes}
    dist = {}
    for i, j in itertools.combinations(range(len(all_genes)), 2):
        dist[(i, j)] = jaccard_distance(gene_to_orgset[all_genes[i]], gene_to_orgset[all_genes[j]])
    gene_order_clustered = upgma_leaf_order(all_genes, dist)

    output = {
        "genes_clustered_order": gene_order_clustered,
        "gene_stats": gene_records,
        "organisms": [
            {"name": o, "status": org_status[o]["status"], "source": org_status[o]["source"]}
            for o in organisms
        ],
        "presence": {o: sorted(presence[o]) for o in organisms},
        "prediction_confidence_threshold": PREDICTION_CONFIDENCE_THRESHOLD,
    }
    with open(OUTPUT_JSON, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Saved: {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
