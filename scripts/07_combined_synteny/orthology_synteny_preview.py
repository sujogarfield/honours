#!/usr/bin/env python3
"""
orthology_synteny_preview.py
Identical to synteny_analysis.py, but run against the combined gene set
(string match + eggNOG-confirmed/reclassified orthology candidates +
confirmed maf candidates -- see build_combined_gene_order.py for the
strict OrthoFinder+eggNOG-agreement inclusion policy) instead of the
string-matched baseline alone. Deliberately kept as a separate sibling
script (like synteny_analysis_full_contigs.py already is) rather than
folded into synteny_analysis.py itself, so the original string-match-only
result stays independently reproducible. Compare its Mantel r/p and
per-organism agreement scores against synteny_analysis.json's to see what
the orthology-recovered genes changed.

Run AFTER build_combined_gene_order.py.
"""

import json
import itertools
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from ete3 import Tree

DEFAULT_GENE_ORDER_JSON = "gene_order/gene_order_combined.json"
NAMED_NWK = "images_and_newick/thesis_b/gtdb_ref_tree.nwk"
DEFAULT_OUTPUT_JSON = "gene_order/synteny_analysis_combined.json"

MIN_SHARED_GENES = 4
N_PERMUTATIONS = 999
RNG_SEED = 0


def load_gene_order(path: str) -> dict:
    """Returns {organism_name: [gene_records]} from gene_order.json."""
    with open(path) as f:
        raw = json.load(f)
    by_organism = {}
    for data in raw.values():
        by_organism.setdefault(data["organism"], []).extend(data["genes"])
    return by_organism


def primary_contig_order(records: list[dict]) -> tuple[list[str], int]:
    """Ordered gene list for the contig with the most flagellar genes.
    Duplicated gene tags on that contig are dropped (ambiguous order mapping).
    Returns (ordered_gene_names, n_genes_on_primary_contig)."""
    by_contig = {}
    for r in records:
        by_contig.setdefault(r["contig"], []).append(r)

    if not by_contig:
        return [], 0

    primary = max(by_contig.values(), key=len)
    counts = Counter(r["gene"].lower() for r in primary)
    ordered = sorted(
        (r for r in primary if counts[r["gene"].lower()] == 1),
        key=lambda r: r["start"],
    )
    return [r["gene"].lower() for r in ordered], len(primary)


def breakpoint_distance(order_a: list[str], order_b: list[str]):
    """Normalised breakpoint distance over genes shared by both orders.
    Direction-agnostic — a whole inverted block still counts as conserved,
    since strand is treated as a display detail rather than part of this
    metric. Returns (distance_or_None, n_shared_genes)."""
    shared = set(order_a) & set(order_b)
    if len(shared) < MIN_SHARED_GENES:
        return None, len(shared)

    sub_a = [g for g in order_a if g in shared]
    sub_b = [g for g in order_b if g in shared]

    adj_a = {frozenset((sub_a[i], sub_a[i + 1])) for i in range(len(sub_a) - 1)}
    adj_b = {frozenset((sub_b[i], sub_b[i + 1])) for i in range(len(sub_b) - 1)}

    possible = len(shared) - 1
    conserved = len(adj_a & adj_b)
    return 1 - conserved / possible, len(shared)


def spearman(x, y) -> float:
    """Spearman rank correlation, implemented by hand (no scipy in the venv)."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    if rx.std() == 0 or ry.std() == 0:
        return 0.0
    return float(np.corrcoef(rx, ry)[0, 1])


def mantel_test(dist_synt, dist_phylo, valid_mask, n_perm=N_PERMUTATIONS, seed=RNG_SEED):
    """Permutation-based Mantel test between two square distance matrices,
    restricted to valid_mask pairs. Returns (r, p_value, n_pairs_used)."""
    n = dist_synt.shape[0]
    iu = np.triu_indices(n, k=1)
    mask = valid_mask[iu]

    a = dist_synt[iu][mask]
    b = dist_phylo[iu][mask]
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return None, None, len(a)

    observed_r = float(np.corrcoef(a, b)[0, 1])

    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(n)
        perm_dist = dist_synt[np.ix_(perm, perm)]
        pb = perm_dist[iu][mask]
        r = float(np.corrcoef(pb, b)[0, 1]) if pb.std() > 0 else 0.0
        if abs(r) >= abs(observed_r):
            count += 1

    p_value = (count + 1) / (n_perm + 1)
    return observed_r, p_value, len(a)


def upgma_tree(organisms: list[str], dist: np.ndarray, valid: np.ndarray) -> dict:
    """Average-linkage hierarchical clustering, built by hand (no scipy).
    Cluster-to-cluster distance is the mean of valid pairwise distances
    between members. Tracked incrementally as (sum, count) per cluster pair:
    when clusters merge, the sum/count of valid original leaf-pairs between
    the merged cluster and any other cluster is just the sum of the two
    parents' sums/counts (their member sets are disjoint, so nothing is
    double-counted) — an O(1) update per pair instead of rescanning every
    member pair from scratch on every merge. A cluster pair with zero valid
    distance between any of their members merges last (treated as +inf)."""
    n = len(organisms)
    total = 2 * n - 1
    sums = np.zeros((total, total))
    counts = np.zeros((total, total), dtype=int)
    for i, j in itertools.combinations(range(n), 2):
        if valid[i, j]:
            sums[i, j] = sums[j, i] = dist[i, j]
            counts[i, j] = counts[j, i] = 1

    nodes = {i: {"name": name, "height": 0.0, "children": []} for i, name in enumerate(organisms)}
    active = list(range(n))
    next_id = n

    while len(active) > 1:
        best = None
        for x, y in itertools.combinations(active, 2):
            c = counts[x, y]
            d = sums[x, y] / c if c > 0 else float("inf")
            if best is None or d < best[0]:
                best = (d, x, y)
        d, a, b = best
        new_id = next_id
        next_id += 1
        height = d if np.isfinite(d) else max(nodes[a]["height"], nodes[b]["height"]) + 1
        nodes[new_id] = {"name": "", "height": height, "children": [nodes[a], nodes[b]]}

        for m in active:
            if m in (a, b):
                continue
            s = sums[a, m] + sums[b, m]
            c = counts[a, m] + counts[b, m]
            sums[new_id, m] = sums[m, new_id] = s
            counts[new_id, m] = counts[m, new_id] = c

        active = [x for x in active if x not in (a, b)] + [new_id]

    def strip(node):
        out = {"name": node["name"], "height": round(node["height"], 4)}
        if node["children"]:
            out["children"] = [strip(c) for c in node["children"]]
        return out

    return strip(nodes[active[0]])


def classical_mds(dist: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Classical MDS / PCoA. Missing pairwise distances are imputed with the
    mean of that organism's valid distances (a simple, standard placeholder —
    classical MDS needs a complete matrix). Returns an (n, 2) coordinate array."""
    n = dist.shape[0]
    complete = dist.copy()
    row_means = np.array([
        dist[i, valid[i]].mean() if valid[i].any() else 0.0
        for i in range(n)
    ])
    for i in range(n):
        for j in range(n):
            if i != j and not valid[i, j]:
                complete[i, j] = (row_means[i] + row_means[j]) / 2

    d2 = complete ** 2
    j_mat = np.eye(n) - np.ones((n, n)) / n
    # All inputs here are verified finite; macOS's Accelerate BLAS backend
    # emits spurious divide-by-zero/overflow warnings for this exact
    # scalar-matmul pattern even though the result is correct.
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        b = -0.5 * j_mat @ d2 @ j_mat
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    top2 = np.clip(eigvals[:2], 0, None)
    return eigvecs[:, :2] * np.sqrt(top2)


def main():
    GENE_ORDER_JSON = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_GENE_ORDER_JSON
    OUTPUT_JSON = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_OUTPUT_JSON

    Path("gene_order").mkdir(exist_ok=True)

    print(f"Loading gene order from {GENE_ORDER_JSON} …")
    gene_order = load_gene_order(GENE_ORDER_JSON)

    print(f"Loading tree from {NAMED_NWK} …")
    tree = Tree(NAMED_NWK, format=1, quoted_node_names=True)
    leaf_nodes = {leaf.name: leaf for leaf in tree.get_leaves()}

    organisms = sorted(name for name in gene_order if name in leaf_nodes)
    skipped = sorted(name for name in gene_order if name not in leaf_nodes)
    print(f"{len(organisms)} organisms present in both gene order data and tree")
    if skipped:
        print(f"{len(skipped)} organisms with gene data skipped (not in pruned tree)")

    orders, primary_sizes, total_genes = {}, {}, {}
    for name in organisms:
        recs = gene_order[name]
        ordered, psize = primary_contig_order(recs)
        orders[name] = ordered
        primary_sizes[name] = psize
        total_genes[name] = len(recs)

    n = len(organisms)
    dist_synt = np.zeros((n, n))
    valid = np.zeros((n, n), dtype=bool)

    for i, j in itertools.combinations(range(n), 2):
        d, _ = breakpoint_distance(orders[organisms[i]], orders[organisms[j]])
        if d is not None:
            dist_synt[i, j] = dist_synt[j, i] = d
            valid[i, j] = valid[j, i] = True

    print("Computing phylogenetic (patristic) distances …")
    dist_phylo = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        d = leaf_nodes[organisms[i]].get_distance(leaf_nodes[organisms[j]])
        dist_phylo[i, j] = dist_phylo[j, i] = d

    print("Running Mantel test …")
    r, p, n_pairs = mantel_test(dist_synt, dist_phylo, valid)
    if r is None:
        print("Not enough valid pairs to run a Mantel test.")
    else:
        print(f"Mantel r = {r:.4f}, p = {p:.4f}  (n = {n_pairs} pairs, {N_PERMUTATIONS} permutations)")

    print("Computing per-organism agreement scores …")
    per_organism = {}
    for i, name in enumerate(organisms):
        neighbours = [j for j in range(n) if j != i and valid[i, j]]
        agreement = None
        if len(neighbours) >= 3:
            agreement = spearman(
                [dist_synt[i, j] for j in neighbours],
                [dist_phylo[i, j] for j in neighbours],
            )
        coverage = round(primary_sizes[name] / total_genes[name], 3) if total_genes[name] else None
        per_organism[name] = {
            "primary_contig_gene_count": primary_sizes[name],
            "total_gene_count": total_genes[name],
            "coverage": coverage,
            "n_valid_comparisons": len(neighbours),
            "agreement_score": round(agreement, 4) if agreement is not None else None,
            "gene_order": orders[name],
        }

    print("Building UPGMA synteny dendrogram …")
    upgma = upgma_tree(organisms, dist_synt, valid)

    print("Computing classical MDS coordinates …")
    coords = classical_mds(dist_synt, valid)
    mds_coords = {
        name: {"x": round(float(coords[i, 0]), 4), "y": round(float(coords[i, 1]), 4)}
        for i, name in enumerate(organisms)
    }

    ranked = sorted(
        ((v["agreement_score"], k) for k, v in per_organism.items() if v["agreement_score"] is not None)
    )

    if ranked:
        print("\nMost discordant (gene-order neighbours diverge from phylogenetic neighbours):")
        for score, name in ranked[:5]:
            print(f"  {score:+.3f}  {name}")
        print("Most concordant (gene order tracks phylogeny):")
        for score, name in ranked[-5:][::-1]:
            print(f"  {score:+.3f}  {name}")

    synteny_distance = [
        [round(dist_synt[i, j], 4) if valid[i, j] or i == j else None for j in range(n)]
        for i in range(n)
    ]
    phylogenetic_distance = [[round(dist_phylo[i, j], 6) for j in range(n)] for i in range(n)]

    output = {
        "mantel": {
            "r": r,
            "p_value": p,
            "n_pairs": n_pairs,
            "n_permutations": N_PERMUTATIONS,
            "n_organisms": n,
            "n_skipped_not_in_tree": len(skipped),
            "min_shared_genes": MIN_SHARED_GENES,
        },
        "organisms": per_organism,
        "matrix": {
            "organisms": organisms,
            "synteny_distance": synteny_distance,
            "phylogenetic_distance": phylogenetic_distance,
        },
        "upgma_tree": upgma,
        "mds_coords": mds_coords,
    }
    with open(OUTPUT_JSON, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
