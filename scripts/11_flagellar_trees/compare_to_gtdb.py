#!/usr/bin/env python3
"""
compare_to_gtdb.py
How closely does the concatenated flagellar tree follow the GTDB species
tree? Run after katana/02_concat_concordance.pbs (it runs this as its last
step), or locally once the Katana outputs have been copied back.

Comparisons (unrooted, on the shared taxon set):
  - shared bipartitions (branches) between the flagellar and GTDB trees,
    also restricted to flagellar branches with UFBoot >= 95;
  - a random-shuffle null for that count: flagellar tip labels permuted
    N_SHUFFLES times (tree shape kept), shared branches recounted -- tells
    you how many shared branches to expect by chance for this tree shape;
  - Robinson-Foulds distance, raw and normalised to [0, 1];
  - Mantel test between the two patristic distance matrices
    (mantel_test from 04_synteny/synteny_analysis.py);
  - per-gene trees: normalised RF against GTDB and against the flagellar tree;
  - gCF / sCF summaries from IQ-TREE's concordance runs on both trees.

Run from the project root:  python scripts/11_flagellar_trees/compare_to_gtdb.py
Output: gene_order/flagellar_trees/comparison.json
"""

import glob
import json
import os
import statistics
import sys

import numpy as np
from ete3 import Tree

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "04_synteny"))
from synteny_analysis import mantel_test

D = "gene_order/flagellar_trees"
FLAG_TREE = f"{D}/concat/concat.treefile"
GTDB_TREE = f"{D}/gtdb_ref_pruned.nwk"
GENE_TREES = f"{D}/gene_trees/*.treefile"
CF_STATS = {"flagellar_tree": [f"{D}/concordance/flagellar_gcf.cf.stat",
                               f"{D}/concordance/flagellar_scfl.cf.stat"],
            "gtdb_tree": [f"{D}/concordance/gtdb_gcf.cf.stat",
                          f"{D}/concordance/gtdb_scfl.cf.stat"]}
OUTPUT_JSON = f"{D}/comparison.json"

N_SHUFFLES = 1000
N_MANTEL_PERM = 9999
UFBOOT_STRONG = 95
RNG_SEED = 0


def splits(tree, index, support_min=None):
    """Non-trivial unrooted bipartitions as bitmasks over `index` (taxon -> bit),
    canonicalised to the side without bit 0. With support_min, only branches
    whose support (UFBoot, read from IQ-TREE's node labels) is >= support_min."""
    full = (1 << len(index)) - 1
    out = set()
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node.add_feature("mask", 1 << index[node.name])
            continue
        node.add_feature("mask", sum(c.mask for c in node.children))
        if node.is_root() or (support_min is not None and node.support < support_min):
            continue
        size = bin(node.mask).count("1")
        if 2 <= size <= len(index) - 2:
            out.add(node.mask ^ full if node.mask & 1 else node.mask)
    return out


def rf(a, b, n):
    d = len(a ^ b)
    return d, d / (2 * (n - 3))


def shuffle_null(tree, index, reference, rng):
    """Shared-branch counts for N_SHUFFLES random relabellings of `tree`'s tips."""
    names = tree.get_leaf_names()
    counts = []
    for _ in range(N_SHUFFLES):
        perm = rng.permutation(names)
        t = tree.copy()
        for leaf, new in zip(t.iter_leaves(), perm):
            leaf.name = new
        counts.append(len(splits(t, index) & reference))
    return np.array(counts)


def patristic(tree, taxa):
    """Leaf-to-leaf path length matrix, ordered as `taxa`."""
    depth, anc = {}, {}
    for node in tree.traverse("preorder"):
        depth[node] = 0.0 if node.is_root() else depth[node.up] + node.dist
    leaves = {l.name: l for l in tree.iter_leaves()}
    for name in taxa:
        anc[name] = set(leaves[name].get_ancestors())
    n = len(taxa)
    m = np.zeros((n, n))
    for i in range(n):
        a = leaves[taxa[i]]
        for j in range(i + 1, n):
            node = leaves[taxa[j]]
            while node not in anc[taxa[i]]:
                node = node.up
            m[i, j] = m[j, i] = depth[a] + depth[leaves[taxa[j]]] - 2 * depth[node]
    return m


def restrict(tree, taxa):
    t = tree.copy()
    t.prune(list(taxa), preserve_branch_length=True)
    return t


def cf_summary(paths):
    """gCF and sCF come from separate IQ-TREE runs (--gcf, --scfl); merge them."""
    out = {}
    for path in paths:
        if not os.path.exists(path):
            continue
        with open(path) as f:
            rows = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#") and l.strip()]
        out.update(cf_columns(rows[0], rows[1:]))
    return out or None


def cf_columns(header, rows):
    col = {h: i for i, h in enumerate(header)}
    out = {"n_branches": len(rows)}
    for key in ("gCF", "sCF"):
        if key not in col:
            continue
        vals = [float(r[col[key]]) for r in rows if r[col[key]] not in ("NA", "")]
        out[key] = {"mean": round(statistics.mean(vals), 2),
                    "median": round(statistics.median(vals), 2),
                    "branches_zero": sum(v == 0 for v in vals),
                    "branches_below_10": sum(v < 10 for v in vals),
                    "branches_50_plus": sum(v >= 50 for v in vals)}
    return out


def main():
    flag = Tree(FLAG_TREE, format=0)
    gtdb = Tree(GTDB_TREE, format=1)
    taxa = sorted(set(flag.get_leaf_names()) & set(gtdb.get_leaf_names()))
    flag, gtdb = restrict(flag, taxa), restrict(gtdb, taxa)
    index = {t: i for i, t in enumerate(taxa)}
    n = len(taxa)

    s_flag, s_gtdb = splits(flag, index), splits(gtdb, index)
    s_flag_strong = splits(flag, index, UFBOOT_STRONG)
    shared = len(s_flag & s_gtdb)
    rf_raw, rf_norm = rf(s_flag, s_gtdb, n)

    rng = np.random.default_rng(RNG_SEED)
    null = shuffle_null(flag, index, s_gtdb, rng)
    p_shared = (int((null >= shared).sum()) + 1) / (N_SHUFFLES + 1)

    print(f"{n} taxa; flagellar tree {len(s_flag)} branches, GTDB {len(s_gtdb)}")
    print(f"Shared branches: {shared}/{len(s_gtdb)} ({100 * shared / len(s_gtdb):.1f}%)")
    print(f"  of which UFBoot >= {UFBOOT_STRONG}: {len(s_flag_strong & s_gtdb)}/{len(s_flag_strong)}")
    print(f"  shuffle null: mean {null.mean():.2f}, max {null.max()}, p = {p_shared:.4f}")
    print(f"RF = {rf_raw} (normalised {rf_norm:.3f})")

    d_flag, d_gtdb = patristic(flag, taxa), patristic(gtdb, taxa)
    r, p_mantel, n_pairs = mantel_test(d_flag, d_gtdb, np.ones((n, n), bool),
                                       n_perm=N_MANTEL_PERM, seed=RNG_SEED)
    print(f"Mantel r = {r:.4f}, p = {p_mantel:.4f} ({n_pairs} pairs, {N_MANTEL_PERM} permutations)")

    per_gene = {}
    for path in sorted(glob.glob(GENE_TREES)):
        gene = os.path.basename(path).removesuffix(".treefile")
        gt = Tree(path, format=0)
        gtaxa = sorted(set(gt.get_leaf_names()) & set(taxa))
        if len(gtaxa) < 4:
            continue
        gi = {t: i for i, t in enumerate(gtaxa)}
        s_gene = splits(restrict(gt, gtaxa), gi)
        s_g = splits(restrict(gtdb, gtaxa), gi)
        s_f = splits(restrict(flag, gtaxa), gi)
        per_gene[gene] = {
            "n_taxa": len(gtaxa),
            "rf_norm_vs_gtdb": round(rf(s_gene, s_g, len(gtaxa))[1], 4),
            "rf_norm_vs_flagellar": round(rf(s_gene, s_f, len(gtaxa))[1], 4),
            "shared_with_gtdb_frac": round(len(s_gene & s_g) / len(s_g), 4),
        }
    if per_gene:
        print(f"\nPer-gene normalised RF vs GTDB ({len(per_gene)} genes, most discordant first):")
        for gene, v in sorted(per_gene.items(), key=lambda kv: -kv[1]["rf_norm_vs_gtdb"]):
            print(f"  {gene:6s} n={v['n_taxa']:3d}  vs GTDB {v['rf_norm_vs_gtdb']:.3f}  "
                  f"vs flagellar {v['rf_norm_vs_flagellar']:.3f}")

    concordance = {k: cf_summary(p) for k, p in CF_STATS.items()}
    for k, v in concordance.items():
        if v:
            print(f"\nConcordance on {k}: {json.dumps(v)}")

    json.dump({
        "n_taxa": n,
        "branches": {"flagellar": len(s_flag), "gtdb": len(s_gtdb),
                     "flagellar_ufboot_95": len(s_flag_strong)},
        "shared_branches": shared,
        "shared_branches_frac_of_gtdb": round(shared / len(s_gtdb), 4),
        "shared_branches_ufboot_95": len(s_flag_strong & s_gtdb),
        "shuffle_null": {"n_shuffles": N_SHUFFLES, "mean": round(float(null.mean()), 3),
                         "sd": round(float(null.std()), 3), "max": int(null.max()),
                         "p_value": round(p_shared, 5)},
        "robinson_foulds": {"raw": rf_raw, "normalised": round(rf_norm, 4)},
        "mantel": {"r": round(r, 4), "p_value": round(p_mantel, 5), "n_pairs": n_pairs,
                   "n_permutations": N_MANTEL_PERM, "distance": "patristic, Pearson"},
        "per_gene": per_gene,
        "concordance": concordance,
    }, open(OUTPUT_JSON, "w"), indent=2)
    print(f"\nWrote {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
