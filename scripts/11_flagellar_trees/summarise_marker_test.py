#!/usr/bin/env python3
"""
summarise_marker_test.py
Collects the species-marker test (katana/06_marker_au.pbs): at each AU
candidate, does the genome backbone (GTDB R232 bac120 markers) agree with GTDB
or with the flagellar genes?

AU tree set per candidate: 1 = unconstrained marker tree, 2 = flagellar
grouping forced, 3 = GTDB branch(es) forced.
  p_flag < 0.05 (Holm across candidates): the markers significantly REJECT the
      flagellar grouping -> the flagellar genes and the genome backbone have
      different histories here -> non-vertical inheritance (HGT/recombination)
      candidate.
  p_flag >= 0.05: the markers cannot reject the flagellar grouping -> the
      species tree is uncertain here; not evidence for HGT.
Also reports whether the unconstrained marker tree contains the GTDB
branch(es) and/or the flagellar grouping, with its UFBoot support.
p-values use the same bp-RELL = 0 rule as summarise_au.py (reported < 1e-4).

Run from the project root:  python scripts/11_flagellar_trees/summarise_marker_test.py
Output: gene_order/flagellar_trees/markers/marker_test.tsv
"""

import os
import sys

from ete3 import Tree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from screen_hgt_candidates import bitmask_splits
from summarise_au import holm, parse_user_trees, N_RELL, PRIMARY

M = "gene_order/flagellar_trees/markers"
ALPHA = 0.05


def p_of(t, i):
    if t[i].get("bp-RELL", 1) == 0:
        return 1 / N_RELL, f"<{1 / N_RELL:g}"
    return t[i]["p-AU"], f"{t[i]['p-AU']:.3g}"


def canon(mask, full):
    return mask ^ full if mask & 1 else mask


def splits_with_support(tree, index):
    full = (1 << len(index)) - 1
    out = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node.add_feature("mask", 1 << index[node.name])
            continue
        node.add_feature("mask", sum(c.mask for c in node.children))
        if not node.is_root():
            out[canon(node.mask, full)] = node.support
    return out


def main():
    free = Tree(f"{M}/free/markers.treefile", format=0)
    taxa = sorted(free.get_leaf_names())
    index = {t: i for i, t in enumerate(taxa)}
    full = (1 << len(taxa)) - 1
    free_splits = splits_with_support(free, index)

    rows = []
    for cand in open(f"{M}/candidates.txt").read().split():
        group = [l.split("\t")[0] for l in open(f"{M}/{cand}/group.txt").read().splitlines()]
        gmask = canon(sum(1 << index[a] for a in group), full)
        gtdb_c = Tree(open(f"{M}/{cand}/gtdb_constraint.nwk").read(), format=9)
        gtdb_masks = list(bitmask_splits(gtdb_c, index))
        has_flag = gmask in free_splits
        has_gtdb = all(m in free_splits for m in gtdb_masks)
        sup_gtdb = min(free_splits[m] for m in gtdb_masks) if has_gtdb else None
        path = f"{M}/{cand}/run/au.iqtree"
        t = parse_user_trees(path) if os.path.exists(path) else None
        if not t or not all(i in t for i in (1, 2, 3)):
            rows.append([cand, has_gtdb, sup_gtdb, has_flag, "NA", "NA", None, "NA", "not finished"])
            continue
        pf, pf_s = p_of(t, 2)
        pg, pg_s = p_of(t, 3)
        rows.append([cand, has_gtdb, sup_gtdb, has_flag, pf_s, round(t[2]["logL"] - t[1]["logL"], 2),
                     pf, pg_s, None])

    done = [r for r in rows if r[6] is not None]
    for r, adj in zip(done, holm([r[6] for r in done])):
        r[6] = round(adj, 5)
        if adj < ALPHA:
            r[8] = "markers reject flagellar grouping -> non-vertical inheritance candidate"
        elif r[3]:
            r[8] = "marker tree itself has the flagellar grouping -> species tree uncertain"
        else:
            r[8] = "markers cannot reject flagellar grouping -> species tree uncertain"

    prim = [r for r in done if r[0] in PRIMARY]
    hp = dict(zip([r[0] for r in prim], holm([float(r[4].lstrip("<")) for r in prim])))
    for r in rows:
        r.append(round(hp[r[0]], 5) if r[0] in hp else "")

    with open(f"{M}/marker_test.tsv", "w") as f:
        f.write("candidate\tmarker_tree_has_gtdb_branch\tmarker_ufboot_gtdb_branch\tmarker_tree_has_flag_group\t"
                "p_flag\tdelta_logL_flag\tp_flag_holm\tp_gtdb\tverdict\tp_flag_holm_primary\n")
        for r in rows:
            f.write("\t".join("" if x is None else str(x) for x in r) + "\n")

    print(f"{'candidate':30s} {'GTDB br in marker tree':>22s} {'flag grp':>8s} {'p_flag':>8s} {'Holm':>8s}  verdict")
    for r in rows:
        g = f"{'yes' if r[1] else 'no'}" + (f" (UFB {r[2]:.0f})" if r[2] is not None else "")
        h = f"{r[6]:.4f}" if isinstance(r[6], float) else "NA"
        print(f"{r[0]:30s} {g:>22s} {'yes' if r[3] else 'no':>8s} {r[4]:>8s} {h:>8s}  {r[8]}")
    print(f"Wrote {M}/marker_test.tsv")


if __name__ == "__main__":
    main()
