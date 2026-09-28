#!/usr/bin/env python3
"""
summarise_single_gene.py
Collects the single-gene screen (katana/08-09_*.pbs). Per case:
  p_gene   -- the gene alignment rejects the marker-tree arrangement
              (tree 2 = marker branch(es) forced)
  p_marker -- the marker alignment rejects the gene's grouping
              (tree 2 = gene grouping forced; shared by cases with the same group)
Both are Benjamini-Hochberg corrected across the tested cases / groups (an
exploratory screen of ~70 cases). A case is a single-gene non-vertical
inheritance candidate only if BOTH are < 0.05 after correction.
Cases whose group contains a genome that fails IQ-TREE's composition test in
>= COMP_FAIL_MIN genes are labelled "composition caveat" (the thermophile
artefact, RESULTS.md section 3). p-values use the bp-RELL = 0 rule of
summarise_au.py.

Run from the project root:  python scripts/11_flagellar_trees/summarise_single_gene.py
Output: gene_order/flagellar_trees/single/single_gene_results.tsv
"""

import glob
import os
import re
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summarise_au import bh, parse_user_trees, N_RELL

D = "gene_order/flagellar_trees"
S = f"{D}/single"
ALPHA = 0.05
COMP_FAIL_MIN = 5


def p2(path):
    t = parse_user_trees(path) if os.path.exists(path) else None
    if not t or 1 not in t or 2 not in t:
        return None
    if t[2].get("bp-RELL", 1) == 0:
        return 1 / N_RELL
    return t[2]["p-AU"]


def composition_flagged():
    fails = Counter()
    for f in glob.glob(f"{D}/gene_trees/*.log"):
        for line in open(f):
            m = re.match(r"\s*\d+\s+(GCF_\S+)\s+[\d.]+%\s+failed", line)
            if m:
                fails[m.group(1)] += 1
    return {a for a, n in fails.items() if n >= COMP_FAIL_MIN}


def main():
    cases = [l.rstrip("\n").split("\t") for l in open(f"{S}/cases.tsv")][1:]
    groups = {r[0]: r for r in (l.rstrip("\n").split("\t") for l in open(f"{S}/groups.tsv")) if r[0] != "group_id"}
    comp = composition_flagged()

    p_gene = {c[0]: p2(f"{S}/run/{c[0]}/au.iqtree") for c in cases}
    p_marker = {g: p2(f"{S}/run/{g}/au.iqtree") for g in groups}
    gc = [c[0] for c in cases if p_gene[c[0]] is not None]
    gg = [g for g in groups if p_marker[g] is not None]
    q_gene = dict(zip(gc, bh([p_gene[c] for c in gc])))
    q_marker = dict(zip(gg, bh([p_marker[g] for g in gg])))

    rows = []
    for c in cases:
        cid, gene, _, gid = c[:4]
        qg, qm = q_gene.get(cid), q_marker.get(gid)
        caveat = "composition caveat" if set(groups[gid][2].split(",")) & comp else ""
        if qg is None or qm is None:
            verdict = "not finished"
        elif qg < ALPHA and qm < ALPHA:
            verdict = "CANDIDATE: gene and markers disagree"
        elif qg < ALPHA:
            verdict = "gene rejects marker tree, markers can't reject gene grouping"
        else:
            verdict = "gene can't reject marker tree"
        fmt = lambda x: "NA" if x is None else (f"<{1 / N_RELL:g}" if x == 1 / N_RELL else f"{x:.3g}")
        rows.append([cid, gene, c[5], gid, fmt(p_gene[cid]), fmt(qg), fmt(p_marker[gid]), fmt(qm),
                     verdict, caveat, groups[gid][3]])

    with open(f"{S}/single_gene_results.tsv", "w") as f:
        f.write("case\tgene\tinformative_sites\tgroup\tp_gene\tq_gene\tp_marker\tq_marker\tverdict\t"
                "caveat\torganisms\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    done = [r for r in rows if r[8] != "not finished"]
    print(f"{len(done)}/{len(rows)} cases finished")
    for v, n in Counter(r[8] for r in done).most_common():
        print(f"  {n:3d}  {v}")
    cand = [r for r in done if r[8].startswith("CANDIDATE")]
    if cand:
        print("\nCandidates (gene q, marker q):")
        for r in sorted(cand, key=lambda r: r[1]):
            orgs = r[10] if len(r[10]) < 90 else r[10][:87] + "..."
            print(f"  {r[1]:5s} {r[5]:>8s} {r[7]:>8s}  {r[9]:18s} {orgs}")
    print(f"Wrote {S}/single_gene_results.tsv")


if __name__ == "__main__":
    main()
