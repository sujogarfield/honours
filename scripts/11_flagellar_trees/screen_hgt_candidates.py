#!/usr/bin/env python3
"""
screen_hgt_candidates.py
First-pass HGT screen: find places where a gene tree (or the concatenated
flagellar tree) STRONGLY contradicts the GTDB species tree, as opposed to
merely failing to resolve it. Run after katana/02_concat_concordance.pbs
(it runs this as its last step), or locally once the outputs are copied back.

A candidate is a branch in a gene tree that
  - has UFBoot >= UFBOOT_MIN, and
  - is incompatible with at least one GTDB branch (the two bipartitions cannot
    both be true on one tree -- not just "missing from GTDB", which would only
    mean one tree is less resolved than the other).
A single moved lineage makes every branch along its new path conflict, so only
the MINIMAL conflicting groups are kept (no smaller conflicting group nested
inside them): that is "these genomes group together in the gene tree but not
in GTDB", usually the moved lineage plus its new sister.

Each candidate is annotated with
  - gtdb_span: size of the smallest GTDB clade containing the group. Group of 3
    whose GTDB span is 4 is a shallow reshuffle; span 90 is a deep conflict;
  - gtdb_extra: the genomes in that GTDB clade that are NOT in the group (listed
    when <= 5). A moved lineage shows up twice: once at its new position
    (moved genome + new sister, large span) and once as the gap it left behind
    (its old relatives, gtdb_extra = the genome that moved out);
  - min_gcf_conflicting: lowest gCF among the GTDB branches it contradicts
    (from concordance/gtdb.cf.*). Low gCF there means the flagellar genes as a
    whole don't support that species branch either;
  - n_sites: trimmed alignment length -- short genes carry little signal;
  - recurrence: how many genes/trees show the same group. The same group in
    several genes suggests a co-transferred operon, not single-gene noise.
Separately lists GTDB branches with gCF <= GCF_LOW (species branches that
almost no flagellar gene supports).

This is a SCREEN, not a test: candidates still need a targeted AU test,
paralogy / model checks and a donor search before being called HGT.

Run from the project root:  python scripts/11_flagellar_trees/screen_hgt_candidates.py
Output:
  gene_order/flagellar_trees/hgt_candidates.tsv   -- one row per candidate, ranked
  gene_order/flagellar_trees/hgt_screen.json      -- candidates + low-gCF GTDB branches
"""

import glob
import json
import os
import re
from collections import defaultdict

from ete3 import Tree

D = "gene_order/flagellar_trees"
GTDB_TREE = f"{D}/gtdb_ref_pruned.nwk"
FLAG_TREE = f"{D}/concat/concat.treefile"
GENE_TREES = f"{D}/gene_trees/*.treefile"
GTDB_CF_BRANCH = f"{D}/concordance/gtdb.cf.branch"
GTDB_CF_STAT = f"{D}/concordance/gtdb.cf.stat"
TAXA_TSV = f"{D}/taxa.tsv"
OUT_TSV = f"{D}/hgt_candidates.tsv"
OUT_JSON = f"{D}/hgt_screen.json"

UFBOOT_MIN = 95
GCF_LOW = 5.0


def bitmask_splits(tree, index, with_support=False):
    """{canonical mask: support} for non-trivial unrooted splits. Canonical =
    side without bit 0 (so identical splits compare equal across trees)."""
    full = (1 << len(index)) - 1
    out = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node.add_feature("mask", 1 << index[node.name])
            continue
        node.add_feature("mask", sum(c.mask for c in node.children))
        if node.is_root():
            continue
        size = bin(node.mask).count("1")
        if 2 <= size <= len(index) - 2:
            key = node.mask ^ full if node.mask & 1 else node.mask
            out[key] = node.support if with_support else node.name
    return out


def compatible(a, b, full):
    """Two bipartitions are compatible iff some pair of their sides is disjoint."""
    ac, bc = full ^ a, full ^ b
    return not (a & b) or not (a & bc) or not (ac & b) or not (ac & bc)


def smaller_side(mask, n):
    full = (1 << n) - 1
    other = full ^ mask
    return mask if bin(mask).count("1") <= bin(other).count("1") else other


def members(mask, taxa):
    return [taxa[i] for i in range(len(taxa)) if mask >> i & 1]


def gtdb_clade(gtdb, names):
    return gtdb.get_common_ancestor(names).get_leaf_names() if len(names) > 1 else list(names)


def load_gtdb_gcf(taxa, index):
    """GTDB split mask -> gCF, by matching IQ-TREE's branch IDs (internal node
    labels in gtdb.cf.branch) to the ID column of gtdb.cf.stat."""
    if not (os.path.exists(GTDB_CF_BRANCH) and os.path.exists(GTDB_CF_STAT)):
        return {}
    with open(GTDB_CF_STAT) as f:
        rows = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#") and l.strip()]
    col = {h: i for i, h in enumerate(rows[0])}
    gcf_by_id = {r[col["ID"]]: float(r[col["gCF"]]) for r in rows[1:] if r[col["gCF"]] != "NA"}
    t = Tree(GTDB_CF_BRANCH, format=1)
    t.prune([x for x in t.get_leaf_names() if x in index], preserve_branch_length=True)
    return {m: gcf_by_id[bid] for m, bid in bitmask_splits(t, index).items() if bid in gcf_by_id}


def n_sites(treefile):
    iq = treefile.removesuffix(".treefile") + ".iqtree"
    if os.path.exists(iq):
        m = re.search(r"Input data: \d+ sequences with (\d+)", open(iq).read())
        if m:
            return int(m.group(1))
    return None


def main():
    organism = {}
    if os.path.exists(TAXA_TSV):
        for line in open(TAXA_TSV).read().splitlines()[1:]:
            acc, org = line.split("\t")
            organism[acc] = org

    gtdb_full = Tree(GTDB_TREE, format=1)
    all_taxa = sorted(gtdb_full.get_leaf_names())
    all_index = {t: i for i, t in enumerate(all_taxa)}
    gcf_of = load_gtdb_gcf(all_taxa, all_index)

    trees = {os.path.basename(p).removesuffix(".treefile"): p for p in sorted(glob.glob(GENE_TREES))}
    if os.path.exists(FLAG_TREE):
        trees["CONCAT"] = FLAG_TREE

    candidates = []
    for name, path in trees.items():
        gt = Tree(path, format=0)
        taxa = sorted(set(gt.get_leaf_names()) & set(all_taxa))
        if len(taxa) < 5:
            continue
        index = {t: i for i, t in enumerate(taxa)}
        full = (1 << len(taxa)) - 1
        gt.prune(taxa, preserve_branch_length=True)
        g = gtdb_full.copy()
        g.prune(taxa, preserve_branch_length=True)

        gene_splits = bitmask_splits(gt, index, with_support=True)
        gtdb_splits = list(bitmask_splits(g, index))
        conflicts = {}
        for s, sup in gene_splits.items():
            if sup < UFBOOT_MIN:
                continue
            bad = [x for x in gtdb_splits if not compatible(s, x, full)]
            if bad:
                conflicts[smaller_side(s, len(taxa))] = (sup, bad)

        # keep minimal groups: no other conflicting group strictly inside
        for side, (sup, bad) in conflicts.items():
            if any(o != side and o & side == o for o in conflicts):
                continue
            names = members(side, taxa)
            clade = gtdb_clade(g, names)
            extra = sorted(set(clade) - set(names))
            to_full = lambda m: sum(1 << all_index[t] for t in members(m, taxa))
            gcfs = []
            for x in bad:
                fm = to_full(x)
                key = fm ^ ((1 << len(all_taxa)) - 1) if fm & 1 else fm
                if key in gcf_of:
                    gcfs.append(gcf_of[key])
            candidates.append({
                "tree": name,
                "n_sites": n_sites(path),
                "ufboot": sup,
                "group_size": len(names),
                "gtdb_span": len(clade),
                "gtdb_extra": [organism.get(t, t) for t in extra] if len(extra) <= 5 else f"{len(extra)} genomes",
                "n_gtdb_branches_contradicted": len(bad),
                "min_gcf_conflicting": min(gcfs) if gcfs else None,
                "group": names,
                "organisms": [organism.get(t, t) for t in names],
            })

    recur = defaultdict(set)
    for c in candidates:
        recur[frozenset(c["group"])].add(c["tree"])
    for c in candidates:
        c["recurrence"] = len(recur[frozenset(c["group"])])
        c["also_in"] = sorted(recur[frozenset(c["group"])] - {c["tree"]})

    # deep, recurrent, strongly supported conflicts first
    candidates.sort(key=lambda c: (-c["recurrence"], -(c["gtdb_span"] - c["group_size"]), -c["ufboot"]))

    with open(OUT_TSV, "w") as f:
        f.write("tree\tn_sites\tufboot\tgroup_size\tgtdb_span\tn_gtdb_branches_contradicted\t"
                "min_gcf_conflicting\trecurrence\talso_in\torganisms\tgtdb_extra\taccessions\n")
        for c in candidates:
            f.write("\t".join(str(x) for x in [
                c["tree"], c["n_sites"], c["ufboot"], c["group_size"], c["gtdb_span"],
                c["n_gtdb_branches_contradicted"], c["min_gcf_conflicting"], c["recurrence"],
                ",".join(c["also_in"]), "; ".join(c["organisms"]),
                c["gtdb_extra"] if isinstance(c["gtdb_extra"], str) else "; ".join(c["gtdb_extra"]),
                ",".join(c["group"])]) + "\n")

    low_gcf = []
    if gcf_of:
        for m, v in gcf_of.items():
            if v <= GCF_LOW:
                side = smaller_side(m, len(all_taxa))
                names = members(side, all_taxa)
                low_gcf.append({"gcf": v, "clade_size": len(names),
                                "organisms": [organism.get(t, t) for t in names]})
        low_gcf.sort(key=lambda x: (x["gcf"], -x["clade_size"]))

    json.dump({"ufboot_min": UFBOOT_MIN, "gcf_low": GCF_LOW,
               "n_trees_screened": len(trees), "candidates": candidates,
               "gtdb_branches_low_gcf": low_gcf}, open(OUT_JSON, "w"), indent=2)

    print(f"Screened {len(trees)} trees: {len(candidates)} strongly supported conflicts "
          f"(UFBoot >= {UFBOOT_MIN}) in {len({c['tree'] for c in candidates})} trees")
    for c in candidates[:15]:
        orgs = "; ".join(c["organisms"][:4]) + (" ..." if len(c["organisms"]) > 4 else "")
        print(f"  {c['tree']:6s} UFBoot {c['ufboot']:5.0f}  group {c['group_size']:3d} "
              f"(GTDB span {c['gtdb_span']:3d})  in {c['recurrence']} tree(s)  {orgs}")
    if gcf_of:
        print(f"{len(low_gcf)} GTDB branches with gCF <= {GCF_LOW}")
    print(f"Wrote {OUT_TSV} and {OUT_JSON}")


if __name__ == "__main__":
    main()
