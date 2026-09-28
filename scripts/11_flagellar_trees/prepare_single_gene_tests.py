#!/usr/bin/env python3
"""
prepare_single_gene_tests.py
Inputs for the single-gene screen (katana/08-10_*.pbs): genes that may have
moved on their own, rather than with the rest of the flagellar gene set.

1. Screen the 37 gene trees against the ML bac120 MARKER tree (not GTDB --
   the marker test showed GTDB is wrong or uncertain at several spots) with
   screen_hgt_candidates.py's logic: branches with UFBoot >= 95 incompatible
   with the marker tree, reduced to minimal groups.
2. Keep single-gene cases worth testing:
     - not also in the concatenated tree, and in <= MAX_GENES gene trees
       (recurring conflicts were handled by the multi-gene AU tests);
     - gene has >= MIN_INFORMATIVE parsimony-informative sites (short genes
       lack power, section 2 of RESULTS.md);
     - the group's smallest enclosing marker-tree clade has >= MIN_EXTRA more
       genomes than the group (not a tip reshuffle).
3. Per case, two tests (run on Katana):
     gene test   -- the gene alignment: unconstrained gene tree vs best tree
                    forced to keep the marker-tree branch(es) the group
                    contradicts (genes/<case>.constraint.nwk);
     marker test -- the marker alignment: unconstrained marker tree vs best
                    tree forcing the gene's grouping (one per distinct group,
                    marker/<group>.constraint.nwk).
   A case is supported only if both reject.

Run from the project root (needs markers/free/markers.treefile copied back
from Katana):  python scripts/11_flagellar_trees/prepare_single_gene_tests.py
Output (gene_order/flagellar_trees/single/): cases.tsv, groups.tsv,
genes/<case>.constraint.nwk, marker/<group>.constraint.nwk
"""

import contextlib
import io
import json
import os
import re
import sys

from ete3 import Tree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import screen_hgt_candidates as screen
from prepare_au_tests import aln_taxa, constraint_tree

D = "gene_order/flagellar_trees"
OUT = f"{D}/single"
MARKER_TREE = f"{D}/markers/free/markers.treefile"

MAX_GENES = 3
MIN_INFORMATIVE = 150
MIN_EXTRA = 3


def main():
    os.makedirs(f"{OUT}/genes", exist_ok=True)
    os.makedirs(f"{OUT}/marker", exist_ok=True)
    marker = Tree(MARKER_TREE, format=0)
    # IQ-TREE roots arbitrarily; root like GTDB (Desulfurellia) so the
    # clade-span filter measures the same thing as in the GTDB screen
    org = {l.split("\t")[0]: l.split("\t")[1] for l in open(f"{D}/taxa.tsv").read().splitlines()[1:]}
    outgroup = [t for t in marker.get_leaf_names() if org[t].split()[0] in ("Hippea", "Desulfurella")]
    anc = marker.get_common_ancestor(outgroup)
    marker.set_outgroup(anc if anc is not marker else marker.get_midpoint_outgroup())
    marker.write(outfile=f"{OUT}/marker_tree_topology.nwk", format=5)
    taxa = marker.get_leaf_names()

    # 1. screen gene trees against the marker tree
    screen.GTDB_TREE = f"{OUT}/marker_tree_topology.nwk"
    screen.GTDB_CF_BRANCH = screen.GTDB_CF_STAT = f"{OUT}/none"
    screen.OUT_TSV, screen.OUT_JSON = f"{OUT}/screen_vs_markers.tsv", f"{OUT}/screen_vs_markers.json"
    with contextlib.redirect_stdout(io.StringIO()):
        screen.main()
    conflicts = json.load(open(screen.OUT_JSON))["candidates"]

    # 2. filter
    informative, model = {}, {}
    for g in open(f"{D}/genes.txt").read().split():
        txt = open(f"{D}/gene_trees/{g}.iqtree").read()
        informative[g] = int(re.search(r"Number of parsimony informative sites: (\d+)", txt).group(1))
        model[g] = re.search(r"Best-fit model according to BIC: (\S+)", txt).group(1)

    cases, groups = [], {}
    for c in conflicts:
        gene = c["tree"]
        trees = set([gene] + c["also_in"])
        if gene == "CONCAT" or "CONCAT" in trees or len(trees) > MAX_GENES:
            continue
        if informative[gene] < MIN_INFORMATIVE or c["gtdb_span"] - c["group_size"] < MIN_EXTRA:
            continue
        key = tuple(sorted(c["group"]))
        if key not in groups:
            groups[key] = f"g{len(groups) + 1:03d}"
        cases.append((gene, key, c))

    # 3. constraints
    with open(f"{OUT}/groups.tsv", "w") as f:
        f.write("group_id\tn_genomes\taccessions\torganisms\n")
        for key, gid in groups.items():
            rest = sorted(set(taxa) - set(key))
            with open(f"{OUT}/marker/{gid}.constraint.nwk", "w") as out:
                out.write(f"(({','.join(key)}),{','.join(rest)});\n")
            orgs = next(c["organisms"] for g, k, c in cases if k == key)
            f.write(f"{gid}\t{len(key)}\t{','.join(key)}\t{'; '.join(orgs)}\n")

    n = 0
    with open(f"{OUT}/cases.tsv", "w") as f:
        f.write("case_id\tgene\tmodel\tgroup_id\tufboot\tinformative_sites\tn_gene_trees\tmarker_span\t"
                "forced_branches\torganisms\n")
        for gene, key, c in cases:
            gene_taxa = aln_taxa(f"{D}/trimmed/{gene}.faa")
            nwk, forced = constraint_tree(Tree(f"{OUT}/marker_tree_topology.nwk", format=1), gene_taxa, list(key))
            if forced == 0:
                continue
            n += 1
            cid = f"c{n:03d}"
            with open(f"{OUT}/genes/{cid}.constraint.nwk", "w") as out:
                out.write(nwk + "\n")
            f.write(f"{cid}\t{gene}\t{model[gene]}\t{groups[key]}\t{c['ufboot']}\t{informative[gene]}\t"
                    f"{1 + len(c['also_in'])}\t{c['gtdb_span']}\t{forced}\t{'; '.join(c['organisms'])}\n")

    print(f"{len(conflicts)} conflicts with the marker tree; {n} single-gene cases in "
          f"{len({g for g, _, _ in cases})} genes, {len(groups)} distinct groups")
    print(f"Wrote {OUT}/")


if __name__ == "__main__":
    main()
