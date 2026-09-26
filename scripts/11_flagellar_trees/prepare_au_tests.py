#!/usr/bin/env python3
"""
prepare_au_tests.py
Builds the inputs for targeted AU tests (katana/03_au_tests.pbs) on the HGT
shortlist from RESULTS.md. Run locally after the Katana tree outputs have been
copied back (needs concat/concat.treefile, gene_trees/, trimmed/).

For each candidate -- a group of genomes that the flagellar genes place
together against GTDB -- the constraint tree forces ONLY the GTDB branches
that the group contradicts (the GTDB tree with every other internal branch
collapsed). Everything else stays free, so the test asks one question: does
the flagellar alignment significantly reject the GTDB arrangement at this
spot? (A whole-GTDB-tree constraint would be rejected by almost any gene for
unrelated reasons, which says nothing about the candidate.)

Per candidate the AU test is run on the concatenated alignment and on every
gene alignment where the constraint is still informative (enough of the
relevant genomes present, and the GTDB branch still contradicts the group).

Run from the project root:  python scripts/11_flagellar_trees/prepare_au_tests.py
Output (gene_order/flagellar_trees/au/):
  candidates.txt                  candidate names, one per line (PBS array index)
  <cand>/group.txt                accessions + organisms in the group
  <cand>/constraint.nwk           constraint for the concatenated alignment
  <cand>/genes.tsv                gene, best-fit model, for informative genes
  <cand>/genes/<gene>.constraint.nwk
"""

import os
import re
import sys

from ete3 import Tree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from screen_hgt_candidates import bitmask_splits, compatible

D = "gene_order/flagellar_trees"
OUT = f"{D}/au"

# shortlist (RESULTS.md section 4): name -> organisms the flagellar genes group together
CANDIDATES = {
    "poseidonibacter_roscoffensis": ["Poseidonibacter lekithochrous", "Arcobacter roscoffensis"],
    "hydrogenimonas": ["Hydrogenimonas urashimensis", "Hydrogenimonas cancrithermarum"],
    "helicobacter_saguini": ["Helicobacter saguini", "Helicobacter didelphidarum"],
    "campylobacter_pinnipediorum": ["Campylobacter pinnipediorum", "Campylobacter gastrosuis",
                                    "Campylobacter majalis"],
    "campylobacter_mucosalis": ["Campylobacter mucosalis", "Campylobacter suis"],
    "nitrosophilus": ["Nitrosophilus labii", "Nitrosophilus alvini", "Nitrosophilus kaiyonis"],
    "sulfurimonas": ["Sulfurimonas denitrificans", "Sulfurimonas crateris", "Sulfurimonas baltica",
                     "Sulfurimonas marisnigri"],
    # GTDB puts S. tamanense sister to Campylobacter; flagellar genes keep the genus together
    "sulfurospirillum_tamanense": "Sulfurospirillum",
}


def mask_of(names, index):
    return sum(1 << index[n] for n in names)


def constraint_tree(gtdb, taxa, group):
    """GTDB pruned to `taxa`, collapsed to the branches incompatible with `group`.
    Returns (newick, n_forced_branches)."""
    t = gtdb.copy()
    t.prune(list(taxa), preserve_branch_length=False)
    index = {x: i for i, x in enumerate(sorted(taxa))}
    full = (1 << len(index)) - 1
    g = mask_of([x for x in group if x in index], index)
    bitmask_splits(t, index)
    keep = 0
    for node in list(t.traverse("postorder")):
        if node.is_leaf() or node.is_root():
            continue
        size = bin(node.mask).count("1")
        informative = 2 <= size <= len(index) - 2
        if informative and not compatible(g, node.mask, full):
            keep += 1
        else:
            node.delete(prevent_nondicotomic=False)
    return t.write(format=9), keep


def group_is_split(tree, group):
    """Is `group` one side of a branch in the (unrooted) tree?"""
    taxa = tree.get_leaf_names()
    index = {x: i for i, x in enumerate(sorted(taxa))}
    full = (1 << len(index)) - 1
    g = mask_of(group, index)
    splits = bitmask_splits(tree, index)
    return any(s in (g, full ^ g) for s in splits)


def aln_taxa(path):
    return [l[1:].split()[0] for l in open(path) if l.startswith(">")]


def main():
    org = {l.split("\t")[0]: l.split("\t")[1] for l in open(f"{D}/taxa.tsv").read().splitlines()[1:]}
    gtdb = Tree(f"{D}/gtdb_ref_pruned.nwk", format=1)
    concat = Tree(f"{D}/concat/concat.treefile", format=0)
    genes = open(f"{D}/genes.txt").read().split()
    os.makedirs(OUT, exist_ok=True)

    names = []
    for cand, spec in CANDIDATES.items():
        if isinstance(spec, str):
            group = sorted(a for a, o in org.items() if o.split()[0] == spec)
        else:
            by_name = {o: a for a, o in org.items()}
            group = sorted(by_name[o] for o in spec)
        d = f"{OUT}/{cand}/genes"
        os.makedirs(d, exist_ok=True)

        nwk, forced = constraint_tree(gtdb, concat.get_leaf_names(), group)
        in_concat = group_is_split(concat, group)
        open(f"{OUT}/{cand}/constraint.nwk", "w").write(nwk + "\n")
        with open(f"{OUT}/{cand}/group.txt", "w") as f:
            for a in group:
                f.write(f"{a}\t{org[a]}\n")

        rows = []
        for gene in genes:
            taxa = aln_taxa(f"{D}/trimmed/{gene}.faa")
            present = [a for a in group if a in taxa]
            if len(present) < 2:
                continue
            gnwk, gforced = constraint_tree(gtdb, taxa, present)
            if gforced == 0:
                continue  # GTDB and the group no longer conflict with this gene's taxa
            model = re.search(r"Best-fit model according to BIC: (\S+)",
                              open(f"{D}/gene_trees/{gene}.iqtree").read()).group(1)
            open(f"{d}/{gene}.constraint.nwk", "w").write(gnwk + "\n")
            rows.append((gene, model))
        with open(f"{OUT}/{cand}/genes.tsv", "w") as f:
            for gene, model in rows:
                f.write(f"{gene}\t{model}\n")

        names.append(cand)
        print(f"{cand}: {len(group)} genomes, {forced} GTDB branch(es) forced, "
              f"group in concat tree: {in_concat}, {len(rows)} informative genes")

    open(f"{OUT}/candidates.txt", "w").write("\n".join(names) + "\n")
    print(f"Wrote {OUT}/")


if __name__ == "__main__":
    main()
