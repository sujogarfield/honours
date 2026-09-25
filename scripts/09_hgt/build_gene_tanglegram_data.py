#!/usr/bin/env python3
"""
build_gene_tanglegram_data.py
Pairs a specific orthogroup's OrthoFinder gene tree (Resolved_Gene_Trees/)
against the GTDB species tree, for the HGT tanglegram: crossing connectors
between the two trees mean a genome's placement by GENE sequence disagrees
with its placement by the SPECIES tree -- phylogenetic incongruence, the
classic HGT signature (see hgt_detection_notes.txt for the full reasoning
and why this is a more direct test than synteny/breakpoint distance).

Handles the leaf-count mismatch between gene trees and the species tree
head-on (a gene tree has one leaf per SEQUENCE, not per genome -- see
findings summary / conversation for real examples: flgI is 154 leaves for
154 genomes, clean 1:1; flaab is 209 leaves for 92 genomes, real
duplication): for a genome with multiple copies in the gene tree, EVERY
copy gets its own leaf in the output gene tree (labelled with its protein
id), and the species tree gets ONE leaf per genome duplicated out to match
-- so a genome with 2 paralogs shows up twice on the species-tree side too,
each connected to its own gene-tree copy. This keeps the tanglegram
honest about duplication rather than arbitrarily picking one copy.

Usage: build_gene_tanglegram_data.py <orthogroup_id> [<orthogroup_id> ...]
Output: gene_order/tanglegram_<orthogroup_id>.json per orthogroup, plus
        gene_order/tanglegram_index.json listing what's available (for the
        dashboard's dropdown).
"""

import json
import re
import sys
from pathlib import Path

from ete3 import Tree

RESOLVED_TREES = "campy_orthologs/orthofinder_results/Results_Sep23/Resolved_Gene_Trees/Resolved_Gene_Trees.txt"
NAMED_NWK = "images_and_newick/thesis_b/gtdb_ref_tree.nwk"
GENE_ORDER_JSON = "gene_order/gene_order_combined.json"
OUTPUT_DIR = Path("gene_order")
INDEX_JSON = OUTPUT_DIR / "tanglegram_index.json"


def load_gcf_prefix_map():
    """{orthofinder-mangled-prefix: (gcf_id, organism)} for all 166 genomes."""
    with open(GENE_ORDER_JSON) as f:
        combined = json.load(f)
    out = {}
    for gcf_id, data in combined.items():
        mangled = gcf_id.replace(".", "_")
        out[mangled] = (gcf_id, data["organism"])
    return out


def load_orthogroup_tree(og_id: str) -> str:
    with open(RESOLVED_TREES) as f:
        for line in f:
            this_og, tree = line.split(":", 1)
            if this_og == og_id:
                return tree.strip()
    raise ValueError(f"{og_id} not found in {RESOLVED_TREES}")


def tree_to_dict(node, leaf_display: dict) -> dict:
    """Nested {name, match, children}. "match" is the unique matching key
    used to pair leaves across the two trees (needed because paralogs mean
    the same organism can appear more than once -- matching by organism
    name alone would collide); "name" is what's actually displayed."""
    if node.is_leaf():
        return {"name": leaf_display.get(node.name, node.name), "match": node.name}
    d = {"name": "", "match": ""}
    children = [tree_to_dict(c, leaf_display) for c in node.children]
    if children:
        d["children"] = children
    return d


def build_one(og_id: str, prefix_map: dict, species_tree: Tree):
    newick = load_orthogroup_tree(og_id)
    gtree = Tree(newick, format=1)

    # Resolve each leaf to (gcf_id, organism, protein_id); drop anything
    # that doesn't match one of our 166 genomes (shouldn't happen, but be
    # defensive rather than crash on an unexpected leaf format).
    leaf_display = {}
    leaf_organism_of = {}  # leaf.name -> organism (for grouping duplicates)
    unmatched = []
    sorted_prefixes = sorted(prefix_map, key=len, reverse=True)
    for leaf in gtree.get_leaves():
        matched = None
        for prefix in sorted_prefixes:
            if leaf.name.startswith(prefix + "_"):
                matched = prefix
                break
        if matched is None:
            unmatched.append(leaf.name)
            continue
        gcf_id, organism = prefix_map[matched]
        protein_id = leaf.name[len(matched) + 1:]
        leaf_organism_of[leaf.name] = organism
        leaf_display[leaf.name] = f"{organism} [{protein_id}]"

    if unmatched:
        print(f"  {og_id}: {len(unmatched)} unmatched leaves (not one of our 166 genomes), dropped")

    gene_tree_dict = tree_to_dict(gtree, leaf_display)

    # Species tree side: one leaf per genome by default, but duplicated out
    # per extra paralog copy so genomes with N gene-tree copies get N
    # species-tree leaves too. Each duplicate is paired to one SPECIFIC
    # gene-tree copy via a shared "match" key (that copy's own gene-tree
    # leaf name) -- matching by organism name alone would collide for any
    # genome with paralogs, connecting to an arbitrary copy instead of the
    # right one.
    from collections import defaultdict
    gene_leaves_by_organism = defaultdict(list)
    for gene_leaf_name, organism in leaf_organism_of.items():
        gene_leaves_by_organism[organism].append(gene_leaf_name)
    organisms_present = set(gene_leaves_by_organism)

    pruned = species_tree.copy()
    pruned.prune([leaf.name for leaf in pruned.get_leaves() if leaf.name in organisms_present],
                 preserve_branch_length=True)

    def species_tree_to_dict(node):
        if node.is_leaf():
            copies = gene_leaves_by_organism.get(node.name, [node.name])
            if len(copies) <= 1:
                return {"name": node.name, "match": copies[0]}
            # duplicate this leaf once per paralog copy under a zero-length
            # spacer node, each carrying that copy's own gene-tree leaf name
            # as its match key, so the tanglegram draws N independent,
            # correctly-paired connector lines instead of one ambiguous one.
            return {"name": "", "match": "", "children": [
                {"name": f"{node.name} (copy {i+1}/{len(copies)})", "match": c}
                for i, c in enumerate(copies)
            ]}
        children = [species_tree_to_dict(c) for c in node.children]
        return {"name": "", "match": "", "children": children} if children else {"name": node.name, "match": node.name}

    species_dict = species_tree_to_dict(pruned)

    n_leaves_gene = len(gtree.get_leaves())
    n_genomes = len(organisms_present)
    print(f"  {og_id}: {n_leaves_gene} gene-tree leaves, {n_genomes} genomes "
          f"({'clean single-copy' if n_leaves_gene == n_genomes else 'has paralogs'})")

    return {
        "orthogroup": og_id,
        "n_leaves": n_leaves_gene,
        "n_genomes": n_genomes,
        "gene_tree": gene_tree_dict,
        "species_tree": species_dict,
    }


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    og_ids = sys.argv[1:]

    print("Loading GCF->organism map ...")
    prefix_map = load_gcf_prefix_map()
    print(f"Loading species tree from {NAMED_NWK} ...")
    species_tree = Tree(NAMED_NWK, format=1, quoted_node_names=True)

    OUTPUT_DIR.mkdir(exist_ok=True)
    index = []
    if INDEX_JSON.exists():
        index = json.load(open(INDEX_JSON))
        index = [e for e in index if e["orthogroup"] not in og_ids]  # replace if rebuilding

    for og_id in og_ids:
        print(f"Building tanglegram data for {og_id} ...")
        data = build_one(og_id, prefix_map, species_tree)
        out_path = OUTPUT_DIR / f"tanglegram_{og_id}.json"
        with open(out_path, "w") as f:
            json.dump(data, f)
        index.append({"orthogroup": og_id, "n_leaves": data["n_leaves"], "n_genomes": data["n_genomes"]})
        print(f"  Saved {out_path}")

    with open(INDEX_JSON, "w") as f:
        json.dump(index, f, indent=2)
    print(f"\nIndex updated: {INDEX_JSON} ({len(index)} orthogroups available)")


if __name__ == "__main__":
    main()
