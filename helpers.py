#!/usr/bin/env python3

from Bio import SeqIO
import matplotlib.pyplot as plt
import os
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, CircleFace
from pathlib import Path
import random

def get_organism_from_fasta(fasta_path):
    record = next(SeqIO.parse(fasta_path, "fasta"))
    description = record.description
    parts = description.split(" ", 1)

    if len(parts) > 1:
        name_part = parts[1].split(",")
        just_name = name_part[0].split(" ")
        genus = just_name[0]
        species = just_name[1].lower()

        if genus == "Candidatus" and len(just_name) > 2:
            genus = species.capitalize()
            species = just_name[2]
            return f"{genus} {species}"

        return f"{genus} {species}"
    return "Unknown"

def colour_tree_by_species(tree):
    for leaf in tree.iter_leaves():
        nstyle = NodeStyle()
        nstyle["size"] = 0
        if leaf.name.startswith("Helicobacter"):
            nstyle["bgcolor"] = "#ffd2d0"
        elif leaf.name.startswith("Campylobacter"):
            nstyle["bgcolor"] = "#ddffdd"
        elif leaf.name.startswith("Sulfurimonas"):
            nstyle["bgcolor"] = "#fffad1"
        elif leaf.name.startswith("Wolinella"):
            nstyle["bgcolor"] = "#edd9ff"
        elif leaf.name.startswith("Hippea"):
            nstyle["bgcolor"] = "#d0e9ff"
        elif leaf.name.startswith("Sulfuricurvum"):
            nstyle["bgcolor"] = "#ecd5ff"
        elif leaf.name.startswith("Candidatus"):
            nstyle["bgcolor"] = "#d7ffeb"
        elif leaf.name.startswith("Sulfurospirillum"):
            nstyle["bgcolor"] = "#feffc9"
        elif leaf.name.startswith("Hydrogenimonas"):
            nstyle["bgcolor"] = "#d6fff8"
        elif leaf.name.startswith("Nitrosophilus"):
            nstyle["bgcolor"] = "#e2dfff"
        elif leaf.name.startswith("Desulfurella"):
            nstyle["bgcolor"] = "#dbf4ff"
        elif leaf.name.startswith("Caminibacter"):
            nstyle["bgcolor"] = "#ffdefd"
        elif leaf.name.startswith("Nautilia"):
            nstyle["bgcolor"] = "#ffdbdb"
        elif leaf.name.startswith("Lebetimonas"):
            nstyle["bgcolor"] = "#e4ffdf"
        elif leaf.name.startswith("Nitratiruptor"):
            nstyle["bgcolor"] = "#faffd8"
        elif leaf.name.startswith("Sulfurovum"):
            nstyle["bgcolor"] = "#fff4da"
        elif leaf.name.startswith("Nitratifractor"):
            nstyle["bgcolor"] = "#ffdbdb"
        elif leaf.name.startswith("Arcobacter"):
            nstyle["bgcolor"] = "#deffe5"
        elif leaf.name.startswith("Aliarcobacter"):
            nstyle["bgcolor"] = "#e5e8ff"
        elif leaf.name.startswith("Malaciobacter"):
            nstyle["bgcolor"] = "#ffdcf9"
        elif leaf.name.startswith("Halarcobacter"):
            nstyle["bgcolor"] = "#e6fcff"
        elif leaf.name.startswith("Poseidonibacter"):
            nstyle["bgcolor"] = "#ffe1fc"
        leaf.set_style(nstyle)

"""
Verify that pruning hasn't distorted phylogenetic relationships.
Checks that pairwise distances between taxa are preserved.
"""
def validate_pruned_tree_integrity():
    large_tree = Tree("bac120.tree", format=1, quoted_node_names=True)
    small_tree = Tree(
        "phylo_trees/gtdb_ref_tree_og.nwk",
        format=1,
        quoted_node_names=True
    )
    
    small_leaves = small_tree.get_leaf_names()

    # 1: all taxa exist
    # no new random taxa are added
    large_leaves = set(large_tree.get_leaf_names())
    missing = set(small_leaves) - large_leaves
    
    if missing:
        print(f"{len(missing)} taxa not in GTDB")
        return False
    
    print(f"all {len(small_leaves)} taxa found in GTDB")
    
    # 2: topology matches when re-pruned
    # proves the branching pattern between the 166 species
    # is identical to what GTDB has for those same species
    # i.e. pruning didn't accidentally swap any relationships
    fresh_pruned = large_tree.copy()
    fresh_pruned.prune(small_leaves, preserve_branch_length=True)
    
    rf, max_rf = fresh_pruned.robinson_foulds(small_tree, unrooted_trees=True)[:2]
    
    if rf != 0:
        print(f"topology mismatch (RF distance = {rf}/{max_rf})")
        return False
    
    print("topology preserved correctly")
    
    # 3: pairwise distances preserved 
    # proves the evolutionary distances between species pairs
    # are numerically identical to GTDB branch lengths
    # were correctly transferred when intermediate nodes collapsed
    sample_size = min(20, len(small_leaves) // 2)
    sampled_pairs = []
    
    for _ in range(sample_size):
        pair = random.sample(small_leaves, 2)
        sampled_pairs.append(pair)
    
    distance_errors = []
    for leaf1, leaf2 in sampled_pairs:
        node1_large = large_tree.search_nodes(name=leaf1)[0]
        node2_large = large_tree.search_nodes(name=leaf2)[0]
        dist_large = node1_large.get_distance(node2_large)
        
        node1_small = small_tree.search_nodes(name=leaf1)[0]
        node2_small = small_tree.search_nodes(name=leaf2)[0]
        dist_small = node1_small.get_distance(node2_small)

        print(f"{dist_large} - {dist_small} = {abs(dist_large - dist_small)}")
        
        if abs(dist_large - dist_small) > 1e-6:
            distance_errors.append({
                'pair': (leaf1.split('_')[0], leaf2.split('_')[0]),
                'gtdb': dist_large,
                'pruned': dist_small,
                'diff': abs(dist_large - dist_small)
            })
    
    if distance_errors:
        print(f"{len(distance_errors)} distance mismatches detected")
        return False
    
    print("pruned tree did not distort original relationships of GTDB.")
    return True

def validate_bootstrap_support(
    small_nwk="phylo_trees/gtdb_ref_tree_og.nwk",
    threshold=0.70,
    output_png="phylo_trees/gtdb_ref_tree_support.png"
):
    """
    Read GTDB SH support values preserved in subtree and render a
    colour-coded tree. Green >= 0.95, orange >= 0.70, red < 0.70.

    Returns dict with support summary statistics.
    """
    print("\n── Bootstrap / SH support analysis ──")

    tree = Tree(small_nwk, format=1, quoted_node_names=True)
    internal = [n for n in tree.traverse() if not n.is_leaf()]
    values = [n.support for n in internal]

    if not any(v > 0 for v in values):
        print("⚠️  No support values found in tree — may have been lost during pruning")
        return None

    strong   = sum(1 for v in values if v >= 0.95)
    moderate = sum(1 for v in values if 0.70 <= v < 0.95)
    weak     = sum(1 for v in values if v < 0.70)
    total    = len(values)

    print(f"Internal nodes : {total}")
    print(f"Strong  ≥ 0.95 : {strong}  ({100*strong/total:.1f}%)")
    print(f"Moderate 0.70–0.95 : {moderate}  ({100*moderate/total:.1f}%)")
    print(f"Weak    < 0.70 : {weak}  ({100*weak/total:.1f}%)")
    print(f"Mean support   : {sum(values)/total:.3f}")

    weak_nodes = [n for n in internal if n.support < threshold and n.support > 0]
    if weak_nodes:
        print(f"\nWeakly supported clades (< {threshold}):")
        for node in weak_nodes:
            leaves = node.get_leaf_names()
            print(f"  SH={node.support:.2f}  covers: {', '.join(leaves[:3])}{'...' if len(leaves) > 3 else ''}")

    # ── Render support-coloured tree ──────────────────────────────────────
    def support_color(v):
        if v >= 0.95: return "#2ecc71"
        if v >= 0.70: return "#f39c12"
        return "#e74c3c"

    for node in tree.traverse():
        ns = NodeStyle()
        if node.is_leaf():
            ns["size"] = 0
        else:
            ns["size"] = 8
            ns["fgcolor"] = support_color(node.support)
            if node.support > 0:
                face = TextFace(f" {node.support:.2f}", fsize=6,
                                fgcolor=support_color(node.support))
                node.add_face(face, column=0, position="branch-top")
        node.set_style(ns)

    fasta_dir = Path("campy_fas")
    for leaf in tree.iter_leaves():
        bare_acc = leaf.name.replace("RS_", "").replace("GB_", "")
        matching_files = list(fasta_dir.glob(f"{bare_acc}*.fna"))
        if matching_files:
            fasta_file = matching_files[0]
            if fasta_file.exists():
                organism_name = get_organism_from_fasta(fasta_file)
                display_name = organism_name
            else:
                display_name = leaf.name
            leaf.name = display_name
            face = TextFace(leaf.name, fstyle="italic", fsize=9)
            leaf.add_face(face, column=0, position="branch-right")

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_scale = False

    for label, color in [("≥ 0.95 strong", "#2ecc71"),
                          ("0.70–0.95 moderate", "#f39c12"),
                          ("<0.70 weak", "#e74c3c")]:
        ts.legend.add_face(CircleFace(6, color), column=0)
        ts.legend.add_face(TextFace(f" {label}", fsize=8), column=1)

    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    tree.render(output_png, tree_style=ts, w=800, h=4000)

    return {"strong": strong, "moderate": moderate,
            "weak": weak, "mean": sum(values)/total}
    

def plot_species(species_name, genes, output_dir="gene_order"):
    os.makedirs(output_dir, exist_ok=True)
    genes = sorted(genes, key=lambda x: (x["start"]))

    contigs = sorted(set(g["contig"] for g in genes))
    y_positions = {contig: i for i, contig in enumerate(contigs)}

    fig_height = max(2, len(contigs) * 0.6)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    min_start = min(g['start'] for g in genes)
    max_end = max(g['end'] for g in genes)
    padding = max((max_end - min_start) * 0.02, 50)
    ax.set_xlim(min_start - padding, max_end + padding)

    prev_end = {}
    overlaps_info = []

    for i, gene in enumerate(genes):
        y = y_positions[gene["contig"]]
        start = gene['start']
        end = gene['end']
        strand = gene['strand']
        width = max(end - start, 50)

        key = gene['contig']
        if key in prev_end and start < prev_end[key]:
            color = "green"
        else:
            color = "blue" if strand == "+" else "red"
        prev_end[key] = end

        ax.add_patch(plt.Rectangle((start, y), width, 0.6, color=color))
        ax.text((start + end)/2, y + 0.6, gene['gene'], ha='center', fontsize=6)

        overlapping_genes = [g2['gene'] for j,g2 in enumerate(genes)
                             if j != i and g2['contig']==key and not (g2['end'] <= start or g2['start'] >= end)]
        overlaps_info.append({
            "gene": gene['gene'],
            "contig": key,
            "start": start,
            "end": end,
            "strand": strand,
            "overlaps": overlapping_genes
        })

    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels(contigs, fontsize=6)

    ax.set_xlabel("Genome position")
    ax.set_title(species_name)
    plt.subplots_adjust(left=0.25, right=0.95)

    safe_name = species_name.replace(" ", "_")
    plt.savefig(os.path.join(output_dir, f"{safe_name}.png"), dpi=300, bbox_inches='tight')
    plt.close()

    txt_path = os.path.join(output_dir, f"{safe_name}_genes.txt")
    with open(txt_path, "w") as f:
        f.write("gene\tcontig\tstart\tend\tstrand\toverlaps\n")
        for g in overlaps_info:
            f.write(f"{g['gene']}\t{g['contig']}\t{g['start']}\t{g['end']}\t{g['strand']}\t{','.join(g['overlaps'])}\n")

    print(f"see plot and overlaps: {safe_name}.png and gene info: {safe_name}_genes.txt")
