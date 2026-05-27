#!/usr/bin/env python3

from Bio import SeqIO
import matplotlib.pyplot as plt
import os
from ete3 import Tree
import random

def get_organism_from_fasta(fasta_path):
    record = next(SeqIO.parse(fasta_path, "fasta"))
    description = record.description

    parts = description.split(" ", 1)

    if len(parts) > 1:
        name_part = parts[1].split(",")

        # name_part = name_part.split(" chromosome")[0]
        # name_part = name_part.split(" complete genome")[0]

        # return name_part.replace(" ", "_")
        just_name = name_part[0].split(" ")
        genus = just_name[0]
        species = just_name[1].lower()

        if genus == "Candidatus" and len(just_name) > 2:
            genus = species.capitalize()
            species = just_name[2]
            print(f"{genus} {species}")
            return f"{genus} {species}"

        return f"{genus} {species}"
        # return name_part[0].replace(" ", "_")

    return "Unknown"
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
    large_leaves = set(large_tree.get_leaf_names())
    missing = set(small_leaves) - large_leaves
    
    if missing:
        print(f"{len(missing)} taxa not in GTDB")
        return False
    
    print(f"all {len(small_leaves)} taxa found in GTDB")
    
    # 2: topology matches when re-pruned
    fresh_pruned = large_tree.copy()
    fresh_pruned.prune(small_leaves, preserve_branch_length=True)
    
    rf, max_rf = fresh_pruned.robinson_foulds(small_tree, unrooted_trees=True)[:2]
    
    if rf != 0:
        print(f"topology mismatch (RF distance = {rf}/{max_rf})")
        return False
    
    print("topology preserved correctly")
    
    # 3: pairwise distances preserved 
    # sample random pairs to verify distances match
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
