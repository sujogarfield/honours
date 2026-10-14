#!/usr/bin/env python3
"""
synteny_analysis_full_contigs.py
Sensitivity check for synteny_analysis.py: reruns the exact same pipeline
(Mantel test, per-organism agreement, UPGMA dendrogram, MDS) restricted to
organisms whose flagellar gene set is entirely captured on a single contig
(coverage == 1.0, i.e. primary_contig_gene_count == total_gene_count).

This is deliberately NOT "genome has exactly one sequence in its .fna" —
that proxy is wrong. A finished bacterial genome routinely has 2+ sequences
(chromosome + plasmid(s)) and is not fragmented; conversely a genuinely
draft genome can still have its whole flagellar operon intact on its
largest contig even while the rest of the genome is scattered across many
small contigs. coverage == 1.0 directly answers "is this genome's ordering
of THIS gene set fully known", which is the only thing that matters for
the breakpoint-distance metric below — so it's used instead of raw contig
count.

Rationale: draft (multi-contig) assemblies can split flagellar genes across
contigs whose relative order is unknown. synteny_analysis.py already works
around this by using only each genome's primary contig, but genomes with
low primary-contig coverage still contribute noisy, low-N breakpoint
distances to the same matrix as fully-covered genomes (see
gene_order/synteny_analysis.json: primary_contig_gene_count / total_gene_count
per organism). Restricting to coverage == 1.0 removes that contig-
fragmentation confound entirely, at the cost of dropping every genome with
any flagellar gene split off onto a second contig.

Reuses every function from synteny_analysis.py unchanged — only the
organism list fed into the pipeline differs.

Usage:
    python3 synteny_analysis_full_contigs.py [gene_order_json] [output_json]

Defaults to gene_order/gene_order.json -> gene_order/synteny_analysis_full_contigs.json.
Pass gene_order/gene_order_scaffolded.json (see lift_ragtag_coordinates.py)
as gene_order_json to also credit the RagTag-recovered genomes.

--closed: use the stricter reliable set shared with interactive_tree.py
(scripts/common/reliability.py) -- every flagellar gene on one CLOSED
chromosome, natively complete or a RagTag scaffold -- instead of
coverage == 1.0 (which also admits single-scaffold draft assemblies), and
treat each gene order as circular: the last and first genes on the
chromosome are neighbours too (circular_breakpoint_distance). Pass
gene_order/gene_order_combined_scaffolded.json as gene_order_json so the
RagTag genomes are in scaffold coordinates. With --closed the distance is
the strand- and overlap-aware signed-adjacency distance
(scripts/common/synteny_distance.py, circular); add --breakpoint for the
unsigned circular breakpoint distance instead.

Run AFTER gene_order.py, gtdb_tree.py, and synteny_analysis.py.
"""

import itertools
import json
import os
import sys
from pathlib import Path

import numpy as np
from ete3 import Tree

from synteny_analysis import (
    GENE_ORDER_JSON as DEFAULT_GENE_ORDER_JSON,
    MIN_SHARED_GENES,
    N_PERMUTATIONS,
    NAMED_NWK,
    breakpoint_distance,
    classical_mds,
    load_gene_order,
    mantel_test,
    primary_contig_order,
    spearman,
    upgma_tree,
)

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from reliability import FASTA_DIR, RAGTAG_DIR, load_complete_chromosomes, is_reliable
from synteny_distance import primary_contig_records, signed_overlap_distance


def circular_breakpoint_distance(order_a: list[str], order_b: list[str]):
    """breakpoint_distance() for circular chromosomes: same normalised,
    direction-agnostic adjacency comparison, but the last and first shared
    genes are adjacent too, so there are len(shared) possible adjacencies
    rather than len(shared) - 1. Rotation-invariant: where the assembly's
    coordinate origin falls no longer matters."""
    shared = set(order_a) & set(order_b)
    if len(shared) < MIN_SHARED_GENES:
        return None, len(shared)

    sub_a = [g for g in order_a if g in shared]
    sub_b = [g for g in order_b if g in shared]
    n = len(shared)
    adj_a = {frozenset((sub_a[i], sub_a[(i + 1) % n])) for i in range(n)}
    adj_b = {frozenset((sub_b[i], sub_b[(i + 1) % n])) for i in range(n)}
    return 1 - len(adj_a & adj_b) / n, n

DEFAULT_OUTPUT_JSON = "campy_orthologs/archive_json/synteny_analysis_full_contigs.json"


def main():
    closed = "--closed" in sys.argv
    use_breakpoint = "--breakpoint" in sys.argv
    args = [a for a in sys.argv[1:] if a not in ("--closed", "--breakpoint")]
    GENE_ORDER_JSON = args[0] if len(args) > 0 else DEFAULT_GENE_ORDER_JSON
    OUTPUT_JSON = args[1] if len(args) > 1 else DEFAULT_OUTPUT_JSON
    signed_mode = closed and not use_breakpoint
    if not closed:
        distance = breakpoint_distance
    elif use_breakpoint:
        distance = circular_breakpoint_distance
    else:
        distance = lambda a, b: signed_overlap_distance(a, b, circular=True,
                                                        min_shared=MIN_SHARED_GENES)

    Path("gene_order").mkdir(exist_ok=True)

    print(f"Loading gene order from {GENE_ORDER_JSON} …")
    gene_order = load_gene_order(GENE_ORDER_JSON)

    print(f"Loading tree from {NAMED_NWK} …")
    tree = Tree(NAMED_NWK, format=1, quoted_node_names=True)
    leaf_nodes = {leaf.name: leaf for leaf in tree.get_leaves()}

    all_organisms = sorted(name for name in gene_order if name in leaf_nodes)
    skipped_not_in_tree = sorted(name for name in gene_order if name not in leaf_nodes)

    # Compute primary-contig order/coverage for every organism first (same as
    # synteny_analysis.py), then filter to coverage == 1.0.
    orders, signed, primary_sizes, total_genes = {}, {}, {}, {}
    for name in all_organisms:
        recs = gene_order[name]
        ordered, psize = primary_contig_order(recs)
        orders[name] = ordered
        if signed_mode:
            signed[name] = primary_contig_records(recs)[0]
        primary_sizes[name] = psize
        total_genes[name] = len(recs)

    # Genomes with fewer than MIN_SHARED_GENES genes reach coverage == 1.0
    # trivially (a single gene is always "on one contig") but can never form
    # a valid pair in the Mantel test, so they'd inflate n_organisms without
    # contributing anything -- drop them too.
    too_few = [name for name in all_organisms if total_genes[name] < MIN_SHARED_GENES]
    if closed:
        complete = load_complete_chromosomes(FASTA_DIR, RAGTAG_DIR)
        keep = lambda name: is_reliable(gene_order[name], complete)
        why = "not all on one closed chromosome"
    else:
        keep = lambda name: primary_sizes[name] == total_genes[name]
        why = "genes split across >1 contig"
    organisms = [
        name for name in all_organisms
        if total_genes[name] >= MIN_SHARED_GENES and keep(name)
    ]
    n_dropped_too_few = len(too_few)
    n_dropped = len(all_organisms) - len(organisms) - n_dropped_too_few
    print(f"{len(organisms)} organisms kept "
          f"({'closed chromosome, circular distance' if closed else 'coverage == 1.0'}), "
          f"{n_dropped} dropped ({why}), "
          f"{n_dropped_too_few} dropped (< {MIN_SHARED_GENES} genes)")

    n = len(organisms)
    dist_synt = np.zeros((n, n))
    valid = np.zeros((n, n), dtype=bool)

    for i, j in itertools.combinations(range(n), 2):
        pair = signed if signed_mode else orders
        d, _ = distance(pair[organisms[i]], pair[organisms[j]])
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
        per_organism[name] = {
            "primary_contig_gene_count": primary_sizes[name],
            "total_gene_count": total_genes[name],
            "coverage": 1.0,
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
        "filter": ((f"closed_chromosome_only (native complete assembly or RagTag scaffold; "
                    f"circular {'signed-adjacency + overlap' if signed_mode else 'breakpoint'} "
                    f"distance), >= {MIN_SHARED_GENES} genes") if closed else
                   (f"full_flagellar_coverage_only (coverage == 1.0, not raw contig count), "
                    f">= {MIN_SHARED_GENES} genes")),
        "mantel": {
            "r": r,
            "p_value": p,
            "n_pairs": n_pairs,
            "n_permutations": N_PERMUTATIONS,
            "n_organisms": n,
            "n_dropped_incomplete_coverage": n_dropped,
            "n_dropped_too_few_genes": n_dropped_too_few,
            "n_skipped_not_in_tree": len(skipped_not_in_tree),
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
