#!/usr/bin/env python3
"""
build_concat_tree.py
Concatenated ("supermatrix") tree across multiple flagellar/chemotaxis
orthogroups -- the whole-operon-level complement to the single-orthogroup
tanglegrams. Where one gene tree can be wrong from single-locus noise (no
species has enough signal in ~150-600bp of one protein to fully resolve
deep splits, and any one gene can carry its own HGT/recombination history),
concatenating many genes into one alignment and building ONE tree averages
that noise out and recovers the dominant "vertical inheritance" signal --
this is the standard phylogenomics approach, and what was suggested as a
complementary check to the tanglegrams.

Any orthogroup whose tree still strongly disagrees with the concat tree in
a single-gene tanglegram, despite good bootstrap support on both sides, is
much more likely to be real HGT rather than reconstruction noise.

Gene selection: single-copy (no genome has 2+ copies -- required, a
concatenated alignment needs exactly one sequence per genome per gene) AND
well-populated (>=100 of the 166 genomes) confirmed genes from the combined
dataset. This gives a core backbone of basal-body/motor/chemotaxis genes
spread across a large majority of genomes, while still allowing genomes
that lack a handful of these genes (gap-filled in the alignment -- standard
supermatrix practice; requiring the full intersection would discard most
genomes over one or two missing genes).

Usage: python3 build_concat_tree.py
Output:
  gene_order/concat_tree/concat_alignment.fasta   -- supermatrix, gap-filled
  gene_order/concat_tree/partitions.txt            -- per-gene partition spec for IQ-TREE
  gene_order/concat_tree/gene_coverage.json        -- per-gene genome coverage report
  gene_order/concat_tree/concat_tree.treefile      -- final ML tree (after IQ-TREE run)
"""

import csv
import glob
import json
import os
from collections import defaultdict

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "06_orthology"))
from helpers import get_organism_from_fasta
from orthology_crossref import build_gff_index, GFF_DIR

GENE_ORDER_JSON = "gene_order/gene_order_combined.json"
ORTHOGROUPS_TSV = "campy_orthologs/orthofinder_results/Results_Sep23/Orthogroups/Orthogroups.tsv"
MSA_DIR = "campy_orthologs/orthofinder_results/Results_Sep23/MultipleSequenceAlignments"
OUTPUT_DIR = "gene_order/concat_tree"

MIN_GENOME_COVERAGE = 100  # out of 166 -- keeps the backbone well-populated
MIN_GENES_PER_GENOME = 9   # out of 15 selected genes (60%) -- drops near-all-gap taxa
                            # that would otherwise be placed by noise, not signal


def load_gene_to_proteins():
    """gene -> {gcf_id: protein_id}, single-copy genes only (dedup safety net;
    the real single-copy filtering/reporting happens in main())."""
    combined = json.load(open(GENE_ORDER_JSON))
    gff_bases = {os.path.basename(p)[:-4] for p in glob.glob(f"{GFF_DIR}/*.gff")}
    locus_idx = {}
    for base in sorted(gff_bases):
        l2p, _ = build_gff_index(f"{GFF_DIR}/{base}.gff")
        locus_idx[base] = l2p

    gene_to_proteins = defaultdict(lambda: defaultdict(list))
    for gcf, data in combined.items():
        for g in data["genes"]:
            pid = g.get("protein_id")
            if not pid:
                locus = (g["contig"], g["start"], g["end"], g["strand"])
                pid = locus_idx.get(gcf, {}).get(locus)
            if pid:
                gene_to_proteins[g["gene"]][gcf].append(pid)
    return gene_to_proteins


def load_protein_to_og():
    """protein_id -> orthogroup, from Orthogroups.tsv (needed for genes whose
    candidate records don't carry an explicit 'orthogroup' field -- i.e.
    purely string-matched genes)."""
    protein_to_og = {}
    with open(ORTHOGROUPS_TSV) as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)  # header
        for row in reader:
            og_id = row[0]
            for cell in row[1:]:
                for p in cell.split(","):
                    p = p.strip()
                    if p:
                        protein_to_og[p] = og_id
    return protein_to_og


def load_alignment(og_id):
    """{protein_id: aligned_seq} from OrthoFinder's existing MSA for this OG."""
    path = os.path.join(MSA_DIR, f"{og_id}.fa")
    seqs = {}
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            seqs[name] = "".join(chunks)
    return seqs


def og_leaf_name_for(gcf_id, protein_id, sample_names):
    """OrthoFinder mangles GCF prefixes (dots -> underscores) in its own
    sequence headers, e.g. GCF_000009085.1 -> GCF_000009085_1. Reconstruct
    that mangled prefix and confirm against a real header from the file
    rather than assuming the exact mangling rule always holds."""
    mangled = gcf_id.replace(".", "_")
    candidate = f"{mangled}_{protein_id}"
    if candidate in sample_names:
        return candidate
    # fallback: protein_id is unique enough to find by suffix match
    for n in sample_names:
        if n.endswith(protein_id):
            return n
    return None


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Loading combined gene->protein map ...")
    gene_to_proteins = load_gene_to_proteins()

    print("Loading Orthogroups.tsv protein->OG map ...")
    protein_to_og = load_protein_to_og()

    combined = json.load(open(GENE_ORDER_JSON))
    all_gcfs = sorted(combined.keys())
    organism_of = {gcf: combined[gcf]["organism"] for gcf in all_gcfs}

    # Select single-copy, well-populated genes
    selected = []
    for gene, gp in sorted(gene_to_proteins.items()):
        if not all(len(v) == 1 for v in gp.values()):
            continue  # paralog somewhere -- can't concatenate cleanly
        if len(gp) < MIN_GENOME_COVERAGE:
            continue
        selected.append(gene)
    print(f"Selected {len(selected)} single-copy genes with >= {MIN_GENOME_COVERAGE} genomes: {selected}")

    # Resolve each selected gene's orthogroup + per-genome aligned sequence
    coverage_report = {}
    gene_blocks = {}  # gene -> {gcf_id: aligned_seq}
    gene_lengths = {}  # gene -> alignment column count
    for gene in selected:
        gp = gene_to_proteins[gene]
        # Determine OG id: majority vote over protein_to_og lookups
        og_votes = defaultdict(int)
        for gcf, pids in gp.items():
            og = protein_to_og.get(pids[0])
            if og:
                og_votes[og] += 1
        if not og_votes:
            print(f"  {gene}: no orthogroup resolvable, skipping")
            continue
        og_id = max(og_votes, key=og_votes.get)
        agree = og_votes[og_id]
        print(f"  {gene}: orthogroup {og_id} ({agree}/{len(gp)} proteins agree)")

        seqs = load_alignment(og_id)
        sample_names = set(seqs.keys())
        block = {}
        matched, missing = 0, 0
        for gcf, pids in gp.items():
            leaf = og_leaf_name_for(gcf, pids[0], sample_names)
            if leaf and leaf in seqs:
                block[gcf] = seqs[leaf]
                matched += 1
            else:
                missing += 1
        if not block:
            print(f"    no sequences matched in {og_id}.fa, skipping gene")
            continue
        aln_len = len(next(iter(block.values())))
        gene_blocks[gene] = block
        gene_lengths[gene] = aln_len
        coverage_report[gene] = {
            "orthogroup": og_id, "alignment_length": aln_len,
            "genomes_matched": matched, "genomes_unmatched_in_msa": missing,
        }
        print(f"    {matched} genomes matched to alignment ({aln_len} columns)")

    kept_genes = list(gene_blocks.keys())
    print(f"\n{len(kept_genes)} genes kept for concatenation: {kept_genes}")

    # Genomes included: union of all genomes appearing in ANY kept gene block,
    # then drop near-all-gap taxa (too few loci to place reliably).
    genomes_union = sorted(set().union(*(set(b.keys()) for b in gene_blocks.values())))
    genes_per_genome = {
        gcf: sum(1 for g in kept_genes if gcf in gene_blocks[g]) for gcf in genomes_union
    }
    genomes_included = sorted(g for g in genomes_union if genes_per_genome[g] >= MIN_GENES_PER_GENOME)
    dropped = sorted(set(genomes_union) - set(genomes_included))
    print(f"{len(genomes_union)} genomes have >=1 gene; keeping {len(genomes_included)} with "
          f">= {MIN_GENES_PER_GENOME}/{len(kept_genes)} genes")
    if dropped:
        print(f"  dropped {len(dropped)} sparse genomes (too few loci to place reliably): {dropped}")

    # Build concatenated, gap-filled alignment
    concat_path = os.path.join(OUTPUT_DIR, "concat_alignment.fasta")
    with open(concat_path, "w") as out:
        for gcf in genomes_included:
            seq_parts = []
            for gene in kept_genes:
                block = gene_blocks[gene]
                seq_parts.append(block.get(gcf, "-" * gene_lengths[gene]))
            seq = "".join(seq_parts)
            header = gcf.replace(".", "_")  # keep tree leaf names filesystem/newick-safe
            out.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")
    total_len = sum(gene_lengths[g] for g in kept_genes)
    print(f"Wrote {concat_path} ({len(genomes_included)} sequences x {total_len} columns)")

    # Partition file for IQ-TREE (-p), RAxML-style, one block per gene,
    # protein data -- "AUTO" lets ModelFinder pick each partition's own
    # substitution model rather than forcing one model on all 15 genes.
    part_path = os.path.join(OUTPUT_DIR, "partitions.txt")
    with open(part_path, "w") as out:
        pos = 1
        for gene in kept_genes:
            length = gene_lengths[gene]
            out.write(f"AUTO, {gene} = {pos}-{pos + length - 1}\n")
            pos += length
    print(f"Wrote {part_path}")

    # Coverage + genome->organism lookup for downstream tree labeling
    report_path = os.path.join(OUTPUT_DIR, "gene_coverage.json")
    json.dump({
        "genes": coverage_report,
        "n_genomes_included": len(genomes_included),
        "genomes_included": {gcf.replace(".", "_"): organism_of[gcf] for gcf in genomes_included},
        "genomes_dropped_sparse": {gcf.replace(".", "_"): organism_of[gcf] for gcf in dropped},
        "min_genes_per_genome": MIN_GENES_PER_GENOME,
        "total_alignment_columns": total_len,
    }, open(report_path, "w"), indent=2)
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()
