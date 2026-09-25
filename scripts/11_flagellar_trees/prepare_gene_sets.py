#!/usr/bin/env python3
"""
prepare_gene_sets.py
Builds the per-gene protein FASTA inputs for the Katana IQ-TREE run
(katana/*.pbs), from the combined OrthoFinder + eggNOG gene set
(gene_order/gene_order_combined.json).

Supersedes the 15-gene selection in 09_hgt/build_concat_tree.py, which
required a gene to be strictly single-copy in every genome. With the
orthology-expanded dataset almost every well-populated gene has a paralog
in a handful of genomes, so that rule now throws away most of the signal.
Rule used here instead:

  - gene present in >= MIN_GENOME_COVERAGE genomes, and
  - multi-copy in <= MAX_MULTICOPY_FRAC of the genomes that carry it;
  - genomes where that gene IS multi-copy are left out of that gene only
    (which copy is the ortholog is ambiguous), not dropped from the study;
  - genomes carrying < MIN_GENE_FRAC of the selected genes are dropped from
    every gene (too few loci to place reliably -- same rationale as the
    MIN_GENES_PER_GENOME cut in build_concat_tree.py).

Sequences are read from campy_fetched/campy_prot (committed), so Katana does
not need the 6 GB campy_orthologs/ folder. Alignment and trimming happen on
Katana (katana/01_gene_trees.pbs), not here.

Also writes the GTDB reference tree pruned to the same taxa and relabelled
to bare accessions, so IQ-TREE can compute concordance factors on the
species-tree branches too.

Run from the project root:  python scripts/11_flagellar_trees/prepare_gene_sets.py
Output (gene_order/flagellar_trees/):
  genes/<gene>.faa          -- unaligned proteins, one per genome, header = accession
  gene_selection.json       -- selection report (per-gene counts, dropped genomes)
  taxa.tsv                  -- accession -> organism, for labelling trees
  gtdb_ref_pruned.nwk       -- GTDB bac120 tree pruned to the included taxa
  genes.txt                 -- selected gene names, one per line (PBS array index)
"""

import glob
import json
import os
import sys
from collections import Counter, defaultdict

from ete3 import Tree

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "06_orthology"))
from orthology_crossref import build_gff_index, GFF_DIR

GENE_ORDER_JSON = "gene_order/gene_order_combined.json"
PROT_DIR = "campy_fetched/campy_prot"
GTDB_TREE = "images_and_newick/thesis_b/gtdb_ref_tree_og.nwk"
OUTPUT_DIR = "gene_order/flagellar_trees"

MIN_GENOME_COVERAGE = 100   # of 164 genomes with flagellar genes
MAX_MULTICOPY_FRAC = 0.10   # excludes flik, fliw, flge (25-94% multi-copy)
MIN_GENE_FRAC = 0.60        # genome must carry >= 60% of the selected genes


def accession(genome_key):
    """GCF_000012965.1_ASM1296v1 -> GCF_000012965.1"""
    return "_".join(genome_key.split("_")[:2])


def read_fasta(path):
    seqs, name, chunks = {}, None, []
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


def main():
    os.makedirs(os.path.join(OUTPUT_DIR, "genes"), exist_ok=True)
    combined = json.load(open(GENE_ORDER_JSON))

    # gene -> genome -> [protein_id, ...]
    gene_to_proteins = defaultdict(lambda: defaultdict(list))
    unresolved = Counter()
    for genome, data in combined.items():
        gff = f"{GFF_DIR}/{genome}.gff"
        locus_idx = build_gff_index(gff)[0] if os.path.exists(gff) else {}
        for g in data["genes"]:
            pid = g.get("protein_id") or locus_idx.get((g["contig"], g["start"], g["end"], g["strand"]))
            if pid:
                gene_to_proteins[g["gene"]][genome].append(pid)
            else:
                unresolved[g["gene"]] += 1

    selected, report = [], {}
    for gene, gp in sorted(gene_to_proteins.items()):
        n_multi = sum(1 for pids in gp.values() if len(set(pids)) > 1)
        if len(gp) < MIN_GENOME_COVERAGE or n_multi / len(gp) > MAX_MULTICOPY_FRAC:
            continue
        selected.append(gene)
        report[gene] = {"genomes_with_gene": len(gp), "genomes_multicopy_excluded": n_multi}
    print(f"Selected {len(selected)} genes: {selected}")

    # single-copy calls only
    calls = {gene: {gn: pids[0] for gn, pids in gene_to_proteins[gene].items() if len(set(pids)) == 1}
             for gene in selected}
    genes_per_genome = Counter(gn for gene in selected for gn in calls[gene])
    min_genes = int(round(MIN_GENE_FRAC * len(selected)))

    gtdb = Tree(GTDB_TREE, format=1)
    for leaf in gtdb:
        leaf.name = leaf.name.replace("RS_", "", 1).replace("GB_", "", 1)
    gtdb_taxa = set(gtdb.get_leaf_names())

    included = sorted(gn for gn, n in genes_per_genome.items()
                      if n >= min_genes and accession(gn) in gtdb_taxa)
    dropped = sorted(set(combined) - set(included))
    print(f"Keeping {len(included)} genomes with >= {min_genes}/{len(selected)} genes; dropping {len(dropped)}")

    prot_cache = {}
    missing_seq = Counter()
    for gene in selected:
        n_written = 0
        with open(os.path.join(OUTPUT_DIR, "genes", f"{gene}.faa"), "w") as out:
            for gn in included:
                pid = calls[gene].get(gn)
                if not pid:
                    continue
                if gn not in prot_cache:
                    prot_cache[gn] = read_fasta(f"{PROT_DIR}/{gn}.faa")
                seq = prot_cache[gn].get(pid)
                if not seq:
                    missing_seq[gene] += 1
                    continue
                out.write(f">{accession(gn)}\n{seq.rstrip('*')}\n")
                n_written += 1
        report[gene]["sequences_written"] = n_written
        report[gene]["protein_missing_from_faa"] = missing_seq[gene]
        print(f"  {gene}: {n_written} sequences")

    with open(os.path.join(OUTPUT_DIR, "genes.txt"), "w") as f:
        f.write("\n".join(selected) + "\n")

    with open(os.path.join(OUTPUT_DIR, "taxa.tsv"), "w") as f:
        f.write("accession\torganism\n")
        for gn in included:
            f.write(f"{accession(gn)}\t{combined[gn]['organism']}\n")

    gtdb.prune([accession(gn) for gn in included], preserve_branch_length=True)
    gtdb.write(outfile=os.path.join(OUTPUT_DIR, "gtdb_ref_pruned.nwk"), format=5)

    json.dump({
        "rule": {"min_genome_coverage": MIN_GENOME_COVERAGE,
                 "max_multicopy_frac": MAX_MULTICOPY_FRAC,
                 "min_genes_per_genome": min_genes},
        "genes": report,
        "genes_not_selected": {g: {"genomes": len(gp),
                                   "multicopy": sum(1 for p in gp.values() if len(set(p)) > 1)}
                               for g, gp in sorted(gene_to_proteins.items()) if g not in selected},
        "n_genomes_included": len(included),
        "genomes_dropped": {accession(gn): {"organism": combined[gn]["organism"],
                                            "selected_genes_carried": genes_per_genome.get(gn, 0)}
                            for gn in dropped},
        "gene_calls_without_protein_id": dict(unresolved),
    }, open(os.path.join(OUTPUT_DIR, "gene_selection.json"), "w"), indent=2)
    print(f"Wrote {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
