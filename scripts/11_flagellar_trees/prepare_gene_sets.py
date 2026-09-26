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
  - copies are resolved with OrthoFinder: each gene's main orthogroup is the
    one most of its proteins fall in. A genome keeps its copy if exactly one
    is in the main orthogroup (a second, differently-clustered copy is a
    paralog or mis-named protein); a lone copy OUTSIDE the main orthogroup is
    dropped (named like the gene, but OrthoFinder says it is a different
    family -- e.g. the 8 thermophile "flgG" proteins in their own OG);
  - genomes with 2+ copies IN the main orthogroup are left out of that gene
    only (which copy is the ortholog is ambiguous), not dropped from the study;
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

import csv
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
ORTHOGROUPS_TSV = "campy_orthologs/orthofinder_results/Results_Sep23/Orthogroups/Orthogroups.tsv"
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


def load_protein_to_og():
    p2og = {}
    with open(ORTHOGROUPS_TSV) as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for row in reader:
            for cell in row[1:]:
                for p in cell.split(","):
                    if p.strip():
                        p2og[p.strip()] = row[0]
    return p2og


def resolve_copies(gene_to_proteins, p2og):
    """gene -> (main OG, {genome: protein | None}, {genome: reason}) where None
    means 2+ copies in the main OG (ambiguous) and reasons record every
    genome whose call was changed or dropped."""
    out = {}
    for gene, gp in gene_to_proteins.items():
        og_counts = Counter(p2og.get(p) for pids in gp.values() for p in set(pids))
        og_counts.pop(None, None)
        main_og = og_counts.most_common(1)[0][0] if og_counts else None
        calls, notes = {}, {}
        for genome, pids in gp.items():
            pids = sorted(set(pids))
            in_main = [p for p in pids if p2og.get(p) == main_og]
            if len(in_main) == 1:
                calls[genome] = in_main[0]
                if len(pids) > 1:
                    notes[genome] = "multi-copy, resolved to main-OG copy"
            elif len(in_main) > 1:
                calls[genome] = None
                notes[genome] = "2+ copies in main OG (ambiguous), excluded"
            else:
                notes[genome] = f"only copy outside main OG ({p2og.get(pids[0])}), dropped"
        out[gene] = (main_og, calls, notes)
    return out


def main():
    os.makedirs(os.path.join(OUTPUT_DIR, "genes"), exist_ok=True)
    for old in glob.glob(os.path.join(OUTPUT_DIR, "genes", "*.faa")):
        os.remove(old)  # a gene that drops out of the selection must not leave a stale input
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

    resolved = resolve_copies(gene_to_proteins, load_protein_to_og())

    selected, report = [], {}
    for gene, (main_og, calls_g, notes) in sorted(resolved.items()):
        n_multi = sum(1 for v in calls_g.values() if v is None)
        if len(calls_g) < MIN_GENOME_COVERAGE or n_multi / len(calls_g) > MAX_MULTICOPY_FRAC:
            continue
        selected.append(gene)
        report[gene] = {
            "main_orthogroup": main_og,
            "genomes_with_gene": len(calls_g),
            "genomes_multicopy_excluded": n_multi,
            "copy_resolution": {accession(g): r for g, r in sorted(notes.items())},
        }
    print(f"Selected {len(selected)} genes: {selected}")

    # single-copy (or OG-resolved) calls only
    calls = {gene: {gn: pid for gn, pid in resolved[gene][1].items() if pid} for gene in selected}
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
        "genes_not_selected": {g: {"genomes": len(r[1]),
                                   "multicopy_in_main_og": sum(1 for v in r[1].values() if v is None)}
                               for g, r in sorted(resolved.items()) if g not in selected},
        "n_genomes_included": len(included),
        "genomes_dropped": {accession(gn): {"organism": combined[gn]["organism"],
                                            "selected_genes_carried": genes_per_genome.get(gn, 0)}
                            for gn in dropped},
        "gene_calls_without_protein_id": dict(unresolved),
    }, open(os.path.join(OUTPUT_DIR, "gene_selection.json"), "w"), indent=2)
    print(f"Wrote {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
