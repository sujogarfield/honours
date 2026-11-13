#!/usr/bin/env python3
"""
synphy.py
The thesis flagellar-gene pipeline for any GTDB taxon and gene list.

    python synphy/synphy.py <step> <run_dir>

Login node (needs internet):
    init          interactive setup -> run_dir/config.json
    genomes       GTDB TSV -> genomes.tsv (filters, caps, outgroup)
    fetch         NCBI genome FASTA, GFF, proteome (+ seed reference proteome)
    tree          GTDB tree pruned to the genomes
    check         how many genomes each gene name is found in
    string-match  gene= tag matches

Katana jobs (katana/submit.sh runs these for you):
    crossref      OrthoFinder orthogroups -> candidate genes (+ eggNOG input)
    confirm       eggNOG agreement -> gene_order_combined.json
    select        copy rules -> per-gene protein FASTA for IQ-TREE
    synteny       gene order vs phylogeny, Mantel tests
    report        interactive tree, synteny scatter, tanglegram
    hgt-screen    gene-tree branches (UFBoot >= 95) that contradict GTDB
    hgt-shortlist thesis shortlist rule -> trees/hgt_shortlist.tsv

See synphy/README.md for the full walkthrough.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

STEPS = {
    "init": ("synphy_lib.genomes", "init"),
    "genomes": ("synphy_lib.genomes", "genomes"),
    "fetch": ("synphy_lib.genomes", "fetch"),
    "tree": ("synphy_lib.genomes", "tree"),
    "check": ("synphy_lib.orthology", "check"),
    "string-match": ("synphy_lib.orthology", "string_match"),
    "crossref": ("synphy_lib.orthology", "crossref"),
    "confirm": ("synphy_lib.orthology", "confirm"),
    "select": ("synphy_lib.select_genes", "select"),
    "synteny": ("synphy_lib.synteny", "synteny"),
    "report": ("synphy_lib.report", "report"),
    "hgt-screen": ("synphy_lib.hgt", "screen"),
    "hgt-shortlist": ("synphy_lib.hgt", "shortlist"),
}


def main():
    if len(sys.argv) < 3 or sys.argv[1] not in STEPS:
        print(__doc__)
        sys.exit(1)
    module, func = STEPS[sys.argv[1]]
    getattr(__import__(module, fromlist=[func]), func)(sys.argv[2], *sys.argv[3:])


if __name__ == "__main__":
    main()
