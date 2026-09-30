#!/usr/bin/env python3
"""
prepare_core_backbone.py
A second, GTDB-independent genome backbone for the candidate tests: the
OrthoFinder single-copy orthologues (present exactly once in all 166
genomes), minus anything flagellar/chemotaxis, for the 149 tree genomes.

Why: the species-marker test (M11) uses GTDB's own bac120 alignment. This
backbone comes from our own proteomes and gene calls, a different gene set
(overlap with bac120 is partial -- both are mostly universal single-copy
genes such as ribosomal proteins), and its own alignment and trimming. If both
backbones reject a flagellar grouping, the discordance does not hinge on
GTDB's marker set.

Orthogroups are EXCLUDED if any member protein is a flagellar/chemotaxis call
in gene_order_combined.json, or if any member's RefSeq product mentions
flagell-/chemotaxis/motility. Members are mapped to genomes through
Orthogroups.tsv (identical RefSeq WP_ proteins share one ID across genomes, so
FASTA headers alone are ambiguous) and sequences read from each genome's own
proteome.

Run from the project root (needs campy_orthologs/, local only):
  python scripts/11_flagellar_trees/prepare_core_backbone.py
Output (gene_order/flagellar_trees/core/): genes/<OG>.faa (unaligned, header =
accession), core_genes.tsv, candidates.txt + <cand>/{flag,gtdb}_constraint.nwk
(copied from markers/, same 149 taxa)
"""

import csv
import glob
import json
import os
import re
import shutil
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "06_orthology"))
from orthology_crossref import build_gff_index

D = "gene_order/flagellar_trees"
OUT = f"{D}/core"
OF = "campy_orthologs/orthofinder_results/Results_Sep23"
EXCLUDE_PRODUCT = re.compile(r"flagell|chemotaxis|motility|\bmot[AB]\b", re.I)


def accession(genome_key):
    return "_".join(genome_key.split("_")[:2])


def main():
    taxa = [l.split("\t")[0] for l in open(f"{D}/taxa.tsv").read().splitlines()[1:]]
    sc_ogs = sorted(os.path.basename(p)[:-3] for p in glob.glob(f"{OF}/Single_Copy_Orthologue_Sequences/*.fa"))

    # OG -> {genome key: protein}
    members = {}
    with open(f"{OF}/Orthogroups/Orthogroups.tsv") as f:
        reader = csv.reader(f, delimiter="\t")
        genomes = next(reader)[1:]
        wanted = set(sc_ogs)
        for row in reader:
            if row[0] in wanted:
                members[row[0]] = {g: c.strip() for g, c in zip(genomes, row[1:]) if c.strip()}

    # flagellar protein ids
    combined = json.load(open("gene_order/gene_order_combined.json"))
    flag_pids = set()
    for gk, d in combined.items():
        idx = build_gff_index(f"campy_fetched/campy_ann/{gk}.gff")[0]
        for x in d["genes"]:
            pid = x.get("protein_id") or idx.get((x["contig"], x["start"], x["end"], x["strand"]))
            if pid:
                flag_pids.add(pid)

    # products
    product = {}
    for gff in glob.glob("campy_fetched/campy_ann/*.gff"):
        for line in open(gff):
            p = line.split("\t")
            if len(p) > 8 and p[2] == "CDS":
                pid = re.search(r"GenBank:([^,;]+)", p[8])
                pr = re.search(r"product=([^;]+)", p[8])
                if pid and pr:
                    product[pid.group(1)] = pr.group(1)

    kept, excluded = [], []
    for og in sc_ogs:
        m = members[og]
        prods = Counter(product.get(p, "?") for p in m.values())
        if any(p in flag_pids for p in m.values()) or any(EXCLUDE_PRODUCT.search(x) for x in prods):
            excluded.append((og, prods.most_common(1)[0][0]))
        else:
            kept.append((og, prods.most_common(1)[0][0]))

    os.makedirs(f"{OUT}/genes", exist_ok=True)
    for old in glob.glob(f"{OUT}/genes/*.faa"):
        os.remove(old)
    by_acc = {accession(g): g for g in genomes}
    prot_cache = {}

    def seqs_of(gk):
        if gk not in prot_cache:
            s, name = {}, None
            for line in open(f"campy_fetched/campy_prot/{gk}.faa"):
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    s[name] = []
                else:
                    s[name].append(line.strip())
            prot_cache[gk] = {k: "".join(v).rstrip("*") for k, v in s.items()}
        return prot_cache[gk]

    total = 0
    for og, _ in kept:
        with open(f"{OUT}/genes/{og}.faa", "w") as out:
            for acc in taxa:
                gk = by_acc[acc]
                seq = seqs_of(gk)[members[og][gk]]
                out.write(f">{acc}\n{seq}\n")
                total += len(seq)

    with open(f"{OUT}/core_genes.tsv", "w") as f:
        f.write("orthogroup\tstatus\tcommonest_product\n")
        for og, pr in kept:
            f.write(f"{og}\tkept\t{pr}\n")
        for og, pr in excluded:
            f.write(f"{og}\texcluded_flagellar\t{pr}\n")

    # candidate constraints: identical taxon set to markers/, so reuse them
    shutil.copy(f"{D}/markers/candidates.txt", f"{OUT}/candidates.txt")
    for cand in open(f"{D}/markers/candidates.txt").read().split():
        os.makedirs(f"{OUT}/{cand}", exist_ok=True)
        for fn in ("flag_constraint.nwk", "gtdb_constraint.nwk", "group.txt"):
            shutil.copy(f"{D}/markers/{cand}/{fn}", f"{OUT}/{cand}/{fn}")

    ribo = sum(1 for _, p in kept if "ribosom" in p.lower())
    print(f"{len(sc_ogs)} single-copy orthogroups; {len(excluded)} excluded as flagellar/chemotaxis; "
          f"{len(kept)} kept ({ribo} ribosomal-protein), {len(taxa)} genomes, ~{total // len(taxa)} aa per genome")
    for og, pr in excluded:
        print(f"  excluded {og}: {pr}")
    print(f"Wrote {OUT}/")


if __name__ == "__main__":
    main()
