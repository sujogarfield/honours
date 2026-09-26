#!/usr/bin/env python3
"""
verify_gene_calls.py
Independent check of every protein that goes into the flagellar trees
(gene_order/flagellar_trees/genes/*.faa), by two methods:

1. Multi-reference best hit. Each protein is searched (DIAMOND blastp,
   very-sensitive, hit must cover >= 50 % of the query) against the COMPLETE proteomes of five reference genomes
   from different lineages. Its best hit's own NCBI annotation (gene= tag, or
   a gene symbol in the product) should name the same gene. Searching whole
   proteomes, not just flagellar genes, is what catches mix-ups between
   homologous families (FlgE/FlgF/FlgG/FlgK rod-hook proteins, FliM/FliN/FliY
   switch proteins).
2. eggNOG-mapper (if --eggnog is given): the protein's Preferred_name /
   Description should name the same gene (same rule as
   06_orthology/orthology_eggnog_confirm.py).

A protein is FLAGGED if, across the reference hits that carry a gene name,
more disagree than agree, or eggNOG names a different specific gene.

Run from the project root (needs `diamond` on PATH, e.g. the `orthology` env):
  python scripts/11_flagellar_trees/verify_gene_calls.py --write-queries   # writes queries.faa
  emapper.py -i gene_order/flagellar_trees/verify/queries.faa --itype proteins \
      -m diamond --data_dir <eggnog_data> --output tree_inputs \
      --output_dir gene_order/flagellar_trees/verify --cpu 8
  python scripts/11_flagellar_trees/verify_gene_calls.py \
      --eggnog gene_order/flagellar_trees/verify/tree_inputs.emapper.annotations
Output (gene_order/flagellar_trees/verify/):
  verification.tsv   one row per protein
  summary.tsv        one row per gene
"""

import argparse
import glob
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "06_orthology"))
from orthology_eggnog_confirm import classify, load_annotations

D = "gene_order/flagellar_trees"
OUT = f"{D}/verify"
PROT_DIR = "campy_fetched/campy_prot"
GFF_DIR = "campy_fetched/campy_ann"
REFERENCES = {
    "Campylobacter jejuni": "GCF_001457695.1_NCTC11351",
    "Helicobacter pylori": "GCF_900478295.1_35377_D02",
    "Aliarcobacter butzleri": "GCF_900187115.1_51699_F01",
    "Sulfurimonas denitrificans": "GCF_000012965.1_ASM1296v1",
    "Nautilia profundicola": "GCF_000021725.1_ASM2172v1",
}
MIN_QUERY_COVER = 50  # % of the query the best hit must span; below this it is a shared
                      # domain (e.g. the receiver domain of every response regulator), not the gene
SYMBOL = re.compile(r"\b(fl[a-z]{2}\d?|che[a-z]|pfl[a-z]|mot[a-z]|maf\d?)\b")


def write_queries():
    os.makedirs(OUT, exist_ok=True)
    n = 0
    with open(f"{OUT}/queries.faa", "w") as out:
        for path in sorted(glob.glob(f"{D}/genes/*.faa")):
            gene = os.path.basename(path)[:-4]
            name = None
            for line in open(path):
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    out.write(f">{gene}|{name}\n")
                    n += 1
                else:
                    out.write(line)
    print(f"Wrote {OUT}/queries.faa ({n} proteins)")


def reference_names(genome):
    """protein_id -> set of lowercase gene symbols from the GFF gene= tag and product."""
    names = defaultdict(set)
    for line in open(f"{GFF_DIR}/{genome}.gff"):
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9 or parts[2] != "CDS":
            continue
        attrs = parts[8]
        pid = re.search(r"Dbxref=[^;]*GenBank:([^,;]+)", attrs)
        if not pid:
            continue
        gene = re.search(r"(?:^|;)gene=([^;]+)", attrs)
        product = re.search(r"product=([^;]+)", attrs)
        if gene:
            names[pid.group(1)].add(gene.group(1).lower())
        if product:
            names[pid.group(1)].update(SYMBOL.findall(product.group(1).lower()))
    return names


def reference_hits():
    """(gene|accession) -> {reference: (hit protein, pident, names)}"""
    hits = defaultdict(dict)
    for ref, genome in REFERENCES.items():
        db = f"{OUT}/ref_{genome}"
        if not os.path.exists(db + ".dmnd"):
            subprocess.run(["diamond", "makedb", "--in", f"{PROT_DIR}/{genome}.faa", "-d", db, "--quiet"],
                           check=True)
        tsv = f"{OUT}/hits_{genome}.tsv"
        if not os.path.exists(tsv):
            subprocess.run(["diamond", "blastp", "-q", f"{OUT}/queries.faa", "-d", db, "-o", tsv,
                            "--very-sensitive", "--max-target-seqs", "1", "--evalue", "1e-5", "--quiet",
                            "--query-cover", str(MIN_QUERY_COVER),
                            "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue", "bitscore"], check=True)
        names = reference_names(genome)
        ref_acc = "_".join(genome.split("_")[:2])
        for line in open(tsv):
            q, s, pident = line.split("\t")[:3]
            if q in hits and ref in hits[q]:
                continue  # keep the first (best) HSP only
            if q.split("|")[1] == ref_acc:
                continue  # the reference hitting itself says nothing
            hits[q][ref] = (s, float(pident), names.get(s, set()))
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-queries", action="store_true")
    ap.add_argument("--eggnog", help="emapper .annotations file for queries.faa")
    args = ap.parse_args()
    if args.write_queries or not os.path.exists(f"{OUT}/queries.faa"):
        write_queries()
        if args.write_queries:
            return

    queries = [l[1:].strip() for l in open(f"{OUT}/queries.faa") if l.startswith(">")]
    hits = reference_hits()
    ann = load_annotations(args.eggnog) if args.eggnog else {}

    rows, per_gene = [], defaultdict(Counter)
    for q in queries:
        gene, acc = q.split("|")
        agree, disagree, unnamed = 0, [], 0
        for ref, (s, pident, names) in hits.get(q, {}).items():
            if gene in names:
                agree += 1
            elif names:
                disagree.append(f"{ref.split()[0][0]}.{ref.split()[1]}:{'/'.join(sorted(names))}")
            else:
                unnamed += 1
        n_ref = len(hits.get(q, {}))
        eg = ann.get(q)
        eg_status = classify(gene, eg["preferred_name"], eg["description"]) if eg else ("no_hit" if ann else "")
        flagged = len(disagree) > agree or eg_status == "reclassified"
        rows.append([gene, acc, n_ref, agree, len(disagree), unnamed, "; ".join(disagree),
                     eg_status, eg["preferred_name"] if eg else "", "FLAG" if flagged else ""])
        c = per_gene[gene]
        c["proteins"] += 1
        c["flagged"] += flagged
        c["ref_majority_agree"] += agree > len(disagree)
        c[f"eggnog_{eg_status or 'not_run'}"] += 1

    with open(f"{OUT}/verification.tsv", "w") as f:
        f.write("gene\taccession\tref_hits\tref_agree\tref_disagree\tref_unnamed\tdisagreeing_refs\t"
                "eggnog_status\teggnog_name\tflag\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    with open(f"{OUT}/summary.tsv", "w") as f:
        f.write("gene\tproteins\tref_majority_agree\teggnog_confirmed\teggnog_reclassified\t"
                "eggnog_unconfirmed\teggnog_no_hit\tflagged\n")
        for gene, c in sorted(per_gene.items()):
            f.write(f"{gene}\t{c['proteins']}\t{c['ref_majority_agree']}\t{c['eggnog_confirmed']}\t"
                    f"{c['eggnog_reclassified']}\t{c['eggnog_unconfirmed']}\t{c['eggnog_no_hit']}\t"
                    f"{c['flagged']}\n")

    total = len(rows)
    n_flag = sum(r[-1] == "FLAG" for r in rows)
    n_agree = sum(r[3] > r[4] for r in rows)
    print(f"{total} proteins; reference majority agrees for {n_agree} ({100 * n_agree / total:.1f}%); "
          f"{n_flag} flagged")
    if ann:
        st = Counter(r[7] for r in rows)
        print("eggNOG:", dict(st))
    for gene, c in sorted(per_gene.items(), key=lambda kv: -kv[1]["flagged"]):
        if c["flagged"]:
            print(f"  {gene:5s} {c['flagged']:3d} flagged of {c['proteins']}")
    print(f"Wrote {OUT}/verification.tsv and {OUT}/summary.tsv")


if __name__ == "__main__":
    main()
