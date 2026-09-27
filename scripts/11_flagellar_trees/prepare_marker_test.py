#!/usr/bin/env python3
"""
prepare_marker_test.py
Inputs for the species-marker test of the AU candidates (katana/05_marker_*.pbs):
does the genome backbone -- GTDB's own 120 bac120 marker genes -- agree with
GTDB or with the flagellar genes at each candidate?

The marker alignment is GTDB R232's masked concatenated bac120 MSA for
representative genomes (the data the GTDB tree was built from), reduced to the
149 genomes in the flagellar trees. It is committed
(markers/bac120_msa_r232.faa), so Katana needs no download. To regenerate:
  curl -s https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/bac120_msa_reps_r232.faa.gz \\
    | gunzip -c | <keep the rows whose accession, minus RS_/GB_, is in taxa.tsv>

Per candidate (same list as prepare_au_tests.py), two constraints on the
marker alignment:
  flag_constraint.nwk  the flagellar grouping forced (group | everything else)
  gtdb_constraint.nwk  the GTDB branch(es) the grouping contradicts forced
Katana then builds the unconstrained marker tree, the best tree under each
constraint (3 seeds), and runs AU on {free, flagellar-forced, GTDB-forced}.

Run from the project root:  python scripts/11_flagellar_trees/prepare_marker_test.py
Output: gene_order/flagellar_trees/markers/{candidates.txt, <cand>/*.nwk, <cand>/group.txt}
"""

import os
import sys

from ete3 import Tree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prepare_au_tests import CANDIDATES, constraint_tree

D = "gene_order/flagellar_trees"
M = f"{D}/markers"
MSA = f"{M}/bac120_msa_r232.faa"


def main():
    org = {l.split("\t")[0]: l.split("\t")[1] for l in open(f"{D}/taxa.tsv").read().splitlines()[1:]}
    taxa = [l[1:].strip() for l in open(MSA) if l.startswith(">")]
    missing = set(org) - set(taxa)
    if missing:
        sys.exit(f"{len(missing)} genomes missing from {MSA}: {sorted(missing)[:5]}")
    gtdb = Tree(f"{D}/gtdb_ref_pruned.nwk", format=1)
    by_name = {o: a for a, o in org.items()}

    names = []
    for cand, spec in CANDIDATES.items():
        if isinstance(spec, str):
            group = sorted(a for a, o in org.items() if o.split()[0] == spec)
        else:
            group = sorted(by_name[o] for o in spec)
        rest = sorted(set(taxa) - set(group))
        os.makedirs(f"{M}/{cand}", exist_ok=True)
        with open(f"{M}/{cand}/flag_constraint.nwk", "w") as f:
            f.write(f"(({','.join(group)}),{','.join(rest)});\n")
        gnwk, forced = constraint_tree(gtdb, taxa, group)
        with open(f"{M}/{cand}/gtdb_constraint.nwk", "w") as f:
            f.write(gnwk + "\n")
        with open(f"{M}/{cand}/group.txt", "w") as f:
            for a in group:
                f.write(f"{a}\t{org[a]}\n")
        names.append(cand)
        print(f"{cand}: {len(group)} genomes; GTDB constraint forces {forced} branch(es)")

    open(f"{M}/candidates.txt", "w").write("\n".join(names) + "\n")
    print(f"Wrote {M}/ ({len(taxa)} genomes x {len(open(MSA).read().splitlines()[1])} marker sites)")


if __name__ == "__main__":
    main()
