#!/usr/bin/env python3
"""
summarise_hgt_spikein.py
Scores the planted-transfer positive control (see gene_order/hgt_spikein/README.md).
Run after the gene trees are back and the screen has been run on them:
  FLAG_DIR=gene_order/hgt_spikein python scripts/11_flagellar_trees/screen_hgt_candidates.py
  python scripts/11_flagellar_trees/summarise_hgt_spikein.py

For each test gene tree, with R = recipient and D = donor:
  donor lineage  = the GTDB clade on D's side of the R-D split (the child of
                   MRCA(R, D) that contains D), minus D itself;
  home lineage   = the GTDB clade on R's side, minus R;
  sister         = the smaller of the two subtrees R joins in the unrooted
                   gene tree (R's nearest relatives);
  placement      = 'donor' if every sister is in the donor lineage, 'home' if
                   every sister is in the home lineage, else 'elsewhere';
  sister_ufboot  = UFBoot of the R + sister branch;
  screen_hit     = the HGT screen (UFBoot >= 95, incompatible with GTDB) has a
                   candidate group made of R plus only donor-lineage genomes.
A transfer counts as DETECTED when the spike tree has a screen hit and its
control tree does not. Writes gene_order/hgt_spikein/results.tsv.
"""
import csv
import json
from pathlib import Path

from ete3 import Tree

D = Path("gene_order/hgt_spikein")


def sister_of(tree, leaf_name):
    """(sister leaf set, UFBoot of the leaf+sister branch) in the unrooted tree."""
    t = tree.copy()
    leaf = t & leaf_name
    # root on a leaf far from R so R's parent node is a proper (unrooted) junction
    far = max(t.get_leaves(), key=lambda l: t.get_distance(leaf, l, topology_only=True))
    t.set_outgroup(far)
    leaf = t & leaf_name
    parent = leaf.up
    sibs = [c for c in parent.children if c is not leaf]
    sister = set().union(*(set(c.get_leaf_names()) for c in sibs))
    rest = set(t.get_leaf_names()) - sister - {leaf_name}
    if len(rest) < len(sister):  # take the smaller side as the sister group
        sister = rest
    return sister, parent.support


def lineages(gtdb, rec, don, present):
    mrca = gtdb.get_common_ancestor(rec, don)
    side = lambda acc: next(set(c.get_leaf_names()) for c in mrca.children if acc in c.get_leaf_names())
    donor = (side(don) - {don}) & present
    home = (side(rec) - {rec}) & present
    return donor, home


def main():
    design = list(csv.DictReader(open(D / "design.tsv"), delimiter="\t"))
    taxa = {l.split("\t")[0]: l.split("\t")[1] for l in (D / "taxa.tsv").read_text().splitlines()[1:]}
    gtdb = Tree(str(D / "gtdb_ref_pruned.nwk"), format=1)
    screen = json.load(open(D / "hgt_screen.json"))
    out = []
    for r in design:
        tid, rec, don = r["test_id"], r["recipient_acc"], r["donor_acc"]
        tf = D / "gene_trees" / f"{tid}.treefile"
        if not tf.exists():
            print(f"missing {tf}")
            continue
        gt = Tree(str(tf), format=0)
        present = set(gt.get_leaf_names())
        donor, home = lineages(gtdb, rec, don, present)
        sister, sup = sister_of(gt, rec)
        placement = "donor" if sister <= donor else "home" if sister <= home else "elsewhere"
        hits = [c for c in screen["candidates"] if c["tree"] == tid and rec in c["group"]
                and set(c["group"]) - {rec} and set(c["group"]) - {rec} <= donor]
        out.append({**r, "placement": placement, "sister_ufboot": sup,
                    "sister": "; ".join(sorted(taxa.get(s, s) for s in sister)) if len(sister) <= 4
                    else f"{len(sister)} genomes, e.g. " + "; ".join(sorted({taxa.get(s, s) for s in sister})[:3]),
                    "screen_hit": bool(hits),
                    "screen_hit_ufboot": max((c["ufboot"] for c in hits), default="")})

    by = {(o["gene"], o["transfer_distance"], o["condition"]): o for o in out}
    for o in out:
        if o["condition"] == "spike":
            ctrl = by.get((o["gene"], o["transfer_distance"], "control"))
            o["detected"] = o["screen_hit"] and ctrl is not None and not ctrl["screen_hit"]
        else:
            o["detected"] = ""
    cols = list(design[0].keys()) + ["placement", "sister_ufboot", "sister", "screen_hit",
                                     "screen_hit_ufboot", "detected"]
    with open(D / "results.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for o in out:
            f.write("\t".join(str(o[c]) for c in cols) + "\n")

    genes = list(dict.fromkeys(o["gene"] for o in out))
    dists = list(dict.fromkeys(o["transfer_distance"] for o in out))
    print("Spike trees: placement of the recipient / UFBoot / screen detection (control flagged?)")
    print(f"{'':16}" + "".join(f"{g:>34}" for g in genes))
    for d in dists:
        cells = []
        for g in genes:
            s, c = by.get((g, d, "spike")), by.get((g, d, "control"))
            if not s:
                cells.append(f"{'-':>34}")
                continue
            ctrl = "ctrl FLAGGED" if c and c["screen_hit"] else "ctrl clean"
            cells.append(f"{s['placement']:>9} {s['sister_ufboot']:>4.0f} "
                         f"{'DETECTED' if s['detected'] else 'missed':>8} ({ctrl})")
        print(f"{d:16}" + "".join(f"{x:>34}" for x in cells))
    n = [o for o in out if o["condition"] == "spike"]
    print(f"\nDetected {sum(bool(o['detected']) for o in n)}/{len(n)}; "
          f"placed with donor lineage {sum(o['placement'] == 'donor' for o in n)}/{len(n)}; "
          f"controls flagged {sum(o['screen_hit'] for o in out if o['condition'] == 'control')}")
    print(f"Wrote {D}/results.tsv")


if __name__ == "__main__":
    main()
