#!/usr/bin/env python3
"""
summarise_hgt_spikein.py
Scores the planted-transfer positive control (see gene_order/hgt_spikein/README.md).
Run after the gene trees are back and the screen has been run on them:
  FLAG_DIR=gene_order/hgt_spikein python scripts/11_flagellar_trees/screen_hgt_candidates.py
  python scripts/11_flagellar_trees/summarise_hgt_spikein.py

For each test gene tree, with R = recipient and D = donor:
  donor side      = the GTDB clade on D's side of the R-D split (the child of
                    MRCA(R, D) containing D), minus D;
  donor relatives = D's closest GTDB group with >= 2 genomes left once D is
                    removed (design.tsv, from prepare_hgt_spikein.py);
  sister          = the smaller of the two subtrees R joins in the unrooted
                    gene tree (R's nearest relatives);
  placement       = 'donor' if the sister lies within the donor side and
                    includes at least one donor relative, 'home' if the sister
                    lies within R's own side, else 'elsewhere';
  sister_ufboot   = UFBoot of the R + sister branch;
  screen_hit      = the HGT screen (UFBoot >= 95, incompatible with GTDB) has a
                    candidate group made of R plus donor-side genomes only,
                    including at least one donor relative;
  rank            = position of that candidate among all of the tree's
                    candidates, ordered by depth (GTDB span - group size), then
                    UFBoot -- 1 means it is the deepest conflict in the tree.
A transfer counts as DETECTED when the spike tree has a screen hit and its
control tree does not. Writes gene_order/hgt_spikein/results.tsv.
"""
import csv
import json
from math import comb
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
    sister = set().union(*(set(c.get_leaf_names()) for c in parent.children if c is not leaf))
    rest = set(t.get_leaf_names()) - sister - {leaf_name}
    if len(rest) < len(sister):  # take the smaller side as the sister group
        sister = rest
    return sister, parent.support


def sides(gtdb, rec, don, present):
    mrca = gtdb.get_common_ancestor(rec, don)
    side = lambda acc: next(set(c.get_leaf_names()) for c in mrca.children if acc in c.get_leaf_names())
    return (side(don) - {don}) & present, (side(rec) - {rec}) & present


def clopper_pearson_low(k, n, alpha=0.05):
    """One-sided-per-tail exact binomial lower bound (95% two-sided)."""
    if k == 0:
        return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(60):
        p = (lo + hi) / 2
        tail = sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))
        lo, hi = (p, hi) if tail < alpha / 2 else (lo, p)
    return lo


def main():
    design = list(csv.DictReader(open(D / "design.tsv"), delimiter="\t"))
    taxa = {l.split("\t")[0]: l.split("\t")[1] for l in (D / "taxa.tsv").read_text().splitlines()[1:]}
    gtdb = Tree(str(D / "gtdb_ref_pruned.nwk"), format=1)
    screen = json.load(open(D / "hgt_screen.json"))
    rank_key = lambda c: (-(c["gtdb_span"] - c["group_size"]), -c["ufboot"])
    out = []
    for r in design:
        tid, rec, don = r["test_id"], r["recipient_acc"], r["donor_acc"]
        tf = D / "gene_trees" / f"{tid}.treefile"
        if not tf.exists():
            print(f"missing {tf}")
            continue
        gt = Tree(str(tf), format=0)
        present = set(gt.get_leaf_names())
        donor_side, home_side = sides(gtdb, rec, don, present)
        relatives = set(r["donor_relatives"].split(","))
        sister, sup = sister_of(gt, rec)
        placement = ("donor" if sister <= donor_side and sister & relatives
                     else "home" if sister <= home_side else "elsewhere")
        cands = sorted((c for c in screen["candidates"] if c["tree"] == tid), key=rank_key)
        hit = next((i for i, c in enumerate(cands) if rec in c["group"]
                    and set(c["group"]) - {rec} <= donor_side and set(c["group"]) & relatives), None)
        out.append({**r, "placement": placement, "sister_ufboot": sup,
                    "sister": "; ".join(sorted(taxa.get(s, s) for s in sister)) if len(sister) <= 4
                    else f"{len(sister)} genomes, e.g. " + "; ".join(sorted({taxa.get(s, s) for s in sister})[:3]),
                    "screen_hit": hit is not None,
                    "screen_hit_ufboot": cands[hit]["ufboot"] if hit is not None else "",
                    "rank": f"{hit + 1}/{len(cands)}" if hit is not None else "",
                    "n_conflicts": len(cands)})

    by = {(o["gene"], o["recipient_code"], o["donor_level"], o["condition"]): o for o in out}
    for o in out:
        if o["condition"] == "spike":
            ctrl = by.get((o["gene"], o["recipient_code"], o["donor_level"], "control"))
            o["detected"] = o["screen_hit"] and ctrl is not None and not ctrl["screen_hit"]
        else:
            o["detected"] = ""
    cols = list(design[0].keys()) + ["placement", "sister_ufboot", "sister", "screen_hit",
                                     "screen_hit_ufboot", "rank", "n_conflicts", "detected"]
    with open(D / "results.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for o in out:
            f.write("\t".join(str(o[c]) for c in cols) + "\n")

    genes = list(dict.fromkeys(o["gene"] for o in out))
    levels = list(dict.fromkeys(o["donor_level"] for o in out))
    for code in dict.fromkeys(o["recipient_code"] for o in out):
        any_s = next(o for o in out if o["recipient_code"] == code)
        print(f"\n=== recipient {any_s['recipient']}: cell = detected/missed, rank among the tree's conflicts")
        print(f"{'donor level':15}{'donor (GTDB dist)':38}" + "".join(f"{g:>16}" for g in genes))
        for lv in levels:
            s0 = by.get((genes[0], code, lv, "spike"))
            if not s0:
                continue
            cells = []
            for g in genes:
                s, c = by.get((g, code, lv, "spike")), by.get((g, code, lv, "control"))
                if not s:
                    cells.append("-")
                    continue
                txt = (f"DET {s['rank']}" if s["detected"] else
                       "ctrl flagged" if c and c["screen_hit"] else f"missed ({s['placement']})")
                cells.append(txt)
            print(f"{lv:15}{s0['donor'] + ' (' + s0['gtdb_distance'] + ')':38}" + "".join(f"{x:>16}" for x in cells))

    spikes = [o for o in out if o["condition"] == "spike"]
    k = sum(bool(o["detected"]) for o in spikes)
    top = sum(o["detected"] and o["rank"].startswith("1/") for o in spikes)
    print(f"\nDetected {k}/{len(spikes)} (95% CI lower bound {100 * clopper_pearson_low(k, len(spikes)):.0f}%); "
          f"deepest conflict in its tree for {top}/{k}; "
          f"placed with donor's relatives {sum(o['placement'] == 'donor' for o in spikes)}/{len(spikes)}; "
          f"controls flagged {sum(o['screen_hit'] for o in out if o['condition'] == 'control')}/"
          f"{sum(o['condition'] == 'control' for o in out)}")
    for lv in levels:
        sub = [o for o in spikes if o["donor_level"] == lv]
        print(f"  {lv:14} {sum(bool(o['detected']) for o in sub)}/{len(sub)}")
    print(f"Wrote {D}/results.tsv")


if __name__ == "__main__":
    main()
