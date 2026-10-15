#!/usr/bin/env python3
"""
plot_hgt_spikein.py
Tanglegrams for the planted-transfer positive control (gene_order/hgt_spikein,
see its README): GTDB species tree (pruned to the genomes in the gene file) vs
each test gene tree, with the recipient's connecting line highlighted. Shows
what a known horizontal transfer looks like in the same kind of figure as
fig1_tanglegram.png.

Run from the project root after summarise_hgt_spikein.py:
  python scripts/11_flagellar_trees/plot_hgt_spikein.py
Output: gene_order/hgt_spikein/figures/<test_id>.png (planted transfers), plus
        gene_order/hgt_spikein/figures/overview_<gene>_<recipient>.png (control vs spike, 4 donor levels)
"""
import csv
import os
import sys

os.environ["FLAG_DIR"] = "gene_order/hgt_spikein"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from ete3 import Tree

from plot_results import (TEXT, TEXT2, GRID, ACCENT, GROUPS, group_color, load_taxa,
                          phylogram_layout, draw_tree)
from screen_hgt_candidates import bitmask_splits

D = "gene_order/hgt_spikein"
OUT = f"{D}/figures"
OUTGROUP_GENERA = {"Hippea", "Desulfurella"}


def rows(tree):
    return {name: -i for i, name in enumerate(tree.get_leaf_names())}


def untangle(gtdb, gene, taxa):
    """Rotate child order (never topology) to reduce crossings, as in fig1."""
    def rotate_to(tree, target):
        for nd in tree.traverse("postorder"):
            if not nd.is_leaf():
                nd.children.sort(key=lambda c: -np.mean([target[l] for l in c.get_leaf_names()]))

    def crossings(a, b):
        seq = [-b[t] for t in sorted(taxa, key=lambda t: -a[t])]
        return sum(1 for i in range(len(seq)) for j in range(i + 1, len(seq)) if seq[i] > seq[j])

    gtdb.ladderize()
    best = None
    for _ in range(8):
        rotate_to(gene, rows(gtdb))
        c = crossings(rows(gtdb), rows(gene))
        if best is None or c < best[0]:
            best = (c, gtdb.write(format=1), gene.write(format=0))
        rotate_to(gtdb, rows(gene))
        c = crossings(rows(gtdb), rows(gene))
        if c < best[0]:
            best = (c, gtdb.write(format=1), gene.write(format=0))
    return Tree(best[1], format=1), Tree(best[2], format=0)


def load_pair(r, org):
    gene = Tree(f"{D}/gene_trees/{r['test_id']}.treefile", format=0)
    taxa = sorted(gene.get_leaf_names())
    gtdb = Tree(f"{D}/gtdb_ref_pruned.nwk", format=1)
    gtdb.prune(taxa, preserve_branch_length=True)
    out = [t for t in taxa if org[t].split()[0] in OUTGROUP_GENERA]
    anc = gene.get_common_ancestor(out) if len(out) > 1 else (gene & out[0] if out else None)
    gene.set_outgroup(anc if anc is not None and anc is not gene else gene.get_midpoint_outgroup())
    src = {name: set(bitmask_splits(t.copy(), {x: i for i, x in enumerate(taxa)}))
           for name, t in (("gtdb", gtdb), ("gene", gene))}
    gtdb, gene = untangle(gtdb, gene, taxa)
    index = {x: i for i, x in enumerate(taxa)}
    if set(bitmask_splits(gtdb, index)) != src["gtdb"] or set(bitmask_splits(gene, index)) != src["gene"]:
        raise SystemExit(f"{r['test_id']}: drawn topology differs from source")
    return gtdb, gene, taxa


def draw(ax, r, org, compact=False):
    gtdb, gene, taxa = load_pair(r, org)
    rec, don = r["recipient_acc"], r["donor_acc"]
    y_g, y_f = rows(gtdb), rows(gene)
    pos_g, _ = phylogram_layout(gtdb, y_g)
    pos_f, _ = phylogram_layout(gene, y_f)
    tree_w, lab_w, gap = (3.0, 0.0, 2.0) if compact else (3.0, 2.6, 2.2)
    xl_tip, xr_tip = tree_w, tree_w + 2 * lab_w + gap

    rec_path = {gene & rec}
    style_g = lambda n: (TEXT2, 0.6, "-") if n.is_leaf() else (TEXT, 0.8, "-")
    style_f = lambda n: (ACCENT, 2.0, "-") if n in rec_path else \
        ((TEXT2, 0.6, "-") if n.is_leaf() else (TEXT, 0.8, "-"))
    draw_tree(ax, gtdb, pos_g, 0, tree_w, False, style_g)
    draw_tree(ax, gene, pos_f, xr_tip, tree_w, True, style_f)
    for nd in gene.traverse():
        if not nd.is_leaf() and not nd.is_root() and nd.support < 95:
            x, y = pos_f[nd]
            ax.plot(xr_tip + (1 - x) * tree_w, y, "o", ms=2.4 if compact else 3.0,
                    mfc="white", mec=TEXT, mew=0.6, zorder=5)

    for t in taxa:
        yl, yr = y_g[t], y_f[t]
        hi = t == rec
        if compact and hi:  # one short label, outside the gene tree
            note = f" ({org[don].split()[0][0]}. {org[don].split()[1]} gene)" if r["condition"] == "spike" else ""
            ax.text(xr_tip + tree_w + 0.08, yr, f"{short(org[rec])}{note}", va="center", ha="left",
                    fontsize=8, color=ACCENT, weight="bold")
        elif not compact:
            fs = 5.4
            lab = org[t] + (f"  ← carries {org[don]}'s gene" if hi and r["condition"] == "spike" else "")
            ax.text(xl_tip + 0.06, yl, org[t], va="center", ha="left", fontsize=fs,
                    color=ACCENT if hi else TEXT, weight="bold" if hi else "normal")
            ax.text(xr_tip - 0.06, yr, lab, va="center", ha="right", fontsize=fs,
                    color=ACCENT if hi else TEXT, weight="bold" if hi else "normal")
        ax.plot([xl_tip + lab_w, xr_tip - lab_w], [yl, yr],
                color=ACCENT if hi else group_color(org[t]),
                lw=2.4 if hi else 0.6, alpha=1 if hi else 0.55, zorder=6 if hi else 2)
    n = len(taxa)
    ax.set_xlim(-0.1, xr_tip + tree_w + (2.6 if compact else 0.1))
    ax.set_ylim(-n - 2, 4)
    ax.axis("off")
    return xr_tip, tree_w, n


def short(name):
    return f"{name.split()[0][0]}. {name.split()[1]}"


def title_of(r):
    level = r["donor_level"].replace("_", " ")
    if r["condition"] == "spike":
        return (f"{r['gene']}  ·  planted transfer ({level}, GTDB distance {r['gtdb_distance']}): "
                f"{r['donor']} → {r['recipient']}")
    return f"{r['gene']}  ·  control ({r['donor']} removed, no transfer)"


def legend(fig, recipient):
    handles = [plt.Line2D([], [], color=c, lw=3) for _, c, _ in GROUPS]
    handles.append(plt.Line2D([], [], color=ACCENT, lw=3))
    labels = [g for g, _, _ in GROUPS] + [f"recipient ({recipient})"]
    fig.legend(handles, labels, loc="lower center", ncol=5, frameon=False, fontsize=8)


def main():
    os.makedirs(OUT, exist_ok=True)
    org = load_taxa()
    design = list(csv.DictReader(open(f"{D}/design.tsv"), delimiter="\t"))

    # full-size, fully labelled tanglegram for every planted transfer
    for r in design:
        if r["condition"] != "spike":
            continue
        fig, ax = plt.subplots(figsize=(14, 20))
        xr_tip, tree_w, n = draw(ax, r, org)
        ax.text(tree_w / 2, 2.6, "GTDB species tree", ha="center", fontsize=12, weight="bold", color=TEXT)
        ax.text(xr_tip + tree_w / 2, 2.6, f"{r['gene']} gene tree", ha="center", fontsize=12,
                weight="bold", color=TEXT)
        fig.suptitle(title_of(r), fontsize=12, y=0.995)
        legend(fig, r["recipient"])
        fig.savefig(f"{OUT}/{r['test_id']}.png", dpi=180, bbox_inches="tight", facecolor="white")
        plt.close(fig)

    # overview per gene x recipient: control (top) vs spike (bottom) for each donor level
    by = {(r["gene"], r["recipient_code"], r["donor_level"], r["condition"]): r for r in design}
    levels = list(dict.fromkeys(r["donor_level"] for r in design))
    for gene, code in dict.fromkeys((r["gene"], r["recipient_code"]) for r in design):
        fig, axes = plt.subplots(2, len(levels), figsize=(6 * len(levels), 18))
        for j, lv in enumerate(levels):
            for i, cond in enumerate(("control", "spike")):
                r = by[(gene, code, lv, cond)]
                draw(axes[i, j], r, org, compact=True)
                level = f"{lv.replace('_', ' ')}, GTDB distance {r['gtdb_distance']}"
                head = (f"PLANTED: {r['donor']} → {short(r['recipient'])}\n({level})" if cond == "spike"
                        else f"CONTROL ({level})\n{r['donor']} removed, no transfer")
                axes[i, j].set_title(head, fontsize=10, color=ACCENT if cond == "spike" else TEXT)
        rec = by[(gene, code, levels[0], "spike")]["recipient"]
        fig.suptitle(f"{gene}: GTDB (left of each panel) vs gene tree (right); red line = {rec}", fontsize=13)
        legend(fig, rec)
        fig.tight_layout(rect=(0, 0.03, 1, 0.97))
        fig.savefig(f"{OUT}/overview_{gene}_{code}.png", dpi=150, facecolor="white")
        plt.close(fig)
    print(f"Wrote planted-transfer tanglegrams + overviews to {OUT}/")


if __name__ == "__main__":
    main()
