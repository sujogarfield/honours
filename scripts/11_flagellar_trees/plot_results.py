#!/usr/bin/env python3
"""
plot_results.py
Figures for the Katana flagellar tree run. Run after the outputs are copied
back from Katana (needs gene_trees/, concat/, concordance/, comparison.json,
hgt_screen.json in gene_order/flagellar_trees/).

Run from the project root:  python scripts/11_flagellar_trees/plot_results.py
Output (gene_order/flagellar_trees/figures/):
  fig1_tanglegram.png          GTDB (branches shaded by gCF) vs concatenated flagellar tree
  fig2_signal_vs_length.png    per-gene discordance with GTDB vs informative sites
  fig3_concordance_gtdb.png    gCF vs sCF on every GTDB branch
  fig4_recurrent_conflicts.png groups placed differently from GTDB in many gene trees
"""

import json
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ete3 import Tree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from screen_hgt_candidates import bitmask_splits, smaller_side

D = "gene_order/flagellar_trees"
OUT = f"{D}/figures"

TEXT, TEXT2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
ACCENT = "#e34948"
# clade groups -> categorical slots in fixed order (reference palette, light mode)
GROUPS = [
    ("Campylobacter", "#2a78d6", {"Campylobacter"}),
    ("Helicobacter + Wolinella", "#eb6834", {"Helicobacter", "Wolinella"}),
    ("Sulfurospirillum", "#1baf7a", {"Sulfurospirillum"}),
    ("Arcobacteraceae", "#eda100", {"Arcobacter", "Aliarcobacter", "Halarcobacter",
                                     "Malaciobacter", "Poseidonibacter"}),
    ("Sulfurimonas + Sulfuricurvum", "#e87ba4", {"Sulfurimonas", "Sulfuricurvum"}),
    ("Nautiliales", "#008300", {"Nautilia", "Caminibacter", "Lebetimonas"}),
    ("Nitratiruptoraceae", "#4a3aa7", {"Nitratiruptor", "Nitrosophilus", "Hydrogenimonas"}),
    ("Desulfurellia (GTDB root)", "#52514e", {"Hippea", "Desulfurella"}),
]

plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.edgecolor": TEXT2,
                     "axes.labelcolor": TEXT, "xtick.color": TEXT2, "ytick.color": TEXT2,
                     "axes.spines.top": False, "axes.spines.right": False})


def load_taxa():
    return {l.split("\t")[0]: l.split("\t")[1] for l in open(f"{D}/taxa.tsv").read().splitlines()[1:]}


def group_color(organism):
    genus = organism.split()[0]
    for _, color, genera in GROUPS:
        if genus in genera:
            return color
    return TEXT2


def cladogram_layout(tree, y_of_leaf):
    """x: leaves aligned at 1, internal nodes at 1 - height/max_height; y: mean of children."""
    height = {}
    for n in tree.traverse("postorder"):
        height[n] = 0 if n.is_leaf() else 1 + max(height[c] for c in n.children)
    top = height[tree]
    pos = {}
    for n in tree.traverse("postorder"):
        y = y_of_leaf[n.name] if n.is_leaf() else np.mean([pos[c][1] for c in n.children])
        pos[n] = (1 - height[n] / top, y)
    return pos


def draw_tree(ax, tree, pos, x0, width, mirror, branch_style):
    for n in tree.traverse():
        if n.is_root():
            continue
        (xp, yp), (x, y) = pos[n.up], pos[n]
        X = lambda v: x0 + (1 - v) * width if mirror else x0 + v * width
        color, lw, ls = branch_style(n)
        ax.plot([X(xp), X(x)], [y, y], color=color, lw=lw, ls=ls, solid_capstyle="butt")
        ax.plot([X(xp), X(xp)], [yp, y], color=color, lw=lw, ls=ls, solid_capstyle="butt")


def fig1_tanglegram(org):
    gtdb = Tree(f"{D}/gtdb_ref_pruned.nwk", format=1)
    flag = Tree(f"{D}/concat/concat.treefile", format=0)
    taxa = sorted(gtdb.get_leaf_names())
    index = {t: i for i, t in enumerate(taxa)}

    # gCF per GTDB branch, matched through bipartitions (cf.branch is unrooted)
    rows = [l.rstrip("\n").split("\t") for l in open(f"{D}/concordance/gtdb_gcf.cf.stat")
            if not l.startswith("#") and l.strip()]
    col = {h: i for i, h in enumerate(rows[0])}
    gcf_by_id = {r[col["ID"]]: float(r[col["gCF"]]) for r in rows[1:] if r[col["gCF"]] != "NA"}
    cfb = Tree(f"{D}/concordance/gtdb_gcf.cf.branch", format=1)
    gcf_of_split = {m: gcf_by_id.get(bid) for m, bid in bitmask_splits(cfb, index).items()}
    full = (1 << len(taxa)) - 1
    canon = lambda m: m ^ full if m & 1 else m

    gtdb.ladderize()
    y_gtdb = {name: -i for i, name in enumerate(gtdb.get_leaf_names())}

    # root the flagellar tree like GTDB where possible, then rotate to minimise crossings
    outgroup = [t for t in taxa if org[t].split()[0] in {"Hippea", "Desulfurella"}]
    anc = flag.get_common_ancestor(outgroup)
    flag.set_outgroup(anc if anc is not flag else flag.get_midpoint_outgroup())
    for n in flag.traverse("postorder"):
        if not n.is_leaf():
            n.children.sort(key=lambda c: -np.mean([y_gtdb[l] for l in c.get_leaf_names()]))
    y_flag = {name: -i for i, name in enumerate(flag.get_leaf_names())}

    bitmask_splits(gtdb, index)  # annotates node.mask
    pos_g = cladogram_layout(gtdb, y_gtdb)
    pos_f = cladogram_layout(flag, y_flag)

    fig, ax = plt.subplots(figsize=(15, 21))
    tree_w, lab_w, gap = 3.0, 2.6, 2.2
    xl_tip = tree_w
    xr_tip = tree_w + 2 * lab_w + gap

    def gtdb_style(n):
        if n.is_leaf():
            return TEXT2, 0.7, "-"
        g = gcf_of_split.get(canon(n.mask))
        if g is None:
            return GRID, 0.9, "-"
        step = SEQ[min(int(g / 100 * len(SEQ)), len(SEQ) - 1)]
        return (ACCENT, 2.2, "-") if g <= 5 else (step, 1.6, "-")

    def flag_style(n):
        if n.is_leaf():
            return TEXT2, 0.7, "-"
        return (TEXT, 1.0, "-") if n.support >= 95 else ("#9a9893", 0.9, (0, (2, 1.5)))

    draw_tree(ax, gtdb, pos_g, 0, tree_w, False, gtdb_style)
    draw_tree(ax, flag, pos_f, xr_tip, tree_w, True, flag_style)

    for t in taxa:
        c = group_color(org[t])
        yl, yr = y_gtdb[t], y_flag[t]
        ax.text(xl_tip + 0.06, yl, org[t], va="center", ha="left", fontsize=5.6, color=TEXT)
        ax.text(xr_tip - 0.06, yr, org[t], va="center", ha="right", fontsize=5.6, color=TEXT)
        moved = abs(yl - yr) > 3
        ax.plot([xl_tip + lab_w, xr_tip - lab_w], [yl, yr], color=c,
                lw=1.1 if moved else 0.6, alpha=0.95 if moved else 0.45)

    n = len(taxa)
    ax.text(tree_w / 2, 2.2, "GTDB species tree", ha="center", fontsize=12, weight="bold", color=TEXT)
    ax.text(tree_w / 2, 0.9, "branch colour = gCF (share of the 34 flagellar gene trees containing it)",
            ha="center", fontsize=7.5, color=TEXT2)
    ax.text(xr_tip + tree_w / 2, 2.2, "Concatenated flagellar tree (34 genes)", ha="center",
            fontsize=12, weight="bold", color=TEXT)
    ax.text(xr_tip + tree_w / 2, 0.9, "solid = UFBoot ≥ 95, dashed = UFBoot < 95; rooted on Desulfurellia",
            ha="center", fontsize=7.5, color=TEXT2)

    # legends: gCF ramp + clade groups
    lx, ly = 0.05, -n - 3
    ax.text(lx, ly, "gCF", fontsize=8, color=TEXT, weight="bold", va="center")
    for i, c in enumerate(SEQ):
        ax.plot([lx + 0.35 + i * 0.28, lx + 0.6 + i * 0.28], [ly, ly], color=c, lw=5, solid_capstyle="butt")
    ax.text(lx + 0.35, ly - 1.6, "0%", fontsize=7, color=TEXT2)
    ax.text(lx + 0.35 + 7 * 0.28, ly - 1.6, "100%", fontsize=7, color=TEXT2, ha="right")
    ax.plot([lx + 2.5, lx + 2.8], [ly, ly], color=ACCENT, lw=5, solid_capstyle="butt")
    ax.text(lx + 2.9, ly, "gCF ≤ 5%", fontsize=7.5, color=TEXT, va="center")
    gx = xl_tip + 1.4
    for i, (name, c, _) in enumerate(GROUPS):
        x = gx + (i % 4) * 2.3
        y = ly - (i // 4) * 2.0
        ax.plot([x, x + 0.3], [y, y], color=c, lw=3, solid_capstyle="butt")
        ax.text(x + 0.4, y, name, fontsize=7.5, color=TEXT, va="center")
    ax.text(gx, ly - 4.3, "Lines join each genome in the two trees; bold lines = genome sits in a different "
            "place (>3 rows apart).", fontsize=7.5, color=TEXT2)

    ax.set_xlim(-0.1, xr_tip + tree_w + 0.1)
    ax.set_ylim(-n - 8.5, 3)
    ax.axis("off")
    fig.savefig(f"{OUT}/fig1_tanglegram.png", dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def per_gene_table():
    comp = json.load(open(f"{D}/comparison.json"))["per_gene"]
    rows = []
    for gene, v in comp.items():
        log = open(f"{D}/gene_trees/{gene}.iqtree").read()
        sites = int(re.search(r"Input data: \d+ sequences with (\d+)", log).group(1))
        pi = int(re.search(r"Number of parsimony informative sites: (\d+)", log).group(1))
        t = Tree(f"{D}/gene_trees/{gene}.treefile", format=0)
        sup = [n.support for n in t.traverse() if not n.is_leaf() and not n.is_root()]
        rows.append({"gene": gene, "sites": sites, "informative": pi,
                     "rf": v["rf_norm_vs_gtdb"], "strong": float(np.mean(np.array(sup) >= 95))})
    return rows


def spearman(a, b):
    ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


def fig2_signal_vs_length(rows):
    x = np.array([r["informative"] for r in rows])
    y = np.array([r["rf"] for r in rows])
    rho = spearman(x, y)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.scatter(x, y, s=42, color="#2a78d6", edgecolor="white", linewidth=1.5, zorder=3)
    for r in rows:
        if r["informative"] < 110 or r["rf"] > 0.5 or r["informative"] > 450:
            ax.annotate(r["gene"], (r["informative"], r["rf"]), xytext=(5, 3), textcoords="offset points",
                        fontsize=8, color=TEXT)
    ax.set_xscale("log")
    ax.set_xlabel("Parsimony-informative sites in trimmed alignment (log scale)")
    ax.set_ylabel("Normalised RF distance to GTDB\n(0 = identical, 1 = no shared branches)")
    ax.set_title("Short genes disagree most with the species tree", loc="left", fontsize=11,
                 color=TEXT, weight="bold")
    ax.text(0.99, 0.97, f"Spearman ρ = {rho:.2f}, n = {len(rows)} genes", transform=ax.transAxes,
            ha="right", va="top", fontsize=8.5, color=TEXT2)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig2_signal_vs_length.png", dpi=220, facecolor="white")
    plt.close(fig)
    return rho


def fig3_concordance(org):
    def read(p):
        rows = [l.rstrip("\n").split("\t") for l in open(p) if not l.startswith("#") and l.strip()]
        return {h: i for i, h in enumerate(rows[0])}, rows[1:]
    hg, rg = read(f"{D}/concordance/gtdb_gcf.cf.stat")
    hs, rs = read(f"{D}/concordance/gtdb_scfl.cf.stat")
    scf = {r[hs["ID"]]: float(r[hs["sCF"]]) for r in rs if r[hs["sCF"]] != "NA"}
    pts = [(r[hg["ID"]], float(r[hg["gCF"]]), scf[r[hg["ID"]]]) for r in rg
           if r[hg["gCF"]] != "NA" and r[hg["ID"]] in scf]

    cfb = Tree(f"{D}/concordance/gtdb_gcf.cf.branch", format=1)
    ntaxa = len(cfb)
    label = {}
    for n in cfb.traverse():
        if n.name and not n.is_leaf():
            leaves = n.get_leaf_names()
            side = leaves if len(leaves) <= ntaxa / 2 else sorted(set(cfb.get_leaf_names()) - set(leaves))
            genera = sorted({org[t].split()[0] for t in side})
            label[n.name] = (", ".join(genera[:2]) + (" …" if len(genera) > 2 else "")) + f" ({len(side)})"

    low = sorted([p for p in pts if p[1] <= 5], key=lambda p: -p[2])
    rest = [p for p in pts if p[1] > 5]
    fig, (ax, key) = plt.subplots(1, 2, figsize=(10.5, 5.8), gridspec_kw={"width_ratios": [2.3, 1]})
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.axhline(100 / 3, color=TEXT2, lw=0.8, ls=(0, (3, 2)), label="sCF expected with no signal (⅓)")
    ax.scatter([p[1] for p in rest], [p[2] for p in rest], s=30, color="#2a78d6", edgecolor="white",
               linewidth=1.2, label="GTDB branch", zorder=3)
    ax.scatter([p[1] for p in low], [p[2] for p in low], s=42, color="#eb6834", edgecolor="white",
               linewidth=1.2, label="gCF ≤ 5%: no flagellar gene tree recovers it", zorder=4)
    # numbered labels in a column left of x=0, joined to their points by leader lines
    ys = np.linspace(82, 26, len(low))
    for i, ((bid, g, s), y) in enumerate(zip(low, ys), 1):
        ax.plot([-9.5, g], [y, s], color=TEXT2, lw=0.5, zorder=2)
        ax.text(-10, y, str(i), fontsize=7.5, color=TEXT, ha="right", va="center")
    ax.set_xlim(-15, 102)
    ax.set_ylim(15, 100)
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_xlabel("gCF: % of flagellar gene trees containing the branch")
    ax.set_ylabel("sCF: % of informative sites supporting the branch")
    ax.set_title("Concordance of the 34 flagellar genes with each GTDB branch", loc="left",
                 fontsize=11, color=TEXT, weight="bold")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False, fontsize=8)

    key.axis("off")
    key.text(0, 1.0, "Branches with gCF ≤ 5%", fontsize=9, weight="bold", color=TEXT, va="top",
             transform=key.transAxes)
    key.text(0, 0.95, "clade (genomes)   gCF / sCF", fontsize=7.5, color=TEXT2, va="top",
             transform=key.transAxes)
    for i, (bid, g, s) in enumerate(low, 1):
        key.text(0, 0.89 - (i - 1) * 0.075, f"{i:>2}  {label.get(bid, bid)}\n     {g:.0f}% / {s:.0f}%",
                 fontsize=7.5, color=TEXT, va="top", transform=key.transAxes)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig3_concordance_gtdb.png", dpi=220, facecolor="white")
    plt.close(fig)


def fig4_recurrent(n_show=12):
    cands = json.load(open(f"{D}/hgt_screen.json"))["candidates"]
    seen, groups = set(), []
    for c in cands:
        key = tuple(c["group"])
        if key in seen:
            continue
        seen.add(key)
        genes = [t for t in [c["tree"]] + c["also_in"] if t != "CONCAT"]
        in_concat = "CONCAT" in [c["tree"]] + c["also_in"]
        groups.append((len(genes), in_concat, c["gtdb_span"], c["organisms"], c["gtdb_extra"]))
    groups.sort(key=lambda g: (-g[0], -g[2]))
    top = groups[:n_show][::-1]

    def name(g):
        orgs = g[3]
        short = [f"{o.split()[0][0]}. {o.split()[1]}" for o in orgs[:3]]
        s = " + ".join(short) + (f" + {len(orgs) - 3} more" if len(orgs) > 3 else "")
        return s

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ys = np.arange(len(top))
    ax.barh(ys, [g[0] for g in top], height=0.62, color="#2a78d6", edgecolor="white", linewidth=2)
    for y, g in zip(ys, top):
        ax.text(g[0] + 0.3, y, f"{g[0]}" + ("  + concat tree" if g[1] else "") + f"   (GTDB span {g[2]})",
                va="center", fontsize=7.5, color=TEXT2)
    ax.set_yticks(ys)
    ax.set_yticklabels([name(g) for g in top], fontsize=8, color=TEXT)
    ax.set_xlabel("Gene trees grouping these genomes together (UFBoot ≥ 95), against GTDB")
    ax.set_xlim(0, 34)
    ax.xaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_title("Groupings repeated across many flagellar genes", loc="left", fontsize=11,
                 color=TEXT, weight="bold", pad=22)
    ax.text(0, 1.02, "GTDB span = size of the smallest GTDB clade containing the group",
            transform=ax.transAxes, fontsize=8, color=TEXT2)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig4_recurrent_conflicts.png", dpi=220, facecolor="white")
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    org = load_taxa()
    fig1_tanglegram(org)
    rows = per_gene_table()
    rho = fig2_signal_vs_length(rows)
    fig3_concordance(org)
    fig4_recurrent()
    print(f"Spearman rho (informative sites vs RF to GTDB): {rho:.3f}")
    print(f"Wrote {OUT}/")


if __name__ == "__main__":
    main()
