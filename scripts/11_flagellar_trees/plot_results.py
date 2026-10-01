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
  fig5_poseidonibacter.png     zoom on the Poseidonibacter/Arcobacter candidate
  fig6_branch_lengths.png      flagellar vs genome distance for every genome pair
"""

import glob
import json
from collections import Counter
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
N_GENES = len(open(f"{D}/genes.txt").read().split())

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


def phylogram_layout(tree, y_of_leaf):
    """x: root-to-node path length / maximum root-to-tip length (true branch
    lengths); y: leaves from y_of_leaf, internal nodes at the mean of children."""
    depth = {}
    for n in tree.traverse("preorder"):
        depth[n] = 0.0 if n.is_root() else depth[n.up] + n.dist
    top = max(depth[l] for l in tree.iter_leaves())
    pos = {}
    for n in tree.traverse("postorder"):
        y = y_of_leaf[n.name] if n.is_leaf() else np.mean([pos[c][1] for c in n.children])
        pos[n] = (depth[n] / top, y)
    return pos, top


def scale_bar(ax, x0, y, width, top, mirror, label_len=None):
    """Scale bar for a tree drawn `width` wide whose deepest tip is `top`
    substitutions from the root; length is a round number near top/8."""
    if label_len is None:
        label_len = min((0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0), key=lambda v: abs(v - top / 8))
    L = width * label_len / top
    xs = [x0, x0 + L] if not mirror else [x0 - L, x0]
    ax.plot(xs, [y, y], color=TEXT, lw=1.2)
    ax.text(np.mean(xs), y - 1.3, f"{label_len} subst./site", ha="center", fontsize=7, color=TEXT2)


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

    # root the flagellar tree like GTDB (Desulfurellia) where possible
    outgroup = [t for t in taxa if org[t].split()[0] in {"Hippea", "Desulfurella"}]
    anc = flag.get_common_ancestor(outgroup)
    flag.set_outgroup(anc if anc is not flag else flag.get_midpoint_outgroup())

    # Rotate child order (never topology) to minimise crossing lines: alternately
    # order each tree's children by the mean row of their tips in the other tree,
    # and keep the arrangement with the fewest crossings.
    def rows(tree):
        return {name: -i for i, name in enumerate(tree.get_leaf_names())}

    def rotate_to(tree, target):
        for nd in tree.traverse("postorder"):
            if not nd.is_leaf():
                nd.children.sort(key=lambda c: -np.mean([target[l] for l in c.get_leaf_names()]))

    def crossings(a, b):
        order = sorted(taxa, key=lambda t: -a[t])
        seq = [-b[t] for t in order]
        return sum(1 for i in range(len(seq)) for j in range(i + 1, len(seq)) if seq[i] > seq[j])

    gtdb.ladderize()
    best = None
    for _ in range(8):
        rotate_to(flag, rows(gtdb))
        c = crossings(rows(gtdb), rows(flag))
        if best is None or c < best[0]:
            best = (c, gtdb.write(format=1), flag.write(format=0))
        rotate_to(gtdb, rows(flag))
        c = crossings(rows(gtdb), rows(flag))
        if c < best[0]:
            best = (c, gtdb.write(format=1), flag.write(format=0))
    gtdb, flag = Tree(best[1], format=1), Tree(best[2], format=0)
    print(f"tanglegram: {best[0]} crossing line pairs after rotation")

    # Faithfulness check: the drawn trees must have exactly the source trees' branches
    src_g = Tree(f"{D}/gtdb_ref_pruned.nwk", format=1)
    src_f = Tree(f"{D}/concat/concat.treefile", format=0)
    for drawn, src, name in ((gtdb, src_g, "GTDB"), (flag, src_f, "flagellar")):
        if set(bitmask_splits(drawn, index)) != set(bitmask_splits(src, index)):
            raise SystemExit(f"{name} tree drawn with a different topology from its source file")
    y_gtdb, y_flag = rows(gtdb), rows(flag)

    bitmask_splits(gtdb, index)  # annotates node.mask
    pos_g, top_g = phylogram_layout(gtdb, y_gtdb)
    pos_f, top_f = phylogram_layout(flag, y_flag)

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
        return (TEXT2, 0.7, "-") if n.is_leaf() else (TEXT, 1.0, "-")

    draw_tree(ax, gtdb, pos_g, 0, tree_w, False, gtdb_style)
    draw_tree(ax, flag, pos_f, xr_tip, tree_w, True, flag_style)
    for nd in flag.traverse():
        if not nd.is_leaf() and not nd.is_root() and nd.support < 95:
            x, y = pos_f[nd]
            ax.plot(xr_tip + (1 - x) * tree_w, y, "o", ms=3.2, mfc="white", mec=TEXT, mew=0.8, zorder=5)
    # dotted leaders from each tip to the aligned labels
    for leaf in gtdb.iter_leaves():
        x, y = pos_g[leaf]
        if x < 0.995:
            ax.plot([x * tree_w, xl_tip], [y, y], color=GRID, lw=0.5, ls=(0, (1, 1.5)))
    for leaf in flag.iter_leaves():
        x, y = pos_f[leaf]
        if x < 0.995:
            ax.plot([xr_tip, xr_tip + (1 - x) * tree_w], [y, y], color=GRID, lw=0.5, ls=(0, (1, 1.5)))

    for t in taxa:
        yl, yr = y_gtdb[t], y_flag[t]
        ax.text(xl_tip + 0.06, yl, org[t], va="center", ha="left", fontsize=5.6, color=TEXT)
        ax.text(xr_tip - 0.06, yr, org[t], va="center", ha="right", fontsize=5.6, color=TEXT)
        ax.plot([xl_tip + lab_w, xr_tip - lab_w], [yl, yr], color=group_color(org[t]), lw=0.8, alpha=0.75)

    n = len(taxa)
    ax.text(tree_w / 2, 2.2, "GTDB species tree", ha="center", fontsize=12, weight="bold", color=TEXT)
    ax.text(tree_w / 2, 0.9, f"branch colour = gCF (share of the {N_GENES} flagellar gene trees containing it)",
            ha="center", fontsize=7.5, color=TEXT2)
    ax.text(xr_tip + tree_w / 2, 2.2, f"Concatenated flagellar tree ({N_GENES} genes)", ha="center",
            fontsize=12, weight="bold", color=TEXT)
    ax.text(xr_tip + tree_w / 2, 0.9, "branch lengths to scale; ○ = node with UFBoot < 95; rooted on Desulfurellia",
            ha="center", fontsize=7.5, color=TEXT2)

    scale_bar(ax, 0.05, -n - 0.5, tree_w, top_g, False)
    scale_bar(ax, xr_tip + tree_w - 0.05, -n - 0.5, tree_w, top_f, True)

    # legends: gCF ramp + clade groups
    lx, ly = 0.05, -n - 4.5
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
    ax.text(gx, ly - 4.3, "Lines join each genome in the two trees. Branch order is rotated only to reduce "
            "crossings (topology unchanged); remaining crossings are real differences.",
            fontsize=7.5, color=TEXT2)

    ax.set_xlim(-0.1, xr_tip + tree_w + 0.1)
    ax.set_ylim(-n - 10, 3)
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
    ax.set_title(f"Concordance of the {N_GENES} flagellar genes with each GTDB branch", loc="left",
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
    ax.set_xlim(0, N_GENES)
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


def _clade_support(t, members):
    """UFBoot of the branch separating `members` from the rest (unrooted), else None."""
    ms, allx = set(members), set(t.get_leaf_names())
    for nd in t.traverse():
        if nd.is_leaf() or nd.is_root():
            continue
        side = set(nd.get_leaf_names())
        if side == ms or allx - side == ms:
            return nd.support
    return None


def _pair_identity(a, b):
    m = t = 0
    for x, y in zip(a, b):
        if x != "-" and y != "-":
            t += 1
            m += x == y
    return m, t


def _read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def fig5_case(org):
    """Zoom on the Poseidonibacter/Arcobacter candidate: marker vs flagellar
    tree (true branch lengths), per-gene support along P. lekithochrous's
    flagellar loci, and the raw-identity check."""
    from Bio import Align
    from Bio.Align import substitution_matrices
    inv = {v: k for k, v in org.items()}
    L, R, P, A = (inv[x] for x in ("Poseidonibacter lekithochrous", "Arcobacter roscoffensis",
                                   "Poseidonibacter parvus", "Poseidonibacter antarcticus"))
    outg = [inv[x] for x in ("Aliarcobacter butzleri", "Aliarcobacter skirrowii", "Aliarcobacter cryaerophilus")]
    keep = [L, R, P, A] + outg
    short = lambda a: org[a].replace("Poseidonibacter", "P.").replace("Aliarcobacter", "Al.").replace("Arcobacter", "A.")

    fig = plt.figure(figsize=(11, 7.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.3, 1], hspace=0.22, wspace=0.35)
    for k, (path, title) in enumerate([(f"{D}/markers/free/markers.treefile", "Genome: GTDB bac120 markers"),
                                       (f"{D}/concat/concat.treefile", f"Flagellar genes ({N_GENES}, concatenated)")]):
        t = Tree(path, format=0)
        t.prune(keep, preserve_branch_length=True)
        og = t.get_common_ancestor(outg)
        t.set_outgroup(og if og is not t else outg[0])
        t.ladderize()
        y = {n: -i for i, n in enumerate(t.get_leaf_names())}
        pos, top = phylogram_layout(t, y)
        ax = fig.add_subplot(gs[0, k])
        for nd in t.traverse():
            if nd.is_root():
                continue
            (xp, yp), (x, yy) = pos[nd.up], pos[nd]
            ax.plot([xp, x], [yy, yy], color=TEXT, lw=1.3)
            ax.plot([xp, xp], [yp, yy], color=TEXT, lw=1.3)
            if not nd.is_leaf():
                ax.text(x - 0.01, yy + 0.12, f"{nd.support:.0f}", fontsize=7, color=TEXT2, ha="right")
        for leaf in t.iter_leaves():
            x, yy = pos[leaf]
            hl = leaf.name in (L, R)
            ax.text(x + 0.02, yy, short(leaf.name), va="center", fontsize=8.5,
                    color=ACCENT if hl else TEXT, weight="bold" if hl else "normal")
        bar = 0.02 if "markers" in path else 0.1
        ax.plot([0, bar / top], [-len(keep) + 0.2] * 2, color=TEXT, lw=1.2)
        ax.text(bar / top / 2, -len(keep) - 0.35, f"{bar} subst./site", ha="center", fontsize=7, color=TEXT2)
        ax.set_xlim(-0.02, 1.75)
        ax.set_ylim(-len(keep) - 0.8, 0.8)
        ax.axis("off")
        ax.set_title(title, loc="left", fontsize=10.5, weight="bold", color=TEXT)

    # per-gene support strip in P. lekithochrous genome order
    comb = json.load(open("gene_order/gene_order_combined.json"))
    gk = next(k for k in comb if k.startswith(L))
    pos_l = {}
    for x in comb[gk]["genes"]:
        pos_l.setdefault(x["gene"], (x["contig"], x["start"]))
    genes = [g for g in open(f"{D}/genes.txt").read().split() if g in pos_l]
    genes.sort(key=lambda g: pos_l[g])
    ax = fig.add_subplot(gs[1, :])
    cols = {"flag": ACCENT, "genome": "#2a78d6", "neither": "#c9c8c3", "missing": "white"}
    counts = Counter()
    for i, g in enumerate(genes):
        t = Tree(f"{D}/gene_trees/{g}.treefile", format=0)
        if not {L, R, P, A} <= set(t.get_leaf_names()):
            kind, sup = "missing", None
        elif (sup := _clade_support(t, [L, R])) is not None:
            kind = "flag"
        elif (sup := _clade_support(t, [L, P, A])) is not None:
            kind = "genome"
        else:
            kind, sup = "neither", None
        counts[kind] += 1
        ax.add_patch(plt.Rectangle((i, 0), 0.9, 1, facecolor=cols[kind], edgecolor=TEXT2, lw=0.6))
        ax.text(i + 0.45, -0.15, g[:3] + g[3:].upper(), rotation=90, ha="center", va="top", fontsize=7.5, color=TEXT)
        if sup is not None:
            ax.text(i + 0.45, 1.12, f"{sup:.0f}", ha="center", fontsize=6.5, color=TEXT2)
    starts = [pos_l[g][1] for g in genes]
    for lo, hi, lab in [(2.9e6, 3.0e6, "locus 1 (≈2.94–2.97 Mb)"), (3.2e6, 3.3e6, "locus 2 (≈3.24–3.25 Mb)")]:
        idx = [i for i, s0 in enumerate(starts) if lo <= s0 <= hi]
        if idx:
            ax.plot([idx[0], idx[-1] + 0.9], [1.55, 1.55], color=TEXT2, lw=1)
            ax.text((idx[0] + idx[-1] + 0.9) / 2, 1.65, lab, ha="center", fontsize=8, color=TEXT2)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=cols[k], edgecolor=TEXT2, lw=0.6) for k in ("flag", "genome", "neither", "missing")]
    ax.legend(handles, [f"groups P. lekithochrous with A. roscoffensis ({counts['flag']})",
                        f"groups P. lekithochrous with Poseidonibacter ({counts['genome']})",
                        f"neither / unresolved ({counts['neither']})",
                        f"gene missing in one of the four ({counts['missing']})"],
              loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=2, frameon=False, fontsize=8)
    ax.set_xlim(-0.3, len(genes) + 0.2)
    ax.set_ylim(-1.0, 1.9)
    ax.axis("off")
    ax.set_title("Each flagellar gene tree, in P. lekithochrous genome order (numbers = UFBoot)",
                 loc="left", fontsize=10.5, weight="bold", color=TEXT)

    # raw identity check (tree-free)
    def pooled(alns, x, y):
        M = T = 0
        for al in alns:
            if x in al and y in al:
                m, t2 = _pair_identity(al[x], al[y])
                M += m
                T += t2
        return 100 * M / T
    flag_alns = [_read_fasta(p2) for p2 in glob.glob(f"{D}/trimmed/*.faa")]
    mark = [_read_fasta(f"{D}/markers/bac120_msa_r232.faa")]
    al = Align.PairwiseAligner(mode="global", open_gap_score=-11, extend_gap_score=-1, end_gap_score=0)
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    core = {}
    for x in (R, P):
        M = T = 0
        for p2 in glob.glob(f"{D}/core/genes/*.faa"):
            sq = _read_fasta(p2)
            a2 = al.align(sq[L], sq[x])[0]
            for (s1, e1), (s2, e2) in zip(*a2.aligned):
                for i, j in zip(range(s1, e1), range(s2, e2)):
                    T += 1
                    M += sq[L][i] == sq[x][j]
        core[x] = 100 * M / T
    txt = ("% protein identity of P. lekithochrous to A. roscoffensis vs to P. parvus:   "
           f"flagellar genes {pooled(flag_alns, L, R):.1f} vs {pooled(flag_alns, L, P):.1f}   ·   "
           f"core genes {core[R]:.1f} vs {core[P]:.1f}   ·   "
           f"GTDB markers {pooled(mark, L, R):.1f} vs {pooled(mark, L, P):.1f}\n"
           "The genome leans the same way more weakly, so this is a candidate pending the core-gene backbone test.")
    fig.text(0.06, -0.04, txt, fontsize=8, color=TEXT2)
    fig.suptitle("Candidate non-vertical inheritance: Poseidonibacter lekithochrous / Arcobacter roscoffensis",
                 x=0.06, ha="left", fontsize=12, weight="bold", color=TEXT)
    fig.savefig(f"{OUT}/fig5_poseidonibacter.png", dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _patristic(path):
    t = Tree(path, format=0)
    depth = {}
    for nd in t.traverse("preorder"):
        depth[nd] = 0.0 if nd.is_root() else depth[nd.up] + nd.dist
    leaves = {l.name: l for l in t.iter_leaves()}
    names = sorted(leaves)
    anc = {x: set(leaves[x].get_ancestors()) for x in names}
    out = {}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            nd = leaves[b]
            while nd not in anc[a]:
                nd = nd.up
            out[(a, b)] = depth[leaves[a]] + depth[leaves[b]] - 2 * depth[nd]
    return out


def fig6_branch_lengths(org):
    """Flagellar vs genome (bac120 marker ML tree) patristic distance for every
    genome pair, log-log, with the power-law fit; pairs involving
    Nitrosophilus/Nitratiruptor highlighted."""
    F = _patristic(f"{D}/concat/concat.treefile")
    G = _patristic(f"{D}/markers/free/markers.treefile")
    pairs = sorted(F)
    f = np.array([F[p] for p in pairs])
    g = np.array([G[p] for p in pairs])
    b, a = np.polyfit(np.log(g), np.log(f), 1)
    therm = {t for t in org if org[t].split()[0] in ("Nitrosophilus", "Nitratiruptor")}
    hit = np.array([bool(set(p) & therm) for p in pairs])
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.scatter(g[~hit], f[~hit], s=6, color="#2a78d6", alpha=0.35, lw=0, label="genome pair")
    ax.scatter(g[hit], f[hit], s=10, color="#eb6834", alpha=0.8, lw=0,
               label="pair involving Nitrosophilus / Nitratiruptor")
    xs = np.linspace(g.min(), g.max(), 100)
    ax.plot(xs, np.exp(a) * xs ** b, color=TEXT, lw=1.2, label=f"fit: flagellar = {np.exp(a):.1f} × genome$^{{{b:.2f}}}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Genome distance (bac120 marker tree, subst./site)")
    ax.set_ylabel("Flagellar distance (concatenated tree, subst./site)")
    ax.set_title("Flagellar vs genome divergence for every pair of genomes", loc="left",
                 fontsize=11, color=TEXT, weight="bold")
    r = np.corrcoef(np.log(f), np.log(g))[0, 1]
    ax.text(0.99, 0.03, f"{len(pairs):,} pairs, log-log r = {r:.2f}", transform=ax.transAxes,
            ha="right", fontsize=8.5, color=TEXT2)
    ax.legend(loc="upper left", frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig6_branch_lengths.png", dpi=220, facecolor="white")
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    org = load_taxa()
    fig1_tanglegram(org)
    rows = per_gene_table()
    rho = fig2_signal_vs_length(rows)
    fig3_concordance(org)
    fig4_recurrent()
    fig5_case(org)
    fig6_branch_lengths(org)
    print(f"Spearman rho (informative sites vs RF to GTDB): {rho:.3f}")
    print(f"Wrote {OUT}/")


if __name__ == "__main__":
    main()
