#!/usr/bin/env python3
"""
collinear_blocks.py
Which flagellar gene clusters are collinear (same order, same orientation,
physically contiguous) across the reliable closed-chromosome genomes, and how
they are arranged in each lineage.

1. Contiguity: two consecutive flagellar genes are physically adjacent when
   the intergenic gap is <= GAP_BP (300 bp -- too short to hide another gene;
   35 % of flagellar neighbour gaps are <= 50 bp, then the distribution
   flattens).
2. Signed adjacencies (the gene ends facing each other, as in
   synteny_distance.py) are recorded between physically adjacent genes.
3. Lineages are the GTDB groups of plot_results.GROUPS; an adjacency is
   "held" by a lineage when >= LINEAGE_FRAC of its genomes have it (only
   lineages with >= MIN_LINEAGE genomes vote).
4. Core collinear blocks: adjacencies held by >= 2 lineages, chained greedily
   from most to least frequent overall (each gene end used once, no cycles),
   so each block is one linear chain, read with its first gene on '+'.
   Lineage-specific adjacencies (held by exactly one lineage) are reported
   separately: they show how that lineage has re-joined or broken the blocks.
5. Lineage consensus loci: the same greedy chaining per lineage, over every
   adjacency that lineage holds -- the collinear clusters a typical member of
   the lineage carries.
6. Per genome x core block: intact (all genes present, every block adjacency
   present), partial (>= 1 block adjacency), dispersed (>= 2 block genes
   present, none adjacent as in the block), absent (< 2 genes present).

Run from the project root:
  python3 scripts/07_combined_synteny/collinear_blocks.py
Output: gene_order/collinear_blocks/{core_blocks.tsv, genome_blocks.tsv,
        group_summary.tsv, lineage_specific.tsv, lineage_loci.tsv,
        fig_collinear_blocks.png, fig_lineage_loci.png}
"""
import json
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
sys.path.insert(0, os.path.join(HERE, "..", "11_flagellar_trees"))

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Patch, Rectangle
from ete3 import Tree

import synteny_distance as sd
from plot_results import TEXT, TEXT2, GROUPS, group_color, phylogram_layout, draw_tree

RELIABLE_JSON = "gene_order/synteny_analysis_combined_reliable.json"
GENE_ORDER_JSON = "gene_order/gene_order_combined_scaffolded.json"
NAMED_NWK = "images_and_newick/thesis_b/gtdb_ref_tree.nwk"
OUT = "gene_order/collinear_blocks"
GAP_BP = 300
LINEAGE_FRAC = 0.5
MIN_LINEAGE = 3
STATUS = ["intact", "partial", "dispersed", "absent"]
STATUS_COLOR = {"intact": "#1f5f8b", "partial": "#8fbcd9", "dispersed": "#f0b46a", "absent": "#f4f3f0"}
BLOCK_COLORS = ["#2a78d6", "#e34948", "#1baf7a", "#eda100", "#8e5cc4", "#e87ba4", "#008300",
                "#4a3aa7", "#c2702b", "#3aa6b9", "#9c9a2f", "#b84f86"]
OTHER_GENE = "#a9a7a2"


def pretty(g):
    return g[:3] + g[3:].upper() if len(g) == 4 and g[:3] in ("fli", "flg", "flh", "che", "mot", "fla") else g


def group_of(org):
    genus = org.split()[0]
    return next((name for name, _, genera in GROUPS if genus in genera), "other")


def short_group(g):
    return g.split(" (")[0]


def physical_adjacencies(recs):
    """{signed adjacency: gap} for consecutive genes with gap <= GAP_BP."""
    out = {}
    for a, b in zip(recs, recs[1:]):
        gap = b["start"] - a["end"] - 1
        if gap <= GAP_BP:
            out[frozenset((sd._ends(a)[1], sd._ends(b)[0]))] = gap
    return out


def adj_of(a, sa, b, sb):
    return frozenset((sd._ends({"gene": a, "strand": sa})[1], sd._ends({"gene": b, "strand": sb})[0]))


def chain(adjs_ranked):
    """Greedy linear chains from (adjacency, weight) pairs, best first."""
    parent, used, acc = {}, set(), {}

    def find(x):
        while parent.get(x, x) != x:
            x = parent[x]
        return x

    for adj, w in adjs_ranked:
        e1, e2 = tuple(adj)
        if e1 in used or e2 in used or find(e1[0]) == find(e2[0]):
            continue
        parent[find(e1[0])] = find(e2[0])
        used |= {e1, e2}
        acc[adj] = w
    nbr = defaultdict(dict)
    for adj in acc:
        e1, e2 = tuple(adj)
        nbr[e1[0]][e1[1]] = e2
        nbr[e2[0]][e2[1]] = e1
    chains, seen = [], set()
    for g in sorted(nbr):
        if g in seen:
            continue
        comp, stack = {g}, [g]
        while stack:
            for e in nbr[stack.pop()].values():
                if e[0] not in comp:
                    comp.add(e[0])
                    stack.append(e[0])
        term = sorted(x for x in comp if len(nbr[x]) == 1)[0]
        enter = "tail" if "tail" not in nbr[term] else "head"
        genes, x = [], term
        while True:
            genes.append((x, "+" if enter == "tail" else "-"))
            seen.add(x)
            out = "head" if enter == "tail" else "tail"
            if out not in nbr[x]:
                break
            x, enter = nbr[x][out]
        if sum(s == "-" for _, s in genes) > len(genes) / 2 or (genes[0][1] == "-" and
                                                               sum(s == "-" for _, s in genes) == len(genes) / 2):
            genes = [(gg, "+" if s == "-" else "-") for gg, s in reversed(genes)]
        adjs = [adj_of(a, sa, b, sb) for (a, sa), (b, sb) in zip(genes, genes[1:])]
        chains.append({"genes": genes, "adjs": adjs, "support": [acc[a] for a in adjs]})
    chains.sort(key=lambda c: (-len(c["genes"]), -sum(c["support"])))
    return chains


def status(block, genes_present, adjs):
    present = sum(g in genes_present for g, _ in block["genes"])
    kept = sum(a in adjs for a in block["adjs"])
    if present < 2:
        return "absent"
    if present == len(block["genes"]) and kept == len(block["adjs"]):
        return "intact"
    return "partial" if kept else "dispersed"


def label(c):
    return " ".join(pretty(g) for g, _ in c["genes"])


def signed_label(c):
    return " ".join(pretty(g) + s for g, s in c["genes"])


def adj_label(adj):
    (g1, e1), (g2, e2) = sorted(adj)
    return f"{pretty(g1)}({e1[0]})–{pretty(g2)}({e2[0]})"


def main():
    os.makedirs(OUT, exist_ok=True)
    orgs = json.load(open(RELIABLE_JSON))["matrix"]["organisms"]
    go = json.load(open(GENE_ORDER_JSON))
    acc_of = {v["organism"]: a for a, v in go.items()}
    recs = {o: sd.primary_contig_records(go[acc_of[o]]["genes"])[0] for o in orgs}
    padj = {o: physical_adjacencies(r) for o, r in recs.items()}
    counts = Counter(a for s in padj.values() for a in s)

    members = defaultdict(list)
    for o in orgs:
        members[group_of(o)].append(o)
    groups = [g for g, _, _ in GROUPS if g in members]
    voting = [g for g in groups if len(members[g]) >= MIN_LINEAGE]
    held = {a: [g for g in voting if sum(a in padj[o] for o in members[g]) >= LINEAGE_FRAC * len(members[g])]
            for a in counts}

    # ---- core blocks and lineage-specific adjacencies
    core = chain(sorted(((a, counts[a]) for a in counts if len(held[a]) >= 2), key=lambda x: -x[1]))
    for i, b in enumerate(core, 1):
        b["id"] = f"C{i}"
    color_of = {g: BLOCK_COLORS[i % len(BLOCK_COLORS)] for i, b in enumerate(core) for g, _ in b["genes"]}
    block_of = {g: b["id"] for b in core for g, _ in b["genes"]}
    specific = sorted((a for a in counts if len(held[a]) == 1), key=lambda a: (held[a][0], -counts[a]))

    # ---- lineage consensus loci
    lineage_loci = {}
    for g in voting:
        mine = [(a, sum(a in padj[o] for o in members[g])) for a in counts if g in held[a]]
        lineage_loci[g] = chain(sorted(mine, key=lambda x: -x[1]))

    # ---- tree order and per-genome status
    tree = Tree(NAMED_NWK, format=1, quoted_node_names=True)
    tree.prune(orgs, preserve_branch_length=True)
    tree.ladderize()
    order = tree.get_leaf_names()
    st = {o: {b["id"]: status(b, {r["gene"] for r in recs[o]}, padj[o]) for b in core} for o in order}
    sp = {o: [a in padj[o] for a in specific] for o in order}

    # ---- tables
    with open(f"{OUT}/core_blocks.tsv", "w") as f:
        f.write("block\tgenes_signed\tadjacency_support\tlineages_holding_every_adjacency\t"
                + "\t".join(f"n_{s}" for s in STATUS) + "\n")
        for b in core:
            c = Counter(st[o][b["id"]] for o in order)
            hold_all = [g for g in voting if all(g in held[a] for a in b["adjs"])]
            f.write(f"{b['id']}\t{signed_label(b)}\t{','.join(map(str, b['support']))}\t"
                    f"{'; '.join(map(short_group, hold_all))}\t" + "\t".join(str(c[s]) for s in STATUS) + "\n")
    with open(f"{OUT}/genome_blocks.tsv", "w") as f:
        f.write("organism\tgroup\t" + "\t".join(b["id"] for b in core) + "\t"
                + "\t".join(adj_label(a) for a in specific) + "\n")
        for o in order:
            f.write(f"{o}\t{group_of(o)}\t" + "\t".join(st[o][b["id"]] for b in core) + "\t"
                    + "\t".join("yes" if x else "no" for x in sp[o]) + "\n")
    summary = {}
    with open(f"{OUT}/group_summary.tsv", "w") as f:
        f.write("block\tgenes\t" + "\t".join(groups) + "\n")
        for b in core:
            row = []
            for g in groups:
                k = sum(st[o][b["id"]] == "intact" for o in members[g])
                summary[(b["id"], g)] = (k, len(members[g]))
                row.append(f"{k}/{len(members[g])}")
            f.write(f"{b['id']}\t{label(b)}\t" + "\t".join(row) + "\n")
    with open(f"{OUT}/lineage_specific.tsv", "w") as f:
        f.write("adjacency\tlineage\tn_in_lineage\tlineage_size\tn_total\tjoins_core_blocks\n")
        for a in specific:
            g = held[a][0]
            (g1, _), (g2, _) = tuple(a)
            f.write(f"{adj_label(a)}\t{g}\t{sum(a in padj[o] for o in members[g])}\t{len(members[g])}\t"
                    f"{counts[a]}\t{block_of.get(g1, '-')}–{block_of.get(g2, '-')}\n")
    with open(f"{OUT}/lineage_loci.tsv", "w") as f:
        f.write("lineage\tlocus\tgenes_signed\tadjacency_support\tlineage_size\n")
        for g in voting:
            for k, c in enumerate(lineage_loci[g], 1):
                f.write(f"{g}\tL{k}\t{signed_label(c)}\t{','.join(map(str, c['support']))}\t{len(members[g])}\n")

    # ---- console
    print(f"{len(order)} reliable genomes; adjacent = gap <= {GAP_BP} bp; lineage holds an adjacency "
          f"if >= {LINEAGE_FRAC:.0%} of its members have it (lineages with >= {MIN_LINEAGE} genomes)")
    print("\nCore collinear blocks (every adjacency held by >= 2 lineages):")
    for b in core:
        c = Counter(st[o][b["id"]] for o in order)
        print(f"  {b['id']:4} {signed_label(b):34} " + "  ".join(f"{s} {c[s]:3}" for s in STATUS))
    print("\nIntact, by lineage:")
    print(f"{'':6}" + "".join(f"{g.split()[0][:12]:>14}" for g in groups))
    for b in core:
        print(f"{b['id']:6}" + "".join(f"{'%d/%d' % summary[(b['id'], g)]:>14}" for g in groups))
    print("\nLineage-specific adjacencies:")
    for a in specific:
        g = held[a][0]
        (g1, _), (g2, _) = tuple(a)
        print(f"  {short_group(g):28} {adj_label(a):30} {sum(a in padj[o] for o in members[g])}/{len(members[g])}"
              f"  joins {block_of.get(g1, '-')}–{block_of.get(g2, '-')}")
    print("\nLineage consensus loci:")
    for g in voting:
        print(f"  {g} ({len(members[g])}):  " + "  |  ".join(label(c) for c in lineage_loci[g]))

    # ---- figure 1: tree + core-block status + lineage-specific adjacencies
    n, nb, ns = len(order), len(core), len(specific)
    fig = plt.figure(figsize=(17, 0.16 * n + 3.0))
    ax = fig.add_axes([0.01, 0.06, 0.98, 0.91])
    tree_w, lab_w, strip_w, cell, gap = 3.0, 3.6, 0.35, 0.5, 0.6
    y = {name: -i for i, name in enumerate(order)}
    pos, _ = phylogram_layout(tree, y)
    draw_tree(ax, tree, pos, 0, tree_w, False, lambda nd: (TEXT2, 0.6, "-"))
    x_lab, x_strip = tree_w + 0.08, tree_w + 0.08 + lab_w
    x_blk = x_strip + strip_w + 0.25
    x_sp = x_blk + nb * cell + gap
    for o in order:
        ax.text(x_lab, y[o], o, fontsize=5.2, va="center", ha="left", style="italic", color=TEXT)
        ax.add_patch(Rectangle((x_strip, y[o] - 0.5), strip_w, 1, color=group_color(o) or TEXT2, lw=0))
        for k, b in enumerate(core):
            ax.add_patch(Rectangle((x_blk + k * cell, y[o] - 0.5), cell, 1,
                                   facecolor=STATUS_COLOR[st[o][b["id"]]], edgecolor="white", lw=0.3))
        for k, a in enumerate(specific):
            col = [c for gname, c, _ in GROUPS if gname == held[a][0]][0]
            ax.add_patch(Rectangle((x_sp + k * cell, y[o] - 0.5), cell, 1,
                                   facecolor=col if sp[o][k] else "#f4f3f0", edgecolor="white", lw=0.3))
    for k, b in enumerate(core):
        ax.text(x_blk + (k + 0.5) * cell, 1.0, f"{b['id']}  {label(b)}", rotation=60, fontsize=6.5,
                ha="left", va="bottom", rotation_mode="anchor", color=TEXT)
    for k, a in enumerate(specific):
        ax.text(x_sp + (k + 0.5) * cell, 1.0, adj_label(a), rotation=60, fontsize=6,
                ha="left", va="bottom", rotation_mode="anchor", color=TEXT2)
    ax.text(x_blk, 9.0, "Core collinear blocks", fontsize=10, weight="bold", color=TEXT)
    ax.text(x_sp, 9.0, "Lineage-specific adjacencies (coloured by the lineage that holds it)",
            fontsize=10, weight="bold", color=TEXT)
    ax.set_xlim(-0.1, x_sp + max(ns, 1) * cell + 3.0)
    ax.set_ylim(-n - 0.5, 10.5)
    ax.axis("off")
    handles = [Patch(color=STATUS_COLOR[s]) for s in STATUS] + [Patch(color=c) for g, c, _ in GROUPS if g in groups]
    labels = ["intact (same order + orientation, contiguous)", "partial (some block adjacencies kept)",
              "dispersed (genes present, never adjacent as in block)", "absent (< 2 block genes)"] + groups
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, fontsize=7.5)
    fig.suptitle(f"Collinear flagellar gene clusters across {n} closed-chromosome genomes "
                 f"(GTDB order; adjacent = gap ≤ {GAP_BP} bp)", fontsize=12, y=0.999)
    fig.savefig(f"{OUT}/fig_collinear_blocks.png", dpi=200, facecolor="white")
    plt.close(fig)

    # ---- figure 2: lineage consensus loci, genes coloured by core block
    gw, gap_l = 0.62, 0.45
    rows = []
    for g in voting:
        loci = [c for c in lineage_loci[g] if len(c["genes"]) >= 2]
        rows.append((g, loci))
    width = max(sum(len(c["genes"]) * gw + gap_l for c in loci) for _, loci in rows)
    fig, ax = plt.subplots(figsize=(3.2 + width * 1.05, 1.15 * len(rows) + 1.6))
    for i, (g, loci) in enumerate(rows):
        yy = -i * 1.15
        col = [c for gname, c, _ in GROUPS if gname == g][0]
        ax.text(-0.15, yy, f"{short_group(g)}\n(n = {len(members[g])})", ha="right", va="center",
                fontsize=8.5, color=col, weight="bold")
        x = 0.0
        for c in loci:
            for k, (gene, s) in enumerate(c["genes"]):
                x0, x1 = x + k * gw + 0.02, x + (k + 1) * gw - 0.02
                start, length = (x0, x1 - x0) if s == "+" else (x1, -(x1 - x0))
                ax.add_patch(FancyArrow(start, yy, length, 0, width=0.36, head_width=0.36, head_length=0.16,
                                        length_includes_head=True, color=color_of.get(gene, OTHER_GENE), lw=0))
                ax.text((x0 + x1) / 2, yy - 0.3, pretty(gene), ha="center", va="top", fontsize=5.6,
                        rotation=90, color=TEXT, style="italic")
            for k, sup in enumerate(c["support"]):
                ax.text(x + (k + 1) * gw, yy + 0.22, str(sup), ha="center", va="bottom", fontsize=4.8, color=TEXT2)
            x += len(c["genes"]) * gw + gap_l
    handles = [Patch(color=BLOCK_COLORS[i % len(BLOCK_COLORS)]) for i in range(nb)] + [Patch(color=OTHER_GENE)]
    labels = [f"{b['id']} {label(b)}" for b in core] + ["not in a core block"]
    ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False, fontsize=7,
              title="gene colour = core block", title_fontsize=7.5)
    ax.set_xlim(-0.2, width + 0.2)
    ax.set_ylim(-(len(rows) - 1) * 1.15 - 1.0, 0.7)
    ax.axis("off")
    ax.set_title("Consensus collinear loci per lineage (arrow = strand; small number = genomes in the "
                 "lineage with that adjacency; loci separated by gaps)", fontsize=9, color=TEXT, loc="left")
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_lineage_loci.png", dpi=200, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"\nWrote {OUT}/")


if __name__ == "__main__":
    main()
