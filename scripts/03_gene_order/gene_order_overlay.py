#!/usr/bin/env python3
"""
gene_order_overlay.py
Overlays flagellar gene order (genomic position) on the GTDB reference tree.

Within each genome, genes are drawn in chromosomal position order.
For draft assemblies with genes on multiple contigs, the contig containing the
most flagellar genes is drawn first; subsequent contigs are separated by '//'.
The relative order of contigs is unknown in draft assemblies — '//' makes that
explicit.

Run AFTER gtdb_tree.py has produced images_and_newick/thesis_b/gtdb_ref_tree.nwk
"""

import json
import sys
import tempfile
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
from ete3 import Tree, TextFace, NodeStyle, TreeStyle, faces

# ── Config ────────────────────────────────────────────────────────────────────
# Usage: gene_order_overlay.py [gene_order_json] [output_png]
DEFAULT_GENE_ORDER_JSON = "campy_orthologs/archive_json/gene_order.json"
NAMED_NWK               = "images_and_newick/thesis_b/gtdb_ref_tree.nwk"
DEFAULT_OUTPUT_PNG      = "images_and_newick/thesis_c/gtdb_ref_tree_gene_overlay.png"

BLOCK_H   = 14
BLOCK_W   = 22
GAP       = 3
ARROW_TIP = 6
LABEL_H   = 13   # px reserved below arrows for axis + position labels

FONT_PATH = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from gene_colors import generate_gene_colors, text_color_for

# Populated in main() from whatever gene names actually appear in the loaded
# data -- randomised (min-distance threshold, see gene_colors.py) rather than
# hand-picked, so a new gene name never needs a manually-chosen hex code.
GENE_COLORS: dict[str, str] = {}

# ── Helpers ───────────────────────────────────────────────────────────────────
def _load_font(size: int):
    try:
        return ImageFont.truetype(FONT_PATH, size)
    except Exception:
        return ImageFont.load_default()


def _save_img(img) -> faces.ImgFace:
    tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
    img.save(tmp.name)
    tmp.close()
    return faces.ImgFace(tmp.name)


def _fmt_pos(bp: int) -> str:
    if bp >= 1_000_000:
        return f"{bp / 1_000_000:.1f}M"
    if bp >= 1_000:
        return f"{bp // 1000}k"
    return str(bp)


# ── Data loading ──────────────────────────────────────────────────────────────
def load_gene_order(path: str) -> dict:
    """Returns {organism_name: [gene_records]} from gene_order.json."""
    with open(path) as f:
        raw = json.load(f)
    by_organism = {}
    for data in raw.values():
        by_organism.setdefault(data["organism"], []).extend(data["genes"])
    return by_organism


# ── Gene arrow row with position track ────────────────────────────────────────
def make_gene_face(records: list[dict]):
    """Draw genes in chromosomal position order with a coordinate track below.

    Contigs are ordered by flagellar gene count (most first) so the primary
    cluster appears on the left. '//' separates contig boundaries.
    """
    if not records:
        return None

    sorted_recs = sorted(records, key=lambda r: r["start"])

    # Genes appearing more than once — all occurrences get a red dot
    from collections import Counter
    gene_counts = Counter(r["gene"].lower() for r in sorted_recs)
    duplicated  = {g for g, c in gene_counts.items() if c > 1}

    total_w = len(sorted_recs) * (BLOCK_W + GAP) + GAP
    total_h   = BLOCK_H + 6 + LABEL_H
    ARROW_TOP = 3
    ARROW_BOT = ARROW_TOP + BLOCK_H
    AXIS_Y    = ARROW_BOT + 3
    TICK_Y    = ARROW_BOT + 1
    LABEL_Y   = AXIS_Y + 2
    DOT_R     = 3   # radius of duplication dot

    img  = Image.new("RGBA", (total_w, total_h), (255, 255, 255, 0))
    draw = ImageDraw.Draw(img)
    font = _load_font(7)

    arrow_centres = []

    x = GAP
    for rec in sorted_recs:
        gene   = rec["gene"].lower()
        color  = GENE_COLORS.get(gene, "#AAAAAA")
        strand = rec.get("strand", "+")

        x0, x1 = x, x + BLOCK_W
        mid_y  = (ARROW_TOP + ARROW_BOT) // 2
        cx     = x0 + BLOCK_W // 2

        if strand == "+":
            pts = [(x0, ARROW_TOP), (x1 - ARROW_TIP, ARROW_TOP),
                   (x1, mid_y),
                   (x1 - ARROW_TIP, ARROW_BOT), (x0, ARROW_BOT)]
        else:
            pts = [(x0 + ARROW_TIP, ARROW_TOP), (x1, ARROW_TOP),
                   (x1, ARROW_BOT),
                   (x0 + ARROW_TIP, ARROW_BOT), (x0, mid_y)]

        draw.polygon(pts, fill=color)
        draw.polygon(pts, outline="#333333")

        label = gene[:4]
        bbox  = draw.textbbox((0, 0), label, font=font)
        tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
        draw.text(
            (x0 + (BLOCK_W - tw) // 2, ARROW_TOP + (BLOCK_H - th) // 2),
            label, fill=text_color_for(color), font=font,
        )

        # Red dot in top-right corner for duplicated genes
        if gene in duplicated:
            dx = x1 - DOT_R - 1
            dy = ARROW_TOP + DOT_R + 1
            draw.ellipse([dx - DOT_R, dy - DOT_R, dx + DOT_R, dy + DOT_R],
                         fill="#FF0000", outline="#AA0000")

        arrow_centres.append((cx, rec["start"]))
        x += BLOCK_W + GAP

    # Axis line across full width
    draw.line([(GAP, AXIS_Y), (total_w - GAP, AXIS_Y)], fill="#888888", width=1)

    # Tick marks and position labels
    for cx, pos in arrow_centres:
        draw.line([(cx, TICK_Y), (cx, AXIS_Y)], fill="#888888", width=1)
        lbl  = _fmt_pos(pos)
        bbox = draw.textbbox((0, 0), lbl, font=font)
        tw   = bbox[2] - bbox[0]
        draw.text((cx - tw // 2, LABEL_Y), lbl, fill="#444444", font=font)

    return _save_img(img)


# ── Legend ─────────────────────────────────────────────────────────────────────
def build_legend(genes: list[str], cols: int = 8) -> str:
    cell_w, cell_h, pad = 80, 20, 4
    rows = (len(genes) + cols - 1) // cols
    img  = Image.new("RGB", (cols * cell_w + pad, rows * cell_h + pad * 2), (255, 255, 255))
    draw = ImageDraw.Draw(img)
    font = _load_font(9)
    for idx, gene in enumerate(genes):
        x = pad + (idx % cols) * cell_w
        y = pad + (idx // cols) * cell_h
        draw.rectangle([x, y + 2, x + 14, y + 14],
                       fill=GENE_COLORS.get(gene, "#AAAAAA"), outline="#333")
        draw.text((x + 18, y + 2), gene, fill="black", font=font)
    path = tempfile.NamedTemporaryFile(suffix=".png", delete=False).name
    img.save(path)
    return path


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    gene_order_json = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_GENE_ORDER_JSON
    output_png = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_OUTPUT_PNG

    Path("images_and_newick/thesis_c").mkdir(exist_ok=True)

    print(f"Loading gene order data from {gene_order_json} …")
    gene_order = load_gene_order(gene_order_json)

    print(f"Loading tree from {NAMED_NWK} …")
    tree = Tree(NAMED_NWK, format=1, quoted_node_names=True)

    for node in tree.traverse():
        node.dist = 1

    nstyle = NodeStyle()
    nstyle["size"] = 0
    for node in tree.traverse():
        node.set_style(nstyle)

    try:
        from helpers import colour_tree_by_species
        colour_tree_by_species(tree)
    except ImportError:
        print("helpers.colour_tree_by_species not found — skipping bg colours")

    all_genes = {g["gene"].lower() for recs in gene_order.values() for g in recs}
    GENE_COLORS.update(generate_gene_colors(all_genes))
    print(f"Generated {len(GENE_COLORS)} distinct gene colours (min-distance randomised)")

    matched = unmatched = 0
    for leaf in tree.iter_leaves():
        name_face = TextFace(leaf.name, fstyle="italic", fsize=9)
        name_face.margin_right = 8
        leaf.add_face(name_face, column=0, position="branch-right")

        records = gene_order.get(leaf.name, [])
        if records:
            gf = make_gene_face(records)
            if gf:
                leaf.add_face(gf, column=1, position="aligned")
            matched += 1
        else:
            unmatched += 1

    print(f"Gene blocks: {matched} matched, {unmatched} unmatched")

    used_genes = sorted({g["gene"].lower() for recs in gene_order.values() for g in recs})
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_scale     = False
    ts.title.add_face(
        TextFace("Flagellar gene synteny — genomic position order", fsize=14, bold=True),
        column=0,
    )
    ts.legend.add_face(faces.ImgFace(build_legend(used_genes)), column=0)
    ts.legend_position = 4

    print(f"Rendering → {output_png}")
    tree.render(output_png, w=3200, units="px", tree_style=ts)
    print("Done.")


if __name__ == "__main__":
    main()
