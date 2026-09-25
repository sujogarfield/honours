#!/usr/bin/env python3
"""
orthology_overlay_preview.py
TEMPORARY review artifact: renders the tree with BOTH the existing
string-matched flagellar genes (gene_order/gene_order.json, solid arrows)
and the new orthology-recovered candidates
(gene_order/gene_order_orthology_candidates.json, hatched arrows) so the
new calls can be eyeballed before touching the real pipeline
(gene_order_overlay.py / synteny_analysis.py / dashboards).

Orthology candidates are drawn with a diagonal hatch and a dashed-look
outline; nothing here is folded into gene_order.json.

Run AFTER orthology_crossref.py.
"""

import json
import tempfile
from pathlib import Path
from collections import Counter
from PIL import Image, ImageDraw, ImageFont
from ete3 import Tree, TextFace, NodeStyle, TreeStyle, faces

GENE_ORDER_JSON       = "campy_orthologs/archive_json/gene_order.json"
CANDIDATES_JSON       = "campy_orthologs/archive_json/gene_order_orthology_candidates.json"
MAF_CANDIDATES_JSON   = "campy_orthologs/archive_json/gene_order_maf_candidates.json"
NAMED_NWK             = "images_and_newick/thesis_b/gtdb_ref_tree.nwk"
OUTPUT_PNG            = "images_and_newick/thesis_c/_tmp_orthology_overlay.png"

BLOCK_H   = 14
BLOCK_W   = 22
GAP       = 3
ARROW_TIP = 6
LABEL_H   = 13

FONT_PATH = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from gene_colors import generate_gene_colors, text_color_for

# Populated in main() from whatever gene names actually appear in the loaded
# data -- randomised (with a minimum colour-distance threshold) rather than
# hand-picked, so a new gene name never needs a manually-chosen hex code.
GENE_COLORS: dict[str, str] = {}


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


def load_gene_order(path: str) -> dict:
    with open(path) as f:
        raw = json.load(f)
    by_organism = {}
    for data in raw.values():
        by_organism.setdefault(data["organism"], []).extend(data["genes"])
    return by_organism


def merge_sources(known: dict, candidates: dict) -> dict:
    """Tags known records source='string_match' and merges in candidates
    (already tagged source='orthology' by orthology_crossref.py)."""
    merged = {}
    for organism, recs in known.items():
        for r in recs:
            r.setdefault("source", "string_match")
        merged.setdefault(organism, []).extend(recs)
    for organism, recs in candidates.items():
        merged.setdefault(organism, []).extend(recs)
    return merged


STATUS_OUTLINE = {
    "confirmed":    "#FF00FF",  # magenta -- eggNOG independently agrees
    "reclassified": "#FF8C00",  # orange  -- eggNOG names a specific, different gene
    "unconfirmed":  "#888888",  # grey    -- eggNOG hit, but no name/description support
    "no_hit":       "#BBBBBB",  # light grey -- no eggNOG hit at all
}


def _draw_hatch(img, pts, x0, x1, top, bot, colour=(20, 20, 20, 160)):
    """Diagonal hatch clipped to the arrow polygon, marking a candidate
    gene that orthology found but the annotation string match missed."""
    w, h = x1 - x0, bot - top
    if w <= 0 or h <= 0:
        return
    mask = Image.new("L", (w, h), 0)
    ImageDraw.Draw(mask).polygon([(px - x0, py - top) for px, py in pts], fill=255)
    stripes = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    sdraw = ImageDraw.Draw(stripes)
    for sx in range(-h, w, 4):
        sdraw.line([(sx, 0), (sx + h, h)], fill=colour, width=1)
    img.paste(stripes, (x0, top), mask=mask)


def make_gene_face(records: list[dict]):
    if not records:
        return None

    sorted_recs = sorted(records, key=lambda r: r["start"])

    gene_counts = Counter(r["gene"].lower() for r in sorted_recs)
    duplicated  = {g for g, c in gene_counts.items() if c > 1}

    total_w = len(sorted_recs) * (BLOCK_W + GAP) + GAP
    total_h   = BLOCK_H + 6 + LABEL_H
    ARROW_TOP = 3
    ARROW_BOT = ARROW_TOP + BLOCK_H
    AXIS_Y    = ARROW_BOT + 3
    TICK_Y    = ARROW_BOT + 1
    LABEL_Y   = AXIS_Y + 2
    DOT_R     = 3

    img  = Image.new("RGBA", (total_w, total_h), (255, 255, 255, 0))
    draw = ImageDraw.Draw(img)
    font = _load_font(7)

    arrow_centres = []

    x = GAP
    for rec in sorted_recs:
        gene   = rec["gene"].lower()
        color  = GENE_COLORS.get(gene, "#AAAAAA")
        strand = rec.get("strand", "+")
        is_new = rec.get("source") in ("orthology", "orthology_external_seed")

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
        if is_new:
            status = rec.get("eggnog_status", "no_hit")
            outline = STATUS_OUTLINE.get(status, "#888888")
            _draw_hatch(img, pts, x0, x1, ARROW_TOP, ARROW_BOT)
            draw.polygon(pts, outline=outline, width=2)
        else:
            draw.polygon(pts, outline="#333333")

        label = gene[:4]
        bbox  = draw.textbbox((0, 0), label, font=font)
        tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
        draw.text(
            (x0 + (BLOCK_W - tw) // 2, ARROW_TOP + (BLOCK_H - th) // 2),
            label, fill=text_color_for(color), font=font,
        )

        if gene in duplicated:
            dx = x1 - DOT_R - 1
            dy = ARROW_TOP + DOT_R + 1
            draw.ellipse([dx - DOT_R, dy - DOT_R, dx + DOT_R, dy + DOT_R],
                         fill="#FF0000", outline="#AA0000")

        arrow_centres.append((cx, rec["start"]))
        x += BLOCK_W + GAP

    draw.line([(GAP, AXIS_Y), (total_w - GAP, AXIS_Y)], fill="#888888", width=1)

    for cx, pos in arrow_centres:
        draw.line([(cx, TICK_Y), (cx, AXIS_Y)], fill="#888888", width=1)
        lbl  = _fmt_pos(pos)
        bbox = draw.textbbox((0, 0), lbl, font=font)
        tw   = bbox[2] - bbox[0]
        draw.text((cx - tw // 2, LABEL_Y), lbl, fill="#444444", font=font)

    return _save_img(img)


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


def build_key() -> str:
    """Solid vs hatched-by-status key, separate from the gene-colour legend."""
    rows = [
        (None, "string match (existing gene= tag)"),
        ("confirmed", "orthology + eggNOG independently agree"),
        ("reclassified", "eggNOG names a different specific gene (inspect)"),
        ("unconfirmed", "orthology only, eggNOG hit gave no name support"),
        ("no_hit", "orthology only, no eggNOG hit at all"),
    ]
    w, h = 480, 22 * len(rows) + 10
    img = Image.new("RGB", (w, h), (255, 255, 255))
    draw = ImageDraw.Draw(img)
    font = _load_font(11)

    for i, (status, label) in enumerate(rows):
        y = 8 + i * 22
        if status is None:
            draw.rectangle([10, y, 34, y + 16], fill="#3498DB", outline="#333")
        else:
            demo = Image.new("RGBA", (24, 16), (255, 255, 255, 255))
            hd = ImageDraw.Draw(demo)
            hd.rectangle([0, 0, 23, 15], fill="#3498DB",
                         outline=STATUS_OUTLINE[status], width=2)
            for sx in range(-16, 24, 4):
                hd.line([(sx, 0), (sx + 16, 16)], fill=(20, 20, 20, 160), width=1)
            img.paste(demo, (10, y))
        draw.text((40, y + 1), label, fill="black", font=font)

    path = tempfile.NamedTemporaryFile(suffix=".png", delete=False).name
    img.save(path)
    return path


def main():
    Path("images_and_newick/thesis_c").mkdir(exist_ok=True)

    print("Loading gene order + orthology candidates ...")
    known          = load_gene_order(GENE_ORDER_JSON)
    candidates     = load_gene_order(CANDIDATES_JSON)
    maf_candidates = load_gene_order(MAF_CANDIDATES_JSON)
    for organism, recs in maf_candidates.items():
        candidates.setdefault(organism, []).extend(recs)
    gene_order = merge_sources(known, candidates)

    print(f"Loading tree from {NAMED_NWK} ...")
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
    n_new_total = 0
    for leaf in tree.iter_leaves():
        name_face = TextFace(leaf.name, fstyle="italic", fsize=9)
        name_face.margin_right = 8
        leaf.add_face(name_face, column=0, position="branch-right")

        records = gene_order.get(leaf.name, [])
        n_new_total += sum(1 for r in records if r.get("source") in ("orthology", "orthology_external_seed"))
        if records:
            gf = make_gene_face(records)
            if gf:
                leaf.add_face(gf, column=1, position="aligned")
            matched += 1
        else:
            unmatched += 1

    print(f"Gene blocks: {matched} matched, {unmatched} unmatched")
    print(f"Orthology-recovered candidate genes drawn: {n_new_total}")

    used_genes = sorted({g["gene"].lower() for recs in gene_order.values() for g in recs})
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_scale     = False
    ts.title.add_face(
        TextFace("TEMP PREVIEW — string match + orthology-recovered candidates (unconfirmed)",
                  fsize=14, bold=True),
        column=0,
    )
    ts.legend.add_face(faces.ImgFace(build_key()), column=0)
    ts.legend.add_face(faces.ImgFace(build_legend(used_genes)), column=0)
    ts.legend_position = 4

    print(f"Rendering → {OUTPUT_PNG}")
    tree.render(OUTPUT_PNG, w=3200, units="px", tree_style=ts)
    print("Done.")


if __name__ == "__main__":
    main()
