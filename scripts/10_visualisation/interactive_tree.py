#!/usr/bin/env python3
"""Generate a standalone interactive HTML tree with flagellar gene-order
glyphs and gene-order/phylogeny agreement scores per leaf.

Usage:
    python3 interactive_tree.py [named_nwk] [og_nwk] [output_html] [gene_order_json] [synteny_json] [motility_json]

Defaults:
    named_nwk       = images_and_newick/thesis_b/gtdb_ref_tree.nwk    (organism names)
    og_nwk          = images_and_newick/thesis_b/gtdb_ref_tree_og.nwk (GCF IDs)
    output          = frontend/interactive_tree.html
    gene_order_json = gene_order/gene_order_combined.json
    synteny_json    = gene_order/synteny_analysis_combined_merged.json (optional, from build_merged_synteny_combined.py)
    motility_json   = gene_order/motility.json (optional, from fetch_motility.py)

Hovering over a leaf shows the GCF accession ID; hovering over a gene mark
shows its position/strand/contig; the small dot next to each leaf name shows
its gene-order/phylogeny agreement score (see synteny_analysis.py); the dot
after that shows its BacDive motility phenotype (see fetch_motility.py). A
checklist lets you filter which genes are drawn, and the layout can switch
between rectangular and circular (radial) dendrogram styles.
"""

import base64
import json
import sys
from collections import Counter
from pathlib import Path
from ete3 import Tree

MOTILITY_ICON_DIR = Path("frontend/motility_icons")
MOTILITY_ICON_FILES = {
    "motile": "motile.png",
    "nonmotile": "nonmotile.png",
    "predicted_motile": "predicted_motile.png",
    "predicted_nonmotile": "predicted_nonmotile.png",
}


def load_motility_icons() -> dict:
    """{key: 'data:image/png;base64,...'} for each motility icon, embedded
    inline so the generated HTML stays a single self-contained file (same
    reasoning as embedding the tree/gene data directly rather than fetching
    it). Empty dict (icons silently skipped, falls back to the plain grey
    dot) if the icon files aren't present -- keeps this script runnable even
    without them rather than hard failing."""
    icons = {}
    for key, filename in MOTILITY_ICON_FILES.items():
        path = MOTILITY_ICON_DIR / filename
        if not path.exists():
            print(f"Warning: {path} not found — motility icon '{key}' will fall back to a plain dot")
            continue
        b64 = base64.b64encode(path.read_bytes()).decode("ascii")
        icons[key] = f"data:image/png;base64,{b64}"
    return icons

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from gene_colors import generate_gene_colors

# Populated in main() from whatever gene names actually appear in the loaded
# gene_order data -- randomised (min-distance threshold, see gene_colors.py)
# rather than hand-picked, so a new gene name never needs a manually-chosen
# hex code. Same dict object is reused by generate_html()'s JS embed so the
# browser draws the same colours as this run's ete3-adjacent Python side.
GENE_COLORS: dict[str, str] = {}


def load_gene_order(path: str) -> dict:
    """Returns {organism_name: [gene_records]} from gene_order.json."""
    if not Path(path).exists():
        print(f"Warning: {path} not found — gene arrows will be omitted")
        return {}
    with open(path) as f:
        raw = json.load(f)
    by_organism = {}
    for data in raw.values():
        by_organism.setdefault(data["organism"], []).extend(data["genes"])
    return by_organism


def load_synteny(path: str) -> dict:
    """Returns {organism_name: {agreement_score, coverage, ...}} from
    synteny_analysis.json, or {} if it hasn't been generated yet."""
    if not Path(path).exists():
        print(f"Warning: {path} not found — agreement scores will be omitted")
        return {}
    with open(path) as f:
        data = json.load(f)
    return data.get("organisms", {})


def load_motility(path: str) -> dict:
    """Returns {organism_name: {motility, flagellum_arrangement, ...}} from
    motility.json (see fetch_motility.py), or {} if it hasn't been generated yet."""
    if not Path(path).exists():
        print(f"Warning: {path} not found — motility markers will be omitted")
        return {}
    with open(path) as f:
        return json.load(f)


def order_records_for_display(records: list[dict]) -> list[dict]:
    """Order gene records for genome-order display: contigs sorted by
    flagellar-gene count (most first, i.e. the primary contig leads), genes
    within a contig sorted by position. A {'separator': True} entry marks
    the boundary between contigs, since their relative order is unknown in
    draft assemblies (mirrors the '//' convention in gene_order_overlay.py).
    Genes duplicated anywhere in the genome are flagged for a dot marker."""
    if not records:
        return []

    by_contig: dict[str, list[dict]] = {}
    for r in records:
        by_contig.setdefault(r["contig"], []).append(r)

    contig_groups = sorted(by_contig.values(), key=len, reverse=True)
    gene_counts = Counter(r["gene"].lower() for r in records)

    flat = []
    for gi, group in enumerate(contig_groups):
        if gi > 0:
            flat.append({"separator": True})
        for r in sorted(group, key=lambda r: r["start"]):
            gene = r["gene"].lower()
            flat.append({
                "gene": gene,
                "start": r["start"],
                "end": r["end"],
                "strand": r.get("strand", "+"),
                "contig": r["contig"],
                "duplicated": gene_counts[gene] > 1,
            })
    return flat


def tree_to_dict(node, gcf_map, gene_order, synteny, motility):
    d = {
        "name": node.name or "",
        "length": round(node.dist, 6),
    }
    if node.is_leaf():
        d["gcf"] = gcf_map.get(node.name, "")
        d["genes"] = order_records_for_display(gene_order.get(node.name, []))
        syn = synteny.get(node.name)
        d["agreement"] = syn["agreement_score"] if syn else None
        d["coverage"] = syn["coverage"] if syn else None
        d["reliable"] = bool(syn["reliable"]) if syn and "reliable" in syn else None
        d["scaffolded"] = bool(syn["scaffolded"]) if syn and "scaffolded" in syn else False
        mot = motility.get(node.name)
        d["motility"] = mot["motility"] if mot else None
        d["flagellum_arrangement"] = mot["flagellum_arrangement"] if mot else None
        d["motility_madin"] = mot.get("motility_madin") if mot else None
        d["predicted_motility"] = mot.get("predicted_motility") if mot else None
    if node.children:
        d["children"] = [tree_to_dict(c, gcf_map, gene_order, synteny, motility) for c in node.children]
    return d


def build_gcf_map(t_named, t_og):
    """Pair leaves from both trees (same topology) to build name → GCF map."""
    named_leaves = t_named.get_leaves()
    og_leaves = t_og.get_leaves()
    if len(named_leaves) != len(og_leaves):
        print(f"Warning: leaf count mismatch ({len(named_leaves)} vs {len(og_leaves)})")
    return {l.name: r.name for l, r in zip(named_leaves, og_leaves)}


def generate_html(tree_json, leaf_count, output_path):
    LEAF_SPACING = 26
    tree_data_js = json.dumps(tree_json)
    gene_colors_js = json.dumps(GENE_COLORS)
    motility_icons = load_motility_icons()
    motility_icons_js = json.dumps(motility_icons)
    # .get(..., "") for the legend <img> tags below: if an icon is missing,
    # falls back to a blank src (broken-image icon) rather than a Python
    # KeyError, matching the JS side's own fallback-to-grey-dot behaviour.
    icon = lambda key: motility_icons.get(key, "")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Interactive Phylogenetic Tree — Flagellar Gene Order</title>
<style>
  html, body {{
    margin: 0; height: 100%;
    font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
    background: #fcfcfb; color: #0b0b0b;
  }}
  #controls {{
    position: sticky; top: 0; z-index: 90;
    padding: 8px 16px; background: #fff; border-bottom: 1px solid #e1e0d9; display: flex; gap: 14px; align-items: center; flex-wrap: wrap;
  }}
  #controls label {{ font-size: 13px; color: #52514e; }}
  #controls button {{
    font-size: 13px; padding: 4px 10px; cursor: pointer;
    border: 1px solid #e1e0d9; background: #f9f9f7; border-radius: 4px; color: #0b0b0b;
  }}
  #modeGroup {{ display: inline-flex; border: 1px solid #e1e0d9; border-radius: 4px; overflow: hidden; }}
  #modeGroup button {{ border: none; border-radius: 0; }}
  #modeGroup button.active {{ background: #2a78d6; color: #fff; }}
  .hint {{ color: #898781; font-size: 12px; }}
  #chart-wrap {{ position: relative; width: 100%; height: calc(100vh - 96px); overflow: hidden; background: #fcfcfb; cursor: grab; }}
  #chart-wrap:active {{ cursor: grabbing; }}
  #chart {{ width: 100%; height: 100%; }}
  svg {{ display: block; width: 100%; height: 100%; }}
  .node circle.node-marker {{ fill: #52514e; stroke: #898781; stroke-width: 1px; }}
  .node.leaf circle.node-marker {{ fill: #2a78d6; }}
  .node text {{ font-size: 12px; fill: #0b0b0b; }}
  .link {{ fill: none; stroke: #c3c2b7; stroke-width: 1.2px; }}
  #tooltip {{
    position: fixed;
    background: rgba(11,11,11,0.9);
    color: #fff;
    padding: 6px 10px;
    border-radius: 4px;
    font-size: 12px;
    pointer-events: none;
    opacity: 0;
    transition: opacity 0.1s;
    max-width: 320px;
    z-index: 100;
  }}
  .leaf-label {{ cursor: default; fill: #0b0b0b; }}
  .leaf-label:hover {{ fill: #2a78d6; font-weight: 600; }}
  .gene-arrow {{ cursor: default; }}
  .gene-arrow:hover polygon {{ stroke: #0b0b0b; stroke-width: 1.4px; }}
  .gene-tick:hover {{ stroke: #0b0b0b !important; }}
  /* Collapsible popup rather than a static footer -- kept out of normal
     document flow (position: fixed) so it never pushes the page taller than
     the viewport and never needs scrolling to find; #legendToggle in the
     header opens/closes it. */
  #legend {{
    display: none;
    position: fixed; right: 16px; bottom: 16px; z-index: 95;
    width: min(760px, calc(100vw - 32px)); max-height: 65vh; overflow-y: auto;
    padding: 12px 16px; background: #fff; border: 1px solid #e1e0d9; border-radius: 8px;
    box-shadow: 0 6px 24px rgba(11,11,11,0.22);
    font-size: 11px; color: #52514e;
  }}
  #legend.open {{ display: block; }}
  #legendCloseRow {{ display: flex; justify-content: flex-end; margin-bottom: 4px; }}
  #legendClose {{ font-size: 12px; padding: 2px 8px; cursor: pointer; border: 1px solid #e1e0d9; background: #f9f9f7; border-radius: 4px; }}
  #legend .legend-head {{ display: flex; align-items: center; gap: 10px; }}
  #legend .legend-head button {{ font-size: 11px; padding: 2px 8px; cursor: pointer; border: 1px solid #e1e0d9; background: #f9f9f7; border-radius: 4px; }}
  #legend .swatches {{ display: flex; flex-wrap: wrap; gap: 6px 14px; margin-top: 6px; }}
  #legend .swatch {{ display: inline-flex; align-items: center; gap: 4px; cursor: pointer; user-select: none; }}
  /* flex-shrink: 0 -- without it, a swatch/icon inside the flex-row .swatch
     or .legend-head squashes out of its aspect ratio (an oval instead of a
     circle) whenever the row runs tight on width, e.g. next to a long label
     like the cluster-break description -- the text should wrap, not the icon. */
  #legend .chip {{ width: 12px; height: 12px; border: 1px solid #333; display: inline-block; flex-shrink: 0; }}
  #legend .motility-icon {{ width: 16px; height: 16px; border-radius: 50%; vertical-align: middle; flex-shrink: 0; }}
  #legend .agreement-scale .bar {{ width: 140px; height: 10px; border: 1px solid #c3c2b7;
    background: linear-gradient(to right, #e15759, #f2e394, #59a14f); }}
</style>
</head>
<body>
<div id="controls">
  <button id="zoomIn">Zoom +</button>
  <button id="zoomOut">Zoom −</button>
  <button id="zoomFit">Fit</button>
  <button id="legendToggle">Legend</button>
  <button id="saveImage">Save Image</button>
  <label><input type="checkbox" id="includeLegendInSave" checked> Include legend in saved image</label>
  <span id="modeGroup">
    <button id="modeRect">Rectangular</button>
    <button id="modeCircular" class="active">Circular</button>
  </span>
  <label><input type="checkbox" id="showBranchLen"> Show branch lengths</label>
  <label><input type="checkbox" id="showAgreement" checked> Show agreement dot</label>
  <label><input type="checkbox" id="showMotility" checked> Show motility dot</label>
  <label title="Keep only genomes whose flagellar genes are entirely on one contig (directly assembled, or recovered by RagTag scaffolding against a real same-species reference — see the ⚡ marker). Drops every genome where gene order is still ambiguous.">
    <input type="checkbox" id="reliableOnly"> Reliable gene order only (<span id="reliableCount"></span>)
  </label>
  <label title="Draws a cluster break (|) between two consecutive genes on the same contig whenever the genomic distance between them is at least this many bp — a rough way to see which genes plausibly share one operon vs. sit at separate loci. Purely visual: doesn't affect coverage, agreement scores, or any other computed statistic. Set to 0 to disable.">
    Cluster gap ≥
    <input type="range" id="gapThreshold" min="0" max="5000" step="50" value="100" style="vertical-align: middle;">
    <span id="gapThresholdValue">100 bp</span>
  </label>
  <span class="hint">Scroll to zoom, drag to pan. Hover a leaf for its GCF accession, a gene mark for position/strand, the agreement dot for its score, or the motility dot for its BacDive phenotype.</span>
</div>
<div id="tooltip"></div>
<div id="chart-wrap"><div id="chart"></div></div>
<div id="legend">
  <div id="legendCloseRow"><button id="legendClose">Close ✕</button></div>
  <div class="legend-head">
    <b>Flagellar gene colours</b>
    <button id="genesAll">All</button>
    <button id="genesNone">None</button>
    <span class="hint">Untick a gene to hide it from every genome's track.</span>
  </div>
  <div class="swatches" id="geneSwatches"></div>
  <div class="legend-head agreement-scale" style="margin-top: 8px;">
    <b>Agreement (gene order vs. phylogeny):</b>
    <span>discordant</span><div class="bar"></div><span>concordant</span>
    <span class="swatch"><span class="chip" style="background:#bbb"></span>insufficient data</span>
  </div>
  <div class="legend-head" style="margin-top: 8px;">
    <b>Motility (BacDive curated, or Madin et al. 2020 where BacDive has no record)</b>
    <span class="swatch"><img class="motility-icon" src="{icon('motile')}">motile</span>
    <span class="swatch"><img class="motility-icon" src="{icon('nonmotile')}">non-motile</span>
    <span class="swatch"><img class="motility-icon" src="{icon('predicted_motile')}">predicted motile &gt;60% (ML, not confirmed)</span>
    <span class="swatch"><img class="motility-icon" src="{icon('predicted_nonmotile')}">predicted non-motile &gt;60% (ML, not confirmed)</span>
    <span class="swatch"><span class="chip" style="background:#bbb"></span>no record (BacDive, Madin, or confident prediction)</span>
  </div>
  <div class="legend-head" style="margin-top: 8px;">
    <b>Gene track markers</b>
    <span class="swatch">// contig break (genes on either side are on separate contigs/scaffolds — true order unknown)</span>
    <span class="swatch" style="color:#c17a2a;">| cluster break (same contig, but genomic gap ≥ the slider threshold — hover it for the exact distance)</span>
    <span class="swatch"><span class="chip" style="background:#FF0000; border-radius:50%;"></span>gene also appears elsewhere in this genome (duplicated)</span>
    <span class="swatch">⚡ gene order includes contigs joined by RagTag scaffolding against a real same-species reference (not directly assembled — hover the leaf name for details)</span>
  </div>
</div>

<script src="https://d3js.org/d3.v7.min.js"></script>
<script>
const RAW = {tree_data_js};
const GENE_COLORS = {gene_colors_js};
const MOTILITY_ICONS = {motility_icons_js};

function textColorFor(hex) {{
  // Same perceived-luminance rule as gene_colors.py's text_color_for(),
  // so a dark gene fill gets a white label instead of unreadable black-on-dark.
  const r = parseInt(hex.slice(1, 3), 16), g = parseInt(hex.slice(3, 5), 16), b = parseInt(hex.slice(5, 7), 16);
  const luminance = 0.299 * r + 0.587 * g + 0.114 * b;
  return luminance < 128 ? "white" : "black";
}}

const BLOCK_W = 22, BLOCK_H = 14, GAP = 3, ARROW_TIP = 6, SEP_W = 14;
const LEAF_SPACING = {LEAF_SPACING};
const MARGIN = {{ top: 30, right: 60, bottom: 30, left: 20 }};
// Depth levels get a fixed pixel allotment each (like a static ete3 render
// just growing the canvas for a deeper tree) rather than squeezing every
// level into one fixed total width — cramming this particular tree's 18
// levels into a fixed 260px left branch points only ~14px apart, close
// enough to look like they overlap. Computed per-render in layoutTree()
// from the tree's actual depth; this is just the minimum floor.
const LEVEL_SPACING = 24;

// Circular-mode geometry: inner radius leaves room for the root; the tree
// spans at least RADIAL_TREE_SPAN_MIN, grown as needed (see renderCircular)
// so the tightest-packed leaves are never closer than MIN_LABEL_ARC pixels
// apart — otherwise a big tree with lots of tightly nested clades can pack
// adjacent leaf labels closer than a line of italic text needs, however big
// the circle is drawn. Leaf labels, the agreement dot, and the gene-tick
// "barcode" stack outward from the tree in that order — the dot and
// tick-ring radii are computed per-render from actual label widths (see
// renderCircular), not fixed constants, since species names vary a lot in length.
const RADIAL_INNER = 40;
const RADIAL_TREE_SPAN_MIN = 260;
const MIN_LABEL_ARC = 15;
const RADIAL_LABEL_GAP = 8;
const TICK_LEN = 6, TICK_PITCH = 4;
const SEAM = 0.035; // fraction of the circle left as a gap at the seam

let showLen = false;
let showAgreement = true;
let showMotility = true;
let reliableOnly = false;
let gapThreshold = 100; // bp; genes this far apart (same contig) are drawn as separate clusters
let layoutMode = "circular";
let activeGenes = new Set();

// Recursively drops any leaf whose flagellar gene order isn't fully known
// (see synteny_analysis_full_contigs.py / build_merged_synteny.py — a
// leaf is "reliable" if its genes are entirely on one contig, either as
// directly assembled or after RagTag-recovering a real same-species
// reference; see the ⚡ marker for which). Internal nodes left with no
// surviving children are dropped too, so pruning never leaves an empty stub.
// A node left with exactly ONE surviving child (its other side was entirely
// pruned away) is collapsed by returning that child directly, rather than
// keeping a single-child pass-through wrapper — the server-side pruning this
// mirrors (ete3's Tree.prune(), see gtdb_tree.py) already does this, and
// skipping it here left ~40 spurious single-child nodes in the reliable-only
// tree, each drawing an unnecessary extra elbow joint. Rectangular mode's
// right-angle connectors make that much more visible than circular mode's
// radial lines, which is why only rectangular looked "tangled".
function filterReliable(node) {{
  if (!node.children || node.children.length === 0) {{
    return node.reliable ? node : null;
  }}
  const kept = node.children.map(filterReliable).filter(c => c !== null);
  if (kept.length === 0) return null;
  if (kept.length === 1) return kept[0];
  return {{ ...node, children: kept }};
}}

function countLeaves(node) {{
  if (!node.children || node.children.length === 0) return 1;
  return node.children.reduce((sum, c) => sum + countLeaves(c), 0);
}}

function currentTreeData() {{
  if (!reliableOnly) return RAW;
  const filtered = filterReliable(RAW);
  return filtered || RAW; // guard against filtering everything out
}}
let currentTransform = d3.zoomIdentity;
const zoomBehavior = d3.zoom().scaleExtent([0.02, 8]).on("zoom", (event) => {{
  currentTransform = event.transform;
  d3.select("#chart svg g.zoom-layer").attr("transform", event.transform);
}});

function fmtPos(bp) {{
  if (bp >= 1000000) return (bp / 1000000).toFixed(1) + "M";
  if (bp >= 1000) return Math.floor(bp / 1000) + "k";
  return String(bp);
}}

function agreementColor(v) {{
  if (v === null || v === undefined) return "#bbb";
  const t = Math.max(-1, Math.min(1, v));
  return t < 0 ? d3.interpolateRgb("#e15759", "#f2e394")(t + 1)
               : d3.interpolateRgb("#f2e394", "#59a14f")(t);
}}

const PREDICTION_CONFIDENCE_THRESHOLD = 60;

// Distinct from agreementColor's red so the two dots aren't confusable —
// motile reuses the same green (concordant) as agreement, since both mean
// "expected/positive", but non-motile is orange rather than red. Yellow and
// violet are two further, deliberately different states: BacDive has no
// curated/observed record for this species, but its genome-based ML
// prediction says motile (yellow) or non-motile (violet) with better than
// PREDICTION_CONFIDENCE_THRESHOLD% confidence — a weaker signal than the
// empirical green/orange, so both stay visually distinct from the empirical
// colours AND from each other rather than being folded into "motile"/
// "non-motile" (a different, independent model on the same pages has been
// seen to disagree with this one for at least one species, so treat this as
// suggestive, not confirmed).
function motilityIconKey(motility, predicted, madin) {{
  if (motility === "motile") return "motile";
  if (motility === "non-motile") return "nonmotile";
  // Madin et al. 2020 is treated as trust-equivalent to BacDive curated
  // (15/15 agreement where both have data, see fetch_madin_traits.py) --
  // same icon, not a separate "predicted"-style marker.
  if (madin && madin.motility === "motile") return "motile";
  if (madin && madin.motility === "non-motile") return "nonmotile";
  if (predicted && predicted.confidence > PREDICTION_CONFIDENCE_THRESHOLD) {{
    if (predicted.answer === "yes") return "predicted_motile";
    if (predicted.answer === "no") return "predicted_nonmotile";
  }}
  return null; // no data / no confident prediction -> falls back to a plain grey dot
}}

// Draws the BacDive motility icon (pixel-art flagellated/non-flagellated
// cell for the two confirmed states, same icon plus a "?" for the two
// ML-predicted-only states) instead of a flat colour dot, since MOTILITY_ICONS
// now carries the actual artwork. Falls back to a plain grey circle when
// there's no icon for this leaf's state (no data at all) or the icons
// weren't embedded (see load_motility_icons in interactive_tree.py).
// parentSel: the per-leaf selection to append into. transformAttr: passed
// straight to .attr("transform", ...) so this works identically whether the
// caller wants a simple translate (rectangular mode) or a radial one
// (circular mode).
function drawMotilityMarker(parentSel, transformAttr, size) {{
  const g = parentSel.append("g")
    .attr("transform", transformAttr)
    .on("mousemove", (event, d) => showTooltip(event, `<b>${{d.data.name}}</b><br>${{motilityTooltipText(d)}}`))
    .on("mouseleave", hideTooltip);
  g.each(function(d) {{
    const key = motilityIconKey(d.data.motility, d.data.predicted_motility, d.data.motility_madin);
    const sel = d3.select(this);
    if (key && MOTILITY_ICONS[key]) {{
      sel.append("image")
        .attr("href", MOTILITY_ICONS[key])
        .attr("x", -size / 2).attr("y", -size / 2)
        .attr("width", size).attr("height", size);
    }} else {{
      sel.append("circle").attr("r", size / 2 - 2)
        .attr("fill", "#bbb").attr("stroke", "#333").attr("stroke-width", 0.6);
    }}
  }});
}}

function motilityTooltipText(d) {{
  const m = d.data.motility;
  const flag = d.data.flagellum_arrangement;
  if (m) return `Motility (BacDive): ${{m}}${{flag ? `<br>Flagellum arrangement: ${{flag}}` : ""}}`;
  const madin = d.data.motility_madin;
  if (madin && (madin.motility === "motile" || madin.motility === "non-motile")) {{
    return `No curated BacDive record — Madin et al. 2020: ${{madin.motility}} ` +
           `(${{madin.n_strains}} strain record${{madin.n_strains === 1 ? "" : "s"}}, ` +
           `${{madin.agreement_pct}}% agreement)`;
  }}
  const pred = d.data.predicted_motility;
  if (pred && pred.confidence > PREDICTION_CONFIDENCE_THRESHOLD) {{
    if (pred.answer === "yes") {{
      return `No curated BacDive or Madin record — genome-based prediction: motile (${{pred.confidence}}% confidence)<br>` +
             `Note: an ML prediction from genome content, not an observed/confirmed result.`;
    }}
    if (pred.answer === "no") {{
      return `No curated BacDive or Madin record — genome-based prediction: non-motile (${{pred.confidence}}% confidence)<br>` +
             `Note: an ML prediction from genome content, not an observed/confirmed result.`;
    }}
  }}
  return "No BacDive or Madin et al. record for this species";
}}

// Drop genes the checklist has unticked, then collapse any contig separator
// left with nothing next to it (leading, trailing, or back-to-back) so
// hiding every gene on a contig doesn't leave a stray "//". Then, if the gap
// slider is above zero, insert an additional "cluster break" between any two
// consecutive (same-contig, still-visible) genes whose genomic distance is
// at least gapThreshold -- distinct from a "//" contig break (order unknown)
// since this is a real, measured distance on one ordered contig, just a
// candidate operon/cluster boundary rather than an assembly limitation.
// Computed on the already-checklist-filtered list, so hiding a gene here
// correctly merges its neighbours' gap rather than using stale distances.
function filterGenes(genes) {{
  if (!genes || !genes.length) return [];
  const kept = genes.filter(g => g.separator || activeGenes.has(g.gene));
  const cleaned = [];
  kept.forEach(g => {{
    if (g.separator) {{
      if (cleaned.length === 0 || cleaned[cleaned.length - 1].separator) return;
      cleaned.push(g);
    }} else {{
      cleaned.push(g);
    }}
  }});
  while (cleaned.length && cleaned[cleaned.length - 1].separator) cleaned.pop();

  if (gapThreshold <= 0) return cleaned;
  const withGaps = [];
  cleaned.forEach((g, i) => {{
    if (i > 0 && !g.separator) {{
      const prev = cleaned[i - 1];
      if (!prev.separator && prev.contig === g.contig) {{
        const gapSize = g.start - prev.end;
        if (gapSize >= gapThreshold) {{
          withGaps.push({{ separator: true, gapBreak: true, gapSize }});
        }}
      }}
    }}
    withGaps.push(g);
  }});
  return withGaps;
}}

function trackWidth(genes) {{
  let w = GAP;
  filterGenes(genes).forEach(g => {{ w += (g.separator ? SEP_W : BLOCK_W) + GAP; }});
  return w;
}}

const tooltip = document.getElementById("tooltip");
function showTooltip(event, html) {{
  tooltip.innerHTML = html;
  tooltip.style.opacity = 1;
  tooltip.style.left = (event.clientX + 14) + "px";
  tooltip.style.top = (event.clientY - 10) + "px";
}}
function hideTooltip() {{ tooltip.style.opacity = 0; }}

function leafLabelTooltip(d) {{
  let html = `<b>${{d.data.name}}</b><br>GCF: ${{d.data.gcf || "—"}}`;
  if (d.data.scaffolded) {{
    html += `<br>⚡ Gene order includes contigs joined by RagTag scaffolding ` +
            `against a real same-species reference genome. Inferred, not ` +
            `directly assembled from this genome's own reads.`;
  }}
  return html;
}}

function drawGeneTrack(sel, genes) {{
  let x = GAP;
  filterGenes(genes).forEach(g => {{
    if (g.separator) {{
      const txt = sel.append("text")
        .attr("x", x + 1).attr("y", BLOCK_H / 2 + 4)
        .attr("font-size", 11)
        .attr("fill", g.gapBreak ? "#c17a2a" : "#898781")
        .text(g.gapBreak ? "|" : "//");
      if (g.gapBreak) {{
        txt.style("cursor", "default")
           .on("mousemove", (event) => showTooltip(event,
             `Cluster break: ${{fmtPos(g.gapSize)}} gap to next gene ` +
             `(≥ ${{fmtPos(gapThreshold)}} threshold)`))
           .on("mouseleave", hideTooltip);
      }}
      x += SEP_W + GAP;
      return;
    }}
    const color = GENE_COLORS[g.gene] || "#AAAAAA";
    const x0 = x, x1 = x + BLOCK_W, midY = BLOCK_H / 2;
    const pts = g.strand === "-"
      ? [[x0 + ARROW_TIP, 0], [x1, 0], [x1, BLOCK_H], [x0 + ARROW_TIP, BLOCK_H], [x0, midY]]
      : [[x0, 0], [x1 - ARROW_TIP, 0], [x1, midY], [x1 - ARROW_TIP, BLOCK_H], [x0, BLOCK_H]];

    const grp = sel.append("g").attr("class", "gene-arrow");
    grp.append("polygon")
      .attr("points", pts.map(p => p.join(",")).join(" "))
      .attr("fill", color).attr("stroke", "#333").attr("stroke-width", 0.6);
    grp.append("text")
      .attr("x", x0 + BLOCK_W / 2).attr("y", BLOCK_H / 2 + 3)
      .attr("text-anchor", "middle").attr("font-size", 7).attr("fill", textColorFor(color))
      .text(g.gene.slice(0, 4));
    if (g.duplicated) {{
      grp.append("circle")
        .attr("cx", x1 - 4).attr("cy", 4).attr("r", 3)
        .attr("fill", "#FF0000").attr("stroke", "#AA0000");
    }}

    grp.on("mousemove", (event) => showTooltip(event,
      `<b>${{g.gene}}</b><br>${{g.contig}}:${{g.start}}-${{g.end}} (${{g.strand}})`
    )).on("mouseleave", hideTooltip);

    x += BLOCK_W + GAP;
  }});
}}

function layoutTree() {{
  // d3.cluster() (dendrogram layout) never collapses two unrelated leaves
  // onto the same y — unlike d3.tree(), which optimises internal-node
  // compactness and, on very ladder-like ("caterpillar") phylogenies like
  // this one, occasionally did exactly that. Its default separation still
  // gives sibling pairs half the gap of cousin pairs, so rescale x until the
  // tightest (sibling) gap equals LEAF_SPACING — enough room for a gene track.
  const root = d3.hierarchy(currentTreeData(), d => d.children);
  const baseHeight = root.leaves().length * LEAF_SPACING;
  const maxDepth = d3.max(root.descendants(), d => d.depth);
  const treeWidth = maxDepth * LEVEL_SPACING;
  d3.cluster().size([baseHeight, treeWidth])(root);

  // d3.cluster() positions internal nodes by how many splits remain BELOW
  // them (their remaining subtree height), not by how far they actually are
  // from the root — so a shallow branch point in one part of the tree can
  // land at the same depth as a deeply-nested one from an unrelated part,
  // making separate branches look connected. Reposition internal nodes by
  // their true depth-from-root instead (leaves stay exactly where cluster()
  // put them — flush at treeWidth — so the gene-track/dot alignment, which
  // assumes every leaf shares that same y, is unaffected).
  root.each(d => {{
    if (d.children) d.y = (d.depth / maxDepth) * treeWidth;
  }});

  const leaves = root.leaves();
  const sorted = leaves.slice().sort((a, b) => a.x - b.x);
  let minGap = Infinity;
  for (let i = 1; i < sorted.length; i++) {{
    const gap = sorted[i].x - sorted[i - 1].x;
    if (gap > 0 && gap < minGap) minGap = gap;
  }}
  if (!isFinite(minGap) || minGap <= 0) minGap = 1;

  const scale = LEAF_SPACING / minGap;
  root.each(d => {{ d.x *= scale; }});

  // d3.cluster() places each internal node at the mean of its children's x
  // — for a large, unevenly-shaped tree, two completely unrelated branch
  // points can coincidentally land within a fraction of a pixel of each
  // other this way, so their horizontal connector lines look like one
  // continuous crossing line even though nothing is actually wrong
  // topologically. Leaves are already spaced apart (above); this sweeps
  // every node (leaves included, harmlessly) and nudges any remaining
  // near-collision apart by a few px — enough to read as separate lines,
  // too small to visibly disturb the layout otherwise.
  const MIN_NODE_GAP = 3;
  const sortedAll = root.descendants().slice().sort((a, b) => a.x - b.x);
  for (let i = 1; i < sortedAll.length; i++) {{
    const prev = sortedAll[i - 1], cur = sortedAll[i];
    if (cur.x - prev.x < MIN_NODE_GAP) cur.x = prev.x + MIN_NODE_GAP;
  }}

  // The nudge above only knows about global x-order, not parent/child
  // relationships — it can push a parent's x past one of its own children's
  // (nudged because of some unrelated node elsewhere in the tree), leaving
  // the parent outside its children's span. The connecting line then has to
  // double back on itself to reach it: a visible notch/kink, previously
  // hidden under the internal-node dot markers. Restore the invariant a
  // clean cluster layout guarantees (every internal node sits at its
  // children's mean x, so it's always within their span) by recomputing
  // every internal node bottom-up from the now-nudged leaves. descendants()
  // is pre-order (parent before its children), so processing it in reverse
  // guarantees every node's children are already finalised before the node
  // itself is recomputed.
  root.descendants().slice().reverse().forEach(d => {{
    if (d.children) d.x = d3.mean(d.children, c => c.x);
  }});

  const contentHeight = d3.max(leaves, d => d.x) + LEAF_SPACING;
  return {{ root, leaves, contentHeight, treeWidth }};
}}

function render() {{
  d3.select("#chart").selectAll("*").remove();

  const wrap = document.getElementById("chart-wrap");
  const viewW = wrap.clientWidth, viewH = wrap.clientHeight;

  const svg = d3.select("#chart").append("svg")
    .attr("viewBox", `0 0 ${{viewW}} ${{viewH}}`);
  const zoomLayer = svg.append("g").attr("class", "zoom-layer")
    .attr("transform", `translate(${{MARGIN.left}},${{MARGIN.top}})`);

  if (layoutMode === "circular") {{
    renderCircular(zoomLayer);
  }} else {{
    renderRect(zoomLayer);
  }}

  svg.call(zoomBehavior).call(zoomBehavior.transform, currentTransform);
}}

function renderRect(zoomLayer) {{
  const {{ root, leaves, treeWidth }} = layoutTree();

  const geneLayer = zoomLayer.append("g").attr("class", "genes");
  const treeLayer = zoomLayer.append("g").attr("class", "tree");

  // Links (elbow connectors)
  treeLayer.selectAll(".link")
    .data(root.links())
    .join("path")
    .attr("class", "link")
    .attr("d", d => `M${{d.source.y}},${{d.source.x}} H${{d.target.y}} V${{d.target.x}}`);

  const node = treeLayer.selectAll(".node")
    .data(root.descendants())
    .join("g")
    .attr("class", d => "node" + (d.children ? "" : " leaf"))
    .attr("transform", d => `translate(${{d.y}},${{d.x}})`);

  node.append("circle").attr("class", "node-marker").attr("r", d => d.children ? 0 : 3.5);

  const leafLabel = node.filter(d => !d.children)
    .append("text")
    .attr("class", "leaf-label")
    .attr("dy", "0.31em")
    .attr("x", 8)
    .style("font-style", "italic")
    .text(d => d.data.name + (d.data.scaffolded ? " ⚡" : ""))
    .on("mousemove", (event, d) => showTooltip(event, leafLabelTooltip(d)))
    .on("mouseleave", hideTooltip);

  if (showLen) {{
    node.filter(d => !d.children)
      .append("text")
      .attr("dy", "-0.5em").attr("x", 8)
      .style("font-size", "9px").style("fill", "#898781")
      .text(d => d.data.length > 0 ? d.data.length.toFixed(4) : "");
  }}

  // Measure label extents so the gene track aligns at a single x for every leaf.
  let maxLabelRight = 0;
  leafLabel.each(function(d) {{
    const right = d.y + 8 + this.getBBox().width;
    if (right > maxLabelRight) maxLabelRight = right;
  }});
  const dotX = maxLabelRight + 18;
  const motilityDotX = dotX + 24;
  const geneTrackX = motilityDotX + 24;
  // dotX/motilityDotX are absolute x (maxLabelRight already has each leaf's
  // own d.y baked in via the measurement above) — but these circles are
  // children of `node`, which is already translated by d.y for that leaf.
  // Using them directly as "cx" would double-count d.y and land the dots
  // deep inside the gene track instead of before it, so subtract it back
  // out; every leaf shares the same d.y (treeWidth) since d3.cluster() keeps
  // all leaves flush at the same depth. Both marker slots are always
  // reserved (so the gene track doesn't jump when a toggle is flipped) —
  // only the circle itself is conditionally drawn.
  const dotCxLocal = dotX - treeWidth;
  const motilityCxLocal = motilityDotX - treeWidth;

  if (showAgreement) {{
    node.filter(d => !d.children).append("circle")
      .attr("cx", dotCxLocal).attr("cy", 0).attr("r", 5)
      .attr("fill", d => agreementColor(d.data.agreement))
      .attr("stroke", "#333").attr("stroke-width", 0.6)
      .on("mousemove", (event, d) => {{
        const a = d.data.agreement;
        const cov = d.data.coverage;
        const txt = a === null || a === undefined
          ? "Not enough shared genes for a score"
          : `Agreement: ${{a.toFixed(3)}}<br>Primary-contig coverage: ${{(cov * 100).toFixed(0)}}%`;
        showTooltip(event, `<b>${{d.data.name}}</b><br>${{txt}}`);
      }})
      .on("mouseleave", hideTooltip);
  }}

  if (showMotility) {{
    drawMotilityMarker(node.filter(d => !d.children), `translate(${{motilityCxLocal}},0)`, 13);
  }}

  // Gene-order tracks, aligned at a single x across all leaves.
  leaves.forEach(d => {{
    if (!d.data.genes || !d.data.genes.length) return;
    const g = geneLayer.append("g")
      .attr("transform", `translate(${{geneTrackX}},${{d.x - BLOCK_H / 2}})`);
    drawGeneTrack(g, d.data.genes);
  }});
}}

// ---------------------------------------------------------------------------
// Circular (radial) layout
// ---------------------------------------------------------------------------
function pointRadial(angle, radius) {{
  const a = angle - Math.PI / 2;
  return [radius * Math.cos(a), radius * Math.sin(a)];
}}

function radialElbow(source, target) {{
  const [x0, y0] = pointRadial(source.angle, source.radius);
  const [xm, ym] = pointRadial(target.angle, source.radius);
  const [x1, y1] = pointRadial(target.angle, target.radius);
  const large = Math.abs(target.angle - source.angle) > Math.PI ? 1 : 0;
  const sweep = target.angle > source.angle ? 1 : 0;
  return `M${{x0}},${{y0}} A${{source.radius}},${{source.radius}} 0 ${{large}},${{sweep}} ${{xm}},${{ym}} L${{x1}},${{y1}}`;
}}

function computeRadialTreeSpan(maxX) {{
  // layoutTree() guarantees the tightest-packed adjacent leaves are
  // LEAF_SPACING apart in the underlying (pre-angle) layout — but converting
  // that straight to an angle can still pack labels closer than a line of
  // italic text needs on a big, lopsided tree, however large the circle is.
  // Grow the radius (never shrink below the default) until that tightest
  // pair is at least MIN_LABEL_ARC pixels apart along the circle.
  const minAngularGap = (LEAF_SPACING / maxX) * (2 * Math.PI) * (1 - SEAM);
  const neededOuterRadius = MIN_LABEL_ARC / minAngularGap;
  return Math.max(RADIAL_TREE_SPAN_MIN, neededOuterRadius - RADIAL_INNER);
}}

function renderCircular(zoomLayer) {{
  const {{ root, leaves, treeWidth }} = layoutTree();
  const maxX = d3.max(leaves, d => d.x) || 1;
  const treeSpan = computeRadialTreeSpan(maxX);

  root.each(d => {{
    d.angle = (d.x / maxX) * (2 * Math.PI) * (1 - SEAM) + (Math.PI * SEAM);
    d.radius = RADIAL_INNER + (d.y / treeWidth) * treeSpan;
  }});

  // No extra offset here — fitToScreen()'s zoom transform already centres
  // the radial origin on the viewport, so this group inherits that directly.
  const g = zoomLayer.append("g").attr("class", "radial");

  const treeLayer = g.append("g").attr("class", "tree");
  const geneLayer = g.append("g").attr("class", "genes");

  treeLayer.selectAll(".link")
    .data(root.links())
    .join("path")
    .attr("class", "link")
    .attr("d", d => radialElbow(d.source, d.target));

  const node = treeLayer.selectAll(".node")
    .data(root.descendants())
    .join("g")
    .attr("class", d => "node" + (d.children ? "" : " leaf"))
    .attr("transform", d => {{ const [x, y] = pointRadial(d.angle, d.radius); return `translate(${{x}},${{y}})`; }});

  node.append("circle").attr("class", "node-marker").attr("r", d => d.children ? 0 : 3.5);

  const outerR = RADIAL_INNER + treeSpan;
  const leafSel = node.filter(d => !d.children);

  const leafLabel = leafSel.append("text")
    .attr("class", "leaf-label")
    .attr("dy", "0.31em")
    .style("font-style", "italic")
    .attr("transform", d => {{
      const [lx, ly] = pointRadial(d.angle, RADIAL_LABEL_GAP);
      const flip = d.angle > Math.PI;
      return `translate(${{lx}},${{ly}}) rotate(${{(d.angle * 180 / Math.PI) - 90 + (flip ? 180 : 0)}})`;
    }})
    .attr("text-anchor", d => (d.angle > Math.PI ? "end" : "start"))
    .text(d => d.data.name + (d.data.scaffolded ? " ⚡" : ""))
    .on("mousemove", (event, d) => showTooltip(event, leafLabelTooltip(d)))
    .on("mouseleave", hideTooltip);

  if (showLen) {{
    leafSel.append("text")
      .attr("transform", d => {{ const [x, y] = pointRadial(d.angle, -10); return `translate(${{x}},${{y}})`; }})
      .style("font-size", "9px").style("fill", "#898781")
      .text(d => d.data.length > 0 ? d.data.length.toFixed(4) : "");
  }}

  // Measure actual label widths (getBBox gives the untransformed local box,
  // i.e. the label's own pixel length) so the agreement dot and gene-tick
  // ring start past every name, however long — a fixed radius here is what
  // let long species names run through the dot/ticks before.
  let maxLabelWidth = 0;
  leafLabel.each(function() {{
    const w = this.getBBox().width;
    if (w > maxLabelWidth) maxLabelWidth = w;
  }});
  const labelClearRadius = outerR + RADIAL_LABEL_GAP + maxLabelWidth + 10;
  const dotRadius = labelClearRadius;
  const motilityRadius = dotRadius + 14;
  const ring0Radius = motilityRadius + 14;
  // Both circles are children of the leaf's own <g>, which is already
  // translated out to radius outerR — so their own transforms only need
  // the remaining distance, not the full absolute radius (that doubled up
  // before: outerR + dotRadius instead of just dotRadius). Both marker
  // slots are always reserved (so the tick ring doesn't jump when a toggle
  // is flipped) — only the circle itself is conditionally drawn.
  const dotOffset = dotRadius - outerR;
  const motilityOffset = motilityRadius - outerR;

  if (showAgreement) {{
    leafSel.append("circle")
      .attr("transform", d => {{ const [x, y] = pointRadial(d.angle, dotOffset); return `translate(${{x}},${{y}})`; }})
      .attr("r", 4)
      .attr("fill", d => agreementColor(d.data.agreement))
      .attr("stroke", "#333").attr("stroke-width", 0.6)
      .on("mousemove", (event, d) => {{
        const a = d.data.agreement;
        const cov = d.data.coverage;
        const txt = a === null || a === undefined
          ? "Not enough shared genes for a score"
          : `Agreement: ${{a.toFixed(3)}}<br>Primary-contig coverage: ${{(cov * 100).toFixed(0)}}%`;
        showTooltip(event, `<b>${{d.data.name}}</b><br>${{txt}}`);
      }})
      .on("mouseleave", hideTooltip);
  }}

  if (showMotility) {{
    drawMotilityMarker(leafSel, d => {{
      const [x, y] = pointRadial(d.angle, motilityOffset);
      return `translate(${{x}},${{y}})`;
    }}, 11);
  }}

  // Gene-order "barcode": a short radial tick per gene, stacked outward at
  // the leaf's angle, starting past the label (and agreement dot, if shown).
  // Arrows don't bend sensibly around a circle, so direction isn't shown
  // here — it's still in the hover tooltip.
  leaves.forEach(d => {{
    const genes = filterGenes(d.data.genes);
    if (!genes.length) return;
    genes.forEach((gRec, i) => {{
      const r0 = ring0Radius + i * TICK_PITCH;
      const r1 = r0 + TICK_LEN;
      if (gRec.separator) return;
      const [x0, y0] = pointRadial(d.angle, r0);
      const [x1, y1] = pointRadial(d.angle, r1);
      geneLayer.append("line")
        .attr("class", "gene-tick")
        .attr("x1", x0).attr("y1", y0).attr("x2", x1).attr("y2", y1)
        .attr("stroke", GENE_COLORS[gRec.gene] || "#AAAAAA")
        .attr("stroke-width", TICK_PITCH * 0.85)
        .on("mousemove", (event) => showTooltip(event,
          `<b>${{gRec.gene}}</b><br>${{gRec.contig}}:${{gRec.start}}-${{gRec.end}} (${{gRec.strand}})`
        ))
        .on("mouseleave", hideTooltip);
    }});
  }});
}}

function fitToScreen() {{
  const wrap = document.getElementById("chart-wrap");
  if (layoutMode === "circular") {{
    // Circular mode benefits from seeing the whole ring up front (unlike the
    // tall rectangular list, where "fit everything" just shrinks labels and
    // genes into illegible slivers) — estimate the full radius including
    // labels/gene ticks and scale it to fit the viewport.
    const {{ leaves }} = layoutTree();
    const maxX = d3.max(leaves, d => d.x) || 1;
    const treeSpan = computeRadialTreeSpan(maxX);
    const maxGenes = d3.max(leaves, d => filterGenes(d.data.genes).length) || 0;
    const maxNameLen = d3.max(leaves, d => (d.data.name || "").length) || 20;
    // No rendered DOM to measure yet, so estimate label width from character
    // count (~6.5px/char for italic 12px) rather than a flat constant —
    // refined automatically once renderCircular() actually measures it.
    const estLabelWidth = maxNameLen * 6.5;
    const fullRadius = RADIAL_INNER + treeSpan + RADIAL_LABEL_GAP + estLabelWidth + 38 + maxGenes * TICK_PITCH + 40;
    const scale = Math.min(wrap.clientWidth, wrap.clientHeight) * 0.48 / fullRadius;
    currentTransform = d3.zoomIdentity.translate(wrap.clientWidth / 2, wrap.clientHeight / 2).scale(Math.max(0.05, scale));
  }} else {{
    // The tree is far taller than any screen once ~150+ leaves each carry a
    // gene track, so "fit" resets to a readable scale at the top of the tree
    // rather than cramming every leaf into view — scrolling/dragging reaches
    // the rest.
    currentTransform = d3.zoomIdentity.translate(10, 10).scale(1);
  }}
  render();
}}

function buildLegend() {{
  const used = new Set();
  (function walk(n) {{
    if (n.genes) n.genes.forEach(g => {{ if (!g.separator) used.add(g.gene); }});
    (n.children || []).forEach(walk);
  }})(RAW);
  const sortedGenes = Array.from(used).sort();
  activeGenes = new Set(sortedGenes);

  const container = d3.select("#geneSwatches");
  sortedGenes.forEach(gene => {{
    const s = container.append("label").attr("class", "swatch");
    s.append("input").attr("type", "checkbox").property("checked", true)
      .on("change", function() {{
        if (this.checked) activeGenes.add(gene); else activeGenes.delete(gene);
        render();
      }});
    s.append("span").attr("class", "chip").style("background", GENE_COLORS[gene] || "#AAAAAA");
    s.append("span").text(gene);
  }});

  document.getElementById("genesAll").addEventListener("click", () => {{
    activeGenes = new Set(sortedGenes);
    container.selectAll("input").property("checked", true);
    render();
  }});
  document.getElementById("genesNone").addEventListener("click", () => {{
    activeGenes = new Set();
    container.selectAll("input").property("checked", false);
    render();
  }});
}}

function downloadBlob(blob, filename) {{
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}}

// A <style> tag embedded in a standalone SVG isn't reliably applied once
// that SVG is loaded via img.src (e.g. .link paths came out filled solid
// black instead of fill:none/stroke-only) — so instead of relying on a
// stylesheet surviving that round-trip, copy each element's actual
// COMPUTED style inline. Slower than embedding a stylesheet, but it's
// reading the literal values the browser is really rendering, so the
// export is guaranteed to match what's on screen.
const INLINE_STYLE_PROPS = [
  "fill", "stroke", "stroke-width", "stroke-opacity", "fill-opacity",
  "font-family", "font-size", "font-style", "font-weight", "text-anchor", "opacity",
];
function inlineComputedStyles(sourceRoot, cloneRoot) {{
  function copy(srcEl, cloneEl) {{
    const cs = getComputedStyle(srcEl);
    let styleStr = "";
    INLINE_STYLE_PROPS.forEach(p => {{
      const v = cs.getPropertyValue(p);
      if (v) styleStr += `${{p}}:${{v}};`;
    }});
    if (styleStr) cloneEl.setAttribute("style", styleStr);
  }}
  copy(sourceRoot, cloneRoot);
  const srcEls = sourceRoot.querySelectorAll("*");
  const cloneEls = cloneRoot.querySelectorAll("*");
  for (let i = 0; i < srcEls.length; i++) copy(srcEls[i], cloneEls[i]);
}}

// Places a row of small items (colour swatches, icons, text) left-to-right,
// wrapping to a new row when the next one would run past maxWidth. Used to
// lay out the exported-image legend without hand-picking coordinates for
// every gene/icon — each item supplies its own width estimate and a
// draw(g, x, y) callback. Returns the y just past the last row, so callers
// can stack sections underneath each other.
function flowItems(g, items, x0, y0, maxWidth, rowGap) {{
  let x = x0, y = y0, rowH = 0;
  items.forEach(item => {{
    if (x + item.w > x0 + maxWidth && x > x0) {{
      x = x0; y += rowH + rowGap; rowH = 0;
    }}
    item.draw(g, x, y);
    x += item.w + 14;
    rowH = Math.max(rowH, item.h);
  }});
  return y + rowH;
}}

// Recreates the on-page HTML legend as SVG so it can be rasterised into the
// same exported PNG as the tree (the live legend is plain HTML/CSS, which a
// canvas export can't capture directly). Not a pixel-identical copy of the
// popup's layout — this is a fresh, export-sized rendering of the same
// content. Returns the y coordinate just past the bottom of what it drew.
function buildLegendSVG(g, x0, y0, width) {{
  let y = y0;
  const heading = (text) => {{
    g.append("text").attr("x", x0).attr("y", y).attr("font-weight", 700).attr("font-size", 13).attr("fill", "#0b0b0b").text(text);
    y += 16;
  }};

  heading("Flagellar gene colours");
  const geneItems = Object.keys(GENE_COLORS).sort().map(gene => ({{
    w: 58, h: 14,
    draw: (gg, x, yy) => {{
      gg.append("rect").attr("x", x).attr("y", yy).attr("width", 12).attr("height", 12)
        .attr("fill", GENE_COLORS[gene]).attr("stroke", "#333").attr("stroke-width", 0.6);
      gg.append("text").attr("x", x + 16).attr("y", yy + 10).attr("font-size", 10).attr("fill", "#0b0b0b").text(gene);
    }},
  }}));
  y = flowItems(g, geneItems, x0, y, width, 6) + 22;

  heading("Agreement (gene order vs. phylogeny)");
  const gradId = "exportAgreementGrad";
  const grad = g.append("defs").append("linearGradient").attr("id", gradId);
  [["0%", "#e15759"], ["50%", "#f2e394"], ["100%", "#59a14f"]].forEach(([offset, stop]) =>
    grad.append("stop").attr("offset", offset).attr("stop-color", stop));
  g.append("rect").attr("x", x0).attr("y", y).attr("width", 140).attr("height", 10)
    .attr("fill", `url(#${{gradId}})`).attr("stroke", "#c3c2b7");
  g.append("text").attr("x", x0 - 6).attr("y", y + 9).attr("text-anchor", "end").attr("font-size", 10).attr("fill", "#52514e").text("discordant");
  g.append("text").attr("x", x0 + 146).attr("y", y + 9).attr("font-size", 10).attr("fill", "#52514e").text("concordant");
  y += 28;

  heading("Motility (BacDive)");
  const motilityItems = [
    ["motile", "motile"],
    ["nonmotile", "non-motile"],
    ["predicted_motile", "predicted motile >60% (ML, not confirmed)"],
    ["predicted_nonmotile", "predicted non-motile >60% (ML, not confirmed)"],
  ].map(([key, label]) => ({{
    w: label.length * 5.6 + 26, h: 18,
    draw: (gg, x, yy) => {{
      if (MOTILITY_ICONS[key]) {{
        gg.append("image").attr("href", MOTILITY_ICONS[key]).attr("x", x).attr("y", yy).attr("width", 16).attr("height", 16);
      }}
      gg.append("text").attr("x", x + 20).attr("y", yy + 12).attr("font-size", 10).attr("fill", "#0b0b0b").text(label);
    }},
  }}));
  y = flowItems(g, motilityItems, x0, y, width, 8) + 22;

  heading("Gene track markers");
  g.append("text").attr("x", x0).attr("y", y).attr("font-size", 10).attr("fill", "#0b0b0b")
    .text("// contig break (order unknown)     | cluster break (gap ≥ slider threshold)     " +
          "⚡ RagTag-scaffolded gene order     ● red = duplicated gene");
  y += 6;

  return y;
}}

function saveTreeImage() {{
  // Exports the full tree (every leaf, at its full rendered extent) rather
  // than just whatever's currently panned/zoomed into view — getBBox() on
  // the zoom-layer reads its own local (pre-transform) extent, so this is
  // unaffected by the current pan/zoom state.
  const zoomLayer = document.querySelector("#chart svg g.zoom-layer");
  const bbox = zoomLayer.getBBox();
  const pad = 20;
  const w = bbox.width + pad * 2;
  const treeH = bbox.height + pad * 2;

  const svgNS = "http://www.w3.org/2000/svg";
  const clone = zoomLayer.cloneNode(true);
  inlineComputedStyles(zoomLayer, clone);
  clone.setAttribute("transform", `translate(${{pad - bbox.x}},${{pad - bbox.y}})`);

  const includeLegend = document.getElementById("includeLegendInSave").checked;

  const exportSvg = document.createElementNS(svgNS, "svg");
  exportSvg.setAttribute("xmlns", svgNS);
  const bg = document.createElementNS(svgNS, "rect");
  exportSvg.appendChild(bg);
  exportSvg.appendChild(clone);

  let h = treeH;
  if (includeLegend) {{
    const legendPad = 16;
    const legendG = d3.select(exportSvg).append("g")
      .attr("transform", `translate(${{legendPad}},${{treeH + legendPad}})`);
    legendG.append("line").attr("x1", 0).attr("x2", w - legendPad * 2).attr("y1", -8).attr("y2", -8)
      .attr("stroke", "#e1e0d9");
    const legendBottom = buildLegendSVG(legendG, 0, 12, w - legendPad * 2);
    h = treeH + legendPad + legendBottom + legendPad;
  }}
  exportSvg.setAttribute("width", w);
  exportSvg.setAttribute("height", h);
  bg.setAttribute("width", w); bg.setAttribute("height", h); bg.setAttribute("fill", "#fcfcfb");

  const svgString = new XMLSerializer().serializeToString(exportSvg);
  const svgUrl = URL.createObjectURL(new Blob([svgString], {{ type: "image/svg+xml;charset=utf-8" }}));

  const scale = 2;
  const img = new Image();
  img.onload = () => {{
    const canvas = document.createElement("canvas");
    canvas.width = w * scale; canvas.height = h * scale;
    const ctx = canvas.getContext("2d");
    ctx.scale(scale, scale);
    ctx.fillStyle = "#fcfcfb";
    ctx.fillRect(0, 0, w, h);
    ctx.drawImage(img, 0, 0, w, h);
    URL.revokeObjectURL(svgUrl);
    canvas.toBlob(blob => downloadBlob(blob, `interactive_tree_${{layoutMode}}.png`), "image/png");
  }};
  img.onerror = (e) => {{ console.error("Save image failed:", e); URL.revokeObjectURL(svgUrl); }};
  img.src = svgUrl;
}}

buildLegend();
fitToScreen();
document.getElementById("reliableCount").textContent =
  `${{countLeaves(filterReliable(RAW) || {{children: []}})}} / ${{countLeaves(RAW)}}`;

window.addEventListener("resize", () => render());

document.getElementById("zoomIn").addEventListener("click", () => {{
  d3.select("#chart svg").transition().call(zoomBehavior.scaleBy, 1.3);
}});
document.getElementById("zoomOut").addEventListener("click", () => {{
  d3.select("#chart svg").transition().call(zoomBehavior.scaleBy, 1 / 1.3);
}});
document.getElementById("zoomFit").addEventListener("click", fitToScreen);
document.getElementById("saveImage").addEventListener("click", saveTreeImage);

document.getElementById("legendToggle").addEventListener("click", () => {{
  document.getElementById("legend").classList.toggle("open");
}});
document.getElementById("legendClose").addEventListener("click", () => {{
  document.getElementById("legend").classList.remove("open");
}});
document.addEventListener("click", (event) => {{
  const legend = document.getElementById("legend");
  const toggle = document.getElementById("legendToggle");
  if (legend.classList.contains("open") && !legend.contains(event.target) && event.target !== toggle) {{
    legend.classList.remove("open");
  }}
}});

document.getElementById("modeRect").addEventListener("click", () => {{
  layoutMode = "rect";
  document.getElementById("modeRect").classList.add("active");
  document.getElementById("modeCircular").classList.remove("active");
  fitToScreen();
}});
document.getElementById("modeCircular").addEventListener("click", () => {{
  layoutMode = "circular";
  document.getElementById("modeCircular").classList.add("active");
  document.getElementById("modeRect").classList.remove("active");
  fitToScreen();
}});

document.getElementById("showBranchLen").addEventListener("change", e => {{
  showLen = e.target.checked;
  render();
}});
document.getElementById("showAgreement").addEventListener("change", e => {{
  showAgreement = e.target.checked;
  render();
}});
document.getElementById("showMotility").addEventListener("change", e => {{
  showMotility = e.target.checked;
  render();
}});
document.getElementById("reliableOnly").addEventListener("change", e => {{
  reliableOnly = e.target.checked;
  fitToScreen(); // leaf count changes a lot (155 -> 108) so re-fit, not just re-render
}});
document.getElementById("gapThreshold").addEventListener("input", e => {{
  gapThreshold = +e.target.value;
  document.getElementById("gapThresholdValue").textContent =
    gapThreshold === 0 ? "off" : `${{gapThreshold}} bp`;
  render(); // gene tracks change width when clusters split, but leaf count doesn't -> no re-fit needed
}});
</script>
</body>
</html>
"""
    output_path.write_text(html)
    print(f"Saved: {output_path}")


def main():
    named_nwk = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("images_and_newick/thesis_b/gtdb_ref_tree.nwk")
    og_nwk = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("images_and_newick/thesis_b/gtdb_ref_tree_og.nwk")
    output = Path(sys.argv[3]) if len(sys.argv) > 3 else Path("frontend/interactive_tree.html")
    gene_order_json = Path(sys.argv[4]) if len(sys.argv) > 4 else Path("gene_order/gene_order_combined.json")
    synteny_json = Path(sys.argv[5]) if len(sys.argv) > 5 else Path("gene_order/synteny_analysis_combined_merged.json")
    motility_json = Path(sys.argv[6]) if len(sys.argv) > 6 else Path("gene_order/motility.json")

    for p in (named_nwk, og_nwk):
        if not p.exists():
            print(f"Not found: {p}")
            sys.exit(1)

    print(f"Loading {named_nwk} ...")
    t_named = Tree(str(named_nwk), format=1, quoted_node_names=True)
    print(f"Loading {og_nwk} ...")
    t_og = Tree(str(og_nwk), format=1, quoted_node_names=True)

    gcf_map = build_gcf_map(t_named, t_og)
    leaf_count = len(t_named.get_leaves())
    print(f"{leaf_count} leaves, {len(gcf_map)} GCF mappings built")

    gene_order = load_gene_order(str(gene_order_json))
    synteny = load_synteny(str(synteny_json))
    motility = load_motility(str(motility_json))

    all_genes = {g["gene"].lower() for recs in gene_order.values() for g in recs}
    GENE_COLORS.update(generate_gene_colors(all_genes))
    print(f"Generated {len(GENE_COLORS)} distinct gene colours (min-distance randomised)")

    tree_json = tree_to_dict(t_named, gcf_map, gene_order, synteny, motility)
    generate_html(tree_json, leaf_count, output)
    print(f"Open in browser: open {output}")


if __name__ == "__main__":
    main()
