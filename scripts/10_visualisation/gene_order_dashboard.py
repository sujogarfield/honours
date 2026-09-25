#!/usr/bin/env python3
"""
gene_order_dashboard.py
Quantitative dashboard for the flagellar gene-order vs. phylogeny analysis:
a synteny-distance vs. phylogenetic-distance scatter (visualises the Mantel
relationship directly), a tree-ordered pairwise distance heatmap, and a
ranked diverging bar chart of per-organism agreement scores.

Run AFTER synteny_analysis.py (reads gene_order/synteny_analysis.json) and
gtdb_tree.py (reads images_and_newick/thesis_b/gtdb_ref_tree.nwk for leaf/genus order).
"""

import json
import sys
from pathlib import Path
from ete3 import Tree

# Usage: gene_order_dashboard.py [synteny_json] [synteny_reliable_json] [output_html]
DEFAULT_SYNTENY_JSON = "gene_order/synteny_analysis_combined.json"
DEFAULT_SYNTENY_RELIABLE_JSON = "gene_order/synteny_analysis_combined_reliable.json"
DEFAULT_MOTILITY_JSON = "gene_order/motility_clustering.json"
NAMED_NWK = "images_and_newick/thesis_b/gtdb_ref_tree.nwk"
DEFAULT_OUTPUT_HTML = "frontend/synteny_dashboard.html"

# Same genus -> clade colour palette as interactive_tree.py / helpers.colour_tree_by_species.
GENUS_COLORS = {
    "Helicobacter": "#ffd2d0",
    "Campylobacter": "#ddffdd",
    "Sulfurimonas": "#fffad1",
    "Wolinella": "#edd9ff",
    "Hippea": "#d0e9ff",
    "Sulfuricurvum": "#ecd5ff",
    "Candidatus": "#d7ffeb",
    "Sulfurospirillum": "#feffc9",
    "Hydrogenimonas": "#d6fff8",
    "Nitrosophilus": "#e2dfff",
    "Desulfurella": "#dbf4ff",
    "Caminibacter": "#ffdefd",
    "Nautilia": "#ffdbdb",
    "Lebetimonas": "#e4ffdf",
    "Nitratiruptor": "#faffd8",
    "Sulfurovum": "#fff4da",
    "Nitratifractor": "#ffdbdb",
    "Arcobacter": "#deffe5",
    "Aliarcobacter": "#e5e8ff",
    "Malaciobacter": "#ffdcf9",
    "Halarcobacter": "#e6fcff",
    "Poseidonibacter": "#ffe1fc",
}


def tree_to_simple_dict(node) -> dict:
    """Bare {name, children} shape — just enough topology for the tanglegram,
    no gene/GCF payload needed on this side."""
    d = {"name": node.name or ""}
    if node.children:
        d["children"] = [tree_to_simple_dict(c) for c in node.children]
    return d


def generate_html(
    data_all: dict, data_reliable: dict,
    tree_order_all: list[str], tree_order_reliable: list[str],
    phylo_subtree_all: dict, phylo_subtree_reliable: dict,
    motility_data=None,
    tanglegram_data=None,
) -> str:
    data_all_js = json.dumps(data_all)
    data_reliable_js = json.dumps(data_reliable)
    tree_order_all_js = json.dumps(tree_order_all)
    tree_order_reliable_js = json.dumps(tree_order_reliable)
    genus_colors_js = json.dumps(GENUS_COLORS)
    phylo_subtree_all_js = json.dumps(phylo_subtree_all)
    phylo_subtree_reliable_js = json.dumps(phylo_subtree_reliable)
    motility_data_js = json.dumps(motility_data or {})
    tanglegram_data_js = json.dumps(tanglegram_data or {})
    n_all = len(data_all["matrix"]["organisms"])
    n_reliable = len(data_reliable["matrix"]["organisms"])

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Synteny Dashboard — Gene Order vs. Phylogeny</title>
<style>
  html, body {{ margin: 0; height: 100%; font-family: sans-serif; background: #fafafa; color: #0b0b0b; }}
  #tabs {{ display: flex; gap: 4px; padding: 8px 16px 0; background: #fff; border-bottom: 1px solid #e1e0d9; align-items: flex-start; }}
  #saveImage {{
    margin-left: auto; font-size: 13px; padding: 8px 16px; cursor: pointer;
    border: 1px solid #e1e0d9; background: #f9f9f7; border-radius: 6px; color: #0b0b0b;
  }}
  #tabs button {{
    font-size: 13px; padding: 8px 16px; cursor: pointer; border: 1px solid #e1e0d9;
    border-bottom: none; background: #f9f9f7; border-radius: 6px 6px 0 0; color: #52514e;
  }}
  #tabs button.active {{ background: #fff; color: #0b0b0b; font-weight: 600; }}
  .panel {{ display: none; padding: 16px; }}
  .panel.active {{ display: block; }}
  .hint {{ color: #898781; font-size: 12px; margin-bottom: 8px; }}
  #tooltip {{
    position: fixed; background: rgba(11,11,11,0.9); color: #fff; padding: 6px 10px;
    border-radius: 4px; font-size: 12px; pointer-events: none; opacity: 0;
    transition: opacity 0.1s; max-width: 320px; z-index: 100;
  }}
  .axis-label {{ font-size: 11px; fill: #52514e; }}
  .axis-line {{ stroke: #c3c2b7; stroke-width: 1px; }}
  .grid-line {{ stroke: #e1e0d9; stroke-width: 1px; }}
  .stat-box {{ font-size: 13px; fill: #0b0b0b; }}
  .legend-tick {{ font-size: 10px; fill: #52514e; }}
  #rankWrap, #tangleWrap {{ height: calc(100vh - 160px); overflow-y: auto; border: 1px solid #e1e0d9; border-radius: 6px; }}
  .rank-name {{ font-size: 11px; font-style: italic; fill: #0b0b0b; }}
  .rank-null {{ font-size: 11px; font-style: italic; fill: #898781; }}
  #heatmapWrap, #scatterWrap, #mdsWrap, #motilityWrap {{ position: relative; }}
  #motilityWrap {{ overflow-x: auto; border: 1px solid #e1e0d9; border-radius: 6px; padding: 8px; }}
  .motility-gene-label {{ font-size: 10px; fill: #0b0b0b; cursor: default; }}
  .motility-group-label {{ font-size: 11px; font-weight: 600; }}
  .motility-toggle {{ font-size: 12px; color: #52514e; margin-bottom: 8px; display: flex; align-items: center; gap: 6px; }}
  #hgtWrap {{ overflow: auto; max-height: calc(100vh - 220px); border: 1px solid #e1e0d9; border-radius: 6px; }}
  #dotPlotModal {{
    display: none; position: fixed; inset: 0; background: rgba(11,11,11,0.45);
    z-index: 200; align-items: center; justify-content: center;
  }}
  #dotPlotModal.open {{ display: flex; }}
  #dotPlotCard {{ background: #fff; border-radius: 8px; padding: 16px; max-width: 90vw; max-height: 90vh; overflow: auto; }}
  #dotPlotHead {{ display: flex; justify-content: space-between; align-items: center; margin-bottom: 4px; gap: 16px; }}
  #dotPlotTitle {{ font-size: 14px; }}
  #dotPlotClose {{ cursor: pointer; border: 1px solid #e1e0d9; background: #f9f9f7; border-radius: 4px; padding: 4px 12px; font-size: 12px; }}
</style>
</head>
<body>
<div id="tabs">
  <button class="tab-btn active" data-tab="scatter">Scatter (synteny vs. phylogeny)</button>
  <button class="tab-btn" data-tab="heatmap">Heatmap (tree-ordered distances)</button>
  <button class="tab-btn" data-tab="ranking">Ranking (agreement score)</button>
  <button class="tab-btn" data-tab="tanglegram">Tanglegram (phylogeny vs. gene order)</button>
  <button class="tab-btn" data-tab="mds">MDS (gene order only)</button>
  <button class="tab-btn" data-tab="motility">Gene content vs. motility</button>
  <button class="tab-btn" data-tab="hgt">HGT tanglegram (gene tree vs. species tree)</button>
  <label id="reliableToggleLabel" style="margin-left: 16px; font-size: 13px; color: #52514e; align-self: center;"
    title="Switch every panel to the genomes whose flagellar genes are entirely on one contig (directly assembled, or recovered by RagTag scaffolding against a real same-species reference). Drops every genome where gene order is still ambiguous. This is a precomputed dataset, not a client-side filter — the Mantel test, UPGMA tree, and MDS coordinates are all recomputed on just this set, not just visually hidden.">
    <input type="checkbox" id="reliableOnly"> Reliable gene order only ({n_reliable} / {n_all})
  </label>
  <button id="saveImage">Save Image</button>
</div>
<div id="tooltip"></div>

<div id="panel-scatter" class="panel active">
  <div class="hint">Each point is one genome pair. Hover for organism names and distances.</div>
  <div id="scatterWrap"><canvas id="scatterCanvas"></canvas><svg id="scatterSvg"></svg></div>
</div>

<div id="panel-heatmap" class="panel">
  <div class="hint">Rows/columns follow tree order, so conserved clades show as light blocks along the diagonal. Grey = no gene-order data. Click a cell to see the two genomes' gene-order dot plot.</div>
  <div id="heatmapWrap"><canvas id="heatmapCanvas"></canvas><svg id="heatmapSvg"></svg></div>
</div>

<div id="panel-ranking" class="panel">
  <div class="hint">Agreement score: how well each genome's gene-order neighbours match its phylogenetic neighbours (-1 to +1). Sorted most concordant → most discordant.</div>
  <div id="rankWrap"><svg id="rankSvg"></svg></div>
</div>

<div id="panel-tanglegram" class="panel">
  <div class="hint">Left: GTDB phylogeny. Right: dendrogram built purely from gene order (UPGMA on synteny distance, no phylogenetic information). A straight, uncrossed connector means a genome's gene-order neighbours match its phylogenetic ones; long diagonal/crossing connectors mark discordant genomes. Hover a connector for the organism name.</div>
  <div id="tangleWrap"><svg id="tangleSvg"></svg></div>
</div>

<div id="panel-mds" class="panel">
  <div class="hint">Classical MDS (PCoA) of the synteny distance matrix — organisms placed using gene order alone, no phylogenetic information. Colored by genus: clusters here mean gene order alone recovers that group. Missing pairwise distances were imputed with each organism's average, so treat exact positions as approximate.</div>
  <div id="mdsWrap"><svg id="mdsSvg"></svg></div>
</div>

<div id="panel-motility" class="panel">
  <div class="hint" id="motilityHint"></div>
  <label class="motility-toggle">
    <input type="checkbox" id="includePredicted">
    Include ML-predicted motility calls as a fallback (genomes with no curated BacDive or Madin et al. record) — off by default, since the prediction model is a secondary, sometimes-disagreeing signal
  </label>
  <div id="motilityWrap"><svg id="motilitySvg"></svg></div>
</div>

<div id="panel-hgt" class="panel">
  <div class="hint">Left: GTDB phylogeny (pruned to genomes with a member of the selected orthogroup). Right: that orthogroup's own OrthoFinder gene tree (built from the actual protein sequences, not gene order). A straight, uncrossed connector means a genome's placement by this gene agrees with its placement genome-wide; long diagonal/crossing connectors mean the gene's evolutionary history disagrees with the species tree — the classic horizontal-gene-transfer signature. Genomes with more than one copy (paralogs) get one connector per copy, each labelled by its protein ID. Each leaf has a genus-coloured dot (same palette as the rest of this dashboard) — a dot sitting inside a differently-coloured cluster on the gene-tree side is the visual tell for HGT, faster to spot than reading every label.</div>
  <div style="margin-bottom:10px;">
    <label style="font-size:13px; color:#52514e;">Gene / orthogroup:
      <select id="hgtSelect"></select>
    </label>
    <span id="hgtMeta" style="font-size:12px; color:#898781; margin-left:12px;"></span>
  </div>
  <div id="hgtWrap"><svg id="hgtSvg"></svg></div>
</div>

<div id="dotPlotModal">
  <div id="dotPlotCard">
    <div id="dotPlotHead">
      <b id="dotPlotTitle"></b>
      <button id="dotPlotClose">Close</button>
    </div>
    <div class="hint">Each point is a gene shared by both genomes. A diagonal line means fully conserved order.</div>
    <svg id="dotPlotSvg"></svg>
  </div>
</div>

<script src="https://d3js.org/d3.v7.min.js"></script>
<script>
const DATA_ALL = {data_all_js};
const DATA_RELIABLE = {data_reliable_js};
const TREE_ORDER_ALL = {tree_order_all_js};
const TREE_ORDER_RELIABLE = {tree_order_reliable_js};
const PHYLO_SUBTREE_ALL = {phylo_subtree_all_js};
const PHYLO_SUBTREE_RELIABLE = {phylo_subtree_reliable_js};
const GENUS_COLORS = {genus_colors_js};
const MOTILITY_DATA = {motility_data_js};
const TANGLEGRAM_DATA = {tanglegram_data_js};

// Mutable "active dataset" pointers -- every render function below reads
// these, never DATA_ALL/DATA_RELIABLE directly, so toggling just swaps
// what they point at and re-renders. Each side is a fully separate
// precomputed result (its own Mantel r/p, UPGMA tree, MDS coordinates) from
// synteny_analysis_full_contigs.py, not a client-side filter of DATA_ALL --
// dropping fragmented genomes changes those statistics, so recomputing them
// properly server-side is what makes this a meaningful toggle rather than
// just hiding points.
let DATA = DATA_ALL;
let TREE_ORDER = TREE_ORDER_ALL;
let PHYLO_SUBTREE = PHYLO_SUBTREE_ALL;

const tooltip = document.getElementById("tooltip");
function showTooltip(event, html) {{
  tooltip.innerHTML = html;
  tooltip.style.opacity = 1;
  tooltip.style.left = (event.clientX + 14) + "px";
  tooltip.style.top = (event.clientY - 10) + "px";
}}
function hideTooltip() {{ tooltip.style.opacity = 0; }}

function genusOf(name) {{ return (name || "").split(" ")[0]; }}

// Same red -> yellow -> green diverging scale as interactive_tree.py's
// agreement dot, reused here for visual consistency across the project.
function agreementColor(v) {{
  if (v === null || v === undefined) return "#bbb";
  const t = Math.max(-1, Math.min(1, v));
  return t < 0 ? d3.interpolateRgb("#e15759", "#f2e394")(t + 1)
               : d3.interpolateRgb("#f2e394", "#59a14f")(t);
}}

// ---------------------------------------------------------------------------
// Tabs
// ---------------------------------------------------------------------------
document.querySelectorAll(".tab-btn").forEach(btn => {{
  btn.addEventListener("click", () => {{
    document.querySelectorAll(".tab-btn").forEach(b => b.classList.remove("active"));
    document.querySelectorAll(".panel").forEach(p => p.classList.remove("active"));
    btn.classList.add("active");
    document.getElementById("panel-" + btn.dataset.tab).classList.add("active");
  }});
}});

// ---------------------------------------------------------------------------
// Scatter: synteny distance (y) vs. phylogenetic distance (x), one dot per pair
// ---------------------------------------------------------------------------
function renderScatter() {{
  const names = DATA.matrix.organisms;
  const synt = DATA.matrix.synteny_distance;
  const phylo = DATA.matrix.phylogenetic_distance;
  const points = [];
  for (let i = 0; i < names.length; i++) {{
    for (let j = i + 1; j < names.length; j++) {{
      if (synt[i][j] === null || synt[i][j] === undefined) continue;
      points.push({{ a: names[i], b: names[j], synt: synt[i][j], phylo: phylo[i][j] }});
    }}
  }}

  const wrap = document.getElementById("scatterWrap");
  const width = Math.max(700, window.innerWidth - 60);
  const height = Math.max(500, window.innerHeight - 200);
  const margin = {{ top: 30, right: 30, bottom: 50, left: 60 }};

  const canvas = document.getElementById("scatterCanvas");
  const dpr = window.devicePixelRatio || 1;
  canvas.width = width * dpr;
  canvas.height = height * dpr;
  canvas.style.width = width + "px";
  canvas.style.height = height + "px";
  canvas.style.position = "absolute";
  const ctx = canvas.getContext("2d");
  ctx.scale(dpr, dpr);

  const svg = d3.select("#scatterSvg")
    .attr("width", width).attr("height", height)
    .style("position", "absolute").style("pointer-events", "none");
  svg.selectAll("*").remove();

  const maxPhylo = d3.max(points, d => d.phylo) || 1;
  const x = d3.scaleLinear().domain([0, maxPhylo]).range([margin.left, width - margin.right]);
  const y = d3.scaleLinear().domain([0, 1]).range([height - margin.bottom, margin.top]);

  // Gridlines + axes (SVG, crisp chrome)
  const xTicks = x.ticks(6), yTicks = y.ticks(6);
  svg.selectAll(".grid-x").data(xTicks).join("line").attr("class", "grid-line")
    .attr("x1", d => x(d)).attr("x2", d => x(d)).attr("y1", margin.top).attr("y2", height - margin.bottom);
  svg.selectAll(".grid-y").data(yTicks).join("line").attr("class", "grid-line")
    .attr("x1", margin.left).attr("x2", width - margin.right).attr("y1", d => y(d)).attr("y2", d => y(d));
  svg.selectAll(".tick-x").data(xTicks).join("text").attr("class", "axis-label")
    .attr("x", d => x(d)).attr("y", height - margin.bottom + 16).attr("text-anchor", "middle").text(d => d.toFixed(2));
  svg.selectAll(".tick-y").data(yTicks).join("text").attr("class", "axis-label")
    .attr("x", margin.left - 8).attr("y", d => y(d) + 3).attr("text-anchor", "end").text(d => d.toFixed(1));
  svg.append("text").attr("class", "axis-label").attr("x", width / 2).attr("y", height - 8)
    .attr("text-anchor", "middle").text("Phylogenetic (patristic) distance");
  svg.append("text").attr("class", "axis-label")
    .attr("transform", `translate(14,${{height / 2}}) rotate(-90)`).attr("text-anchor", "middle")
    .text("Synteny (breakpoint) distance");

  // Points on canvas (~9k+ — cheaper than one DOM node each)
  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = "rgba(42,120,214,0.35)";
  points.forEach(p => {{
    ctx.beginPath();
    ctx.arc(x(p.phylo), y(p.synt), 2.6, 0, 2 * Math.PI);
    ctx.fill();
  }});

  // Least-squares trend line
  const n = points.length;
  const sx = d3.sum(points, d => d.phylo), sy = d3.sum(points, d => d.synt);
  const sxx = d3.sum(points, d => d.phylo * d.phylo), sxy = d3.sum(points, d => d.phylo * d.synt);
  const slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
  const intercept = (sy - slope * sx) / n;
  svg.append("line")
    .attr("x1", x(0)).attr("y1", y(intercept))
    .attr("x2", x(maxPhylo)).attr("y2", y(intercept + slope * maxPhylo))
    .attr("stroke", "#eb6834").attr("stroke-width", 2);

  const m = DATA.mantel;
  svg.append("text").attr("class", "stat-box").attr("x", width - margin.right).attr("y", margin.top)
    .attr("text-anchor", "end")
    .text(`Mantel r = ${{m.r.toFixed(3)}}, p = ${{m.p_value.toFixed(3)}}  (n = ${{m.n_pairs}} pairs)`);

  // Hover: nearest point via quadtree
  const quad = d3.quadtree(points, d => x(d.phylo), d => y(d.synt));
  const hoverLayer = svg.append("g");
  canvas.style.pointerEvents = "auto";
  canvas.onmousemove = (event) => {{
    const rect = canvas.getBoundingClientRect();
    const mx = event.clientX - rect.left, my = event.clientY - rect.top;
    const found = quad.find(mx, my, 20);
    hoverLayer.selectAll("*").remove();
    if (!found) {{ hideTooltip(); return; }}
    hoverLayer.append("circle").attr("cx", x(found.phylo)).attr("cy", y(found.synt))
      .attr("r", 5).attr("fill", "none").attr("stroke", "#0b0b0b").attr("stroke-width", 1.5);
    showTooltip(event,
      `<b>${{found.a}}</b> &times; <b>${{found.b}}</b><br>Synteny distance: ${{found.synt.toFixed(3)}}<br>Phylogenetic distance: ${{found.phylo.toFixed(4)}}`);
  }};
  canvas.onmouseleave = () => {{ hoverLayer.selectAll("*").remove(); hideTooltip(); }};
}}

// ---------------------------------------------------------------------------
// Heatmap: tree-ordered pairwise synteny distance
// ---------------------------------------------------------------------------
function renderHeatmap() {{
  const orgIndex = new Map(DATA.matrix.organisms.map((name, i) => [name, i]));
  const synt = DATA.matrix.synteny_distance;
  const order = TREE_ORDER;
  const n = order.length;

  const CELL = Math.max(3, Math.min(7, Math.floor((window.innerHeight - 260) / n)));
  const STRIP = 6, GAP = 2;
  const margin = {{ top: 30, right: 160, bottom: 30, left: 160 }};
  const gridSize = n * CELL;
  const width = margin.left + gridSize + margin.right;
  const height = margin.top + STRIP + GAP + gridSize + 40;

  const canvas = document.getElementById("heatmapCanvas");
  const dpr = window.devicePixelRatio || 1;
  canvas.width = width * dpr; canvas.height = height * dpr;
  canvas.style.width = width + "px"; canvas.style.height = height + "px";
  canvas.style.position = "absolute";
  const ctx = canvas.getContext("2d");
  ctx.scale(dpr, dpr);
  ctx.clearRect(0, 0, width, height);

  // Sequential blue ramp (light = low distance/conserved, dark = high distance/rearranged)
  const steps = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"];
  const color = d3.scaleLinear().domain(d3.range(steps.length).map(i => i / (steps.length - 1)))
    .range(steps).clamp(true);
  const NO_DATA = "#dedede";

  const gridX0 = margin.left, gridY0 = margin.top + STRIP + GAP;

  // Genus tick strips (top + left)
  for (let i = 0; i < n; i++) {{
    const g = GENUS_COLORS[genusOf(order[i])] || "#fff";
    ctx.fillStyle = g;
    ctx.fillRect(gridX0 + i * CELL, margin.top, CELL, STRIP);
    ctx.fillRect(margin.left - STRIP - GAP, gridY0 + i * CELL, STRIP, CELL);
  }}

  // Grid
  for (let r = 0; r < n; r++) {{
    const gi = orgIndex.get(order[r]);
    for (let c = 0; c < n; c++) {{
      const gj = orgIndex.get(order[c]);
      let val = null;
      if (gi !== undefined && gj !== undefined) val = synt[gi][gj];
      ctx.fillStyle = (val === null || val === undefined) ? NO_DATA : color(val);
      ctx.fillRect(gridX0 + c * CELL, gridY0 + r * CELL, CELL, CELL);
    }}
  }}

  // Sparse genus labels (SVG, one per contiguous run)
  const svg = d3.select("#heatmapSvg").attr("width", width).attr("height", height)
    .style("position", "absolute").style("pointer-events", "none");
  svg.selectAll("*").remove();
  // Skip a label if it would sit closer than one text-line to the previous
  // one — small genera clustered together otherwise collide illegibly.
  let runStart = 0;
  let lastLabelY = -Infinity;
  const MIN_LABEL_GAP = 11;
  for (let i = 1; i <= n; i++) {{
    if (i === n || genusOf(order[i]) !== genusOf(order[runStart])) {{
      const mid = gridY0 + (runStart + i) / 2 * CELL;
      if (mid - lastLabelY >= MIN_LABEL_GAP) {{
        svg.append("text").attr("class", "axis-label")
          .attr("x", margin.left - STRIP - GAP - 4).attr("y", mid + 3).attr("text-anchor", "end")
          .text(genusOf(order[runStart]));
        lastLabelY = mid;
      }}
      runStart = i;
    }}
  }}

  // Colorbar legend
  const legendW = 120, legendX = width - margin.right + 20, legendY = margin.top;
  const grad = svg.append("defs").append("linearGradient").attr("id", "hmGrad");
  steps.forEach((s, i) => grad.append("stop").attr("offset", `${{i / (steps.length - 1) * 100}}%`).attr("stop-color", s));
  svg.append("rect").attr("x", legendX).attr("y", legendY).attr("width", legendW).attr("height", 10)
    .attr("fill", "url(#hmGrad)").attr("stroke", "#c3c2b7");
  [0, 0.5, 1].forEach(v => {{
    svg.append("text").attr("class", "legend-tick").attr("x", legendX + v * legendW)
      .attr("y", legendY + 24).attr("text-anchor", "middle").text(v);
  }});
  svg.append("text").attr("class", "axis-label").attr("x", legendX).attr("y", legendY - 6).text("Synteny distance");
  svg.append("rect").attr("x", legendX).attr("y", legendY + 40).attr("width", 12).attr("height", 12).attr("fill", NO_DATA);
  svg.append("text").attr("class", "legend-tick").attr("x", legendX + 18).attr("y", legendY + 50).text("no gene-order data");

  // Hover: compute cell from pixel position directly (no per-cell DOM)
  canvas.style.pointerEvents = "auto";
  const hoverLayer = svg.append("g");
  canvas.onmousemove = (event) => {{
    const rect = canvas.getBoundingClientRect();
    const mx = event.clientX - rect.left, my = event.clientY - rect.top;
    const c = Math.floor((mx - gridX0) / CELL), r = Math.floor((my - gridY0) / CELL);
    hoverLayer.selectAll("*").remove();
    if (r < 0 || r >= n || c < 0 || c >= n) {{ hideTooltip(); return; }}
    hoverLayer.append("rect").attr("x", gridX0 + c * CELL).attr("y", gridY0 + r * CELL)
      .attr("width", CELL).attr("height", CELL).attr("fill", "none").attr("stroke", "#0b0b0b");
    const gi = orgIndex.get(order[r]), gj = orgIndex.get(order[c]);
    const val = (gi !== undefined && gj !== undefined) ? synt[gi][gj] : null;
    const valText = (val === null || val === undefined) ? "no gene-order data" : val.toFixed(3);
    showTooltip(event, `<b>${{order[r]}}</b><br>&times; <b>${{order[c]}}</b><br>Synteny distance: ${{valText}}`);
  }};
  canvas.onmouseleave = () => {{ hoverLayer.selectAll("*").remove(); hideTooltip(); }};
  canvas.onclick = (event) => {{
    const rect = canvas.getBoundingClientRect();
    const mx = event.clientX - rect.left, my = event.clientY - rect.top;
    const c = Math.floor((mx - gridX0) / CELL), r = Math.floor((my - gridY0) / CELL);
    if (r < 0 || r >= n || c < 0 || c >= n || r === c) return;
    openDotPlot(order[r], order[c]);
  }};
}}

// ---------------------------------------------------------------------------
// Ranking: diverging bar chart of per-organism agreement score
// ---------------------------------------------------------------------------
function renderRanking() {{
  const entries = Object.entries(DATA.organisms)
    .filter(([, v]) => v.agreement_score !== null && v.agreement_score !== undefined)
    .sort((a, b) => b[1].agreement_score - a[1].agreement_score);
  const nullEntries = Object.entries(DATA.organisms)
    .filter(([, v]) => v.agreement_score === null || v.agreement_score === undefined);

  const ROW_H = 15;
  const width = Math.max(700, document.getElementById("rankWrap").clientWidth - 20);
  const margin = {{ top: 30, right: 100, bottom: 10, left: 220 }};
  const plotW = width - margin.left - margin.right;
  const chartH = entries.length * ROW_H;
  const nullH = nullEntries.length * 12 + 30;
  const height = margin.top + chartH + nullH + 20;

  const svg = d3.select("#rankSvg").attr("width", width).attr("height", height);
  svg.selectAll("*").remove();

  const xScale = d3.scaleLinear().domain([-1, 1]).range([margin.left, margin.left + plotW]);
  const zeroX = xScale(0);

  [-1, -0.5, 0, 0.5, 1].forEach(v => {{
    svg.append("line").attr("class", "grid-line")
      .attr("x1", xScale(v)).attr("x2", xScale(v))
      .attr("y1", margin.top - 6).attr("y2", margin.top + chartH);
    svg.append("text").attr("class", "axis-label").attr("x", xScale(v)).attr("y", margin.top - 12)
      .attr("text-anchor", "middle").text(v > 0 ? "+" + v : v);
  }});
  svg.append("text").attr("class", "axis-label").attr("x", margin.left).attr("y", margin.top - 24)
    .text("discordant ←     agreement score     → concordant");

  const rows = svg.selectAll(".rank-row").data(entries).join("g")
    .attr("transform", (d, i) => `translate(0,${{margin.top + i * ROW_H}})`);

  rows.append("text").attr("class", "rank-name")
    .attr("x", margin.left - 8).attr("y", ROW_H / 2 + 4).attr("text-anchor", "end")
    .text(d => d[0]);

  rows.append("rect")
    .attr("y", 2).attr("height", ROW_H - 5)
    .attr("x", d => Math.min(zeroX, xScale(d[1].agreement_score)))
    .attr("width", d => Math.abs(xScale(d[1].agreement_score) - zeroX))
    .attr("fill", d => agreementColor(d[1].agreement_score))
    .attr("stroke", "#333").attr("stroke-width", 0.4)
    .on("mousemove", (event, d) => showTooltip(event,
      `<b>${{d[0]}}</b><br>Agreement: ${{d[1].agreement_score.toFixed(3)}}<br>Primary-contig coverage: ${{(d[1].coverage * 100).toFixed(0)}}%`))
    .on("mouseleave", hideTooltip);

  // Insufficient-data section
  const nullTop = margin.top + chartH + 24;
  svg.append("text").attr("class", "axis-label").attr("x", margin.left).attr("y", nullTop - 8)
    .text(`Insufficient data for a score (${{nullEntries.length}} genomes, <3 valid comparisons):`);
  const nullRows = svg.selectAll(".null-row").data(nullEntries).join("g")
    .attr("transform", (d, i) => `translate(${{margin.left}},${{nullTop + i * 12}})`);
  nullRows.append("text").attr("class", "rank-null").attr("y", 8).text(d => d[0]);
}}

// ---------------------------------------------------------------------------
// Tanglegram: GTDB phylogeny vs. UPGMA gene-order dendrogram, connected
// ---------------------------------------------------------------------------
function clusterLayout(rootData, leafSpacing) {{
  // Same "d3.cluster + rescale to the tightest gap" approach as
  // interactive_tree.py's layoutTree() — d3.tree() occasionally collapsed
  // unrelated leaves onto the same position on this ladder-like phylogeny;
  // d3.cluster() never does, but its default separation still gives sibling
  // pairs a smaller gap than cousin pairs, so rescale to guarantee a minimum.
  const root = d3.hierarchy(rootData, d => d.children);
  const baseHeight = root.leaves().length * leafSpacing;
  d3.cluster().size([baseHeight, 1])(root);

  // d3.cluster() positions internal nodes by how many splits remain BELOW
  // them (remaining subtree height), not by true distance from the root —
  // so a shallow branch point from one part of the tree can land at the
  // same depth as a deeply-nested one from an unrelated part, making
  // separate branches look connected/overlapping. Reposition by true
  // depth-from-root instead (leaves stay at y=1, flush, as cluster() set them).
  const maxDepth = d3.max(root.descendants(), d => d.depth);
  root.each(d => {{
    if (d.children) d.y = d.depth / maxDepth;
  }});

  const leaves = root.leaves();
  const sorted = leaves.slice().sort((a, b) => a.x - b.x);
  let minGap = Infinity;
  for (let i = 1; i < sorted.length; i++) {{
    const gap = sorted[i].x - sorted[i - 1].x;
    if (gap > 0 && gap < minGap) minGap = gap;
  }}
  if (!isFinite(minGap) || minGap <= 0) minGap = 1;
  const scale = leafSpacing / minGap;
  root.each(d => {{ d.x *= scale; }});

  // d3.cluster() places each internal node at the mean of its children's x
  // — for a large, unevenly-shaped tree, two unrelated branch points can
  // coincidentally land within a fraction of a pixel of each other this
  // way, making their connector lines look like one continuous crossing
  // line even though nothing is topologically wrong. Sweep every node
  // (leaves included, harmlessly) and nudge any remaining near-collision
  // apart by a few px.
  const MIN_NODE_GAP = 3;
  const sortedAll = root.descendants().slice().sort((a, b) => a.x - b.x);
  for (let i = 1; i < sortedAll.length; i++) {{
    const prev = sortedAll[i - 1], cur = sortedAll[i];
    if (cur.x - prev.x < MIN_NODE_GAP) cur.x = prev.x + MIN_NODE_GAP;
  }}

  const contentHeight = d3.max(leaves, d => d.x) + leafSpacing;
  return {{ root, leaves, contentHeight, maxDepth }};
}}

function renderTanglegram() {{
  const LEAF_SPACING = 16;
  // Each depth level gets a fixed pixel allotment (like a static ete3
  // render just growing the canvas for a deeper tree) rather than
  // squeezing every level of a possibly-18-deep tree into one fixed total
  // width — that cramming is what made branch points look like they
  // overlap. Both sides share the deeper tree's width so they stay visually
  // symmetric.
  const LEVEL_SPACING = 24;
  const CONNECTOR_W = 140, LABEL_MARGIN = 170;
  const margin = {{ top: 20, right: 20, bottom: 20, left: 20 }};

  const left = clusterLayout(PHYLO_SUBTREE, LEAF_SPACING);
  const right = clusterLayout(DATA.upgma_tree, LEAF_SPACING);
  const TREE_W = Math.max(left.maxDepth, right.maxDepth, 1) * LEVEL_SPACING;
  const contentHeight = Math.max(left.contentHeight, right.contentHeight);
  const width = margin.left + TREE_W + LABEL_MARGIN + CONNECTOR_W + LABEL_MARGIN + TREE_W + margin.right;
  const height = margin.top + contentHeight + margin.bottom;

  const svg = d3.select("#tangleSvg").attr("width", width).attr("height", height);
  svg.selectAll("*").remove();
  const g = svg.append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const leftX0 = 0, leftX1 = TREE_W;
  const leftConnX = leftX1 + LABEL_MARGIN;
  const rightConnX = leftConnX + CONNECTOR_W;
  const rightX0 = rightConnX + LABEL_MARGIN, rightX1 = rightX0 + TREE_W;

  function drawTree(layout, x0, x1, mirrored) {{
    // Right-angle elbow links (horizontal - vertical - horizontal) — the
    // standard rectangular dendrogram look.
    const xOf = d => x0 + (mirrored ? (1 - d.y) : d.y) * (x1 - x0);
    g.append("g").selectAll("path").data(layout.root.links()).join("path")
      .attr("d", d => {{
        const sx = xOf(d.source), tx = xOf(d.target);
        const mid = (sx + tx) / 2;
        return `M${{sx}},${{d.source.x}} H${{mid}} V${{d.target.x}} H${{tx}}`;
      }})
      .attr("fill", "none").attr("stroke", "#c3c2b7").attr("stroke-width", 1);
  }}
  drawTree(left, leftX0, leftX1, false);
  drawTree(right, rightX0, rightX1, true);

  // Leaf labels sit in the reserved margin between each tree's tips and the
  // connector zone — left side reads left-to-right from the tips, right
  // side mirrors (text-anchor end) so it reads correctly next to its tree.
  g.selectAll(".tangle-label-left").data(left.leaves).join("text")
    .attr("x", leftX1 + 4).attr("y", d => d.x)
    .attr("dy", "0.31em").attr("font-size", 9).style("font-style", "italic")
    .attr("fill", "#0b0b0b").text(d => d.data.name);
  g.selectAll(".tangle-label-right").data(right.leaves).join("text")
    .attr("x", rightX0 - 4).attr("y", d => d.x)
    .attr("dy", "0.31em").attr("text-anchor", "end")
    .attr("font-size", 9).style("font-style", "italic")
    .attr("fill", "#0b0b0b").text(d => d.data.name);

  const leftPos = new Map(left.leaves.map(d => [d.data.name, d.x]));
  const rightPos = new Map(right.leaves.map(d => [d.data.name, d.x]));
  const names = [...leftPos.keys()].filter(n => rightPos.has(n));

  const conn = g.append("g").selectAll("path").data(names).join("path")
    .attr("d", n => {{
      const y0 = leftPos.get(n), y1 = rightPos.get(n);
      // Control points share their endpoint's y, so the curve departs/arrives
      // horizontally (flat) at both ends; pushing them toward the midpoint
      // (0.7 of the way across) stretches those flat runs further in and
      // concentrates the actual bend into a shorter stretch in the middle.
      const cx0 = leftConnX + CONNECTOR_W * 0.7, cx1 = rightConnX - CONNECTOR_W * 0.7;
      return `M${{leftConnX}},${{y0}} C${{cx0}},${{y0}} ${{cx1}},${{y1}} ${{rightConnX}},${{y1}}`;
    }})
    .attr("fill", "none").attr("stroke", "#898781").attr("stroke-width", 1.2).attr("stroke-opacity", 0.35)
    .style("cursor", "default")
    .on("mousemove", (event, n) => {{
      d3.select(event.currentTarget).attr("stroke", "#2a78d6").attr("stroke-width", 2.5).attr("stroke-opacity", 1);
      showTooltip(event, `<b>${{n}}</b>`);
    }})
    .on("mouseleave", (event) => {{
      d3.select(event.currentTarget).attr("stroke", "#898781").attr("stroke-width", 1.2).attr("stroke-opacity", 0.35);
      hideTooltip();
    }});

  g.append("text").attr("class", "axis-label").attr("x", leftX1).attr("y", -6).attr("text-anchor", "end").text("Phylogeny (GTDB)");
  g.append("text").attr("class", "axis-label").attr("x", rightX0).attr("y", -6).attr("text-anchor", "start").text("Gene order (UPGMA)");
}}

// ---------------------------------------------------------------------------
// MDS: classical MDS / PCoA of the synteny distance matrix, coloured by genus
// ---------------------------------------------------------------------------
function mdsColor(name) {{
  const base = GENUS_COLORS[genusOf(name)];
  return base ? d3.color(base).darker(1.4) : "#999";
}}

function renderMDS() {{
  const width = Math.max(700, window.innerWidth - 220);
  const height = Math.max(500, window.innerHeight - 220);
  const margin = {{ top: 30, right: 170, bottom: 50, left: 60 }};

  const svg = d3.select("#mdsSvg").attr("width", width).attr("height", height);
  svg.selectAll("*").remove();

  const entries = Object.entries(DATA.mds_coords);
  const xExtent = d3.extent(entries, d => d[1].x);
  const yExtent = d3.extent(entries, d => d[1].y);
  const padX = (xExtent[1] - xExtent[0]) * 0.08 || 0.5;
  const padY = (yExtent[1] - yExtent[0]) * 0.08 || 0.5;
  const x = d3.scaleLinear().domain([xExtent[0] - padX, xExtent[1] + padX]).range([margin.left, width - margin.right]);
  const y = d3.scaleLinear().domain([yExtent[0] - padY, yExtent[1] + padY]).range([height - margin.bottom, margin.top]);

  svg.append("line").attr("class", "axis-line").attr("x1", margin.left).attr("x2", width - margin.right)
    .attr("y1", y(0)).attr("y2", y(0));
  svg.append("line").attr("class", "axis-line").attr("x1", x(0)).attr("x2", x(0))
    .attr("y1", margin.top).attr("y2", height - margin.bottom);

  svg.selectAll(".mds-pt").data(entries).join("circle").attr("class", "mds-pt")
    .attr("cx", d => x(d[1].x)).attr("cy", d => y(d[1].y)).attr("r", 4.5)
    .attr("fill", d => mdsColor(d[0])).attr("stroke", "#333").attr("stroke-width", 0.4)
    .on("mousemove", (event, d) => showTooltip(event, `<b>${{d[0]}}</b><br>${{genusOf(d[0])}}`))
    .on("mouseleave", hideTooltip);

  const genera = [...new Set(entries.map(d => genusOf(d[0])))].sort();
  const legend = svg.append("g").attr("transform", `translate(${{width - margin.right + 20}},${{margin.top}})`);
  genera.forEach((genus, i) => {{
    const row = legend.append("g").attr("transform", `translate(0,${{i * 14}})`);
    row.append("rect").attr("width", 10).attr("height", 10).attr("fill", mdsColor(genus + " x"));
    row.append("text").attr("class", "legend-tick").attr("x", 14).attr("y", 9).text(genus);
  }});

  svg.append("text").attr("class", "axis-label").attr("x", width / 2).attr("y", height - 10)
    .attr("text-anchor", "middle").text("MDS dimension 1");
  svg.append("text").attr("class", "axis-label")
    .attr("transform", `translate(14,${{height / 2}}) rotate(-90)`).attr("text-anchor", "middle")
    .text("MDS dimension 2");
}}

// ---------------------------------------------------------------------------
// Dot-plot modal: gene rank in A vs. rank in B, for one clicked pair
// ---------------------------------------------------------------------------
function openDotPlot(orgA, orgB) {{
  document.getElementById("dotPlotModal").classList.add("open");
  document.getElementById("dotPlotTitle").textContent = `${{orgA}} × ${{orgB}}`;

  const svg = d3.select("#dotPlotSvg");
  svg.selectAll("*").remove();

  const goA = (DATA.organisms[orgA] || {{}}).gene_order;
  const goB = (DATA.organisms[orgB] || {{}}).gene_order;
  const width = 480, height = 480, margin = {{ top: 20, right: 20, bottom: 50, left: 55 }};
  svg.attr("width", width).attr("height", height);

  if (!goA || !goB || !goA.length || !goB.length) {{
    svg.append("text").attr("class", "axis-label").attr("x", 20).attr("y", 30)
      .text("No gene-order data available for one or both genomes.");
    return;
  }}

  const rankA = new Map(goA.map((gene, i) => [gene, i]));
  const rankB = new Map(goB.map((gene, i) => [gene, i]));
  const shared = goA.filter(gene => rankB.has(gene));

  const x = d3.scaleLinear().domain([0, Math.max(1, goA.length - 1)]).range([margin.left, width - margin.right]);
  const y = d3.scaleLinear().domain([0, Math.max(1, goB.length - 1)]).range([height - margin.bottom, margin.top]);

  svg.append("line").attr("x1", x(0)).attr("y1", y(0)).attr("x2", x(goA.length - 1)).attr("y2", y(goA.length - 1))
    .attr("stroke", "#e1e0d9").attr("stroke-dasharray", "4,3");

  svg.selectAll(".dp-pt").data(shared).join("circle").attr("class", "dp-pt")
    .attr("cx", gene => x(rankA.get(gene))).attr("cy", gene => y(rankB.get(gene))).attr("r", 4)
    .attr("fill", "#2a78d6").attr("fill-opacity", 0.75).attr("stroke", "#184f95").attr("stroke-width", 0.4)
    .on("mousemove", (event, gene) => showTooltip(event,
      `<b>${{gene}}</b><br>Rank in ${{orgA}}: ${{rankA.get(gene)}}<br>Rank in ${{orgB}}: ${{rankB.get(gene)}}`))
    .on("mouseleave", hideTooltip);

  svg.append("text").attr("class", "axis-label").attr("x", width / 2).attr("y", height - 12)
    .attr("text-anchor", "middle").text(orgA + " — gene rank");
  svg.append("text").attr("class", "axis-label")
    .attr("transform", `translate(16,${{height / 2}}) rotate(-90)`).attr("text-anchor", "middle")
    .text(orgB + " — gene rank");
  svg.append("text").attr("class", "hint-svg").attr("x", margin.left).attr("y", 12)
    .attr("font-size", 11).attr("fill", "#898781")
    .text(`${{shared.length}} shared genes`);
}}

document.getElementById("dotPlotClose").addEventListener("click", () => {{
  document.getElementById("dotPlotModal").classList.remove("open");
}});
document.getElementById("dotPlotModal").addEventListener("click", (event) => {{
  if (event.target.id === "dotPlotModal") event.currentTarget.classList.remove("open");
}});

function downloadBlob(blob, filename) {{
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}}

// Scatter/heatmap draw their data marks on a <canvas> (cheap for thousands
// of points/cells) with an <svg> overlay for axes/legend/labels — export
// composites both; the other tabs are plain SVG.
const TAB_EXPORT_TARGETS = {{
  scatter: {{ canvasId: "scatterCanvas", svgId: "scatterSvg" }},
  heatmap: {{ canvasId: "heatmapCanvas", svgId: "heatmapSvg" }},
  ranking: {{ canvasId: null, svgId: "rankSvg" }},
  tanglegram: {{ canvasId: null, svgId: "tangleSvg" }},
  mds: {{ canvasId: null, svgId: "mdsSvg" }},
}};

// A <style> tag embedded in a standalone SVG isn't reliably applied once
// that SVG is loaded via img.src (links/marks can come out with the wrong
// fill/stroke) — so instead of relying on a stylesheet surviving that
// round-trip, copy each element's actual COMPUTED style inline. Slower
// than embedding a stylesheet, but it reads the literal values the browser
// is really rendering, so the export is guaranteed to match the screen.
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

function saveActiveTabImage() {{
  const tab = document.querySelector(".tab-btn.active").dataset.tab;
  const {{ canvasId, svgId }} = TAB_EXPORT_TARGETS[tab];
  const svgEl = document.getElementById(svgId);
  const w = parseFloat(svgEl.getAttribute("width"));
  const h = parseFloat(svgEl.getAttribute("height"));

  const svgNS = "http://www.w3.org/2000/svg";
  const clone = svgEl.cloneNode(true);
  inlineComputedStyles(svgEl, clone);
  clone.setAttribute("xmlns", svgNS);

  const svgString = new XMLSerializer().serializeToString(clone);
  const svgUrl = URL.createObjectURL(new Blob([svgString], {{ type: "image/svg+xml;charset=utf-8" }}));

  const scale = 2;
  const exportCanvas = document.createElement("canvas");
  exportCanvas.width = w * scale; exportCanvas.height = h * scale;
  const ctx = exportCanvas.getContext("2d");
  ctx.scale(scale, scale);
  ctx.fillStyle = "#fcfcfb";
  ctx.fillRect(0, 0, w, h);
  if (canvasId) ctx.drawImage(document.getElementById(canvasId), 0, 0, w, h);

  const img = new Image();
  img.onload = () => {{
    ctx.drawImage(img, 0, 0, w, h);
    URL.revokeObjectURL(svgUrl);
    exportCanvas.toBlob(blob => downloadBlob(blob, `synteny_dashboard_${{tab}}.png`), "image/png");
  }};
  img.onerror = (e) => {{ console.error("Save image failed:", e); URL.revokeObjectURL(svgUrl); }};
  img.src = svgUrl;
}}
document.getElementById("saveImage").addEventListener("click", saveActiveTabImage);

// ---------------------------------------------------------------------------
// Motility: gene presence/absence clustered by co-occurrence pattern,
// grouped into motile / non-motile genome blocks, with a per-gene
// motile% vs non-motile% dumbbell alongside the heatmap.
// ---------------------------------------------------------------------------
const MOTILE_COLOR = "#2a78d6";     // categorical slot 1 (blue)
const NONMOTILE_COLOR = "#eb6834";  // categorical slot 2 (orange)
const PRESENT_COLOR = "#52514e";
const ABSENT_COLOR = "#ece9e2";

function renderMotility() {{
  if (!MOTILITY_DATA.genes_clustered_order) return;
  const includePredicted = document.getElementById("includePredicted").checked;
  const statKey = includePredicted ? "extended" : "curated";
  const genes = MOTILITY_DATA.genes_clustered_order;
  const presence = MOTILITY_DATA.presence;
  const stats = MOTILITY_DATA.gene_stats;

  const orgsAllowed = new Set(
    MOTILITY_DATA.organisms
      .filter(o => o.status !== "unknown" && (o.source === "curated" || o.source === "madin" || includePredicted))
      .map(o => o.name)
  );
  const motileOrgs = TREE_ORDER_ALL.filter(n => orgsAllowed.has(n) &&
    MOTILITY_DATA.organisms.find(o => o.name === n).status === "motile");
  const nonmotileOrgs = TREE_ORDER_ALL.filter(n => orgsAllowed.has(n) &&
    MOTILITY_DATA.organisms.find(o => o.name === n).status === "non-motile");
  const orgOrder = [...motileOrgs, ...nonmotileOrgs];
  const nMotile = motileOrgs.length, nNonmotile = nonmotileOrgs.length;

  const orgStatusByName = new Map(MOTILITY_DATA.organisms.map(o => [o.name, o]));

  document.getElementById("motilityHint").innerHTML =
    `Rows: ${{genes.length}} genes, ordered by hierarchical clustering on presence-pattern similarity ` +
    `(genes with similar per-genome presence sit adjacent). Columns: genomes grouped Motile ` +
    `(<b style="color:${{MOTILE_COLOR}}">${{nMotile}}</b>) vs. Non-motile ` +
    `(<b style="color:${{NONMOTILE_COLOR}}">${{nNonmotile}}</b>)` +
    (nNonmotile < 10 ? ` — <b>non-motile n=${{nNonmotile}} is small, treat differences as suggestive, not conclusive</b>.` : ".") +
    ` Dark cell = gene present. Dumbbell: % of each group with the gene present.`;

  const ROW_H = 13;
  const CELL_W = 6;
  const DUMBBELL_W = 130;
  const margin = {{ top: 36, right: 20, bottom: 10, left: 92 }};
  const gapAfterLabel = 8, gapAfterDumbbell = 16, gapBetweenGroups = 8;
  const heatmapX0 = margin.left + DUMBBELL_W + gapAfterLabel + gapAfterDumbbell;
  const gridW = (nMotile + nNonmotile) * CELL_W + gapBetweenGroups;
  const width = heatmapX0 + gridW + margin.right;
  const height = margin.top + genes.length * ROW_H + margin.bottom;

  const svg = d3.select("#motilitySvg").attr("width", width).attr("height", height);
  svg.selectAll("*").remove();

  // Group header bands
  const motileX0 = heatmapX0, motileW = nMotile * CELL_W;
  const nonmotileX0 = heatmapX0 + motileW + gapBetweenGroups, nonmotileW = nNonmotile * CELL_W;
  svg.append("rect").attr("x", motileX0).attr("y", margin.top - 20).attr("width", motileW).attr("height", 6).attr("fill", MOTILE_COLOR);
  svg.append("rect").attr("x", nonmotileX0).attr("y", margin.top - 20).attr("width", nonmotileW).attr("height", 6).attr("fill", NONMOTILE_COLOR);
  svg.append("text").attr("class", "motility-group-label").attr("x", motileX0).attr("y", margin.top - 24)
    .attr("fill", MOTILE_COLOR).text(`Motile (${{nMotile}})`);
  svg.append("text").attr("class", "motility-group-label").attr("x", nonmotileX0).attr("y", margin.top - 24)
    .attr("fill", NONMOTILE_COLOR).text(`Non-motile (${{nNonmotile}})`);

  // Dumbbell axis (0-100%)
  const dx = d3.scaleLinear().domain([0, 100]).range([margin.left, margin.left + DUMBBELL_W]);
  [0, 50, 100].forEach(v => {{
    svg.append("line").attr("class", "grid-line")
      .attr("x1", dx(v)).attr("x2", dx(v)).attr("y1", margin.top - 6).attr("y2", margin.top + genes.length * ROW_H);
    svg.append("text").attr("class", "axis-label").attr("x", dx(v)).attr("y", margin.top - 10)
      .attr("text-anchor", "middle").text(v + "%");
  }});

  const rows = svg.selectAll(".motility-row").data(genes).join("g")
    .attr("transform", (g, i) => `translate(0,${{margin.top + i * ROW_H}})`);

  rows.append("text").attr("class", "motility-gene-label")
    .attr("x", margin.left - 6).attr("y", ROW_H / 2 + 3).attr("text-anchor", "end")
    .text(g => g)
    .on("mousemove", (event, g) => {{
      const s = stats[g][statKey];
      const pm = s.pct_motile === null ? "n/a" : s.pct_motile.toFixed(0) + "%";
      const pn = s.pct_nonmotile === null ? "n/a" : s.pct_nonmotile.toFixed(0) + "%";
      showTooltip(event, `<b>${{g}}</b><br>Motile: ${{pm}} (n=${{s.n_motile}})<br>Non-motile: ${{pn}} (n=${{s.n_nonmotile}})`);
    }})
    .on("mouseleave", hideTooltip);

  // Dumbbell: connecting line + two dots
  rows.each(function(g) {{
    const s = stats[g][statKey];
    const row = d3.select(this);
    if (s.pct_motile === null || s.pct_nonmotile === null) return;
    row.append("line")
      .attr("x1", dx(s.pct_motile)).attr("x2", dx(s.pct_nonmotile))
      .attr("y1", ROW_H / 2).attr("y2", ROW_H / 2)
      .attr("stroke", "#c3c2b7").attr("stroke-width", 1.5);
    row.append("circle").attr("cx", dx(s.pct_motile)).attr("cy", ROW_H / 2).attr("r", 3).attr("fill", MOTILE_COLOR);
    row.append("circle").attr("cx", dx(s.pct_nonmotile)).attr("cy", ROW_H / 2).attr("r", 3).attr("fill", NONMOTILE_COLOR);
  }});

  // Heatmap cells
  rows.each(function(g) {{
    const row = d3.select(this);
    orgOrder.forEach((org, c) => {{
      const present = (presence[org] || []).includes(g);
      const st = orgStatusByName.get(org);
      const x = heatmapX0 + (c < nMotile ? c * CELL_W : (c - nMotile) * CELL_W + motileW + gapBetweenGroups);
      row.append("rect")
        .attr("x", x).attr("y", 1).attr("width", CELL_W - 0.6).attr("height", ROW_H - 2)
        .attr("fill", present ? PRESENT_COLOR : ABSENT_COLOR)
        .on("mousemove", (event) => showTooltip(event,
          `<b>${{g}}</b> in <i>${{org}}</i><br>${{present ? "present" : "absent"}}<br>` +
          `${{st.status}} (${{st.source}})`))
        .on("mouseleave", hideTooltip);
    }});
  }});
}}

document.getElementById("includePredicted").addEventListener("change", renderMotility);

// ---------------------------------------------------------------------------
// HGT tanglegram: GTDB species tree vs. a specific orthogroup's own gene
// tree (OrthoFinder's Resolved_Gene_Trees), connected by genome. Crossing
// connectors = phylogenetic incongruence = HGT signature. Reuses
// clusterLayout() (defined above, in the synteny tanglegram section).
// ---------------------------------------------------------------------------
const hgtSelect = document.getElementById("hgtSelect");
(function populateHgtSelect() {{
  const ogs = Object.keys(TANGLEGRAM_DATA).sort();
  hgtSelect.innerHTML = ogs.map(og => {{
    const d = TANGLEGRAM_DATA[og];
    const label = d.n_leaves === d.n_genomes
      ? `${{og}} (${{d.n_genomes}} genomes, single-copy)`
      : `${{og}} (${{d.n_leaves}} copies / ${{d.n_genomes}} genomes, has paralogs)`;
    return `<option value="${{og}}">${{label}}</option>`;
  }}).join("");
}})();
hgtSelect.addEventListener("change", renderHGT);

function renderHGT() {{
  const og = hgtSelect.value;
  if (!og || !TANGLEGRAM_DATA[og]) {{
    document.getElementById("hgtMeta").textContent = "No orthogroups available -- run build_gene_tanglegram_data.py.";
    return;
  }}
  const data = TANGLEGRAM_DATA[og];
  document.getElementById("hgtMeta").textContent =
    `${{data.n_leaves}} gene-tree leaves across ${{data.n_genomes}} genomes` +
    (data.n_leaves !== data.n_genomes ? " (includes paralogs -- one connector per copy)" : "");

  const LEAF_SPACING = 16, LEVEL_SPACING = 24;
  const CONNECTOR_W = 140, LABEL_MARGIN = 205;
  const margin = {{ top: 30, right: 20, bottom: 20, left: 20 }};

  const left = clusterLayout(data.species_tree, LEAF_SPACING);

  // Rotate the gene tree's branches (children order at every internal node
  // can be swapped without changing what the tree means) to minimise
  // crossings against the species tree's fixed leaf order -- otherwise
  // even two topologically IDENTICAL trees can render looking maximally
  // tangled purely from arbitrary left/right choices at each branch
  // point. Standard barycenter heuristic: each leaf's "position" is its
  // matching left-tree x; each internal node's position is the mean of
  // its children's; sort children by that mean, bottom-up. Crossings that
  // remain after this are the real, unavoidable discordance signal.
  const leftPosForReorder = new Map(left.leaves.map(d => [d.data.match, d.x]));
  function reorderByBarycenter(node) {{
    if (!node.children || node.children.length === 0) {{
      const pos = leftPosForReorder.get(node.match);
      return pos !== undefined ? pos : 0;
    }}
    const barys = node.children.map(reorderByBarycenter);
    const order = node.children.map((c, i) => i).sort((a, b) => barys[a] - barys[b]);
    node.children = order.map(i => node.children[i]);
    return d3.mean(order.map(i => barys[i]));
  }}
  reorderByBarycenter(data.gene_tree);

  const right = clusterLayout(data.gene_tree, LEAF_SPACING);
  const TREE_W = Math.max(left.maxDepth, right.maxDepth, 1) * LEVEL_SPACING;
  const contentHeight = Math.max(left.contentHeight, right.contentHeight);
  const width = margin.left + TREE_W + LABEL_MARGIN + CONNECTOR_W + LABEL_MARGIN + TREE_W + margin.right;
  const height = margin.top + contentHeight + margin.bottom;

  const svg = d3.select("#hgtSvg").attr("width", width).attr("height", height);
  svg.selectAll("*").remove();
  const g = svg.append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const leftX0 = 0, leftX1 = TREE_W;
  const leftConnX = leftX1 + LABEL_MARGIN;
  const rightConnX = leftConnX + CONNECTOR_W;
  const rightX0 = rightConnX + LABEL_MARGIN, rightX1 = rightX0 + TREE_W;

  function drawTree(layout, x0, x1, mirrored) {{
    const xOf = d => x0 + (mirrored ? (1 - d.y) : d.y) * (x1 - x0);
    g.append("g").selectAll("path").data(layout.root.links()).join("path")
      .attr("d", d => {{
        const sx = xOf(d.source), tx = xOf(d.target);
        const mid = (sx + tx) / 2;
        return `M${{sx}},${{d.source.x}} H${{mid}} V${{d.target.x}} H${{tx}}`;
      }})
      .attr("fill", "none").attr("stroke", "#c3c2b7").attr("stroke-width", 1);
  }}
  drawTree(left, leftX0, leftX1, false);
  drawTree(right, rightX0, rightX1, true);

  // Small genus-coloured dot before each label (same palette used
  // elsewhere in this dashboard/interactive_tree.py) -- lets a mismatched
  // genus jump out visually in the gene tree without reading every label:
  // a dot sitting in the "wrong" coloured neighbourhood is the HGT tell.
  const DOT_R = 3;
  // GENUS_COLORS is a pale/pastel palette (designed for large background
  // fills elsewhere in this dashboard) -- too washed out as small point
  // markers, so darken it specifically for the dots. Hue/identity stays
  // the same, just more visible against the light page background.
  const dotColourOf = name => d3.color(GENUS_COLORS[genusOf(name)] || "#ccc").darker(1.4);
  g.selectAll(".hgt-dot-left").data(left.leaves).join("circle")
    .attr("cx", leftX1 + 4 + DOT_R).attr("cy", d => d.x).attr("r", DOT_R)
    .attr("fill", d => dotColourOf(d.data.name))
    .attr("stroke", "#333").attr("stroke-width", 0.4);
  g.selectAll(".hgt-label-left").data(left.leaves).join("text")
    .attr("x", leftX1 + 4 + DOT_R * 2 + 4).attr("y", d => d.x)
    .attr("dy", "0.31em").attr("font-size", 9).style("font-style", "italic")
    .attr("fill", "#0b0b0b").text(d => d.data.name);
  g.selectAll(".hgt-dot-right").data(right.leaves).join("circle")
    .attr("cx", rightX0 - 4 - DOT_R).attr("cy", d => d.x).attr("r", DOT_R)
    .attr("fill", d => dotColourOf(d.data.name))
    .attr("stroke", "#333").attr("stroke-width", 0.4);
  g.selectAll(".hgt-label-right").data(right.leaves).join("text")
    .attr("x", rightX0 - 4 - DOT_R * 2 - 4).attr("y", d => d.x)
    .attr("dy", "0.31em").attr("text-anchor", "end")
    .attr("font-size", 9).style("font-style", "italic")
    .attr("fill", "#0b0b0b").text(d => d.data.name);

  // Match by the "match" key (unique per gene copy), not "name" -- a
  // paralog-duplicated genome has several same-organism leaves that would
  // otherwise collide onto one ambiguous connector. See
  // build_gene_tanglegram_data.py for how this pairing is constructed.
  const leftPos = new Map(left.leaves.map(d => [d.data.match, d.x]));
  const rightPos = new Map(right.leaves.map(d => [d.data.match, d.x]));
  const keys = [...leftPos.keys()].filter(k => k && rightPos.has(k));
  // Line colour = same genus colour as its dots, so one lineage reads as
  // one consistent colour through the tangle instead of a flat grey mass
  // -- darkened the same way as the dots (see dotColourOf above), since
  // GENUS_COLORS is a pale palette meant for large background fills.
  const matchGenus = new Map(left.leaves.map(d => [d.data.match, genusOf(d.data.name)]));
  const colourOf = k => d3.color(GENUS_COLORS[matchGenus.get(k)] || "#898781").darker(1.7);

  // Horizontal-diagonal-horizontal ("stub" at each end, straight run in
  // the middle) instead of a single edge-to-edge diagonal -- the stub
  // makes the line leave/arrive perpendicular to its leaf tip, which
  // reads more cleanly than a diagonal starting right at the tip.
  const STUB_W = 24;
  function connectorPath(y0, y1) {{
    return `M${{leftConnX}},${{y0}} H${{leftConnX + STUB_W}} L${{rightConnX - STUB_W}},${{y1}} H${{rightConnX}}`;
  }}

  g.append("g").selectAll("path").data(keys).join("path")
    .attr("d", k => connectorPath(leftPos.get(k), rightPos.get(k)))
    .attr("fill", "none")
    .attr("stroke", colourOf).attr("stroke-width", 1.5).attr("stroke-opacity", 0.85)
    .on("mousemove", (event, k) => {{
      d3.select(event.currentTarget).attr("stroke-width", 3).attr("stroke-opacity", 1);
      showTooltip(event, `<b>${{k}}</b>`);
    }})
    .on("mouseleave", (event, k) => {{
      d3.select(event.currentTarget).attr("stroke-width", 1.5).attr("stroke-opacity", 0.85);
      hideTooltip();
    }});

  g.append("text").attr("class", "axis-label").attr("x", leftX1).attr("y", -10).attr("text-anchor", "end").text("Species tree (GTDB)");
  g.append("text").attr("class", "axis-label").attr("x", rightX0).attr("y", -10).text("Gene tree (this orthogroup)");
}}

function renderAll() {{
  renderScatter(); renderHeatmap(); renderRanking(); renderTanglegram(); renderMDS(); renderMotility(); renderHGT();
}}

document.getElementById("reliableOnly").addEventListener("change", e => {{
  if (e.target.checked) {{
    DATA = DATA_RELIABLE; TREE_ORDER = TREE_ORDER_RELIABLE; PHYLO_SUBTREE = PHYLO_SUBTREE_RELIABLE;
  }} else {{
    DATA = DATA_ALL; TREE_ORDER = TREE_ORDER_ALL; PHYLO_SUBTREE = PHYLO_SUBTREE_ALL;
  }}
  renderAll();
}});

renderAll();
window.addEventListener("resize", renderAll);
</script>
</body>
</html>
"""
    return html


def build_tree_order_and_subtree(tree: Tree, organisms: list[str]) -> tuple[list[str], dict]:
    """tree_order: full-tree leaf order, filtered to just this dataset's
    organisms (for heatmap row/column order). phylo_subtree: the tree
    pruned to exactly this dataset's organisms (for the tanglegram)."""
    full_order = [leaf.name for leaf in tree.get_leaves()]
    organism_set = set(organisms)
    tree_order = [name for name in full_order if name in organism_set]

    pruned = tree.copy()
    pruned.prune(organisms, preserve_branch_length=True)
    phylo_subtree = tree_to_simple_dict(pruned)
    return tree_order, phylo_subtree


def main():
    synteny_json = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SYNTENY_JSON
    synteny_reliable_json = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_SYNTENY_RELIABLE_JSON
    output_html = sys.argv[3] if len(sys.argv) > 3 else DEFAULT_OUTPUT_HTML
    motility_json = sys.argv[4] if len(sys.argv) > 4 else DEFAULT_MOTILITY_JSON

    with open(synteny_json) as f:
        data_all = json.load(f)
    with open(synteny_reliable_json) as f:
        data_reliable = json.load(f)

    try:
        with open(motility_json) as f:
            motility_data = json.load(f)
    except FileNotFoundError:
        print(f"Note: {motility_json} not found -- 'Gene content vs. motility' tab will be empty "
              f"(run motility_clustering.py first)")
        motility_data = None

    tanglegram_data = {}
    try:
        with open("gene_order/tanglegram_index.json") as f:
            tg_index = json.load(f)
        for entry in tg_index:
            og = entry["orthogroup"]
            with open(f"gene_order/tanglegram_{og}.json") as f:
                tanglegram_data[og] = json.load(f)
        print(f"Loaded {len(tanglegram_data)} HGT tanglegram orthogroup(s)")
    except FileNotFoundError:
        print("Note: gene_order/tanglegram_index.json not found -- 'HGT tanglegram' tab will be empty "
              "(run build_gene_tanglegram_data.py first)")

    print(f"Loading {NAMED_NWK} ...")
    tree = Tree(NAMED_NWK, format=1, quoted_node_names=True)

    print("Pruning tree to organisms with synteny data (for the tanglegram) ...")
    tree_order_all, phylo_subtree_all = build_tree_order_and_subtree(
        tree, data_all["matrix"]["organisms"])
    tree_order_reliable, phylo_subtree_reliable = build_tree_order_and_subtree(
        tree, data_reliable["matrix"]["organisms"])

    html = generate_html(
        data_all, data_reliable,
        tree_order_all, tree_order_reliable,
        phylo_subtree_all, phylo_subtree_reliable,
        motility_data,
        tanglegram_data,
    )
    output_path = Path(output_html)
    output_path.write_text(html)
    print(f"Saved: {output_path}")
    print(f"Open in browser: open {output_path}")


if __name__ == "__main__":
    main()
