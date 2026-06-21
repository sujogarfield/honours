#!/usr/bin/env python3
"""Generate a standalone interactive HTML tree.

Usage:
    python3 interactive_tree.py [named_nwk] [og_nwk] [output_html]

Defaults:
    named_nwk = phylo_trees/gtdb_ref_tree.nwk   (organism names)
    og_nwk    = phylo_trees/gtdb_ref_tree_og.nwk (GCF IDs)
    output    = phylo_trees/interactive_tree.html

Hovering over a leaf shows the GCF accession ID.
"""

import json
import sys
from pathlib import Path
from ete3 import Tree


def tree_to_dict(node, gcf_map):
    d = {
        "name": node.name or "",
        "length": round(node.dist, 6),
    }
    if node.is_leaf():
        d["gcf"] = gcf_map.get(node.name, "")
    if node.children:
        d["children"] = [tree_to_dict(c, gcf_map) for c in node.children]
    return d


def build_gcf_map(t_named, t_og):
    """Pair leaves from both trees (same topology) to build name → GCF map."""
    named_leaves = t_named.get_leaves()
    og_leaves = t_og.get_leaves()
    if len(named_leaves) != len(og_leaves):
        print(f"Warning: leaf count mismatch ({len(named_leaves)} vs {len(og_leaves)})")
    return {l.name: r.name for l, r in zip(named_leaves, og_leaves)}


def generate_html(tree_json, leaf_count, output_path):
    height = max(600, leaf_count * 18)
    tree_data_js = json.dumps(tree_json)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Interactive Phylogenetic Tree</title>
<style>
  body {{ margin: 0; font-family: sans-serif; background: #fafafa; }}
  #controls {{ padding: 8px 16px; background: #fff; border-bottom: 1px solid #ddd; display: flex; gap: 12px; align-items: center; }}
  #controls label {{ font-size: 13px; color: #444; }}
  svg {{ display: block; }}
  .node circle {{ fill: #555; stroke: #999; stroke-width: 1px; }}
  .node.leaf circle {{ fill: #2c7bb6; }}
  .node text {{ font-size: 12px; fill: #222; }}
  .link {{ fill: none; stroke: #bbb; stroke-width: 1.2px; }}
  #tooltip {{
    position: fixed;
    background: rgba(0,0,0,0.82);
    color: #fff;
    padding: 6px 10px;
    border-radius: 4px;
    font-size: 12px;
    pointer-events: none;
    opacity: 0;
    transition: opacity 0.15s;
    max-width: 320px;
    word-break: break-all;
    z-index: 100;
  }}
  .leaf-label {{ cursor: default; }}
  .leaf-label:hover {{ fill: #e44; font-weight: bold; }}
</style>
</head>
<body>
<div id="controls">
  <label>Zoom: <input type="range" id="zoomSlider" min="0.3" max="3" step="0.05" value="1" style="width:120px"></label>
  <label><input type="checkbox" id="showBranchLen"> Show branch lengths</label>
  <span style="color:#888; font-size:12px; margin-left:8px">Hover a leaf to see its GCF accession</span>
</div>
<div id="tooltip"></div>
<div id="chart"></div>

<script src="https://d3js.org/d3.v7.min.js"></script>
<script>
const RAW = {tree_data_js};
const MARGIN = {{ top: 20, right: 280, bottom: 20, left: 60 }};
const HEIGHT = {height};

let showLen = false;
let currentZoom = 1;

function render(zoomScale) {{
  d3.select("#chart").selectAll("*").remove();

  const width = Math.max(900, window.innerWidth - 20);
  const innerW = (width - MARGIN.left - MARGIN.right) * zoomScale;
  const innerH = (HEIGHT - MARGIN.top - MARGIN.bottom) * zoomScale;

  const svg = d3.select("#chart")
    .append("svg")
    .attr("width", width)
    .attr("height", HEIGHT * zoomScale + 40)
    .append("g")
    .attr("transform", `translate(${{MARGIN.left}},${{MARGIN.top}})`);

  const root = d3.hierarchy(RAW, d => d.children);

  // compute max depth for scaling
  let maxDepth = 0;
  root.each(d => {{ if (d.depth > maxDepth) maxDepth = d.depth; }});

  const treeLayout = d3.tree().size([innerH, innerW]);
  treeLayout(root);

  // Elbow connector (right-angle cladogram style)
  const link = svg.selectAll(".link")
    .data(root.links())
    .join("path")
    .attr("class", "link")
    .attr("d", d => {{
      return `M${{d.source.y}},${{d.source.x}}
              H${{d.target.y}}
              V${{d.target.x}}`;
    }});

  const node = svg.selectAll(".node")
    .data(root.descendants())
    .join("g")
    .attr("class", d => "node" + (d.children ? "" : " leaf"))
    .attr("transform", d => `translate(${{d.y}},${{d.x}})`);

  node.append("circle").attr("r", d => d.children ? 2 : 3.5);

  // Leaf labels with tooltip
  const tooltip = document.getElementById("tooltip");

  node.filter(d => !d.children)
    .append("text")
    .attr("class", "leaf-label")
    .attr("dy", "0.31em")
    .attr("x", 8)
    .style("font-style", "italic")
    .text(d => d.data.name)
    .on("mousemove", (event, d) => {{
      const gcf = d.data.gcf || "—";
      tooltip.innerHTML = `<b>${{d.data.name}}</b><br>GCF: ${{gcf}}`;
      tooltip.style.opacity = 1;
      tooltip.style.left = (event.clientX + 14) + "px";
      tooltip.style.top  = (event.clientY - 10) + "px";
    }})
    .on("mouseleave", () => {{ tooltip.style.opacity = 0; }});

  // Internal node labels (bootstrap / branch length)
  if (showLen) {{
    node.filter(d => !d.children)
      .append("text")
      .attr("dy", "-0.5em")
      .attr("x", 8)
      .style("font-size", "9px")
      .style("fill", "#888")
      .text(d => d.data.length > 0 ? d.data.length.toFixed(4) : "");
  }}
}}

render(1);

document.getElementById("zoomSlider").addEventListener("input", e => {{
  render(parseFloat(e.target.value));
}});

document.getElementById("showBranchLen").addEventListener("change", e => {{
  showLen = e.target.checked;
  render(parseFloat(document.getElementById("zoomSlider").value));
}});
</script>
</body>
</html>
"""
    output_path.write_text(html)
    print(f"Saved: {output_path}")


def main():
    named_nwk = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("phylo_trees/gtdb_ref_tree.nwk")
    og_nwk    = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("phylo_trees/gtdb_ref_tree_og.nwk")
    output    = Path(sys.argv[3]) if len(sys.argv) > 3 else Path("phylo_trees/interactive_tree.html")

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

    tree_json = tree_to_dict(t_named, gcf_map)
    generate_html(tree_json, leaf_count, output)
    print(f"Open in browser: open {output}")


if __name__ == "__main__":
    main()
