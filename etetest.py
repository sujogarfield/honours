#!usr/bin/env/ python3

from ete3 import Tree, TreeStyle, TextFace, NodeStyle
from pathlib import Path
import os, sys

path = Path("gene_phylo_trees") / "motb_tree.nwk"

if not path.exists():
    print(f"{path} path not found")
    sys.exit(1)

ete_dir = "ete_visualisation"
os.makedirs(ete_dir, exist_ok=True)

save_path = Path("ete_visualisation") / "motbnormlen.png"


# NORMAL TREE WITH BRANCH LENGTHS
t = Tree(str(path), format=1)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
# t.render(str(save_path), w=300, units="mm", tree_style=ts)
# NORMAL TREE WITH BRANCH LENGTHS

# NORMAL TREE
# t = Tree(str(path), format=1)
# t.render(str(save_path), w=300, units="mm")
# NORMAL TREE

# CIRCULAR TREES
# t = Tree(str(path), format=1)
# circular_style = TreeStyle()
# circular_style.mode = "c"
# circular_style.scale = 200
# t.render(str(save_path), w=800, units="mm", tree_style=circular_style)
# CIRCULAR TREES

# ARC TREES
# t = Tree(str(path), format=1)
# ts = TreeStyle()
# ts.mode = "c"
# ts.arc_start = -180
# ts.arc_span = 180
# ARC TREES

ts.title.add_face(TextFace("motB Phylogenetic Tree of Representative GTDB Campylobacterota Species", fsize=20), column=0)

for n in t.traverse():
   nstyle = NodeStyle()
   nstyle["fgcolor"] = "red"
   nstyle["size"] = 8
   n.set_style(nstyle)

t.img_style["size"] = 15
t.img_style["fgcolor"] = "blue"

for leaf in t.iter_leaves():
    nstyle = NodeStyle()
    nstyle["size"] = 6
    if leaf.name.startswith("Helicobacter"):
        nstyle["bgcolor"] = "#ff5b50"
    elif leaf.name.startswith("Campylobacter"):
        nstyle["bgcolor"] = "#ffd54c"
    elif leaf.name.startswith("Sulfurimonas"):
        nstyle["bgcolor"] = "#7c9fd0"
    elif leaf.name.startswith("Wolinella"):
        nstyle["bgcolor"] = "#deb9ff"
    elif leaf.name.startswith("Hippea"):
        nstyle["bgcolor"] = "#b8deff"
    elif leaf.name.startswith("Sulfuricurvum"):
        nstyle["bgcolor"] = "#ab45ff"
    elif leaf.name.startswith("Candidatus"):
        nstyle["bgcolor"] = "#cdd1cf"
    elif leaf.name.startswith("Sulfurospirillum"):
        nstyle["bgcolor"] = "#ffb5ff"
    elif leaf.name.startswith("Hydrogenimonas"):
        nstyle["bgcolor"] = "#b7fff2"
    elif leaf.name.startswith("Nitrosophilus"):
        nstyle["bgcolor"] = "#8a86c7"
    elif leaf.name.startswith("Desulfurella"):
        nstyle["bgcolor"] = "#86b4c7"
    elif leaf.name.startswith("Caminibacter"):
        nstyle["bgcolor"] = "#c786c3"
    elif leaf.name.startswith("Nautilia"):
        nstyle["bgcolor"] = "#c78686"
    elif leaf.name.startswith("Lebetimonas"):
        nstyle["bgcolor"] = "#90c786"
    elif leaf.name.startswith("Nitratiruptor"):
        nstyle["bgcolor"] = "#bec786"
    elif leaf.name.startswith("Sulfurovum"):
        nstyle["bgcolor"] = "#c7b386"
    elif leaf.name.startswith("Nitratifractor"):
        nstyle["bgcolor"] = "#ffa9a9"
    elif leaf.name.startswith("Arcobacter"):
        nstyle["bgcolor"] = "#14ff47"
    elif leaf.name.startswith("Aliarcobacter"):
        nstyle["bgcolor"] = "#6577ff"
    elif leaf.name.startswith("Malaciobacter"):
        nstyle["bgcolor"] = "#ff43dc"
    elif leaf.name.startswith("Halarcobacter"):
        nstyle["bgcolor"] = "#43e6ff"
    elif leaf.name.startswith("Poseidonibacter"):
        nstyle["bgcolor"] = "#fffa9b"
    leaf.set_style(nstyle)

t.render(str(save_path), w=2500, units="mm", tree_style=ts)

# COLOUR
# ts.title.add_face(TextFace("16S Phylogenetic Tree of Representative GTDB Campylobacterota Species", fsize=20), column=0)
# for n in t.traverse():
#    nstyle = NodeStyle()
#    nstyle["fgcolor"] = "red"
#    nstyle["size"] = 8
#    n.set_style(nstyle)

# t.img_style["size"] = 15
# t.img_style["fgcolor"] = "blue"

# for leaf in t.iter_leaves():
#     nstyle = NodeStyle()
#     nstyle["size"] = 6
#     if leaf.name.startswith("Helicobacter"):
#         nstyle["bgcolor"] = "#ff5b50"
#     elif leaf.name.startswith("Campylobacter"):
#         nstyle["bgcolor"] = "#ffd54c"
#     elif leaf.name.startswith("Sulfurimonas"):
#         nstyle["bgcolor"] = "#7c9fd0"
#     elif leaf.name.startswith("Wolinella"):
#         nstyle["bgcolor"] = "#deb9ff"
#     elif leaf.name.startswith("Hippea"):
#         nstyle["bgcolor"] = "#b8deff"
#     elif leaf.name.startswith("Sulfuricurvum"):
#         nstyle["bgcolor"] = "#ab45ff"
#     elif leaf.name.startswith("Candidatus"):
#         nstyle["bgcolor"] = "#cdd1cf"
#     elif leaf.name.startswith("Sulfurospirillum"):
#         nstyle["bgcolor"] = "#ffb5ff"
#     elif leaf.name.startswith("Hydrogenimonas"):
#         nstyle["bgcolor"] = "#b7fff2"
#     elif leaf.name.startswith("Nitrosophilus"):
#         nstyle["bgcolor"] = "#8a86c7"
#     elif leaf.name.startswith("Desulfurella"):
#         nstyle["bgcolor"] = "#86b4c7"
#     elif leaf.name.startswith("Caminibacter"):
#         nstyle["bgcolor"] = "#c786c3"
#     elif leaf.name.startswith("Nautilia"):
#         nstyle["bgcolor"] = "#c78686"
#     elif leaf.name.startswith("Lebetimonas"):
#         nstyle["bgcolor"] = "#90c786"
#     elif leaf.name.startswith("Nitratiruptor"):
#         nstyle["bgcolor"] = "#bec786"
#     elif leaf.name.startswith("Sulfurovum"):
#         nstyle["bgcolor"] = "#c7b386"
#     elif leaf.name.startswith("Nitratifractor"):
#         nstyle["bgcolor"] = "#ffa9a9"
#     elif leaf.name.startswith("Arcobacter"):
#         nstyle["bgcolor"] = "#14ff47"
#     elif leaf.name.startswith("Aliarcobacter"):
#         nstyle["bgcolor"] = "#6577ff"
#     elif leaf.name.startswith("Malaciobacter"):
#         nstyle["bgcolor"] = "#ff43dc"
#     elif leaf.name.startswith("Halarcobacter"):
#         nstyle["bgcolor"] = "#43e6ff"
#     elif leaf.name.startswith("Poseidonibacter"):
#         nstyle["bgcolor"] = "#fffa9b"
#     leaf.set_style(nstyle)

# IMPORTANT STYLE
# ts.scale =  120 # 120 pixels per branch length unit
# ts.branch_vertical_margin = 10 # 10 pixels between adjacent branches
# ts.rotation = 90 # rotate tree
