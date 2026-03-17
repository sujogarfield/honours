#!usr/bin/env/ python3

from ete3 import Tree, TreeStyle, TextFace, NodeStyle
from pathlib import Path
import os, sys

path = Path("gene_phylo_trees") / "flge_tree.nwk"

if not path.exists():
    print(f"{path} path not found")
    sys.exit(1)

ete_dir = "ete_visualisation"
os.makedirs(ete_dir, exist_ok=True)

save_path = Path("ete_visualisation") / "flgearc.png"


# NORMAL TREE WITH BRANCH LENGTHS
# t = Tree(str(path), format=1)
# ts = TreeStyle()
# ts.show_leaf_name = True
# ts.show_branch_length = True
# ts.show_branch_support = True
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
t = Tree(str(path), format=1)
ts = TreeStyle()
ts.mode = "c"
ts.arc_start = -180
ts.arc_span = 180
ts.title.add_face(TextFace("FlgE Phylogenetic Tree of Representative GTDB Campylobacterota Species", fsize=20), column=0)

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
        nstyle["bgcolor"] = "#f2a5a0"
    elif leaf.name.startswith("Campylobacter"):
        nstyle["bgcolor"] = "#f2dfa0"
    elif leaf.name.startswith("Sulfurimonas"):
        nstyle["bgcolor"] = "#c4f2a0"
    elif leaf.name.startswith("Wolinella"):
        nstyle["bgcolor"] = "#cca0f2"
    elif leaf.name.startswith("Hippea"):
        nstyle["bgcolor"] = "#a0ccf2"
    elif leaf.name.startswith("Sulfuricurvum"):
        nstyle["bgcolor"] = "#76d67b"
    elif leaf.name.startswith("Candidatus"):
        nstyle["bgcolor"] = "#cdd1cf"
    elif leaf.name.startswith("Sulfurospirillum"):
        nstyle["bgcolor"] = "#ab91cb"
    elif leaf.name.startswith("Hydrogenimonas"):
        nstyle["bgcolor"] = "#adb7f0"
    leaf.set_style(nstyle)

t.render(str(save_path), w=1000, units="mm", tree_style=ts)
# ARC TREES

# IMPORTANT STYLE
# ts.scale =  120 # 120 pixels per branch length unit
# ts.branch_vertical_margin = 10 # 10 pixels between adjacent branches
# ts.rotation = 90 # rotate tree
