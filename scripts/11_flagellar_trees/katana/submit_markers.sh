#!/bin/bash
# Species-marker test of the AU candidates, on Katana, from anywhere in the repo:
#   bash scripts/11_flagellar_trees/katana/submit_markers.sh
# Chains: marker tree -> per-candidate constrained searches + AU -> summary.
# Inputs (markers/) are committed; see prepare_marker_test.py.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
K=scripts/11_flagellar_trees/katana
N=$(grep -c . gene_order/flagellar_trees/markers/candidates.txt)
TREE=$(qsub "$K/05_marker_tree.pbs")
echo "marker tree: $TREE"
ARRAY=$(qsub -W "depend=afterok:$TREE" -J "1-$N" "$K/06_marker_au.pbs")
echo "marker AU:   $ARRAY ($N candidates)"
SUMMARY=$(qsub -W "depend=afterany:$ARRAY" "$K/07_marker_summary.pbs")
echo "summary:     $SUMMARY"
echo "Watch with: qstat -u \$USER -t"
