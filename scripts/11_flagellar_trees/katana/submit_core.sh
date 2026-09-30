#!/bin/bash
# Core-gene (OrthoFinder single-copy) backbone test on Katana, from anywhere in the repo:
#   bash scripts/11_flagellar_trees/katana/submit_core.sh
# Chains: align/trim -> core tree -> per-candidate constrained searches + AU -> summary.
# Inputs (core/) are committed; see prepare_core_backbone.py.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
K=scripts/11_flagellar_trees/katana
N=$(grep -c . gene_order/flagellar_trees/core/candidates.txt)
ALIGN=$(qsub "$K/11_core_align.pbs");                               echo "align/trim: $ALIGN"
TREE=$(qsub -W "depend=afterok:$ALIGN" "$K/12_core_tree.pbs");       echo "core tree:  $TREE"
AU=$(qsub -W "depend=afterok:$TREE" -J "1-$N" "$K/13_core_au.pbs");  echo "core AU:    $AU ($N candidates)"
SUM=$(qsub -W "depend=afterany:$AU" "$K/14_core_summary.pbs");       echo "summary:    $SUM"
echo "Watch with: qstat -u \$USER -t"
