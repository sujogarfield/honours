#!/bin/bash
# Submit the AU tests on Katana, from anywhere in the repo:
#   bash scripts/11_flagellar_trees/katana/submit_au.sh
# Needs the tree run (submit.sh) finished and au/ inputs from prepare_au_tests.py (committed).

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
K=scripts/11_flagellar_trees/katana
N=$(grep -c . gene_order/flagellar_trees/au/candidates.txt)

ARRAY=$(qsub -J "1-$N" "$K/03_au_tests.pbs")
echo "AU tests: $ARRAY ($N candidates)"
SUMMARY=$(qsub -W "depend=afterany:$ARRAY" "$K/04_au_summary.pbs")
echo "summary:  $SUMMARY"
echo "Watch with: qstat -u \$USER -t"
