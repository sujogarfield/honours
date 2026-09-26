#!/bin/bash
# Submit the AU tests on Katana, from anywhere in the repo:
#   bash scripts/11_flagellar_trees/katana/submit_au.sh
# Needs the tree run (submit.sh) finished. Builds the AU inputs (au/) from those
# trees first, on the login node (seconds), then submits the jobs.

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
K=scripts/11_flagellar_trees/katana
.katana_env/bin/python scripts/11_flagellar_trees/prepare_au_tests.py
N=$(grep -c . gene_order/flagellar_trees/au/candidates.txt)

ARRAY=$(qsub -J "1-$N" "$K/03_au_tests.pbs")
echo "AU tests: $ARRAY ($N candidates)"
SUMMARY=$(qsub -W "depend=afterany:$ARRAY" "$K/04_au_summary.pbs")
echo "summary:  $SUMMARY"
echo "Watch with: qstat -u \$USER -t"
