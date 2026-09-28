#!/bin/bash
# Single-gene screen on Katana, from anywhere in the repo:
#   bash scripts/11_flagellar_trees/katana/submit_single.sh
# Needs the marker test (submit_markers.sh) finished. Inputs (single/) are committed;
# see prepare_single_gene_tests.py.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
K=scripts/11_flagellar_trees/katana
S=gene_order/flagellar_trees/single
NC=$(($(wc -l < "$S/cases.tsv") - 1))
NG=$(($(wc -l < "$S/groups.tsv") - 1))
GENES=$(qsub -J "1-$NC" "$K/08_single_genes.pbs")
echo "gene tests:   $GENES ($NC cases)"
MARKERS=$(qsub -J "1-$NG" "$K/09_single_markers.pbs")
echo "marker tests: $MARKERS ($NG groups)"
SUMMARY=$(qsub -W "depend=afterany:$GENES:$MARKERS" "$K/10_single_summary.pbs")
echo "summary:      $SUMMARY"
echo "Watch with: qstat -u \$USER -t"
