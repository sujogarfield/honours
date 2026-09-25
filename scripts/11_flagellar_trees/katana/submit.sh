#!/bin/bash
# Submit the whole flagellar tree pipeline on Katana, from anywhere in the repo:
#   bash scripts/11_flagellar_trees/katana/submit.sh
# Chains: [env setup, only if missing] -> per-gene array -> concat + concordance + comparison.
# Override the conda install with CONDA_BASE=/path/to/miniforge3 if needed.

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
K=scripts/11_flagellar_trees/katana
N=$(grep -c . gene_order/flagellar_trees/genes.txt)
QV=(-v "CONDA_BASE=${CONDA_BASE:-/srv/scratch/baker/miniforge3}")

DEP=()
if [ ! -x .katana_env/bin/iqtree2 ]; then
    SETUP=$(qsub "${QV[@]}" "$K/00_setup_env.pbs")
    echo "env setup:  $SETUP"
    DEP=(-W "depend=afterok:$SETUP")
fi

ARRAY=$(qsub "${QV[@]}" "${DEP[@]}" -J "1-$N" "$K/01_gene_trees.pbs")
echo "gene trees: $ARRAY ($N genes)"
CONCAT=$(qsub "${QV[@]}" -W "depend=afterok:$ARRAY" "$K/02_concat_concordance.pbs")
echo "concat:     $CONCAT"
echo "Watch with: qstat -u \$USER -t"
