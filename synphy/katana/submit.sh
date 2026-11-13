#!/bin/bash
# Submit a synphy stage on Katana, from anywhere in the repo:
#   bash synphy/katana/submit.sh <run_dir> orthology   # stage 1
#   bash synphy/katana/submit.sh <run_dir> trees       # stage 2, after checking stage 1
# Watch with: qstat -u $USER -t

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
RUN=$(cd "$1" && pwd); STAGE=$2
K=synphy/katana
source $K/common.sh
QV=(-v "RUN=$RUN,SYNPHY_ENV=$SYNPHY_ENV")
[ -x "$SYNPHY_ENV/bin/iqtree2" ] || { echo "No env at $SYNPHY_ENV -- qsub $K/setup_env.pbs first"; exit 1; }
case "$STAGE" in
  orthology)
    [ -d "$RUN/genomes/faa" ] && [ -s "$RUN/tree/gtdb_acc.nwk" ] || { echo "Run init, genomes, fetch and tree first"; exit 1; }
    J=$(qsub "${QV[@]}" "$K/orthology.pbs"); echo "orthology: $J (log: synphy_ortho.o* in the repo root)" ;;
  trees)
    N=$(grep -c . "$RUN/trees/genes.txt")
    [ "$N" -ge 4 ] || { echo "Only $N genes selected -- too few for a concatenated tree (see gene_selection.json)"; exit 1; }
    A=$(qsub "${QV[@]}" -J "1-$N" "$K/gene_trees.pbs"); echo "gene trees: $A ($N genes)"
    C=$(qsub "${QV[@]}" -W "depend=afterok:$A" "$K/concat.pbs"); echo "concat: $C" ;;
  *) echo "stage must be 'orthology' or 'trees'"; exit 1 ;;
esac
