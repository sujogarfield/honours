#!/bin/bash
# Submit a flagpipe stage on Katana, from anywhere in the repo:
#   bash flagpipe/katana/submit.sh <run_dir> orthology   # stage 1
#   bash flagpipe/katana/submit.sh <run_dir> trees       # stage 2, after checking stage 1
#   bash flagpipe/katana/submit.sh <run_dir> hgt         # stage 3, AU tests on the HGT shortlist
# Watch with: qstat -u $USER -t

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
RUN=$(cd "$1" && pwd); STAGE=$2
K=flagpipe/katana
source $K/common.sh
QV=(-v "RUN=$RUN,FLAGPIPE_ENV=$FLAGPIPE_ENV")
[ -x "$FLAGPIPE_ENV/bin/iqtree2" ] || { echo "No env at $FLAGPIPE_ENV -- qsub $K/setup_env.pbs first"; exit 1; }
case "$STAGE" in
  orthology)
    [ -d "$RUN/genomes/faa" ] && [ -s "$RUN/tree/gtdb_acc.nwk" ] || { echo "Run init, genomes, fetch and tree first"; exit 1; }
    J=$(qsub "${QV[@]}" "$K/orthology.pbs"); echo "orthology: $J (log: fp_orthology.o* in the repo root)" ;;
  trees)
    N=$(grep -c . "$RUN/trees/genes.txt")
    [ "$N" -ge 4 ] || { echo "Only $N genes selected -- too few for a concatenated tree (see gene_selection.json)"; exit 1; }
    A=$(qsub "${QV[@]}" -J "1-$N" "$K/gene_trees.pbs"); echo "gene trees: $A ($N genes)"
    C=$(qsub "${QV[@]}" -W "depend=afterok:$A" "$K/concat.pbs"); echo "concat: $C" ;;
  hgt)
    L=$RUN/trees/au/candidates.txt
    [ -f "$L" ] || { echo "No $L -- stage 2 (trees) must finish first"; exit 1; }
    N=$(grep -c . "$L" || true)
    [ "$N" -ge 1 ] || { echo "No candidate passed the shortlist rule (see $RUN/trees/hgt_shortlist.tsv)"; exit 0; }
    if [ "$N" -eq 1 ]; then   # PBS arrays need >= 2 elements
        A=$(qsub -v "RUN=$RUN,FLAGPIPE_ENV=$FLAGPIPE_ENV,PBS_ARRAY_INDEX=1" "$K/au_tests.pbs")
    else
        A=$(qsub "${QV[@]}" -J "1-$N" "$K/au_tests.pbs")
    fi
    echo "AU tests: $A ($N candidates)"
    S=$(qsub "${QV[@]}" -W "depend=afterany:$A" "$K/au_summary.pbs"); echo "summary: $S" ;;
  *) echo "stage must be 'orthology', 'trees' or 'hgt'"; exit 1 ;;
esac
