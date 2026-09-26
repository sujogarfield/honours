#!/bin/bash
# One AU test: constrained ML search + AU on {unconstrained, constrained} trees.
# Usage: au_one.sh <candidate> <gene|concat> <model|-> <threads>
# Called by 03_au_tests.pbs; each step is skipped if its output exists.
set -euo pipefail
CAND=$1; ALN=$2; MODEL=$3; T=$4
D=gene_order/flagellar_trees
A=$D/au
W=$A/$CAND/run
mkdir -p "$W"

if [ "$ALN" = concat ]; then
    DATA=(-p "$D/concat/concat.best_model.nex")
    CONSTRAINT=$A/$CAND/constraint.nwk
    BEST=$D/concat/concat.treefile
else
    DATA=(-s "$D/trimmed/$ALN.faa" -m "$MODEL")
    CONSTRAINT=$A/$CAND/genes/$ALN.constraint.nwk
    BEST=$D/gene_trees/$ALN.treefile
fi

[ -s "$W/${ALN}_con.treefile" ] || \
    iqtree2 "${DATA[@]}" -g "$CONSTRAINT" -T "$T" --prefix "$W/${ALN}_con" --quiet
cat "$BEST" "$W/${ALN}_con.treefile" > "$W/${ALN}_trees.nwk"
[ -s "$W/${ALN}_au.iqtree" ] || \
    iqtree2 "${DATA[@]}" -z "$W/${ALN}_trees.nwk" -n 0 -zb 10000 -au -T "$T" \
        --prefix "$W/${ALN}_au" --quiet
echo "done $CAND $ALN"
