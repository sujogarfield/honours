#!/bin/bash
# One AU test, as scripts/11_flagellar_trees/katana/au_one.sh but for a flagpipe
# run: constrained ML search (3 seeds, best kept) + AU on {unconstrained best
# tree, best constrained tree}, 10,000 RELL. Each step is skipped if its output exists.
# Usage: au_one.sh <trees dir> <candidate> <gene|concat> <model|-> <threads>
set -euo pipefail
D=$1; CAND=$2; ALN=$3; MODEL=$4; T=$5
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
for SEED in 1 2 3; do
    P=$W/${ALN}_con_s$SEED
    [ -s "$P.treefile" ] || \
        iqtree2 "${DATA[@]}" -g "$CONSTRAINT" -T "$T" --seed "$SEED" --prefix "$P" --quiet
done
TOP=$(for SEED in 1 2 3; do
        printf '%s %s\n' "$(awk '/^Log-likelihood of the tree:/{print $5; exit}' "$W/${ALN}_con_s$SEED.iqtree")" "$SEED"
      done | sort -g -r | head -1 | cut -d' ' -f2)
echo "$TOP" > "$W/${ALN}_con_best_seed.txt"
cat "$BEST" "$W/${ALN}_con_s$TOP.treefile" > "$W/${ALN}_trees.nwk"
[ -s "$W/${ALN}_au.iqtree" ] || \
    iqtree2 "${DATA[@]}" -z "$W/${ALN}_trees.nwk" -n 0 -zb 10000 -au -T "$T" --seed 1 \
        --prefix "$W/${ALN}_au" --quiet
echo "done $CAND $ALN (best constrained replicate: seed $TOP)"
