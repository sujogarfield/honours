#!/bin/bash
# One AU test: constrained ML search (3 replicates, different seeds; the best
# is used) + AU on {unconstrained best tree, best constrained tree}.
# A single constrained search can stall on a poor tree, which would make the
# AU test over-reject GTDB -- the replicates guard against that, and their
# logL spread is reported by summarise_au.py.
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

for SEED in 1 2 3; do
    P=$W/${ALN}_con_s$SEED
    [ -s "$P.treefile" ] || \
        iqtree2 "${DATA[@]}" -g "$CONSTRAINT" -T "$T" --seed "$SEED" --prefix "$P" --quiet
done

# highest log-likelihood replicate
TOP=$(for SEED in 1 2 3; do
        printf '%s %s\n' "$(awk '/^Log-likelihood of the tree:/{print $5; exit}' "$W/${ALN}_con_s$SEED.iqtree")" "$SEED"
      done | sort -g -r | head -1 | cut -d' ' -f2)
echo "$TOP" > "$W/${ALN}_con_best_seed.txt"

cat "$BEST" "$W/${ALN}_con_s$TOP.treefile" > "$W/${ALN}_trees.nwk"
[ -s "$W/${ALN}_au.iqtree" ] || \
    iqtree2 "${DATA[@]}" -z "$W/${ALN}_trees.nwk" -n 0 -zb 10000 -au -T "$T" \
        --prefix "$W/${ALN}_au" --quiet
echo "done $CAND $ALN (best constrained replicate: seed $TOP)"
