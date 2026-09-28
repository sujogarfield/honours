#!/bin/bash
# One single-gene-screen test: 3 seeded constrained searches (best kept) + AU
# on {unconstrained tree, best constrained tree}, 10,000 RELL.
# Usage: single_one.sh gene   <case_id>  <gene> <model> <threads>
#        single_one.sh marker <group_id> -      <model> <threads>
# gene:   gene alignment, marker-tree branch(es) forced (does the gene reject the marker arrangement?)
# marker: marker alignment, the gene's grouping forced (do the markers reject the gene's grouping?)
set -euo pipefail
KIND=$1; ID=$2; GENE=$3; MODEL=$4; T=$5
D=gene_order/flagellar_trees
S=$D/single
W=$S/run/$ID
mkdir -p "$W"
if [ "$KIND" = gene ]; then
    ALN=$D/trimmed/$GENE.faa; CON=$S/genes/$ID.constraint.nwk; BEST=$D/gene_trees/$GENE.treefile
else
    ALN=$D/markers/bac120_msa_r232.faa; CON=$S/marker/$ID.constraint.nwk; BEST=$D/markers/free/markers.treefile
fi
for SEED in 1 2 3; do
    [ -s "$W/con_s$SEED.treefile" ] || \
        iqtree2 -s "$ALN" -m "$MODEL" -g "$CON" -T "$T" --seed "$SEED" --prefix "$W/con_s$SEED" --quiet
done
TOP=$(for SEED in 1 2 3; do
        printf '%s %s\n' "$(awk '/^Log-likelihood of the tree:/{print $5; exit}' "$W/con_s$SEED.iqtree")" "$SEED"
      done | sort -g -r | head -1 | cut -d' ' -f2)
echo "$TOP" > "$W/con_best_seed.txt"
cat "$BEST" "$W/con_s$TOP.treefile" > "$W/trees.nwk"
[ -s "$W/au.iqtree" ] || \
    iqtree2 -s "$ALN" -m "$MODEL" -z "$W/trees.nwk" -n 0 -zb 10000 -au -T "$T" --prefix "$W/au" --quiet
echo "done $KIND $ID"
