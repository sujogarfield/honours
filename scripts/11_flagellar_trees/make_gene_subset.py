#!/usr/bin/env python3
"""
make_gene_subset.py
Builds a gene-subset run directory from an existing run, for sensitivity
analyses of the concatenated tree and concordance factors (e.g. the 42-gene
set without motor and chemotaxis genes).

Every run with KEEP_TAXA fixes the same 149 genomes, so a gene's sequences,
alignment, trimmed alignment and gene tree are identical in any subset that
contains it. Instead of re-running prepare_gene_sets.py and the per-gene
Katana array (01_gene_trees.pbs), this copies those files for the kept
genes; only the concatenated tree and concordance step
(02_concat_concordance.pbs) needs to run.

Copies whatever exists in the source: locally that is the committed inputs
(genes/, genes.txt, taxa.tsv, gtdb_ref_pruned.nwk, gene_selection.json); on
Katana it is also aligned/, trimmed/ and gene_trees/. Run it in both places.

Usage (from the project root):
  python scripts/11_flagellar_trees/make_gene_subset.py DEST [--exclude PREFIX ...] [--src SRC]
  e.g.  python scripts/11_flagellar_trees/make_gene_subset.py \
            gene_order/flagellar_trees_no_mot_che --exclude mot che
Then on Katana:
  qsub -v FLAG_DIR=DEST scripts/11_flagellar_trees/katana/02_concat_concordance.pbs
"""
import argparse
import json
import shutil
from pathlib import Path

SHARED_FILES = ["taxa.tsv", "gtdb_ref_pruned.nwk"]
PER_GENE_DIRS = ["genes", "aligned", "trimmed", "gene_trees"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dest", type=Path)
    ap.add_argument("--src", type=Path, default=Path("gene_order/flagellar_trees_ext"))
    ap.add_argument("--exclude", nargs="+", default=[], metavar="PREFIX",
                    help="drop every gene whose name starts with one of these (e.g. mot che)")
    args = ap.parse_args()

    genes = [g for g in (args.src / "genes.txt").read_text().split() if g]
    kept = [g for g in genes if not any(g.startswith(p.lower()) for p in args.exclude)]
    dropped = [g for g in genes if g not in kept]
    print(f"{args.src}: {len(genes)} genes -> {len(kept)} kept, dropped {', '.join(dropped) or 'none'}")

    args.dest.mkdir(parents=True, exist_ok=True)
    (args.dest / "genes.txt").write_text("\n".join(kept) + "\n")
    for name in SHARED_FILES:
        shutil.copy2(args.src / name, args.dest / name)

    # Selection report: the source's, restricted to the kept genes, with a note
    # saying where the subset came from.
    sel = json.loads((args.src / "gene_selection.json").read_text())
    sel["genes"] = {g: v for g, v in sel["genes"].items() if g in kept}
    sel["subset_of"] = {"source": str(args.src), "excluded_prefixes": args.exclude,
                        "excluded_genes": dropped}
    (args.dest / "gene_selection.json").write_text(json.dumps(sel, indent=2) + "\n")

    for d in PER_GENE_DIRS:
        src_dir = args.src / d
        if not src_dir.is_dir():
            print(f"  {d}/: not in source, skipped")
            continue
        (args.dest / d).mkdir(exist_ok=True)
        n = 0
        for f in src_dir.iterdir():
            # gene files are <gene>.<ext> (e.g. flgb.faa, flgb.treefile)
            if f.name.split(".")[0] in kept:
                shutil.copy2(f, args.dest / d / f.name)
                n += 1
        print(f"  {d}/: {n} files")


if __name__ == "__main__":
    main()
