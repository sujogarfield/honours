"""
synteny.py
Gene order vs phylogeny, as in the thesis (scripts/07_combined_synteny and
scripts/04_synteny/synteny_analysis_full_contigs.py --closed), on the
combined gene set:

  all       every genome, primary contig, linear signed-adjacency + overlap distance
  reliable  every flagellar gene on one natively closed chromosome, circular
            distance (no RagTag step here: only natively complete assemblies count)

Each is compared with GTDB patristic distance by a Mantel test (Pearson r,
label permutations). Per-genome agreement = Spearman correlation between a
genome's synteny and phylogenetic distances to every other genome.
Output: synteny/synteny_all.json, synteny_reliable.json, synteny_tree.json
(per-genome values for the interactive tree).
"""
import itertools
import json

import numpy as np

from synphy_lib.common import Run, add_scripts_path, write_json

add_scripts_path("common")
from synteny_distance import primary_contig_records, signed_overlap_distance  # noqa: E402
from reliability import load_complete_chromosomes, is_reliable  # noqa: E402


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
    if rx.std() == 0 or ry.std() == 0:
        return 0.0
    return float(np.corrcoef(rx, ry)[0, 1])


def mantel_test(dist_synt, dist_phylo, valid, n_perm, seed=0):
    """Same as scripts/04_synteny/synteny_analysis.py:mantel_test."""
    n = dist_synt.shape[0]
    iu = np.triu_indices(n, k=1)
    mask = valid[iu]
    a, b = dist_synt[iu][mask], dist_phylo[iu][mask]
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return None, None, len(a)
    observed = float(np.corrcoef(a, b)[0, 1])
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(n)
        pb = dist_synt[np.ix_(perm, perm)][iu][mask]
        r = float(np.corrcoef(pb, b)[0, 1]) if pb.std() > 0 else 0.0
        count += abs(r) >= abs(observed)
    return observed, (count + 1) / (n_perm + 1), len(a)


def analyse(organisms, signed, leaf, circular, min_shared, n_perm):
    n = len(organisms)
    dist_synt, valid = np.zeros((n, n)), np.zeros((n, n), bool)
    dist_phylo = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        d, _ = signed_overlap_distance(signed[organisms[i]], signed[organisms[j]],
                                       circular=circular, min_shared=min_shared)
        if d is not None:
            dist_synt[i, j] = dist_synt[j, i] = d
            valid[i, j] = valid[j, i] = True
        dist_phylo[i, j] = dist_phylo[j, i] = leaf[organisms[i]].get_distance(leaf[organisms[j]])
    r, p, n_pairs = mantel_test(dist_synt, dist_phylo, valid, n_perm)
    per = {}
    for i, name in enumerate(organisms):
        nb = [j for j in range(n) if j != i and valid[i, j]]
        agree = spearman([dist_synt[i, j] for j in nb], [dist_phylo[i, j] for j in nb]) if len(nb) >= 3 else None
        per[name] = {"n_valid_comparisons": len(nb),
                     "agreement_score": round(agree, 4) if agree is not None else None,
                     "gene_order": [rec["gene"] for rec in signed[name]]}
    return {
        "mantel": {"r": r, "p_value": p, "n_pairs": n_pairs, "n_permutations": n_perm,
                   "n_organisms": n, "min_shared_genes": min_shared,
                   "distance": f"signed adjacency + overlap, {'circular' if circular else 'linear'}"},
        "organisms": per,
        "matrix": {"organisms": organisms,
                   "synteny_distance": [[round(dist_synt[i, j], 4) if valid[i, j] or i == j else None
                                         for j in range(n)] for i in range(n)],
                   "phylogenetic_distance": [[round(dist_phylo[i, j], 6) for j in range(n)]
                                             for i in range(n)]},
    }


def synteny(run_dir):
    from ete3 import Tree
    run = Run(run_dir)
    cfg = run.config()
    min_shared, n_perm = cfg["synteny"]["min_shared_genes"], cfg["synteny"]["permutations"]
    raw = json.loads((run.root / "gene_order_combined.json").read_text())
    gene_order = {d["organism"]: d["genes"] for d in raw.values()}
    tree = Tree(str(run.tree / "gtdb_named.nwk"), format=1, quoted_node_names=True)
    leaf = {l.name: l for l in tree.iter_leaves()}
    organisms = sorted(o for o in gene_order if o in leaf)

    signed, coverage = {}, {}
    for o in organisms:
        recs, n_primary = primary_contig_records(gene_order[o])
        signed[o] = recs
        coverage[o] = round(n_primary / len(gene_order[o]), 3) if gene_order[o] else None

    complete = load_complete_chromosomes(run.fna, run.root / "_no_ragtag")
    reliable = [o for o in organisms
                if len(gene_order[o]) >= min_shared and is_reliable(gene_order[o], complete)]

    results = {}
    for name, orgs, circular in (("all", organisms, False), ("reliable", reliable, True)):
        res = analyse(orgs, signed, leaf, circular, min_shared, n_perm)
        write_json(run.synteny / f"synteny_{name}.json", res)
        m = res["mantel"]
        r = f"{m['r']:.3f}" if m["r"] is not None else "NA"
        print(f"{name:8s}: {m['n_organisms']} genomes, {m['n_pairs']} pairs, Mantel r = {r}, p = {m['p_value']}")
        results[name] = res

    tree_values = {}
    for o in organisms:
        src = results["reliable"]["organisms"].get(o) or results["all"]["organisms"][o]
        tree_values[o] = {"coverage": coverage[o], "agreement_score": src["agreement_score"],
                          "n_valid_comparisons": src["n_valid_comparisons"],
                          "reliable": o in results["reliable"]["organisms"], "scaffolded": False}
    write_json(run.synteny / "synteny_tree.json", {"organisms": tree_values})
    print(f"Wrote {run.synteny}/")
