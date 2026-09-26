#!/usr/bin/env python3
"""
summarise_au.py
Collects the AU test results from katana/03_au_tests.pbs (run automatically by
04_au_summary.pbs, or locally after copying au/ back).

Each test compares tree 1 = unconstrained best tree (flagellar grouping
allowed) with tree 2 = best tree forced to keep the candidate's GTDB branch(es).
  p_AU_gtdb < 0.05  -> the alignment significantly REJECTS the GTDB arrangement
  p_AU_free < 0.05  -> the alignment significantly rejects its own unconstrained
                       tree, i.e. prefers GTDB (rare; means the free search was poor)
Per-gene p-values are Benjamini-Hochberg corrected within each candidate;
the concatenated-alignment p-values are Holm corrected across candidates.
con_rep_spread = logL range across the 3 constrained-search replicates: a
large spread means the constrained search is unstable, so treat that p-value
with caution (the best replicate is what the AU test used).

Run from the project root:  python scripts/11_flagellar_trees/summarise_au.py
Output (gene_order/flagellar_trees/au/):
  au_results.tsv   one row per candidate x alignment
  au_summary.tsv   one row per candidate
"""

import glob
import os

A = "gene_order/flagellar_trees/au"
ALPHA = 0.05


def parse_user_trees(path):
    """{tree number: {column: value}} from the USER TREES table of an .iqtree file."""
    lines = open(path).read().splitlines()
    try:
        start = next(i for i, l in enumerate(lines) if l.startswith("Tree") and "p-AU" in l)
    except StopIteration:
        return None
    header = lines[start].split()[1:]
    out = {}
    for line in lines[start + 2:]:
        tok = [t for t in line.split() if t not in ("+", "-")]
        if not tok or not tok[0].isdigit():
            break
        out[int(tok[0])] = dict(zip(header, (float(t) for t in tok[1:])))
    return out


def bh(pvals):
    """Benjamini-Hochberg q-values, same order as input."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        running = min(running, pvals[i] * n / rank)
        q[i] = running
    return q


def holm(pvals):
    """Holm-adjusted p-values, same order as input."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    adj = [0.0] * n
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (n - rank) * pvals[i]))
        adj[i] = running
    return adj


def replicate_spread(cand, aln):
    lls = []
    for path in glob.glob(f"{A}/{cand}/run/{aln}_con_s*.iqtree"):
        for line in open(path):
            if line.startswith("Log-likelihood of the tree:"):
                lls.append(float(line.split()[4]))
                break
    return round(max(lls) - min(lls), 3) if len(lls) > 1 else ""


def main():
    cands = open(f"{A}/candidates.txt").read().split()
    rows, summary = [], []
    for cand in cands:
        group = [l.split("\t")[1] for l in open(f"{A}/{cand}/group.txt").read().splitlines()]
        genes = [l.split("\t")[0] for l in open(f"{A}/{cand}/genes.tsv").read().splitlines() if l]
        results = {}
        for aln in ["concat"] + genes:
            path = f"{A}/{cand}/run/{aln}_au.iqtree"
            t = parse_user_trees(path) if os.path.exists(path) else None
            if t and 1 in t and 2 in t:
                results[aln] = t
        gene_alns = [g for g in genes if g in results]
        q = dict(zip(gene_alns, bh([results[g][2]["p-AU"] for g in gene_alns]))) if gene_alns else {}

        for aln, t in results.items():
            rows.append([cand, aln, t[1]["logL"], t[2]["logL"], round(t[1]["logL"] - t[2]["logL"], 3),
                         t[2]["p-AU"], q.get(aln, ""), t[1]["p-AU"], replicate_spread(cand, aln)])

        c = results.get("concat")
        n_rej = sum(results[g][2]["p-AU"] < ALPHA for g in gene_alns)
        n_rej_bh = sum(q[g] < ALPHA for g in gene_alns)
        n_pref_gtdb = sum(results[g][1]["p-AU"] < ALPHA for g in gene_alns)
        summary.append([cand, "; ".join(group),
                        c[2]["p-AU"] if c else "NA", None,
                        round(c[1]["logL"] - c[2]["logL"], 3) if c else "NA",
                        f"{len(gene_alns)}/{len(genes)}", n_rej, n_rej_bh, n_pref_gtdb,
                        replicate_spread(cand, "concat")])

    have = [r for r in summary if r[2] != "NA"]
    for r, adj in zip(have, holm([r[2] for r in have])):
        r[3] = round(adj, 5)
    for r in summary:
        if r[3] is None:
            r[3] = "NA"

    with open(f"{A}/au_results.tsv", "w") as f:
        f.write("candidate\talignment\tlogL_free\tlogL_gtdb_constrained\tdelta_logL\t"
                "p_AU_gtdb\tq_BH_gtdb\tp_AU_free\tcon_rep_spread\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    with open(f"{A}/au_summary.tsv", "w") as f:
        f.write("candidate\tgroup\tconcat_p_AU_gtdb\tconcat_p_holm\tconcat_delta_logL\tgenes_done\t"
                "genes_reject_gtdb_p05\tgenes_reject_gtdb_BH05\tgenes_prefer_gtdb_p05\tconcat_con_rep_spread\n")
        for r in summary:
            f.write("\t".join(str(x) for x in r) + "\n")

    print(f"{'candidate':30s} {'concat p-AU':>11s} {'Holm':>8s} {'ΔlogL':>8s} {'genes':>6s} "
          f"{'reject':>6s} {'BH':>4s} {'rep spread':>10s}")
    for r in summary:
        p = f"{r[2]:.4f}" if r[2] != "NA" else "NA"
        h = f"{r[3]:.4f}" if r[3] != "NA" else "NA"
        print(f"{r[0]:30s} {p:>11s} {h:>8s} {str(r[4]):>8s} {r[5]:>6s} {r[6]:>6d} {r[7]:>4d} {str(r[9]):>10s}")
    print(f"\nHolm-corrected concat p-AU < {ALPHA}: the flagellar alignment significantly rejects "
          f"GTDB at that spot.\nA replicate spread of more than a few logL units means the "
          f"constrained search was unstable; check before trusting that p-value.")
    print(f"Wrote {A}/au_results.tsv and {A}/au_summary.tsv")


if __name__ == "__main__":
    main()
