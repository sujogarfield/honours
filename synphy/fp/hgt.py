"""
hgt.py
HGT screen -> shortlist: the thesis method (scripts/11_flagellar_trees/
screen_hgt_candidates.py) with the shortlist rule of
gene_order/flagellar_trees/RESULTS.md section 4, which the thesis applied by
hand, written as code so it runs on any taxon. AU tests are not run here: the
shortlist says which conflicts would be worth testing.

  hgt-screen     every gene/concatenated-tree branch with UFBoot >= 95 that is
                 incompatible with GTDB -> trees/hgt_candidates.tsv, hgt_screen.json
  hgt-shortlist  the thesis rule, per distinct group of genomes:
                   1. the group is a branch in >= min_gene_trees gene trees AND in
                      the concatenated tree;
                   2. one-sided: at the GTDB branch it contradicts, >= 50% of the
                      decisive gene trees support the group's alternative
                      arrangement, and >= 3x as many as the other alternative;
                   3. that GTDB branch is not unusually short (>= 25th percentile
                      of GTDB internal branch lengths).
                 -> trees/hgt_shortlist.tsv (every group, with the numbers)

Criterion 2 is counted directly from the gene trees, as IQ-TREE's gCF/gDF are:
the GTDB branch splits the genomes into four subtrees A,B | C,D; a gene tree is
decisive if it has at least one genome from each; it supports GTDB if it has the
branch A,B | C,D, or one of the two alternatives A,C | B,D or A,D | B,C.
A group that contradicts several GTDB branches passes 2 and 3 if any of them does.
"""
import json
import math
import os
import re
import sys

from fp.common import Run, add_scripts_path, write_json

add_scripts_path("11_flagellar_trees")

DEFAULTS = {"min_gene_trees": "auto",       # auto = 25% of genes, rounded up (10 of 37 in the thesis)
            "one_sided_min_frac": 0.50, "one_sided_ratio": 3.0,
            "min_branch_percentile": 25}


def _settings(run):
    s = dict(DEFAULTS, **run.config().get("hgt", {}))
    n_genes = len((run.trees / "genes.txt").read_text().split())
    if s["min_gene_trees"] == "auto":
        s["min_gene_trees"] = max(2, math.ceil(0.25 * n_genes))
    return s


def screen(run_dir):
    """screen_hgt_candidates.py on this run (it reads FLAG_DIR at import)."""
    run = Run(run_dir)
    os.environ["FLAG_DIR"] = str(run.trees)
    import screen_hgt_candidates
    screen_hgt_candidates.main()


def _bits(names, index):
    return sum(1 << index[n] for n in names if n in index)


def _four_subtrees(node, all_mask, mask_of):
    """The four subtrees around the branch above `node` in an unrooted tree, as
    leaf masks, or None if either end is not bifurcating."""
    if len(node.children) != 2:
        return None
    a, b = (mask_of[c] for c in node.children)
    parent = node.up
    if parent.is_root():
        others = [mask_of[c] for c in parent.children if c is not node]
        if len(others) != 2:
            return None
        c, d = others
    else:
        sibs = [mask_of[c] for c in parent.children if c is not node]
        if len(sibs) != 1:
            return None
        c = sibs[0]
        d = all_mask ^ mask_of[parent]
    return a, b, c, d


def _gtdb_branches(gtdb, index):
    """{canonical split mask: (A, B, C, D, branch length)} for internal GTDB branches."""
    t = gtdb.copy()
    t.unroot()
    full = (1 << len(index)) - 1
    mask_of = {}
    for node in t.traverse("postorder"):
        mask_of[node] = (1 << index[node.name]) if node.is_leaf() else sum(mask_of[c] for c in node.children)
    out = {}
    for node in t.traverse():
        if node.is_leaf() or node.is_root():
            continue
        m = mask_of[node]
        if not 2 <= bin(m).count("1") <= len(index) - 2:
            continue
        four = _four_subtrees(node, full, mask_of)
        if four:
            key = m ^ full if m & 1 else m
            out[key] = (*four, node.dist)
    return out


def _canon(mask, full):
    return mask ^ full if mask & 1 else mask


def _quartet_support(gene_trees, all_taxa, a, b, c, d):
    """(n_decisive, n_gtdb, n_alt_ac, n_alt_ad) over the gene trees."""
    from screen_hgt_candidates import bitmask_splits
    counts = [0, 0, 0, 0]
    for gt in gene_trees:
        taxa = sorted(gt["taxa"])
        index = {x: i for i, x in enumerate(taxa)}
        full = (1 << len(taxa)) - 1
        sub = []
        for m in (a, b, c, d):
            sub.append(_bits([all_taxa[i] for i in range(len(all_taxa)) if m >> i & 1], index))
        if not all(sub):
            continue
        counts[0] += 1
        A, B, C, D = sub
        splits = gt.setdefault("splits", set(bitmask_splits(gt["tree"], index)))
        for k, s in ((1, A | B), (2, A | C), (3, A | D)):
            if _canon(s, full) in splits:
                counts[k] += 1
                break
    return counts


def _alternative(group_mask, a, b, c, d, full):
    """Which NNI alternative the group belongs to: 'ac' (A,C | B,D), 'ad', or None."""
    from screen_hgt_candidates import compatible
    fits = [name for name, s in (("ac", a | c), ("ad", a | d)) if compatible(group_mask, s, full)]
    return fits[0] if len(fits) == 1 else None


def _slug(organisms):
    first = re.sub(r"[^a-z0-9]+", "_", organisms[0].lower()).strip("_")
    return f"{first}_and_{len(organisms) - 1}" if len(organisms) > 1 else first


def shortlist(run_dir):
    from ete3 import Tree
    from screen_hgt_candidates import compatible
    run = Run(run_dir)
    s = _settings(run)
    D = run.trees
    screen_json = D / "hgt_screen.json"
    if not screen_json.exists():
        sys.exit(f"No {screen_json} -- run `synphy.py hgt-screen {run.root}` first.")
    cands = json.loads(screen_json.read_text())["candidates"]
    org = {l.split("\t")[0]: l.split("\t")[1] for l in (D / "taxa.tsv").read_text().splitlines()[1:]}

    gtdb = Tree(str(D / "gtdb_ref_pruned.nwk"), format=1)
    all_taxa = sorted(gtdb.get_leaf_names())
    index = {t: i for i, t in enumerate(all_taxa)}
    full = (1 << len(all_taxa)) - 1
    branches = _gtdb_branches(gtdb, index)
    lengths = sorted(v[4] for v in branches.values())
    min_len = lengths[int(len(lengths) * s["min_branch_percentile"] / 100)] if lengths else 0.0

    gene_trees = []
    for p in sorted((D / "gene_trees").glob("*.treefile")):
        t = Tree(str(p), format=0)
        taxa = set(t.get_leaf_names()) & set(all_taxa)
        t.prune(list(taxa), preserve_branch_length=True)
        gene_trees.append({"gene": p.stem, "tree": t, "taxa": taxa})

    groups = {}
    for c in cands:
        g = groups.setdefault(frozenset(c["group"]), set())
        g.add(c["tree"])

    rows = []
    for group, trees in groups.items():
        n_gene = len(trees - {"CONCAT"})
        in_concat = "CONCAT" in trees
        gmask = _bits(group, index)
        best = None
        for key, (a, b, c, d, length) in branches.items():
            if compatible(gmask, key, full):
                continue
            alt = _alternative(gmask, a, b, c, d, full)
            if alt is None:
                continue
            n_dec, n_gtdb, n_ac, n_ad = _quartet_support(gene_trees, all_taxa, a, b, c, d)
            if n_dec == 0:
                continue
            fav, other = (n_ac, n_ad) if alt == "ac" else (n_ad, n_ac)
            frac = fav / n_dec
            one_sided = frac >= s["one_sided_min_frac"] and fav >= s["one_sided_ratio"] * other
            long_enough = length >= min_len
            rec = {"branch_length": round(length, 5), "decisive": n_dec, "gtdb": n_gtdb,
                   "favoured": fav, "other_alt": other, "favoured_frac": round(frac, 3),
                   "one_sided": one_sided, "long_enough": long_enough}
            rank = (one_sided and long_enough, one_sided, frac)
            if best is None or rank > best[0]:
                best = (rank, rec)
        rec = best[1] if best else {}
        passes = (n_gene >= s["min_gene_trees"] and in_concat
                  and bool(rec) and rec["one_sided"] and rec["long_enough"])
        organisms = sorted(org.get(a, a) for a in group)
        rows.append({"name": _slug(organisms), "passes": passes, "n_gene_trees": n_gene,
                     "in_concat": in_concat, **rec, "organisms": organisms, "group": sorted(group)})

    rows.sort(key=lambda r: (not r["passes"], -r["n_gene_trees"]))
    cols = ["name", "passes", "n_gene_trees", "in_concat", "favoured_frac", "favoured", "other_alt",
            "gtdb", "decisive", "one_sided", "branch_length", "long_enough", "organisms"]
    with open(D / "hgt_shortlist.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join("; ".join(r[c]) if c == "organisms" else str(r.get(c, "")) for c in cols) + "\n")
    write_json(D / "hgt_shortlist.json", {"settings": s, "min_branch_length": min_len, "groups": rows})

    passing = [r for r in rows if r["passes"]]
    print(f"{len(groups)} distinct conflicting groups; {len(passing)} pass the shortlist rule "
          f"(>= {s['min_gene_trees']} gene trees + concat, one-sided >= {s['one_sided_min_frac']:.0%} "
          f"and {s['one_sided_ratio']:g}x, GTDB branch >= {min_len:.4f})")
    print(f"Wrote {D}/hgt_shortlist.tsv")

