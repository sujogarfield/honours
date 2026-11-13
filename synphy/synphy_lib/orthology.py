"""
orthology.py
The thesis gene-finding logic (scripts/03_gene_order + scripts/06_orthology),
for any genome set and gene list:

  check         how many genomes each gene name is found in (gene= tags), plus
                flagellar-looking names in the GFFs that are not on the list
  string-match  gene= tag matches                       -> orthology/string_match.json
  crossref      seeds -> OrthoFinder orthogroups -> candidates
                  internal: string-matched proteins of each name
                  external: seed-reference proteins whose product matches the
                            keywords, DIAMOND best hit -> orthogroup
                                                        -> orthology/candidates.json, .faa
  confirm       eggNOG agreement + combine              -> gene_order_combined.json

Confirmation rules (same as the thesis):
  internal  confirmed if eggNOG's Preferred_name equals the seed name, or the
            name appears as a whole word in the eggNOG description. Reclassified,
            unconfirmed and no-hit candidates are dropped.
  external  orthogroup labelled with eggNOG's dominant Preferred_name if >= 80% of
            its named members agree (orthology_eggnog_confirm_regulon.py). If not,
            but every named member carries a name on the gene list (a paralog
            family of target genes, like flaA/flaB in OG0000010), each protein
            takes its own eggNOG name (orthology_eggnog_confirm_flagellin.py).
            Otherwise dropped.
"""
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict

from synphy_lib.common import (Run, build_gff_index, load_eggnog, load_orthogroups, read_fasta,
                       write_json)


def _gene_tag(attrs):
    for field in attrs.split(";"):
        if field.lower().startswith("gene="):
            return field.split("=", 1)[1].strip().lower()
    return ""


def _gff_genes(gff_path):
    """[(gene_tag, contig, start, end, strand)] for GFF `gene` lines."""
    out = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            tag = _gene_tag(p[8])
            if tag:
                out.append((tag, p[0], int(p[3]), int(p[4]), p[6]))
    return out


# Names that mean different, non-flagellar genes outside Campylobacterota.
# String matching can't tell them apart, so check() warns when they're used.
COLLISIONS = {
    "pfla": "pflA = pyruvate formate-lyase activating enzyme in most bacteria (flagellar only in Campylobacterota)",
    "pflb": "pflB = pyruvate formate-lyase in most bacteria (flagellar only in Campylobacterota)",
    "maf": "maf = Maf septum-formation nucleotide pyrophosphatase in most bacteria "
           "(flagellin glycosylation only in Campylobacterota/Helicobacter)",
    "fliy": "fliY = cystine-binding protein in E. coli (a flagellar switch protein in Bacillus and Campylobacter)",
}


# ── check ───────────────────────────────────────────────────────────────────
def check(run_dir):
    run = Run(run_dir)
    cfg = run.config()
    genomes = run.genomes()
    wanted = set(cfg["genes"])
    found, other = Counter(), Counter()
    for g in genomes:
        tags = {t[0] for t in _gff_genes(run.gff / f"{g['base']}.gff")}
        found.update(tags & wanted)
        other.update(t for t in tags - wanted if re.match(r"^(fl[aghi]|mot|che|hag)", t))
    n = len(genomes)
    print(f"Genomes carrying each gene name (gene= tag, {n} genomes):")
    for gene in cfg["genes"]:
        bar = "#" * round(30 * found[gene] / n) if n else ""
        print(f"  {gene:8s} {found[gene]:4d}  {bar}")
    absent = [g for g in cfg["genes"] if not found[g]]
    if absent:
        print(f"\nNot found in any GFF ({len(absent)}): {', '.join(absent)}")
        print("  These can still be recovered by the external seed reference, but not by internal seeds.")
    risky = sorted({k for k in COLLISIONS for g in cfg["genes"] if found[g] and g.startswith(k)})
    if risky:
        print("\nWARNING -- names on your list that often mean a non-flagellar gene:")
        for k in risky:
            print(f"  {COLLISIONS[k]}")
        print("  Remove them from \"genes\" unless you are sure; the seed reference still finds real homologs.")
    if other:
        print("\nFlagellar/chemotaxis-looking names in the GFFs that are NOT on your list:")
        for t, c in other.most_common(25):
            print(f"  {t:10s} {c:4d}")
        print('Add any you want to "genes" in config.json.')


# ── string match ────────────────────────────────────────────────────────────
def string_match(run_dir):
    run = Run(run_dir)
    cfg = run.config()
    wanted = set(cfg["genes"])
    out = {}
    for g in run.genomes():
        recs = [{"gene": t, "contig": c, "start": s, "end": e, "strand": st,
                 "organism": g["label"], "source": "string_match"}
                for t, c, s, e, st in _gff_genes(run.gff / f"{g['base']}.gff") if t in wanted]
        if recs:
            out[g["base"]] = {"organism": g["label"], "genes": recs}
    write_json(run.orthology / "string_match.json", out)
    print(f"{sum(len(d['genes']) for d in out.values())} string-matched genes in {len(out)} genomes "
          f"-> {run.orthology / 'string_match.json'}")


# ── crossref ────────────────────────────────────────────────────────────────
def _reference_queries(run, keywords):
    """Seed-reference proteins whose product (FASTA header) matches a keyword."""
    faa = next(run.reference.glob("*.faa"), None)
    if faa is None:
        sys.exit(f"No seed reference proteome in {run.reference} -- rerun fetch.")
    seqs, products = read_fasta(faa), {}
    with open(faa) as f:
        for line in f:
            if line.startswith(">"):
                pid, _, rest = line[1:].rstrip("\n").partition(" ")
                products[pid] = rest
    pattern = re.compile("|".join(re.escape(k) for k in keywords), re.IGNORECASE)
    hits = {p: seqs[p] for p, prod in products.items() if pattern.search(prod)}
    return hits, products


def _diamond_best_hits(run, queries, genomes):
    """Best DIAMOND hit of each query among all genomes' proteins."""
    if not shutil.which("diamond"):
        sys.exit("diamond not on PATH (it is in the synphy conda env).")
    work = run.orthology / "diamond"
    work.mkdir(parents=True, exist_ok=True)
    db_faa, seen = work / "all_proteins.faa", set()
    with open(db_faa, "w") as out:
        for g in genomes:
            for pid, seq in read_fasta(run.faa / f"{g['base']}.faa").items():
                if pid not in seen:  # WP_ IDs are shared between genomes
                    seen.add(pid)
                    out.write(f">{pid}\n{seq}\n")
    q_faa = work / "reference_queries.faa"
    q_faa.write_text("".join(f">{p}\n{s}\n" for p, s in queries.items()))
    subprocess.run(["diamond", "makedb", "--in", str(db_faa), "-d", str(work / "all_proteins"),
                    "--quiet"], check=True)
    hits_tsv = work / "reference_hits.tsv"
    subprocess.run(["diamond", "blastp", "-q", str(q_faa), "-d", str(work / "all_proteins"),
                    "-o", str(hits_tsv), "--max-target-seqs", "1", "--evalue", "1e-5", "--quiet",
                    "--ultra-sensitive",  # default (fast) mode misses cross-genus flagellar homologs
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", "evalue"], check=True)
    best = {}
    for line in hits_tsv.read_text().splitlines():
        q, s, pident, qcov, ev = line.split("\t")
        best.setdefault(q, (s, float(pident), float(qcov)))
    return best


def crossref(run_dir):
    run = Run(run_dir)
    cfg = run.config()
    genomes = run.genomes()
    bases = {g["base"] for g in genomes}
    label = {g["base"]: g["label"] for g in genomes}
    sm = __import__("json").loads((run.orthology / "string_match.json").read_text())
    p2og, members = load_orthogroups(run.orthogroups_tsv())
    idx = {b: build_gff_index(run.gff / f"{b}.gff") for b in bases}
    print(f"{len(members)} orthogroups, {len(p2og)} proteins; {len(bases)} genomes indexed")

    # internal seeds: string-matched gene -> protein -> orthogroup
    known = defaultdict(dict)               # base -> {protein: gene}
    gene_to_ogs = defaultdict(set)
    unresolved = 0
    for base, data in sm.items():
        for r in data["genes"]:
            pid = idx[base][0].get((r["contig"], r["start"], r["end"], r["strand"]))
            if pid is None:
                unresolved += 1
                continue
            known[base][pid] = r["gene"]
            if pid in p2og:
                gene_to_ogs[r["gene"]].add(p2og[pid])
    print(f"Internal seeds: {len(gene_to_ogs)} gene names -> "
          f"{len({o for v in gene_to_ogs.values() for o in v})} orthogroups "
          f"({unresolved} string-matched records had no CDS protein_id)")

    candidates = defaultdict(list)

    def add(base, pid, gene, og, source, note=None):
        locus = idx[base][1].get(pid)
        if locus is None:
            return
        rec = {"gene": gene, "contig": locus[0], "start": locus[1], "end": locus[2],
               "strand": locus[3], "organism": label[base], "source": source,
               "orthogroup": og, "protein_id": pid}
        if note:
            rec["seed_note"] = note
        candidates[base].append(rec)

    for gene, ogs in gene_to_ogs.items():
        for og in ogs:
            for base, pids in members[og].items():
                if base in bases:
                    for pid in pids:
                        if pid not in known[base]:
                            add(base, pid, gene, og, "orthology")
    n_internal = sum(len(v) for v in candidates.values())
    print(f"Internal-seed candidates: {n_internal}")

    ref = cfg.get("seed_reference")
    ref_ogs = {}
    if ref:
        queries, products = _reference_queries(run, ref["keywords"])
        print(f"External seed: {len(queries)} reference proteins match {ref['keywords']}")
        best = _diamond_best_hits(run, queries, genomes)
        for q in queries:
            s = best.get(q)
            ref_ogs[q] = {"product": products[q], "best_hit": s[0] if s else None,
                          "pident": s[1] if s else None, "qcov": s[2] if s else None,
                          "orthogroup": p2og.get(s[0]) if s else None}
        ext_ogs = sorted({v["orthogroup"] for v in ref_ogs.values() if v["orthogroup"]})
        print(f"  -> {len(ext_ogs)} orthogroups ({sum(1 for v in ref_ogs.values() if not v['orthogroup'])} "
              f"reference proteins had no hit in an orthogroup)")
        for og in ext_ogs:
            for base, pids in members[og].items():
                if base in bases:
                    for pid in pids:
                        add(base, pid, f"ext_{og}", og, "orthology_external_seed",
                            f"seed reference {ref['accession']}")
        write_json(run.orthology / "reference_seed_orthogroups.json", ref_ogs)
    print(f"External-seed candidates: {sum(len(v) for v in candidates.values()) - n_internal}")

    out = {b: {"organism": label[b], "genes": recs} for b, recs in candidates.items()}
    write_json(run.orthology / "candidates.json", out)

    # proteins for eggNOG-mapper (each protein ID once)
    seqs_needed = defaultdict(set)
    for b, recs in candidates.items():
        for r in recs:
            seqs_needed[b].add(r["protein_id"])
    written = set()
    with open(run.orthology / "candidates.faa", "w") as f:
        for b, pids in seqs_needed.items():
            seqs = read_fasta(run.faa / f"{b}.faa")
            for pid in sorted(pids - written):
                if pid in seqs:
                    f.write(f">{pid}\n{seqs[pid].rstrip('*')}\n")
                    written.add(pid)
    print(f"{len(written)} unique candidate proteins -> {run.orthology / 'candidates.faa'} (eggNOG input)")


# ── confirm + combine ───────────────────────────────────────────────────────
def classify(gene, pref, desc):
    """orthology_eggnog_confirm.py:classify"""
    has_pref = pref and pref not in ("-", "")
    if has_pref and pref.lower() == gene.lower():
        return "confirmed"
    if desc and desc != "-" and re.search(rf"\b{re.escape(gene)}\b", desc, re.IGNORECASE):
        return "confirmed"
    return "reclassified" if has_pref else "unconfirmed"


def confirm(run_dir, annotations=None):
    import json
    run = Run(run_dir)
    cfg = run.config()
    ann_path = annotations or run.orthology / "eggnog" / "candidates.emapper.annotations"
    ann = load_eggnog(ann_path)
    cand = json.loads((run.orthology / "candidates.json").read_text())
    sm = json.loads((run.orthology / "string_match.json").read_text())
    targets = set(cfg["genes"])
    excluded = set(cfg.get("exclude_genes", []))
    print(f"{len(ann)} eggNOG annotations from {ann_path}")

    # external orthogroups: dominant eggNOG name, or per-protein for target paralog families
    votes = defaultdict(Counter)
    for data in cand.values():
        for r in data["genes"]:
            if r["source"] == "orthology_external_seed" and r["protein_id"] in ann:
                name = ann[r["protein_id"]]["preferred_name"]
                if name not in ("-", ""):
                    votes[r["orthogroup"]][name.lower()] += 1
    og_rule = {}
    for og, v in votes.items():
        name, n = v.most_common(1)[0]
        if n / sum(v.values()) >= cfg.get("dominance", 0.80):
            og_rule[og] = ("dominant", name)
        elif set(v) <= targets:
            og_rule[og] = ("per_protein", None)

    status = Counter()
    for data in cand.values():
        for r in data["genes"]:
            a = ann.get(r["protein_id"])
            r["eggnog_preferred_name"] = a["preferred_name"] if a else None
            r["eggnog_description"] = a["description"] if a else None
            if a is None:
                r["eggnog_status"] = "no_hit"
            elif r["source"] == "orthology":
                r["eggnog_status"] = classify(r["gene"], a["preferred_name"], a["description"])
            else:
                rule = og_rule.get(r["orthogroup"])
                own = a["preferred_name"].lower()
                if rule and rule[0] == "dominant":
                    r["gene"], r["eggnog_status"] = rule[1], "confirmed"
                elif rule and own not in ("-", ""):
                    r["gene"], r["eggnog_status"] = own, "confirmed"
                else:
                    r["eggnog_status"] = "unconfirmed"
            status[(r["source"], r["eggnog_status"])] += 1
    write_json(run.orthology / "candidates_confirmed.json", cand)
    for (src, st), n in sorted(status.items()):
        print(f"  {src:24s} {st:12s} {n}")

    # combine: string match + confirmed internal + confirmed external at new loci
    # exclude_genes applies to orthology candidates only, as in the thesis
    # (build_combined_gene_order.py): string-matched genes are always kept.
    combined = {b: {"organism": d["organism"], "genes": [dict(r) for r in d["genes"]]}
                for b, d in sm.items()}
    loci = defaultdict(set)
    for b, d in combined.items():
        for r in d["genes"]:
            loci[b].add((r["contig"], r["start"], r["end"], r["strand"]))
    added = Counter()
    for src in ("orthology", "orthology_external_seed"):
        for b, d in cand.items():
            for r in d["genes"]:
                if r["source"] != src or r["eggnog_status"] != "confirmed" or r["gene"] in excluded:
                    continue
                locus = (r["contig"], r["start"], r["end"], r["strand"])
                if src == "orthology_external_seed" and locus in loci[b]:
                    added["external_duplicate"] += 1
                    continue
                combined.setdefault(b, {"organism": d["organism"], "genes": []})["genes"].append(r)
                loci[b].add(locus)
                added[src] += 1
    write_json(run.root / "gene_order_combined.json", combined)
    n_sm = sum(len(d["genes"]) for d in sm.values())
    print(f"\nCombined: {sum(len(d['genes']) for d in combined.values())} genes in {len(combined)} genomes "
          f"({n_sm} string-matched, {added['orthology']} internal-seed, "
          f"{added['orthology_external_seed']} external-seed; "
          f"{added['external_duplicate']} external already present)")
    print(f"-> {run.root / 'gene_order_combined.json'}")
