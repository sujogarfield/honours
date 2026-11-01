"""
genomes.py
Steps that run on the login node (they need internet):

  init     interactive questions -> config.json
  genomes  GTDB advanced-search TSV -> genomes.tsv (filter, cap, outgroup)
  fetch    NCBI genome FASTA + GFF + proteome for every genome (and the seed
           reference proteome, if one is set)
  tree     GTDB bac120 tree pruned to the genomes, with organism labels
"""
import csv
import gzip
import json
import shutil
import sys
import time
import urllib.request
from collections import defaultdict
from pathlib import Path

from fp.common import FLAGPIPE, Run, accession, clean_label, gtdb_tree_path, write_json

DEFAULT_KEYWORDS = ["flagell", "chemotax", "chemorecept", "methyl-accepting", "motility"]


# ── init ────────────────────────────────────────────────────────────────────
def ask(question, default=None):
    suffix = f" [{default}]" if default not in (None, "") else ""
    answer = input(f"{question}{suffix}: ").strip()
    return answer if answer else default


def ask_yes(question, default=True):
    d = "Y/n" if default else "y/N"
    answer = input(f"{question} [{d}]: ").strip().lower()
    return default if not answer else answer.startswith("y")


def read_gene_list(path):
    return [l.strip().lower() for l in Path(path).read_text().splitlines()
            if l.strip() and not l.startswith("#")]


def choose_genes():
    presets = sorted((FLAGPIPE / "presets").glob("*.txt"))
    print("\nGene list. Presets:")
    for i, p in enumerate(presets, 1):
        print(f"  {i}. {p.stem} ({len(read_gene_list(p))} names)")
    print(f"  {len(presets) + 1}. my own file (one gene name per line)")
    print(f"  {len(presets) + 2}. type names now")
    choice = int(ask("Choose", "1"))
    if choice <= len(presets):
        genes = read_gene_list(presets[choice - 1])
    elif choice == len(presets) + 1:
        genes = read_gene_list(ask("Path to gene list"))
    else:
        genes = [g.strip().lower() for g in ask("Gene names, comma-separated").split(",") if g.strip()]
    extra = ask("Extra names to add (comma-separated, blank for none)", "")
    genes += [g.strip().lower() for g in extra.split(",") if g.strip()]
    remove = ask("Names to remove (comma-separated, blank for none)", "")
    genes = [g for g in genes if g not in {x.strip().lower() for x in remove.split(",")}]
    genes = list(dict.fromkeys(genes))
    print(f"{len(genes)} gene names: {', '.join(genes)}")
    print("(After `fetch`, `flagpipe.py check` shows how many genomes each name is found in.)")
    return genes


def init(run_dir):
    run = Run(run_dir)
    if run.config_path.exists() and not ask_yes(f"{run.config_path} exists. Overwrite?", False):
        return
    print("flagpipe setup -- press Enter to accept the [default].\n")
    tsv = ask("GTDB advanced-search TSV (downloaded from gtdb.ecogenomic.org)")
    if not tsv or not Path(tsv).exists():
        sys.exit(f"TSV not found: {tsv}")
    cfg = {
        "name": ask("Run name", run.root.name),
        "gtdb_tsv": str(Path(tsv).resolve()),
        "taxon_filter": ask("Keep only rows whose GTDB taxonomy contains (e.g. f__Listeriaceae; "
                            "blank = all rows)", ""),
        "representatives_only": ask_yes("Only GTDB species representatives?", True),
        "max_genomes": int(ask("Maximum genomes (OrthoFinder is all-vs-all; ~300 is practical)", "300")),
        "max_per_genus": int(ask("Maximum genomes per genus", "40")),
        "outgroup": ask("Outgroup: 'auto' (sister clade in the GTDB tree), 'none', or "
                        "comma-separated GCF accessions", "auto"),
        "n_outgroup": 3,
        "genes": choose_genes(),
        "exclude_genes": [],
        "seed_reference": None,
        "email": ask("Email for NCBI Entrez (required by NCBI)"),
        "eggnog_data_dir": ask("eggNOG data directory (shared, ~50 GB)",
                               "/srv/scratch/baker/eggnog_data"),
        "eggnog_mode": "diamond",
        "dominance": 0.80,
        "selection": {"min_genome_frac": 0.60, "max_multicopy_frac": 0.10, "min_gene_frac": 0.60},
        "synteny": {"min_shared_genes": 4, "permutations": 999},
        "hgt": {"min_gene_trees": "auto", "one_sided_min_frac": 0.50, "one_sided_ratio": 3.0,
                "min_branch_percentile": 25},
    }
    if ask_yes("\nUse an external seed reference genome (e.g. B. subtilis 168, C. jejuni NCTC 11168)? "
               "Its proteins whose product matches flagellar/chemotaxis keywords seed extra orthogroups",
               True):
        cfg["seed_reference"] = {
            "accession": ask("Reference GCF accession", "GCF_000009045.1"),
            "keywords": [k.strip() for k in ask("Product keywords", ",".join(DEFAULT_KEYWORDS)).split(",")],
        }
    run.root.mkdir(parents=True, exist_ok=True)
    write_json(run.config_path, cfg)
    print(f"\nWrote {run.config_path}. Edit it directly to change anything later.")
    print(f"Next: python flagpipe/flagpipe.py genomes {run.root}")


# ── genomes ─────────────────────────────────────────────────────────────────
def genomes(run_dir):
    from ete3 import Tree
    run = Run(run_dir)
    cfg = run.config()
    with open(cfg["gtdb_tsv"]) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    print(f"{len(rows)} rows in {cfg['gtdb_tsv']}")
    keep = []
    for r in rows:
        if cfg.get("taxon_filter") and cfg["taxon_filter"] not in r["GTDB Taxonomy"]:
            continue
        if cfg.get("representatives_only", True) and str(r["GTDB Representative of Species"]) != "True":
            continue
        if not r["Accession"].startswith("GCF_"):
            continue  # GenBank-only assemblies have no RefSeq annotation
        keep.append(r)
    print(f"{len(keep)} RefSeq genomes after filters")

    # Cap per genus, then overall, keeping every genus where possible:
    # type material first, then accession order (deterministic).
    by_genus = defaultdict(list)
    for r in keep:
        genus = next((t[3:] for t in r["GTDB Taxonomy"].split(";") if t.startswith("g__")), "?")
        by_genus[genus].append(r)
    for g in by_genus:
        by_genus[g].sort(key=lambda r: ("type" not in r.get("GTDB Type Material", "").lower(), r["Accession"]))
        by_genus[g] = by_genus[g][:cfg["max_per_genus"]]
    chosen, i = [], 0
    while len(chosen) < cfg["max_genomes"] and any(len(v) > i for v in by_genus.values()):
        for g in sorted(by_genus):
            if i < len(by_genus[g]) and len(chosen) < cfg["max_genomes"]:
                chosen.append(by_genus[g][i])
        i += 1
    print(f"{len(chosen)} genomes chosen across {len(by_genus)} genera "
          f"(caps: {cfg['max_per_genus']} per genus, {cfg['max_genomes']} total)")

    tree = Tree(str(gtdb_tree_path(cfg)), format=1, quoted_node_names=True)
    leaf_by_acc = {l.name[3:]: l for l in tree.iter_leaves()}  # RS_GCF_... -> GCF_...
    out = [{"accession": r["Accession"], "label": r["Organism Name"],
            "taxonomy": r["GTDB Taxonomy"], "role": "ingroup"} for r in chosen]
    missing = [r["accession"] for r in out if r["accession"] not in leaf_by_acc]
    if missing:
        print(f"Warning: {len(missing)} genomes not in the GTDB tree (left out of tree steps): {missing[:5]}")

    og = cfg.get("outgroup", "auto")
    if og == "auto":
        ingroup = [leaf_by_acc[r["accession"]] for r in out if r["accession"] in leaf_by_acc]
        if not ingroup:
            sys.exit("No chosen genome is in the GTDB tree, so the outgroup can't be found: "
                     "use GTDB species representatives, or set \"outgroup\" to accessions or \"none\".")
        mrca = tree.get_common_ancestor(ingroup) if len(ingroup) > 1 else ingroup[0]
        sisters = [c for c in mrca.up.children if c is not mrca] if mrca.up else []
        pool = [l.name[3:] for s in sisters for l in s.iter_leaves() if l.name.startswith("RS_GCF_")]
        n = min(cfg.get("n_outgroup", 3), len(pool))
        picks = [pool[round(k * (len(pool) - 1) / max(n - 1, 1))] for k in range(n)] if n else []
        print(f"Outgroup (auto, sister clade of {len(pool)} RefSeq genomes): {picks}")
    elif og in ("none", "", None):
        picks = []
    else:
        picks = [a.strip() for a in (og if isinstance(og, list) else og.split(",")) if a.strip()]
    chosen_acc = {r["accession"] for r in out}
    out += [{"accession": a, "label": "", "taxonomy": "", "role": "outgroup"}
            for a in dict.fromkeys(picks) if a not in chosen_acc]

    with open(run.genomes_tsv, "w") as f:
        f.write("accession\tbase\tlabel\ttaxonomy\trole\n")
        for r in out:
            f.write(f"{r['accession']}\t\t{r['label']}\t{r['taxonomy']}\t{r['role']}\n")
    print(f"Wrote {run.genomes_tsv}. Next: python flagpipe/flagpipe.py fetch {run.root}")


# ── fetch ───────────────────────────────────────────────────────────────────
def _download(url, dest):
    """url (.gz) -> unzipped dest, with 3 tries."""
    for attempt in range(3):
        try:
            tmp = Path(str(dest) + ".gz")
            urllib.request.urlretrieve(url, tmp)
            with gzip.open(tmp, "rb") as f_in, open(dest, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            tmp.unlink()
            return True
        except Exception as e:
            print(f"    attempt {attempt + 1} failed: {e}")
            time.sleep(5)
    return False


def _assembly_summary(acc):
    from Bio import Entrez
    rec = Entrez.read(Entrez.esearch(db="assembly", term=acc))
    if not rec["IdList"]:
        return None
    s = Entrez.read(Entrez.esummary(db="assembly", id=rec["IdList"][0], report="full"))
    return s["DocumentSummarySet"]["DocumentSummary"][0]


def fetch_one(acc, dirs):
    """Download the RefSeq files for one accession. dirs: {suffix: folder}.
    Returns (base, organism) or (None, None)."""
    doc = _assembly_summary(acc)
    time.sleep(0.4)  # NCBI allows ~3 requests/s without an API key
    if doc is None or not doc["FtpPath_RefSeq"]:
        print(f"  {acc}: not found in NCBI Assembly / no RefSeq FTP path")
        return None, None
    ftp = doc["FtpPath_RefSeq"].replace("ftp://", "https://")
    base = ftp.rstrip("/").split("/")[-1]
    for suffix, folder in dirs.items():
        folder.mkdir(parents=True, exist_ok=True)
        ext = {"_genomic.fna.gz": ".fna", "_genomic.gff.gz": ".gff", "_protein.faa.gz": ".faa"}[suffix]
        dest = folder / f"{base}{ext}"
        if dest.exists():
            continue
        print(f"  {acc}: {suffix}")
        if not _download(f"{ftp}/{base}{suffix}", dest):
            print(f"  {acc}: FAILED {suffix}")
            return None, None
    return base, str(doc["Organism"])


def fetch(run_dir):
    from Bio import Entrez
    run = Run(run_dir)
    cfg = run.config()
    if not cfg.get("email"):
        sys.exit("Set \"email\" in config.json (NCBI requires it).")
    Entrez.email = cfg["email"]
    lines = run.genomes_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, l.split("\t"))) for l in lines[1:] if l.strip()]
    dirs = {"_genomic.fna.gz": run.fna, "_genomic.gff.gz": run.gff, "_protein.faa.gz": run.faa}
    failed = []
    for i, r in enumerate(rows, 1):
        print(f"[{i}/{len(rows)}] {r['accession']}")
        base, organism = fetch_one(r["accession"], dirs)
        if base is None:
            failed.append(r["accession"])
            r["role"] = "failed"
            continue
        r["base"] = base
        if not r["label"]:
            r["label"] = organism

    # Unique, newick-safe labels (NCBI names like "Bacillus sp." repeat).
    counts = defaultdict(int)
    for r in rows:
        r["label"] = clean_label(r["label"]) or r["accession"]
        counts[r["label"]] += 1
    for r in rows:
        if counts[r["label"]] > 1:
            r["label"] = f"{r['label']} {r['accession']}"
    with open(run.genomes_tsv, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r[h] for h in header) + "\n")

    ref = cfg.get("seed_reference")
    if ref:
        print(f"Seed reference {ref['accession']}")
        base, _ = fetch_one(ref["accession"], {"_genomic.gff.gz": run.reference,
                                               "_protein.faa.gz": run.reference})
        if base is None:
            failed.append(ref["accession"])
    print(f"\nDone: {len(rows) - len([r for r in rows if r['role'] == 'failed'])} genomes downloaded, "
          f"{len(failed)} failed {failed if failed else ''}")
    print(f"Next: python flagpipe/flagpipe.py tree {run.root}  and  flagpipe.py check {run.root}")


# ── tree ────────────────────────────────────────────────────────────────────
def tree(run_dir):
    """GTDB tree pruned to the downloaded genomes: one copy labelled with
    organism names (interactive tree) and one with accessions (IQ-TREE
    concordance and the tanglegram). Same leaf order in both."""
    from ete3 import Tree
    run = Run(run_dir)
    cfg = run.config()
    genomes = run.genomes()
    label = {g["accession"]: g["label"] for g in genomes}
    full = Tree(str(gtdb_tree_path(cfg)), format=1, quoted_node_names=True)
    keep = [l.name for l in full.iter_leaves() if l.name[3:] in label]
    print(f"{len(keep)} of {len(genomes)} genomes are in the GTDB tree")
    full.prune(keep, preserve_branch_length=True)
    for l in full.iter_leaves():
        l.name = l.name[3:]
    run.tree.mkdir(parents=True, exist_ok=True)
    full.write(outfile=str(run.tree / "gtdb_acc.nwk"), format=5)
    for l in full.iter_leaves():
        l.name = label[l.name]
    full.write(outfile=str(run.tree / "gtdb_named.nwk"), format=5)
    back = set(Tree(str(run.tree / "gtdb_named.nwk"), format=1, quoted_node_names=True).get_leaf_names())
    broken = sorted(set(label[a] for a in label if a in {l[3:] for l in keep}) - back)
    if broken:
        sys.exit(f"Labels changed when written to newick: {broken} -- fix them in genomes.tsv")
    print(f"Wrote {run.tree}/gtdb_acc.nwk and gtdb_named.nwk")
