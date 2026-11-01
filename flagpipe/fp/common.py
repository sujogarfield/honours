"""
common.py
Shared helpers for every flagpipe step: the run directory layout, config
loading, FASTA/GFF readers and the GFF CDS <-> protein_id index used to link
string-matched genes (GFF coordinates) to OrthoFinder (protein IDs).

Every step takes the run directory as its only argument; all paths below are
relative to it, so a run can be moved or copied as one folder.
"""
import json
import os
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]          # honours repo root
SCRIPTS = REPO / "scripts"
FLAGPIPE = Path(__file__).resolve().parents[1]

# GTDB release the thesis used. The tree is downloaded once into a shared cache.
GTDB_RELEASE = 232
GTDB_TREE_URL = "https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_r232.tree"


class Run:
    """Paths inside one run directory."""

    def __init__(self, root):
        self.root = Path(root).resolve()
        self.config_path = self.root / "config.json"
        self.genomes_tsv = self.root / "genomes.tsv"
        self.fna = self.root / "genomes" / "fna"
        self.gff = self.root / "genomes" / "gff"
        self.faa = self.root / "genomes" / "faa"
        self.reference = self.root / "reference"
        self.tree = self.root / "tree"
        self.orthology = self.root / "orthology"
        self.trees = self.root / "trees"
        self.synteny = self.root / "synteny"
        self.report = self.root / "report"

    def config(self):
        if not self.config_path.exists():
            sys.exit(f"No config at {self.config_path} -- run `flagpipe.py init {self.root}` first.")
        return json.loads(self.config_path.read_text())

    def genomes(self, roles=("ingroup", "outgroup")):
        """[{accession, base, label, taxonomy, role}] from genomes.tsv."""
        if not self.genomes_tsv.exists():
            sys.exit(f"No {self.genomes_tsv} -- run `flagpipe.py genomes` first.")
        lines = self.genomes_tsv.read_text().splitlines()
        header = lines[0].split("\t")
        rows = [dict(zip(header, l.split("\t"))) for l in lines[1:] if l.strip()]
        return [r for r in rows if r["role"] in roles]

    def orthogroups_tsv(self):
        hits = sorted(self.orthology.glob("orthofinder/Results_*/Orthogroups/Orthogroups.tsv"))
        if not hits:
            sys.exit(f"No Orthogroups.tsv under {self.orthology}/orthofinder -- run OrthoFinder first.")
        return hits[-1]


def accession(base):
    """GCF_000012965.1_ASM1296v1 -> GCF_000012965.1"""
    return "_".join(base.split("_")[:2])


def clean_label(name):
    """Organism name safe to use as a newick leaf label: letters, digits,
    spaces, '.', '_' and '-' only (ete3 rewrites others, e.g. '=' -> '_',
    which would break the label -> tree leaf link)."""
    name = re.sub(r"\s*\([^)]*\)", "", name)  # NCBI suffixes like "(firmicutes)"
    return re.sub(r"\s+", " ", re.sub(r"[^\w .\-]", "", name)).strip()


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


PROTEIN_ID_RE = re.compile(r"Dbxref=[^;]*GenBank:([^,;]+)")


def build_gff_index(gff_path):
    """(locus -> protein_id, protein_id -> locus) for one genome's GFF, keyed
    on CDS (contig, start, end, strand). Same as orthology_crossref.py."""
    locus_to_protein, protein_to_locus = {}, {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            m = PROTEIN_ID_RE.search(parts[8])
            if not m:
                continue
            locus = (parts[0], int(parts[3]), int(parts[4]), parts[6])
            locus_to_protein[locus] = m.group(1)
            protein_to_locus[m.group(1)] = locus
    return locus_to_protein, protein_to_locus


def load_orthogroups(path):
    """(protein -> OG, OG -> {genome_base: [protein, ...]})"""
    import csv
    from collections import defaultdict
    p2og, members = {}, defaultdict(lambda: defaultdict(list))
    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        genomes = next(reader)[1:]
        for row in reader:
            for genome, cell in zip(genomes, row[1:]):
                for p in cell.split(","):
                    p = p.strip()
                    if p:
                        p2og[p] = row[0]
                        members[row[0]][genome].append(p)
    return p2og, members


def load_eggnog(path):
    """protein_id -> {preferred_name, description} from an .emapper.annotations file."""
    ann, header = {}, None
    with open(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            ann[row["query"]] = {"preferred_name": row.get("Preferred_name", "-"),
                                 "description": row.get("Description", "-")}
    return ann


def write_json(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def add_scripts_path(*folders):
    """Make thesis modules (scripts/common, scripts/10_visualisation, ...) importable."""
    for folder in folders:
        p = str(SCRIPTS / folder)
        if p not in sys.path:
            sys.path.insert(0, p)


def gtdb_tree_path(cfg):
    """The cached GTDB bac120 tree, downloaded on first use (login node: needs internet)."""
    p = Path(os.path.expandvars(cfg.get("gtdb_tree", ""))) if cfg.get("gtdb_tree") else None
    if p and p.exists():
        return p
    thesis_copy = SCRIPTS / "02_species_tree" / "bac120.tree"
    if thesis_copy.exists():
        return thesis_copy
    cache = FLAGPIPE / "cache" / "bac120_r232.tree"
    if not cache.exists():
        import urllib.request
        cache.parent.mkdir(parents=True, exist_ok=True)
        print(f"Downloading GTDB R{GTDB_RELEASE} bac120 tree -> {cache}")
        urllib.request.urlretrieve(GTDB_TREE_URL, cache)
    return cache
