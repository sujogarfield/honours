"""
reliability.py
Which genomes have a fully known flagellar gene order, shared by every
script that uses the "reliable" set (interactive_tree.py,
synteny_analysis_full_contigs.py --closed) so the definitions can't drift.

Reliable = every flagellar gene on one closed, circular chromosome: a
natively complete assembly, or a RagTag scaffold against a closed
same-species reference (gene_order_combined_scaffolded.json). Paths are
relative to the project root.
"""
import json
from pathlib import Path

FASTA_DIR = Path("campy_fetched/campy_fas")
RAGTAG_DIR = Path("campy_fetched/ragtag_scaffolding/output")
SCAFFOLDED_GENE_ORDER = Path("gene_order/gene_order_combined_scaffolded.json")


def _sequence_lengths(path: Path, keep) -> dict:
    """{seq_id: length} for the FASTA records whose header line passes keep()."""
    lengths = {}
    current = None
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                current = line[1:].split()[0] if keep(line) else None
                if current:
                    lengths[current] = 0
            elif current:
                lengths[current] += len(line.strip())
    return lengths


def load_complete_chromosomes(fasta_dir: Path, ragtag_dir: Path) -> dict:
    """{contig_id: length} for every chromosome whose flagellar gene order is
    fully known:
    - natively closed: the FASTA header calls it a complete genome/sequence
      (and not a plasmid). Single-sequence WGS scaffolds ("whole genome
      shotgun sequence") are deliberately excluded: one piece, but not
      confirmed closed.
    - RagTag scaffolds (*_RagTag): draft contigs ordered against a closed
      same-species reference, so they span the reference chromosome. Only
      reliable for a genome if all its flagellar genes landed on one of them
      (see is_reliable); flagged ⚡ in the tree."""
    lengths = {}
    if not fasta_dir.exists():
        print(f"Warning: {fasta_dir} not found — no genome will count as reliable")
        return lengths

    def native(header):
        h = header.lower()
        return ("complete genome" in h or "complete sequence" in h) and "plasmid" not in h

    for path in sorted(fasta_dir.glob("*.fna")):
        lengths.update(_sequence_lengths(path, native))
    for path in sorted(ragtag_dir.glob("*/ragtag_out/ragtag.scaffold.fasta")):
        lengths.update(_sequence_lengths(path, lambda h: h.split()[0].endswith("_RagTag")))
    return lengths


def apply_ragtag(gene_order: dict, path: Path) -> set:
    """Swap in the RagTag-scaffolded gene records for the genomes RagTag was
    run on, in place. Returns the set of organisms swapped."""
    if not path.exists():
        print(f"Warning: {path} not found — RagTag-scaffolded gene order omitted")
        return set()
    with open(path) as f:
        raw = json.load(f)
    swapped = set()
    for data in raw.values():
        if data.get("scaffolded") and data["genes"]:
            gene_order[data["organism"]] = data["genes"]
            swapped.add(data["organism"])
    return swapped


def is_reliable(records: list[dict], complete: dict) -> bool:
    """Gene order is fully known only if every flagellar gene sits on one
    closed chromosome: natively assembled, or a RagTag scaffold against a
    same-species reference (see load_complete_chromosomes)."""
    contigs = {r["contig"] for r in records}
    return len(contigs) == 1 and next(iter(contigs)) in complete
