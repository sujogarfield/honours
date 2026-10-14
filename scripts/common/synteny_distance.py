"""
synteny_distance.py
Strand- and overlap-aware gene-order distance, the Thesis C default for
orthology_synteny_preview.py and synteny_analysis_full_contigs.py --closed
(both keep the old unsigned breakpoint distance behind --breakpoint).

Signed adjacencies: each gene has two ends, tail and head, and an adjacency
is the pair of gene ENDS that face each other along the chromosome (gene on
'+' reads tail->head left to right, on '-' head->tail). So
  fliA+ fliM+  and  fliA+ fliM-
share the neighbour pair but not the adjacency (fliM is inverted), while
reversing a whole chromosome (or block) leaves its internal adjacencies
unchanged -- an inverted block breaks only the adjacencies at its two ends,
and the arbitrary orientation of an assembly does not matter.

Overlaps: a conserved adjacency counts fully only if both genomes agree on
whether the two genes overlap (coordinates intersect); if one overlaps and
the other doesn't, it counts as half conserved. Overlap is read from the
annotated coordinates, so it depends on start-codon calls.

distance = 1 - conserved / possible, over genes shared by both genomes
(possible = n - 1 adjacencies for a linear order, n for a circular one).
"""
from collections import Counter

MIN_SHARED_GENES = 4


def primary_contig_records(records: list[dict]) -> tuple[list[dict], int]:
    """Gene records on the contig with the most flagellar genes, sorted by
    position, with duplicated gene tags dropped (ambiguous order mapping) --
    the same selection as primary_contig_order(), but keeping strand and
    coordinates. Returns (records, n_genes_on_primary_contig)."""
    by_contig = {}
    for r in records:
        by_contig.setdefault(r["contig"], []).append(r)
    if not by_contig:
        return [], 0
    primary = max(by_contig.values(), key=len)
    counts = Counter(r["gene"].lower() for r in primary)
    ordered = sorted((r for r in primary if counts[r["gene"].lower()] == 1),
                     key=lambda r: r["start"])
    slim = [{"gene": r["gene"].lower(), "strand": r.get("strand", "+"),
             "start": r["start"], "end": r["end"]} for r in ordered]
    return slim, len(primary)


def _ends(rec: dict) -> tuple:
    """(left end, right end) of a gene as drawn left to right."""
    g = rec["gene"]
    return ((g, "tail"), (g, "head")) if rec["strand"] == "+" else ((g, "head"), (g, "tail"))


def _overlap(a: dict, b: dict) -> bool:
    return a["start"] <= b["end"] and b["start"] <= a["end"]


def _adjacencies(recs: list[dict], circular: bool) -> dict:
    """{frozenset(facing gene ends): genes_overlap} for consecutive genes."""
    n = len(recs)
    pairs = range(n) if circular else range(n - 1)
    out = {}
    for i in pairs:
        a, b = recs[i], recs[(i + 1) % n]
        out[frozenset((_ends(a)[1], _ends(b)[0]))] = _overlap(a, b)
    return out


def signed_overlap_distance(recs_a: list[dict], recs_b: list[dict], circular: bool = False,
                            min_shared: int = MIN_SHARED_GENES):
    """Normalised signed-adjacency distance with overlap half-credit (see
    module docstring). recs_*: primary_contig_records() output. Returns
    (distance_or_None, n_shared_genes)."""
    shared = {r["gene"] for r in recs_a} & {r["gene"] for r in recs_b}
    n = len(shared)
    if n < min_shared:
        return None, n
    adj_a = _adjacencies([r for r in recs_a if r["gene"] in shared], circular)
    adj_b = _adjacencies([r for r in recs_b if r["gene"] in shared], circular)
    conserved = sum(1.0 if adj_a[k] == adj_b[k] else 0.5 for k in adj_a.keys() & adj_b.keys())
    possible = n if circular else n - 1
    return 1 - conserved / possible, n
