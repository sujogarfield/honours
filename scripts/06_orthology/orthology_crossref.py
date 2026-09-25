#!/usr/bin/env python3
"""
orthology_crossref.py
Cross-references the existing string-matched flagellar gene set
(gene_order/gene_order.json) against OrthoFinder orthogroups to find
genes present in a genome's proteome but missed by the gene= string match
in campy_fetched/campy_ann/*.gff (e.g. missing/nonstandard gene tag, PGAP naming drift).

For each of the 46 known flagellar gene names:
  1. Find which protein(s) the string match already identified as that gene,
     across all genomes (via GFF CDS coordinates -> protein_id).
  2. Look up which OrthoFinder orthogroup(s) those proteins fall into.
  3. Pull every protein OrthoFinder placed in that same orthogroup, across
     all 166 genomes.
  4. Anything not already covered by the string match is a candidate --
     map its protein_id back to genomic coordinates via that genome's GFF.

Output: gene_order/gene_order_orthology_candidates.json, same record shape
as gene_order.json plus source/orthogroup/protein_id metadata.

Run AFTER orthofinder has produced orthofinder_results/*/Orthogroups/Orthogroups.tsv
and AFTER gene_order.py has produced gene_order/gene_order.json.
"""

import csv
import glob
import json
import os
import re
import sys
from collections import defaultdict
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))

GFF_DIR = "campy_fetched/campy_ann"
GENE_ORDER_JSON = "campy_orthologs/archive_json/gene_order.json"
OUTPUT_JSON = "campy_orthologs/archive_json/gene_order_orthology_candidates.json"

ORTHOGROUPS_TSV = sorted(
    glob.glob("campy_orthologs/orthofinder_results/Results_*/Orthogroups/Orthogroups.tsv")
)
if not ORTHOGROUPS_TSV:
    sys.exit("No Orthogroups.tsv found under campy_orthologs/orthofinder_results/ -- run OrthoFinder first.")
ORTHOGROUPS_TSV = ORTHOGROUPS_TSV[-1]

PROTEIN_ID_RE = re.compile(r"Dbxref=[^;]*GenBank:([^,;]+)")


# ── GFF locus <-> protein_id index, per genome ──────────────────────────────
def build_gff_index(gff_path: str):
    """Returns (locus_to_protein, protein_to_locus) for one genome's GFF,
    keyed on CDS (contig, start, end, strand) <-> protein_id."""
    locus_to_protein = {}
    protein_to_locus = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            contig, start, end, strand, attrs = (
                parts[0], int(parts[3]), int(parts[4]), parts[6], parts[8]
            )
            m = PROTEIN_ID_RE.search(attrs)
            if not m:
                continue
            protein_id = m.group(1)
            locus = (contig, start, end, strand)
            locus_to_protein[locus] = protein_id
            protein_to_locus[protein_id] = locus
    return locus_to_protein, protein_to_locus


def main():
    print(f"Using orthogroups file: {ORTHOGROUPS_TSV}")

    with open(GENE_ORDER_JSON) as f:
        gene_order = json.load(f)  # gcf_id -> {organism, genes:[...]}

    gff_bases = {
        os.path.basename(p)[:-4] for p in glob.glob(os.path.join(GFF_DIR, "*.gff"))
    }

    # ── Build per-genome GFF indexes (only for genomes present in gene_order.json
    #    plus a fallback for every genome, since candidates can come from any of them) ──
    print("Indexing GFF CDS <-> protein_id for all genomes ...")
    locus_idx = {}   # gcf_id -> {locus: protein_id}
    protein_idx = {} # gcf_id -> {protein_id: locus}
    for base in sorted(gff_bases):
        gff_path = os.path.join(GFF_DIR, f"{base}.gff")
        l2p, p2l = build_gff_index(gff_path)
        locus_idx[base] = l2p
        protein_idx[base] = p2l
    print(f"Indexed {len(gff_bases)} genomes")

    # ── For each string-matched record, resolve its protein_id ────────────
    # known_protein_to_gene[gcf_id][protein_id] = gene_name  (already found)
    known_protein_to_gene = defaultdict(dict)
    gene_name_to_proteins = defaultdict(set)  # gene_name -> {protein_id, ...}
    unresolved = 0
    for gcf_id, data in gene_order.items():
        for rec in data["genes"]:
            locus = (rec["contig"], rec["start"], rec["end"], rec["strand"])
            protein_id = locus_idx.get(gcf_id, {}).get(locus)
            if protein_id is None:
                unresolved += 1
                continue
            gene = rec["gene"].lower()
            known_protein_to_gene[gcf_id][protein_id] = gene
            gene_name_to_proteins[gene].add(protein_id)
    print(f"Resolved protein_id for string-matched genes "
          f"({unresolved} records had no matching CDS/protein_id)")

    # ── Load Orthogroups.tsv: header gives genome column order ────────────
    print("Loading Orthogroups.tsv ...")
    with open(ORTHOGROUPS_TSV) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        species_cols = header[1:]  # column name == gff/faa basename

        protein_to_og = {}          # protein_id -> orthogroup id
        og_members = defaultdict(lambda: defaultdict(list))  # og -> gcf_id -> [protein_id,...]

        for row in reader:
            if not row:
                continue
            og_id = row[0]
            for species, cell in zip(species_cols, row[1:]):
                if not cell.strip():
                    continue
                prot_ids = [p.strip() for p in cell.split(",") if p.strip()]
                og_members[og_id][species].extend(prot_ids)
                for p in prot_ids:
                    protein_to_og[p] = og_id
    print(f"Loaded {len(og_members)} orthogroups, {len(protein_to_og)} proteins")

    # ── gene_name -> set of orthogroups it seeds ───────────────────────────
    gene_to_ogs = defaultdict(set)
    for gene, proteins in gene_name_to_proteins.items():
        for p in proteins:
            og = protein_to_og.get(p)
            if og:
                gene_to_ogs[gene].add(og)

    print("\nSeed orthogroups per gene:")
    for gene in sorted(gene_to_ogs):
        ogs = gene_to_ogs[gene]
        flag = "  <-- split across multiple OGs, inspect" if len(ogs) > 1 else ""
        print(f"  {gene:6s}: {sorted(ogs)}{flag}")

    # ── Find candidates: OG members not already known for that gene in that genome ──
    candidates = defaultdict(list)  # gcf_id -> [records]
    n_candidates = 0
    for gene, ogs in gene_to_ogs.items():
        for og in ogs:
            for gcf_id, prot_ids in og_members[og].items():
                if gcf_id not in gff_bases:
                    continue
                already = known_protein_to_gene.get(gcf_id, {})
                for protein_id in prot_ids:
                    if protein_id in already:
                        continue  # string match already found this exact protein
                    locus = protein_idx.get(gcf_id, {}).get(protein_id)
                    if locus is None:
                        continue  # shouldn't happen, but be defensive
                    contig, start, end, strand = locus
                    candidates[gcf_id].append({
                        "gene": gene,
                        "contig": contig,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "source": "orthology",
                        "orthogroup": og,
                        "protein_id": protein_id,
                    })
                    n_candidates += 1

    print(f"\n{n_candidates} candidate genes found via orthology "
          f"but missed by string matching, across {len(candidates)} genomes")

    # ── Attach organism name (reuse gene_order.json's when available) ─────
    organism_by_gcf = {gcf: d["organism"] for gcf, d in gene_order.items()}
    from helpers import get_organism_from_fasta
    out = {}
    for gcf_id, recs in candidates.items():
        organism = organism_by_gcf.get(gcf_id)
        if organism is None:
            fasta_path = os.path.join("campy_fetched/campy_fas", f"{gcf_id}.fna")
            organism = get_organism_from_fasta(fasta_path) if os.path.exists(fasta_path) else gcf_id
        for r in recs:
            r["organism"] = organism
        out[gcf_id] = {"organism": organism, "genes": recs}

    os.makedirs("gene_order", exist_ok=True)
    with open(OUTPUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved candidates to {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
