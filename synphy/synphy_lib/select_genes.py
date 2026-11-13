"""
select_genes.py
The thesis copy-resolution and gene-selection rules
(scripts/11_flagellar_trees/prepare_gene_sets.py), with genome coverage as a
fraction so it works for any number of genomes:

  - each gene's main orthogroup = the one most of its proteins fall in;
  - 2+ same-named copies: keep the one in the main orthogroup;
  - 2+ copies IN the main orthogroup: ambiguous, genome left out of that gene;
  - a lone copy in another orthogroup: kept only if that orthogroup never
    co-occurs with the main one in any genome (a lineage-specific version),
    otherwise dropped as a likely paralog;
  - a gene is selected if present in >= min_genome_frac of genomes and
    ambiguous in <= max_multicopy_frac of them;
  - genomes carrying < min_gene_frac of the selected genes are dropped.

Output (trees/): genes/<gene>.faa, genes.txt, taxa.tsv, gtdb_ref_pruned.nwk,
gene_selection.json -- the same inputs the thesis Katana tree jobs used.
"""
import json
from collections import Counter, defaultdict

from synphy_lib.common import Run, accession, build_gff_index, load_orthogroups, read_fasta, write_json


def resolve_copies(gene_to_proteins, p2og, og_genomes):
    out = {}
    for gene, gp in gene_to_proteins.items():
        og_counts = Counter(p2og.get(p) for pids in gp.values() for p in set(pids))
        og_counts.pop(None, None)
        main_og = og_counts.most_common(1)[0][0] if og_counts else None
        calls, notes = {}, {}
        for genome, pids in gp.items():
            pids = sorted(set(pids))
            in_main = [p for p in pids if p2og.get(p) == main_og]
            if len(in_main) == 1:
                calls[genome] = in_main[0]
                if len(pids) > 1:
                    notes[genome] = "multi-copy, resolved to main-OG copy"
            elif len(in_main) > 1:
                calls[genome] = None
                notes[genome] = "2+ copies in main OG (ambiguous), excluded"
            elif len(pids) == 1:
                other = p2og.get(pids[0])
                shared = len(og_genomes.get(other, set()) & og_genomes.get(main_og, set()))
                if other and shared == 0:
                    calls[genome] = pids[0]
                    notes[genome] = f"single copy in lineage-specific OG ({other}), kept"
                else:
                    notes[genome] = (f"single copy in paralog OG ({other}, co-occurs with main OG "
                                     f"in {shared} genomes), dropped")
            else:
                calls[genome] = None
                notes[genome] = "2+ copies, none in main OG (ambiguous), excluded"
        out[gene] = (main_og, calls, notes)
    return out


def select(run_dir):
    from ete3 import Tree
    run = Run(run_dir)
    cfg = run.config()
    rules = cfg["selection"]
    out_dir = run.trees
    (out_dir / "genes").mkdir(parents=True, exist_ok=True)
    for old in (out_dir / "genes").glob("*.faa"):
        old.unlink()
    combined = json.loads((run.root / "gene_order_combined.json").read_text())

    gene_to_proteins = defaultdict(lambda: defaultdict(list))
    unresolved = Counter()
    for genome, data in combined.items():
        locus_idx = build_gff_index(run.gff / f"{genome}.gff")[0]
        for g in data["genes"]:
            pid = g.get("protein_id") or locus_idx.get((g["contig"], g["start"], g["end"], g["strand"]))
            if pid:
                gene_to_proteins[g["gene"]][genome].append(pid)
            else:
                unresolved[g["gene"]] += 1

    p2og, members = load_orthogroups(run.orthogroups_tsv())
    og_genomes = {og: set(m) for og, m in members.items()}
    resolved = resolve_copies(gene_to_proteins, p2og, og_genomes)

    n_genomes = len(run.genomes())
    min_cov = int(round(rules["min_genome_frac"] * n_genomes))
    selected, report = [], {}
    for gene, (main_og, calls_g, notes) in sorted(resolved.items()):
        n_multi = sum(1 for v in calls_g.values() if v is None)
        if len(calls_g) < min_cov or n_multi / len(calls_g) > rules["max_multicopy_frac"]:
            continue
        selected.append(gene)
        report[gene] = {"main_orthogroup": main_og, "genomes_with_gene": len(calls_g),
                        "genomes_multicopy_excluded": n_multi,
                        "copy_resolution": {accession(g): r for g, r in sorted(notes.items())}}
    print(f"Selected {len(selected)} genes (in >= {min_cov}/{n_genomes} genomes, "
          f"<= {rules['max_multicopy_frac']:.0%} ambiguous): {selected}")
    if len(selected) < 4:
        print("WARNING: fewer than 4 genes selected -- the concatenated tree and gCF will be weak. "
              "Check `synphy.py check` and the seed reference.")

    calls = {gene: {gn: pid for gn, pid in resolved[gene][1].items() if pid} for gene in selected}
    per_genome = Counter(gn for gene in selected for gn in calls[gene])
    min_genes = int(round(rules["min_gene_frac"] * len(selected)))

    gtdb = Tree(str(run.tree / "gtdb_acc.nwk"), format=1)
    gtdb_taxa = set(gtdb.get_leaf_names())
    included = sorted(gn for gn, n in per_genome.items() if n >= min_genes and accession(gn) in gtdb_taxa)
    dropped = sorted(set(combined) - set(included))
    print(f"Keeping {len(included)} genomes with >= {min_genes}/{len(selected)} genes; dropping {len(dropped)}")

    prot_cache = {}
    for gene in selected:
        n = 0
        with open(out_dir / "genes" / f"{gene}.faa", "w") as out:
            for gn in included:
                pid = calls[gene].get(gn)
                if not pid:
                    continue
                if gn not in prot_cache:
                    prot_cache[gn] = read_fasta(run.faa / f"{gn}.faa")
                seq = prot_cache[gn].get(pid)
                if seq:
                    out.write(f">{accession(gn)}\n{seq.rstrip('*')}\n")
                    n += 1
        report[gene]["sequences_written"] = n
        print(f"  {gene}: {n} sequences")

    (out_dir / "genes.txt").write_text("\n".join(selected) + "\n")
    with open(out_dir / "taxa.tsv", "w") as f:
        f.write("accession\torganism\n")
        for gn in included:
            f.write(f"{accession(gn)}\t{combined[gn]['organism']}\n")
    gtdb.prune([accession(gn) for gn in included], preserve_branch_length=True)
    gtdb.write(outfile=str(out_dir / "gtdb_ref_pruned.nwk"), format=5)
    write_json(out_dir / "gene_selection.json", {
        "rule": {**rules, "min_genome_coverage": min_cov, "min_genes_per_genome": min_genes},
        "genes": report,
        "genes_not_selected": {g: {"genomes": len(r[1]),
                                   "multicopy_in_main_og": sum(1 for v in r[1].values() if v is None)}
                               for g, r in sorted(resolved.items()) if g not in selected},
        "n_genomes_included": len(included),
        "genomes_dropped": {accession(gn): {"organism": combined[gn]["organism"],
                                            "selected_genes_carried": per_genome.get(gn, 0)}
                            for gn in dropped},
        "gene_calls_without_protein_id": dict(unresolved),
    })
    print(f"Wrote {out_dir}/ (genes/, genes.txt, taxa.tsv, gtdb_ref_pruned.nwk, gene_selection.json)")
