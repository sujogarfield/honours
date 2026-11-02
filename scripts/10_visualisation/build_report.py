#!/usr/bin/env python3
"""
build_report.py
The thesis frontend as one page, in the same format as the flagpipe report
(flagpipe/fp/report.py), built from the thesis results:

  Tree                  interactive_tree.py with its thesis inputs (RagTag
                        scaffolds, merged agreement scores, motility dots)
  Synteny vs phylogeny  all genomes, RagTag-scaffolded (synteny_analysis_combined_scaffolded.json)
                        vs reliable set (synteny_analysis_combined_reliable.json):
                        the same two analyses the tree's agreement scores come from
  Flagellar vs GTDB     thesis tanglegram (figures/fig1_tanglegram.png),
                        compare_to_gtdb.py numbers and the AU-tested HGT candidates

Nothing in gene_order/ is changed. Run from the project root:
    python scripts/10_visualisation/build_report.py
Output: frontend/report.html
"""
import csv
import html
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "flagpipe"))
sys.path.insert(0, str(ROOT / "scripts" / "10_visualisation"))
sys.path.insert(0, str(ROOT / "scripts" / "common"))
from fp import report as fpr  # noqa: E402

TITLE = "Campylobacterota"
OUT = Path("frontend/report.html")
ALL_JSON = Path("gene_order/synteny_analysis_combined_scaffolded.json")
RELIABLE_JSON = Path("gene_order/synteny_analysis_combined_reliable.json")
FLAG = Path("gene_order/flagellar_trees")
TANGLEGRAM = FLAG / "figures" / "fig1_tanglegram.png"
AU_SUMMARY = FLAG / "au" / "au_summary.tsv"
# RESULTS.md section 4: the candidates passing the explicit shortlist rule
PRIMARY = {"poseidonibacter_roscoffensis", "helicobacter_saguini", "campylobacter_mucosalis",
           "campylobacter_pinnipediorum"}


def tree_page():
    """interactive_tree.py with its default (thesis) inputs, written to a temp file."""
    import interactive_tree as it
    tmp = OUT.with_name("_tree_tmp.html")
    argv = sys.argv
    sys.argv = ["interactive_tree.py", "images_and_newick/thesis_b/gtdb_ref_tree.nwk",
                "images_and_newick/thesis_b/gtdb_ref_tree_og.nwk", str(tmp),
                "gene_order/gene_order_combined.json", "gene_order/synteny_analysis_combined_merged.json",
                "gene_order/motility.json"]
    try:
        it.main()
    finally:
        sys.argv = argv
    page = tmp.read_text()
    tmp.unlink()
    return page


def au_html():
    e = html.escape
    rows = list(csv.DictReader(open(AU_SUMMARY), delimiter="\t"))
    rows.sort(key=lambda r: (r["candidate"] not in PRIMARY, r["candidate"]))
    head = ("<tr><th>Genomes grouped by the flagellar genes</th><th></th><th>Concat p-AU (GTDB)</th>"
            "<th>Holm</th><th>&Delta;logL</th><th>Genes rejecting GTDB (p &lt; 0.05 / BH)</th></tr>")
    body = "".join(
        f"<tr><td>{e(r['group'])}</td><td>{'primary' if r['candidate'] in PRIMARY else 'exploratory'}</td>"
        f"<td>{e(r['concat_p_AU_gtdb'])}</td><td>{e(r['concat_p_holm'])}</td>"
        f"<td>{float(r['concat_delta_logL']):.0f}</td>"
        f"<td>{r['genes_reject_gtdb_p05']} / {r['genes_reject_gtdb_BH05']} of {r['genes_done'].split('/')[1]}</td></tr>"
        for r in rows)
    return ("<h2>HGT candidates (AU tests)</h2><p class=\"note\">Places where the flagellar genes "
            "group genomes against GTDB. Primary = passes the shortlist rule (in &ge; 10 gene trees and "
            "the concatenated tree; one-sided; GTDB branch not unusually short); exploratory = fails one "
            "criterion, tested for a stated reason. AU test: does the flagellar alignment reject the GTDB "
            "arrangement? p &lt; 0.05 = yes. Candidates were chosen from the same trees, so p-values are "
            "post-selection.</p>"
            f"<table class=\"grid\">{head}{body}</table>")


def main():
    scatter = fpr.scatter_page(
        json.loads(ALL_JSON.read_text()), json.loads(RELIABLE_JSON.read_text()), TITLE,
        "All genomes = main flagellar contig, RagTag-scaffolded where a same-species reference allowed, "
        "compared as linear. Reliable = every flagellar gene on one closed chromosome (complete assembly "
        "or RagTag scaffold), compared as circular.")
    cmp = (fpr.tanglegram_html(TANGLEGRAM.read_bytes(), fpr.TANGLE_NOTE)
           + fpr.comparison_html(FLAG / "comparison.json") + au_html())
    OUT.write_text(fpr.build_page(TITLE, fpr.tree_defaults(tree_page()), scatter, cmp))
    print(f"Wrote {OUT} ({OUT.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
