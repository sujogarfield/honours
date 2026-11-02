"""
report.py
One page, report/report.html, with three tabs:

  Tree                    pruned GTDB tree, gene track and synteny agreement dot
                          per genome (scripts/10_visualisation/interactive_tree.py,
                          pointed at this run)
  Synteny vs phylogeny    synteny vs patristic distance, one point per genome pair,
                          all genomes / reliable toggle, Mantel r
  Flagellar vs GTDB       tanglegram (GTDB tree shaded by gCF vs concatenated tree),
                          comparison numbers and the HGT candidate shortlist
                          (filled in once the Katana tree stage has finished)

The tanglegram is also saved as report/tanglegram.png (for figures). `report`
builds whatever its inputs exist for, so it is run after each stage.
"""
import base64
import html as html_lib
import json
import os
import sys
from collections import Counter

from fp.common import REPO, Run, add_scripts_path


# ── interactive tree ────────────────────────────────────────────────────────
def interactive_tree(run):
    add_scripts_path("common", "10_visualisation")
    import interactive_tree as it
    # Point the thesis script at this run; no RagTag here (missing paths are skipped with a warning).
    it.FASTA_DIR = run.fna
    it.RAGTAG_DIR = run.root / "_no_ragtag"
    it.SCAFFOLDED_GENE_ORDER = run.root / "_no_ragtag.json"
    it.MOTILITY_ICON_DIR = REPO / "frontend" / "motility_icons"
    out = run.report / "_tree.html"
    argv = sys.argv
    sys.argv = ["interactive_tree.py", str(run.tree / "gtdb_named.nwk"), str(run.tree / "gtdb_acc.nwk"),
                str(out), str(run.root / "gene_order_combined.json"),
                str(run.synteny / "synteny_tree.json"), str(run.root / "_no_motility.json")]
    try:
        it.main()
    finally:
        sys.argv = argv
    page = out.read_text()
    out.unlink()
    return page


# ── synteny scatter ─────────────────────────────────────────────────────────
SCATTER_HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Synteny vs Phylogeny</title>
<style>
  /* always white, like the tree and tanglegram, so the panels match */
  :root { color-scheme:light; --bg:#fff; --panel:#fff; --text:#0b0b0b; --muted:#52514e; --grid:#e1e0d9; --dot:rgba(42,120,214,0.35); --trend:#eb6834; }
  body { margin:0; font-family:system-ui,sans-serif; background:var(--bg); color:var(--text); }
  header { padding:12px 16px; background:var(--panel); border-bottom:1px solid var(--grid); display:flex; flex-wrap:wrap; gap:12px; align-items:center; }
  h1 { font-size:16px; margin:0; }
  .hint { color:var(--muted); font-size:12px; padding:8px 16px 0; max-width:900px; }
  #stats { font-size:13px; margin-left:auto; }
  #wrap { position:relative; margin:8px 16px 16px; }
  canvas, svg { position:absolute; left:0; top:0; }
  .axis { font-size:11px; fill:var(--muted); }
  .gridline { stroke:var(--grid); }
  #tip { position:fixed; background:rgba(11,11,11,0.9); color:#fff; padding:6px 10px; border-radius:4px; font-size:12px; pointer-events:none; opacity:0; }
</style>
</head>
<body>
<header>
  <h1>__TITLE__: flagellar gene order vs Phylogeny</h1>
  <label><input type="checkbox" id="reliable" __RELIABLE_CHECKED__> Reliable gene order only (<span id="nrel"></span> / <span id="nall"></span> genomes)</label>
  <span id="stats"></span>
</header>
<div class="hint">Each point is one genome pair. x = patristic distance on the GTDB tree; y = synteny distance (share of neighbouring flagellar-gene pairs not conserved, counting strand, half credit for an overlap mismatch). __RELIABLE_NOTE__ Hover for names.</div>
<div id="wrap"><canvas id="c"></canvas><svg id="s"></svg></div>
<div id="tip"></div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.9.0/d3.min.js"></script>
<script>
const DATA = { all: __ALL__, reliable: __RELIABLE__ };
const tip = document.getElementById("tip");
document.getElementById("nall").textContent = DATA.all.matrix.organisms.length;
document.getElementById("nrel").textContent = DATA.reliable.matrix.organisms.length;
function css(v) { return getComputedStyle(document.documentElement).getPropertyValue(v).trim(); }
function render() {
  const D = document.getElementById("reliable").checked ? DATA.reliable : DATA.all;
  const names = D.matrix.organisms, S = D.matrix.synteny_distance, P = D.matrix.phylogenetic_distance;
  const pts = [];
  for (let i = 0; i < names.length; i++) for (let j = i + 1; j < names.length; j++)
    if (S[i][j] !== null) pts.push({ a: names[i], b: names[j], s: S[i][j], p: P[i][j] });
  const W = Math.max(320, window.innerWidth - 32), H = Math.max(320, Math.min(window.innerHeight - 150, W * 0.75));
  const m = { t: 20, r: 20, b: 46, l: 54 };
  const c = document.getElementById("c"), dpr = window.devicePixelRatio || 1;
  c.width = W * dpr; c.height = H * dpr; c.style.width = W + "px"; c.style.height = H + "px";
  document.getElementById("wrap").style.height = H + "px";
  const ctx = c.getContext("2d"); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
  const x = d3.scaleLinear().domain([0, d3.max(pts, d => d.p) || 1]).nice().range([m.l, W - m.r]);
  const y = d3.scaleLinear().domain([0, 1]).range([H - m.b, m.t]);
  const svg = d3.select("#s").attr("width", W).attr("height", H); svg.selectAll("*").remove();
  svg.selectAll(".gx").data(x.ticks(6)).join("line").attr("class", "gridline").attr("x1", d => x(d)).attr("x2", d => x(d)).attr("y1", m.t).attr("y2", H - m.b);
  svg.selectAll(".gy").data(y.ticks(5)).join("line").attr("class", "gridline").attr("x1", m.l).attr("x2", W - m.r).attr("y1", d => y(d)).attr("y2", d => y(d));
  svg.selectAll(".tx").data(x.ticks(6)).join("text").attr("class", "axis").attr("x", d => x(d)).attr("y", H - m.b + 16).attr("text-anchor", "middle").text(d => d.toFixed(2));
  svg.selectAll(".ty").data(y.ticks(5)).join("text").attr("class", "axis").attr("x", m.l - 8).attr("y", d => y(d) + 3).attr("text-anchor", "end").text(d => d.toFixed(1));
  svg.append("text").attr("class", "axis").attr("x", (m.l + W - m.r) / 2).attr("y", H - 8).attr("text-anchor", "middle").text("Phylogenetic (patristic) distance");
  svg.append("text").attr("class", "axis").attr("transform", `translate(14,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Synteny distance");
  ctx.fillStyle = css("--dot");
  pts.forEach(p => { ctx.beginPath(); ctx.arc(x(p.p), y(p.s), 2.6, 0, 2 * Math.PI); ctx.fill(); });
  if (pts.length > 2) {
    const n = pts.length, sx = d3.sum(pts, d => d.p), sy = d3.sum(pts, d => d.s);
    const sxx = d3.sum(pts, d => d.p * d.p), sxy = d3.sum(pts, d => d.p * d.s);
    const k = (n * sxy - sx * sy) / (n * sxx - sx * sx), b0 = (sy - k * sx) / n, X = x.domain()[1];
    svg.append("line").attr("x1", x(0)).attr("y1", y(b0)).attr("x2", x(X)).attr("y2", y(b0 + k * X)).attr("stroke", css("--trend")).attr("stroke-width", 2);
  }
  const M = D.mantel;
  document.getElementById("stats").textContent = M.r === null ? "Mantel test: not enough pairs"
    : `Mantel r = ${M.r.toFixed(3)}, p = ${M.p_value.toFixed(3)} (${M.n_pairs} pairs, ${M.n_permutations} permutations)`;
  const q = d3.quadtree(pts, d => x(d.p), d => y(d.s)), hover = svg.append("g");
  c.onmousemove = e => {
    const r = c.getBoundingClientRect(), f = q.find(e.clientX - r.left, e.clientY - r.top, 20);
    hover.selectAll("*").remove();
    if (!f) { tip.style.opacity = 0; return; }
    hover.append("circle").attr("cx", x(f.p)).attr("cy", y(f.s)).attr("r", 5).attr("fill", "none").attr("stroke", css("--text"));
    tip.innerHTML = `<b>${f.a}</b> × <b>${f.b}</b><br>synteny ${f.s.toFixed(3)} · phylogenetic ${f.p.toFixed(4)}`;
    tip.style.left = (e.clientX + 14) + "px"; tip.style.top = (e.clientY - 10) + "px"; tip.style.opacity = 1;
  };
  c.onmouseleave = () => { hover.selectAll("*").remove(); tip.style.opacity = 0; };
}
document.getElementById("reliable").addEventListener("change", render);
window.addEventListener("resize", render);
render();
</script>
</body>
</html>
"""


RELIABLE_NOTE = "Reliable = every flagellar gene on one natively closed chromosome, compared as circular."


def scatter_page(all_data, reliable_data, title, reliable_note=RELIABLE_NOTE):
    """Scatter page from two synteny results ({matrix, mantel}): all genomes and reliable."""
    slim = lambda d: {"matrix": d["matrix"], "mantel": d["mantel"]}  # all the page needs
    has_reliable = len(reliable_data["matrix"]["organisms"]) >= 3
    return (SCATTER_HTML.replace("__TITLE__", title).replace("__RELIABLE_NOTE__", reliable_note)
            .replace("__RELIABLE_CHECKED__", "checked" if has_reliable else "")
            .replace("__ALL__", json.dumps(slim(all_data))).replace("__RELIABLE__", json.dumps(slim(reliable_data))))


def synteny_scatter(run, title):
    data = {k: json.loads((run.synteny / f"synteny_{k}.json").read_text()) for k in ("all", "reliable")}
    return scatter_page(data["all"], data["reliable"], title)


# ── tanglegram ──────────────────────────────────────────────────────────────
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#9c6b30"]


def tanglegram(run):
    """GTDB tree (branches shaded by gCF) vs the concatenated tree, genomes
    joined by lines coloured by GTDB family. Child order is rotated to reduce
    crossings; topology is checked to be unchanged before drawing."""
    os.environ["FLAG_DIR"] = str(run.trees)
    add_scripts_path("11_flagellar_trees")
    import numpy as np
    from ete3 import Tree
    import plot_results as pr  # reuse the thesis tree layout/drawing helpers
    from screen_hgt_candidates import bitmask_splits
    plt = pr.plt
    D = run.trees
    n_genes = len((D / "genes.txt").read_text().split())
    taxa_rows = [l.split("\t") for l in (D / "taxa.tsv").read_text().splitlines()[1:]]
    org = {a: o for a, o in taxa_rows}
    meta = {g["accession"]: g for g in run.genomes()}

    def family(acc):
        if meta.get(acc, {}).get("role") == "outgroup":
            return "Outgroup"
        tax = meta.get(acc, {}).get("taxonomy", "")
        return next((t[3:] for t in tax.split(";") if t.startswith("f__")), "?")

    gtdb = Tree(str(D / "gtdb_ref_pruned.nwk"), format=1)
    flag = Tree(str(D / "concat" / "concat.treefile"), format=0)
    taxa = sorted(gtdb.get_leaf_names())
    index = {t: i for i, t in enumerate(taxa)}
    fams = Counter(family(t) for t in taxa if family(t) != "Outgroup")
    groups = [f for f, _ in fams.most_common(len(PALETTE))]
    colour = {f: PALETTE[i] for i, f in enumerate(groups)}

    rows = [l.rstrip("\n").split("\t") for l in open(D / "concordance" / "gtdb_gcf.cf.stat")
            if not l.startswith("#") and l.strip()]
    col = {h: i for i, h in enumerate(rows[0])}
    gcf_by_id = {r[col["ID"]]: float(r[col["gCF"]]) for r in rows[1:] if r[col["gCF"]] != "NA"}
    cfb = Tree(str(D / "concordance" / "gtdb_gcf.cf.branch"), format=1)
    gcf_of_split = {m: gcf_by_id.get(bid) for m, bid in bitmask_splits(cfb, index).items()}
    full = (1 << len(taxa)) - 1
    canon = lambda m: m ^ full if m & 1 else m

    outgroup = [t for t in taxa if family(t) == "Outgroup"]
    if outgroup:
        anc = flag.get_common_ancestor(outgroup) if len(outgroup) > 1 else flag & outgroup[0]
        flag.set_outgroup(anc if anc is not flag else flag.get_midpoint_outgroup())
        g_anc = gtdb.get_common_ancestor(outgroup) if len(outgroup) > 1 else gtdb & outgroup[0]
        if g_anc is not gtdb:
            gtdb.set_outgroup(g_anc)
    else:
        flag.set_outgroup(flag.get_midpoint_outgroup())

    def leaf_rows(tree):
        return {name: -i for i, name in enumerate(tree.get_leaf_names())}

    def rotate_to(tree, target):
        for nd in tree.traverse("postorder"):
            if not nd.is_leaf():
                nd.children.sort(key=lambda c: -np.mean([target[l] for l in c.get_leaf_names()]))

    def crossings(a, b):
        seq = [-b[t] for t in sorted(taxa, key=lambda t: -a[t])]
        return sum(1 for i in range(len(seq)) for j in range(i + 1, len(seq)) if seq[i] > seq[j])

    src_splits = {"GTDB": set(bitmask_splits(gtdb.copy(), index)), "flag": set(bitmask_splits(flag.copy(), index))}
    gtdb.ladderize()
    best = None
    for _ in range(8):
        rotate_to(flag, leaf_rows(gtdb))
        c = crossings(leaf_rows(gtdb), leaf_rows(flag))
        if best is None or c < best[0]:
            best = (c, gtdb.write(format=5), flag.write(format=0))
        rotate_to(gtdb, leaf_rows(flag))
        c = crossings(leaf_rows(gtdb), leaf_rows(flag))
        if c < best[0]:
            best = (c, gtdb.write(format=5), flag.write(format=0))
    gtdb, flag = Tree(best[1], format=1), Tree(best[2], format=0)
    if set(bitmask_splits(gtdb, index)) != src_splits["GTDB"] or set(bitmask_splits(flag, index)) != src_splits["flag"]:
        raise SystemExit("tanglegram: drawn topology differs from the source trees")
    y_g, y_f = leaf_rows(gtdb), leaf_rows(flag)
    bitmask_splits(gtdb, index)  # annotates node.mask
    pos_g, top_g = pr.phylogram_layout(gtdb, y_g)
    pos_f, top_f = pr.phylogram_layout(flag, y_f)

    n = len(taxa)
    fig, ax = plt.subplots(figsize=(15, max(8, 0.13 * n + 4)))
    tree_w, lab_w, gap = 3.0, 2.6, 2.2
    xr_tip = tree_w + 2 * lab_w + gap

    def gtdb_style(nd):
        if nd.is_leaf():
            return pr.TEXT2, 0.7, "-"
        g = gcf_of_split.get(canon(nd.mask))
        if g is None:
            return pr.GRID, 0.9, "-"
        return (pr.ACCENT, 2.2, "-") if g <= 5 else (pr.SEQ[min(int(g / 100 * len(pr.SEQ)), len(pr.SEQ) - 1)], 1.6, "-")

    pr.draw_tree(ax, gtdb, pos_g, 0, tree_w, False, gtdb_style)
    pr.draw_tree(ax, flag, pos_f, xr_tip, tree_w, True,
                 lambda nd: (pr.TEXT2, 0.7, "-") if nd.is_leaf() else (pr.TEXT, 1.0, "-"))
    for nd in flag.traverse():
        if not nd.is_leaf() and not nd.is_root() and nd.support < 95:
            x, y = pos_f[nd]
            ax.plot(xr_tip + (1 - x) * tree_w, y, "o", ms=3.2, mfc="white", mec=pr.TEXT, mew=0.8, zorder=5)
    fs = 5.6 if n > 60 else 7.5
    for t in taxa:
        ax.text(tree_w + 0.06, y_g[t], org.get(t, t), va="center", ha="left", fontsize=fs, color=pr.TEXT)
        ax.text(xr_tip - 0.06, y_f[t], org.get(t, t), va="center", ha="right", fontsize=fs, color=pr.TEXT)
        ax.plot([tree_w + lab_w, xr_tip - lab_w], [y_g[t], y_f[t]],
                color=colour.get(family(t), pr.TEXT2), lw=0.8, alpha=0.75)
    ax.text(tree_w / 2, 2.2, "GTDB species tree", ha="center", fontsize=12, weight="bold", color=pr.TEXT)
    ax.text(tree_w / 2, 0.9, f"branch colour = gCF (share of the {n_genes} gene trees containing it)",
            ha="center", fontsize=7.5, color=pr.TEXT2)
    ax.text(xr_tip + tree_w / 2, 2.2, f"Concatenated tree ({n_genes} genes)", ha="center",
            fontsize=12, weight="bold", color=pr.TEXT)
    ax.text(xr_tip + tree_w / 2, 0.9, "○ = node with UFBoot < 95; " +
            ("rooted on the outgroup" if outgroup else "midpoint-rooted"),
            ha="center", fontsize=7.5, color=pr.TEXT2)
    pr.scale_bar(ax, 0.05, -n - 0.5, tree_w, top_g, False)
    pr.scale_bar(ax, xr_tip + tree_w - 0.05, -n - 0.5, tree_w, top_f, True)
    ly = -n - 4.5
    ax.text(0.05, ly, "gCF", fontsize=8, color=pr.TEXT, weight="bold", va="center")
    for i, c in enumerate(pr.SEQ):
        ax.plot([0.4 + i * 0.28, 0.65 + i * 0.28], [ly, ly], color=c, lw=5, solid_capstyle="butt")
    ax.plot([2.55, 2.85], [ly, ly], color=pr.ACCENT, lw=5, solid_capstyle="butt")
    ax.text(2.95, ly, "gCF ≤ 5%", fontsize=7.5, color=pr.TEXT, va="center")
    legend = groups + ([] if not outgroup else ["Outgroup"])
    for i, name in enumerate(legend + (["other"] if len(fams) > len(groups) else [])):
        x, y = tree_w + 1.4 + (i % 4) * 2.3, ly - (i // 4) * 2.0
        ax.plot([x, x + 0.3], [y, y], color=colour.get(name, pr.TEXT2), lw=3, solid_capstyle="butt")
        ax.text(x + 0.4, y, name, fontsize=7.5, color=pr.TEXT, va="center")
    ax.set_xlim(-0.1, xr_tip + tree_w + 0.1)
    ax.set_ylim(ly - 2.0 * (len(legend) // 4 + 2), 3)
    ax.axis("off")
    out = run.report / "tanglegram.png"
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"tanglegram: {best[0]} crossing line pairs after rotation -> {out}")
    return out


def tanglegram_html(png_bytes, note):
    img = base64.b64encode(png_bytes).decode()
    return (f'<h2>Tanglegram</h2><p class="note">{note}</p>'
            f'<img alt="Tanglegram" src="data:image/png;base64,{img}">')


TANGLE_NOTE = ("Left: GTDB tree, branches shaded by gCF (share of decisive flagellar gene trees that "
               "contain the branch). Right: concatenated flagellar tree. Lines join the same genome, "
               "coloured by GTDB family.")


def comparison_html(comp_path, sl_path=None):
    """Numbers from compare_to_gtdb.py (comparison.json) and, if given, the
    HGT shortlist (hgt_shortlist.json), as HTML."""
    e = html_lib.escape
    parts = []
    if comp_path and comp_path.exists():
        c = json.loads(comp_path.read_text())
        gt = c.get("concordance", {}).get("gtdb_tree", {}).get("gCF", {})
        rows = [("Genomes in both trees", c.get("n_taxa")),
                ("GTDB branches also in the flagellar tree",
                 f"{c['shared_branches']} of {c['branches']['gtdb']} ({c['shared_branches_frac_of_gtdb']:.0%})"),
                ("Shared branches expected by chance (shuffled labels, mean)", c.get("shuffle_null", {}).get("mean")),
                ("Normalised Robinson-Foulds distance", c.get("robinson_foulds", {}).get("normalised")),
                ("Mantel r, patristic distances", c.get("mantel", {}).get("r")),
                ("Mean gCF on GTDB branches", gt.get("mean")),
                ("GTDB branches with gCF < 10", gt.get("branches_below_10"))]
        parts.append("<h2>Flagellar tree vs GTDB</h2><table class=\"kv\">" + "".join(
            f"<tr><th>{e(k)}</th><td>{e(str(v))}</td></tr>" for k, v in rows if v is not None) + "</table>")
    if sl_path and sl_path.exists():
        sl = json.loads(sl_path.read_text())
        st, groups = sl["settings"], sl["groups"]
        passing = [g for g in groups if g["passes"]]
        near = [g for g in groups if not g["passes"] and g["n_gene_trees"] >= 3][:10]
        parts.append(
            "<h2>HGT candidates</h2><p class=\"note\">Groups of genomes that the flagellar genes place "
            "together against GTDB (UFBoot &ge; 95). A candidate passes if the group is in "
            f"&ge; {st['min_gene_trees']} gene trees and the concatenated tree; at the GTDB branch it "
            f"contradicts, &ge; {st['one_sided_min_frac']:.0%} of decisive gene trees support it and "
            f"&ge; {st['one_sided_ratio']:g}&times; the other alternative; and that branch is "
            f"&ge; {sl['min_branch_length']:.4f} substitutions/site (25th percentile). "
            "Candidates are not tested for significance (no AU tests).</p>")

        def table(rows_, caption):
            if not rows_:
                return f"<p class=\"note\">{caption}: none.</p>"
            head = ("<tr><th>Genomes</th><th>Gene trees</th><th>In concat</th>"
                    "<th>Support at GTDB branch</th><th>GTDB branch length</th></tr>")
            body = "".join(
                f"<tr><td>{e('; '.join(g['organisms']))}</td><td>{g['n_gene_trees']}</td>"
                f"<td>{'yes' if g['in_concat'] else 'no'}</td>"
                f"<td>{g.get('favoured', '')}/{g.get('decisive', '')} "
                f"({g.get('favoured_frac', 0):.0%}; other alt. {g.get('other_alt', '')})</td>"
                f"<td>{g.get('branch_length', '')}</td></tr>" for g in rows_)
            return f"<h3>{caption} ({len(rows_)})</h3><table class=\"grid\">{head}{body}</table>"
        parts.append(table(passing, "Candidates passing the rule"))
        parts.append(table(near, "Closest groups that did not pass"))
    return "".join(parts)


PAGE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__ flagellar report</title>
<style>
  :root { --bg:#fafafa; --panel:#fff; --text:#0b0b0b; --muted:#52514e; --line:#e1e0d9; --accent:#2a78d6; }
  @media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
    --bg:#161616; --panel:#1f1f1f; --text:#f2f2f0; --muted:#b5b3ad; --line:#33322f; --accent:#6eaaf0; } }
  :root[data-theme="dark"] { --bg:#161616; --panel:#1f1f1f; --text:#f2f2f0; --muted:#b5b3ad; --line:#33322f; --accent:#6eaaf0; }
  html, body { margin:0; height:100%; background:var(--bg); color:var(--text); font-family:system-ui,sans-serif; }
  body { display:flex; flex-direction:column; }
  header { display:flex; flex-wrap:wrap; align-items:center; gap:8px 16px; padding:10px 16px;
           background:var(--panel); border-bottom:1px solid var(--line); }
  h1 { font-size:16px; margin:0; }
  nav { display:flex; gap:4px; flex-wrap:wrap; }
  nav button { font:inherit; font-size:13px; padding:6px 12px; border:1px solid var(--line); border-radius:6px;
               background:transparent; color:var(--text); cursor:pointer; }
  nav button[aria-selected="true"] { background:var(--accent); border-color:var(--accent); color:#fff; }
  nav button:disabled { opacity:.45; cursor:not-allowed; }
  .tab { flex:1; min-height:0; display:none; }
  .tab.on { display:flex; }
  iframe { flex:1; border:0; width:100%; background:#fff; }
  #cmp { overflow:auto; padding:8px 16px 24px; flex-direction:column; }
  #cmp img { max-width:100%; height:auto; background:#fff; border:1px solid var(--line); border-radius:6px; }
  h2 { font-size:15px; margin:20px 0 8px; } h3 { font-size:13px; margin:16px 0 6px; }
  .note { color:var(--muted); font-size:12px; max-width:900px; }
  table { border-collapse:collapse; font-size:12px; margin-bottom:8px; }
  th, td { text-align:left; padding:4px 10px 4px 0; border-bottom:1px solid var(--line); vertical-align:top; }
  .kv th { font-weight:normal; color:var(--muted); }
  .grid { width:100%; max-width:1100px; }
  .wait { padding:24px 16px; color:var(--muted); }
</style>
</head>
<body>
<header>
  <h1>__TITLE__</h1>
  <nav role="tablist">
    <button data-tab="tree" __TREE_STATE__>Tree</button>
    <button data-tab="scatter" __SCATTER_STATE__>Synteny vs Phylogeny</button>
    <button data-tab="cmp" __CMP_STATE__>Flagellar vs GTDB</button>
  </nav>
</header>
<div id="tree" class="tab"><iframe title="Interactive tree"></iframe></div>
<div id="scatter" class="tab"><iframe title="Synteny vs phylogeny"></iframe></div>
<div id="cmp" class="tab">__CMP__</div>
<script type="application/json" id="page-tree">__TREE_PAGE__</script>
<script type="application/json" id="page-scatter">__SCATTER_PAGE__</script>
<script>
// Each embedded page is loaded into its iframe the first time its tab is shown,
// so it lays itself out at the visible size.
const loaded = {};
function show(id) {
  document.querySelectorAll(".tab").forEach(t => t.classList.toggle("on", t.id === id));
  document.querySelectorAll("nav button").forEach(b => b.setAttribute("aria-selected", b.dataset.tab === id));
  const src = document.getElementById("page-" + id);
  if (src && !loaded[id]) { document.querySelector("#" + id + " iframe").srcdoc = JSON.parse(src.textContent); loaded[id] = true; }
  try { localStorage.setItem("fp-tab", id); } catch (e) {}
}
document.querySelectorAll("nav button").forEach(b => b.addEventListener("click", () => show(b.dataset.tab)));
let start = null;
try { start = localStorage.getItem("fp-tab"); } catch (e) {}
const usable = [...document.querySelectorAll("nav button:not([disabled])")].map(b => b.dataset.tab);
show(usable.includes(start) ? start : usable[0] || "cmp");
</script>
</body>
</html>
"""


def _embed(page):
    """An HTML page as a JSON string that is safe inside <script>."""
    return json.dumps(page).replace("</", "<\\/") if page else "null"


def tree_defaults(page, reliable_only=True, adjacent_bp=0):
    """Start the interactive tree with 'Reliable gene order only' ticked and
    'Adjacent <=' at adjacent_bp (interactive_tree.py itself defaults to off / 20 bp)."""
    swaps = [('<input type="number" id="adjThreshold" min="0" step="1" value="20"',
              f'<input type="number" id="adjThreshold" min="0" step="1" value="{adjacent_bp}"'),
             ("let adjThreshold = 20;", f"let adjThreshold = {adjacent_bp};")]
    if reliable_only:
        swaps += [('<input type="checkbox" id="reliableOnly">', '<input type="checkbox" id="reliableOnly" checked>'),
                  ("let reliableOnly = false;", "let reliableOnly = true;")]
    for old, new in swaps:
        if old not in page:
            raise ValueError(f"interactive tree page changed; can't set default: {old}")
        page = page.replace(old, new, 1)
    return page


def build_page(title, tree_page, scatter_page_html, cmp_html):
    """The one-page report: Tree | Synteny vs phylogeny | Flagellar vs GTDB."""
    state = lambda x: "" if x else "disabled"
    return (PAGE.replace("__TITLE__", html_lib.escape(title))
            .replace("__TREE_STATE__", state(tree_page)).replace("__SCATTER_STATE__", state(scatter_page_html))
            .replace("__CMP_STATE__", "")
            .replace("__TREE_PAGE__", _embed(tree_page)).replace("__SCATTER_PAGE__", _embed(scatter_page_html))
            .replace("__CMP__", cmp_html))


def report(run_dir):
    run = Run(run_dir)
    cfg = run.config()
    title = cfg.get("name", run.root.name)
    run.report.mkdir(parents=True, exist_ok=True)

    tree_page = scatter_html = None
    if (run.root / "gene_order_combined.json").exists():
        tree_vals = run.synteny / "synteny_tree.json"
        n_rel = sum(v.get("reliable", False) for v in json.loads(tree_vals.read_text())["organisms"].values()) \
            if tree_vals.exists() else 0
        tree_page = tree_defaults(interactive_tree(run), reliable_only=n_rel >= 3)
    if (run.synteny / "synteny_all.json").exists():
        scatter_html = synteny_scatter(run, title)

    if (run.trees / "concordance" / "gtdb_gcf.cf.stat").exists() and (run.trees / "concat" / "concat.treefile").exists():
        png = tanglegram(run)
        cmp_html = tanglegram_html(png.read_bytes(), TANGLE_NOTE) + comparison_html(
            run.trees / "comparison.json", run.trees / "hgt_shortlist.json")
    else:
        print("Flagellar vs GTDB tab empty: the tree stage (katana/submit.sh RUN trees) hasn't finished yet.")
        cmp_html = '<p class="wait">Available after the tree stage (katana/submit.sh RUN trees).</p>'

    for old in ("interactive_tree.html", "synteny_scatter.html"):   # replaced by report.html
        (run.report / old).unlink(missing_ok=True)
    out = run.report / "report.html"
    out.write_text(build_page(title, tree_page, scatter_html, cmp_html))
    print(f"Wrote {out} ({out.stat().st_size / 1e6:.1f} MB)")
