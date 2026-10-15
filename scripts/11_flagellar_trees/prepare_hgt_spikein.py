#!/usr/bin/env python3
"""
prepare_hgt_spikein.py
Positive control for the HGT screen: plants known (artificial) horizontal
transfers in copies of real flagellar gene files, so the pipeline's power to
detect a transfer can be measured. Full methodology: gene_order/hgt_spikein/README.md.

For each (gene, recipient, donor) three-way combination two test files are
written, both built from gene_order/flagellar_trees_ext/genes/<gene>.faa (the
source is only read, never modified; md5 checked):
  spike    the recipient's sequence is replaced by the donor's sequence
           (header kept, so the recipient "carries the donor's gene"), and the
           donor itself is removed -- so the transferred copy is not identical
           to any sampled sequence and must group with the donor's relatives,
           as a gene acquired from that lineage some time ago would;
  control  the donor is removed, the recipient is untouched -- isolates the
           effect of the transfer from the effect of dropping the donor.

3 genes (short / medium / long alignment) x 3 recipients (Campylobacter,
Helicobacter, Arcobacteraceae) x 4 donor levels x {spike, control} = 72 tests.
Donor levels are defined by position in the GTDB tree (GTDB distance is
reported alongside, since branch lengths differ a lot between families):
  very_close     nearest GTDB clade that still leaves >= 2 of the donor's
                 relatives once the donor is removed
  within_family  the other half of the recipient's family (deepest split)
  sister_family  the GTDB clade sister to the recipient's family
  distant        a family further away in the GTDB tree
Every donor has >= 2 relatives left in all three gene files, so "grouped with
the donor's relatives" is always testable.

Usage (from the project root):
  python scripts/11_flagellar_trees/prepare_hgt_spikein.py
Then on Katana (gene trees use exactly the per-gene pipeline of the real run):
  qsub -v FLAG_DIR=gene_order/hgt_spikein -J 1-72 scripts/11_flagellar_trees/katana/01_gene_trees.pbs
Then locally, after copying gene_trees/ back:
  FLAG_DIR=gene_order/hgt_spikein python scripts/11_flagellar_trees/screen_hgt_candidates.py
  python scripts/11_flagellar_trees/summarise_hgt_spikein.py
  python scripts/11_flagellar_trees/plot_hgt_spikein.py
"""
import hashlib
import shutil
from pathlib import Path

from ete3 import Tree

SRC = Path("gene_order/flagellar_trees_ext")
DEST = Path("gene_order/hgt_spikein")

# trimmed alignment length in the 42-gene run in brackets
GENES = [("flie", "short (67 sites)"), ("flhb", "medium (297 sites)"), ("flha", "long (556 sites)")]
LEVELS = ["very_close", "within_family", "sister_family", "distant"]
# recipient code, accession, organism, donors in LEVELS order (accession, organism)
RECIPIENTS = [
    ("cjej", "GCF_001457695.1", "Campylobacter jejuni", [
        ("GCF_000254135.1", "Campylobacter coli"),
        (None, "Campylobacter fetus"),
        (None, "Sulfurospirillum deleyianum"),
        ("GCF_000012965.1", "Sulfurimonas denitrificans")]),
    ("hfel", None, "Helicobacter felis", [
        (None, "Helicobacter salomonis"),
        (None, "Helicobacter pullorum"),
        (None, "Sulfurospirillum deleyianum"),
        ("GCF_000012965.1", "Sulfurimonas denitrificans")]),
    ("abut", None, "Aliarcobacter butzleri", [
        (None, "Aliarcobacter vitoriensis"),
        (None, "Malaciobacter marinus"),
        (None, "Nitrosophilus alvini"),
        (None, "Campylobacter fetus")]),
]


def read_fasta(path):
    seqs, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(path, seqs):
    path.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))


def md5(path):
    return hashlib.md5(path.read_bytes()).hexdigest()


def donor_relatives(gtdb, rec, don, present):
    """Smallest GTDB clade around the donor with >= 2 other genomes present,
    not containing the recipient (None if there is none)."""
    n = gtdb & don
    while n.up:
        n = n.up
        clade = set(n.get_leaf_names())
        if rec in clade:
            return None
        if len((clade - {don}) & present) >= 2:
            return (clade - {don}) & present
    return None


def main():
    taxa = {l.split("\t")[0]: l.split("\t")[1] for l in (SRC / "taxa.tsv").read_text().splitlines()[1:]}
    acc_of = {v: k for k, v in taxa.items()}
    gtdb = Tree(str(SRC / "gtdb_ref_pruned.nwk"), format=1)

    out = DEST / "genes"
    if out.exists():
        shutil.rmtree(out)  # regenerated in full; only ever holds this script's output
    out.mkdir(parents=True)
    for name in ["taxa.tsv", "gtdb_ref_pruned.nwk"]:
        shutil.copy2(SRC / name, DEST / name)

    before = {g: md5(SRC / "genes" / f"{g}.faa") for g, _ in GENES}
    design, ids = [], []
    for gene, length_class in GENES:
        seqs = read_fasta(SRC / "genes" / f"{gene}.faa")
        for code, rec_acc, rec_org, donors in RECIPIENTS:
            rec_acc = rec_acc or acc_of[rec_org]
            for level, (don_acc, don_org) in zip(LEVELS, donors):
                don_acc = don_acc or acc_of[don_org]
                assert taxa[don_acc] == don_org and taxa[rec_acc] == rec_org
                for acc in (rec_acc, don_acc):
                    if acc not in seqs:
                        raise SystemExit(f"{gene}: {taxa[acc]} missing from {SRC}/genes/{gene}.faa")
                rel = donor_relatives(gtdb, rec_acc, don_acc, set(seqs))
                if not rel:
                    raise SystemExit(f"{gene}: {don_org} has < 2 relatives apart from {rec_org}")
                dist = gtdb.get_distance(rec_acc, don_acc)
                for condition in ("spike", "control"):
                    test = {k: v for k, v in seqs.items() if k != don_acc}
                    if condition == "spike":
                        test[rec_acc] = seqs[don_acc]
                    tid = f"{gene}_{code}_{level}_{condition}"
                    write_fasta(out / f"{tid}.faa", test)
                    ids.append(tid)
                    design.append([tid, gene, length_class, code, level, condition, rec_acc, rec_org,
                                   don_acc, don_org, f"{dist:.3f}", str(len(test)),
                                   ",".join(sorted(rel))])
    after = {g: md5(SRC / "genes" / f"{g}.faa") for g, _ in GENES}
    assert before == after, "source gene files changed"

    (DEST / "genes.txt").write_text("\n".join(ids) + "\n")
    with open(DEST / "design.tsv", "w") as f:
        f.write("test_id\tgene\tlength_class\trecipient_code\tdonor_level\tcondition\trecipient_acc\t"
                "recipient\tdonor_acc\tdonor\tgtdb_distance\tn_seqs\tdonor_relatives\n")
        for row in design:
            f.write("\t".join(row) + "\n")
    print(f"Wrote {len(ids)} test genes to {out}/ (source files unchanged, md5 checked)")


if __name__ == "__main__":
    main()
