#!/usr/bin/env python3
"""
prepare_hgt_spikein.py
Positive control for the HGT screen: plants known (artificial) horizontal
transfers in copies of real flagellar gene files, so the pipeline's power to
detect a transfer can be measured. Full methodology: gene_order/hgt_spikein/README.md.

For each (gene, donor) pair two test files are written, both built from
gene_order/flagellar_trees_ext/genes/<gene>.faa (the source is only read,
never modified):
  spike    the recipient's sequence is replaced by the donor's sequence
           (header kept, so the recipient "carries the donor's gene"), and the
           donor itself is removed -- so the transferred copy is not identical
           to any sampled sequence and must group with the donor's relatives,
           as a gene acquired from that lineage some time ago would;
  control  the donor is removed, the recipient is untouched -- isolates the
           effect of the transfer from the effect of dropping the donor.

3 genes (short / medium / long alignment) x 3 donors (within genus / close
genus / distant genus) x {spike, control} = 18 test genes.

Usage (from the project root):
  python scripts/11_flagellar_trees/prepare_hgt_spikein.py
Then on Katana (gene trees use exactly the per-gene pipeline of the real run):
  qsub -v FLAG_DIR=gene_order/hgt_spikein -J 1-18 scripts/11_flagellar_trees/katana/01_gene_trees.pbs
Then locally, after copying gene_trees/ back:
  FLAG_DIR=gene_order/hgt_spikein python scripts/11_flagellar_trees/screen_hgt_candidates.py
  python scripts/11_flagellar_trees/summarise_hgt_spikein.py
"""
import hashlib
import shutil
from pathlib import Path

SRC = Path("gene_order/flagellar_trees_ext")
DEST = Path("gene_order/hgt_spikein")

RECIPIENT = ("GCF_001457695.1", "Campylobacter jejuni")
# trimmed alignment length in the 42-gene run in brackets
GENES = [("flie", "short (67 sites)"), ("flhb", "medium (297 sites)"), ("flha", "long (556 sites)")]
DONORS = [
    ("within_genus", "GCF_013372185.1", "Campylobacter lari"),
    ("close_genus", "GCF_900478295.1", "Helicobacter pylori"),
    ("distant_genus", "GCF_000012965.1", "Sulfurimonas denitrificans"),
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


def main():
    (DEST / "genes").mkdir(parents=True, exist_ok=True)
    for name in ["taxa.tsv", "gtdb_ref_pruned.nwk"]:
        shutil.copy2(SRC / name, DEST / name)

    rec_acc, rec_org = RECIPIENT
    before = {g: md5(SRC / "genes" / f"{g}.faa") for g, _ in GENES}
    design, ids = [], []
    for gene, length_class in GENES:
        seqs = read_fasta(SRC / "genes" / f"{gene}.faa")
        for distance, don_acc, don_org in DONORS:
            for acc in (rec_acc, don_acc):
                if acc not in seqs:
                    raise SystemExit(f"{gene}: {acc} missing from {SRC}/genes/{gene}.faa")
            for condition in ("spike", "control"):
                test = {k: v for k, v in seqs.items() if k != don_acc}
                if condition == "spike":
                    test[rec_acc] = seqs[don_acc]
                tid = f"{gene}_{distance}_{condition}"
                write_fasta(DEST / "genes" / f"{tid}.faa", test)
                ids.append(tid)
                design.append([tid, gene, length_class, distance, condition, rec_acc, rec_org,
                               don_acc, don_org, str(len(test))])
    after = {g: md5(SRC / "genes" / f"{g}.faa") for g, _ in GENES}
    assert before == after, "source gene files changed"

    (DEST / "genes.txt").write_text("\n".join(ids) + "\n")
    with open(DEST / "design.tsv", "w") as f:
        f.write("test_id\tgene\tlength_class\ttransfer_distance\tcondition\trecipient_acc\t"
                "recipient\tdonor_acc\tdonor\tn_seqs\n")
        for row in design:
            f.write("\t".join(row) + "\n")
    print(f"Wrote {len(ids)} test genes to {DEST}/genes/ (source files unchanged, md5 checked)")


if __name__ == "__main__":
    main()
