#!/usr/bin/env python3
"""
run_ragtag_scaffold.py
For every organism flagged ragtag_possible in gene_order/ragtag_feasibility.json
(see check_ragtag_feasibility.py), downloads its reference genome and runs
`ragtag.py scaffold` against our draft assembly. Verified deterministic
(reran on the pilot genome, byte-identical AGP/fasta output both times) --
unlike the AMOS/minimus2 self-merge approach this project tried and
dropped, which was not.

Requires minimap2 (native, e.g. `brew install minimap2`) and the `ragtag`
pip package -- both already in this project's .venv.

Skips any organism whose output already exists, so it's safe to re-run
after new species become ragtag_possible.

Run AFTER check_ragtag_feasibility.py.
"""
import gzip
import json
import shutil
import subprocess
import urllib.request
from pathlib import Path

FEASIBILITY_JSON = "gene_order/ragtag_feasibility.json"
CAMPY_FAS = Path("campy_fetched/campy_fas")
OUTPUT_DIR = Path("campy_fetched/ragtag_scaffolding/output")
RAGTAG_BIN = ".venv/bin/ragtag.py"


def download_reference(ftp_path: str, dest: Path) -> None:
    fname = ftp_path.rsplit("/", 1)[-1] + "_genomic.fna.gz"
    url = ftp_path.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov") + "/" + fname
    gz_path = dest.with_suffix(dest.suffix + ".gz")
    urllib.request.urlretrieve(url, gz_path)
    with gzip.open(gz_path, "rb") as f_in, open(dest, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    gz_path.unlink()


def main():
    feasibility = json.load(open(FEASIBILITY_JSON))
    targets = {org: v for org, v in feasibility.items() if v["ragtag_possible"]}
    print(f"{len(targets)} organisms flagged ragtag_possible")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for i, (org, info) in enumerate(sorted(targets.items()), 1):
        acc = info["accession"]
        gdir = OUTPUT_DIR / acc
        result_fasta = gdir / "ragtag_out" / "ragtag.scaffold.fasta"
        if result_fasta.exists():
            print(f"[{i}/{len(targets)}] {org}: already scaffolded, skipping")
            continue

        print(f"[{i}/{len(targets)}] {org} ({acc}) vs {info['reference_accession']} ({info['reference_status']})")
        gdir.mkdir(parents=True, exist_ok=True)

        ref_fna = gdir / "reference.fna"
        if not ref_fna.exists():
            download_reference(info["reference_ftp_path"], ref_fna)

        query_fna = gdir / "query.fna"
        shutil.copyfile(CAMPY_FAS / f"{acc}.fna", query_fna)

        out_dir = gdir / "ragtag_out"
        proc = subprocess.run(
            [RAGTAG_BIN, "scaffold", str(ref_fna), str(query_fna), "-o", str(out_dir), "-t", "4"],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            print(f"  RAGTAG FAILED: {proc.stderr[-500:]}")
            continue
        print(f"  done -> {out_dir}")

    print("\nAll flagged organisms processed. Next: lift_ragtag_coordinates.py")


if __name__ == "__main__":
    main()
