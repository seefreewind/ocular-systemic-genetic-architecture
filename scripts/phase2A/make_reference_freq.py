#!/usr/bin/env python3
"""Create allele-frequency tables from the frozen EUR LD panel."""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "data" / "phase1C" / "SUPERGNOVA_reference" / "bfiles"
OUT = ROOT / "results" / "phase2A" / "reference_freq"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for chrom in range(1, 23):
        prefix = PANEL / f"eur_chr{chrom}_SNPmaf5"
        out_prefix = OUT / f"eur_chr{chrom}"
        if (out_prefix.with_suffix(".afreq")).exists():
            continue
        subprocess.run(
            [
                "plink2",
                "--bfile", str(prefix),
                "--freq",
                "--threads", "4",
                "--memory", "4096",
                "--out", str(out_prefix),
            ],
            check=True,
            cwd=str(ROOT),
        )

    combined = OUT / "eur_reference.afreq.tsv"
    with combined.open("w", newline="") as handle:
        writer = None
        for chrom in range(1, 23):
            path = OUT / f"eur_chr{chrom}.afreq"
            with path.open() as source:
                reader = csv.DictReader(source, delimiter="\t")
                if writer is None:
                    fields = ["CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT"]
                    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                for row in reader:
                    writer.writerow({
                        "CHROM": row.get("CHROM", row.get("#CHROM", "")),
                        "ID": row.get("ID", ""),
                        "REF": row.get("REF", ""),
                        "ALT": row.get("ALT", ""),
                        "ALT_FREQS": row.get("ALT_FREQS", ""),
                        "OBS_CT": row.get("OBS_CT", ""),
                    })


if __name__ == "__main__":
    main()
