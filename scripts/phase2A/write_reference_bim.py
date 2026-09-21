#!/usr/bin/env python3
"""Concatenate the frozen EUR panel BIM coordinates for rsID reconciliation."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "data" / "phase1C" / "SUPERGNOVA_reference" / "bfiles"
OUT = ROOT / "results" / "phase2A" / "reference_freq" / "eur_reference.bim.tsv"


def main() -> None:
    with OUT.open("w") as handle:
        handle.write("CHR\tSNP\tCM\tBP\tREF\tALT\n")
        for chrom in range(1, 23):
            path = PANEL / f"eur_chr{chrom}_SNPmaf5.bim"
            with path.open() as source:
                for line in source:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) >= 6:
                        handle.write("\t".join(fields[:6]) + "\n")


if __name__ == "__main__":
    main()
