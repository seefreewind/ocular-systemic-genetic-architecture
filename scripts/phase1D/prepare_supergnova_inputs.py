#!/usr/bin/env python3
"""Prepare all frozen Phase 1D traits for the frozen SUPERGNOVA runtime."""

from __future__ import annotations

import csv
import gzip
import math
import re
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HARM = ROOT / "data" / "harmonized"
REF_BFILES = ROOT / "data" / "phase1C" / "SUPERGNOVA_reference" / "bfiles"
OUT = ROOT / "data" / "phase1D" / "sumstats"
QC = ROOT / "results" / "phase1D" / "SUPERGNOVA_INPUT_QC.tsv"
TRAITS = ["AMD", "POAG", "CAT", "CAD", "STROKE", "T2D", "CKD", "AD", "PD"]
NUMBER = re.compile(r"^[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def is_number(value: str) -> bool:
    return bool(NUMBER.match(value.strip()))


def load_freeze() -> dict[str, dict[str, str]]:
    with (ROOT / "config" / "phase1D_dataset_freeze.tsv").open(newline="") as handle:
        return {row["trait"]: row for row in csv.DictReader(handle, delimiter="\t")}


def load_reference() -> dict[str, tuple[str, str]]:
    reference: dict[str, tuple[str, str]] = {}
    for chrom in range(1, 23):
        with (REF_BFILES / f"eur_chr{chrom}_SNPmaf5.bim").open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split()
                if len(fields) >= 6:
                    reference.setdefault(fields[1], (fields[4].upper(), fields[5].upper()))
    return reference


def prepare(trait: str, metadata: dict[str, str], reference: dict[str, tuple[str, str]]) -> dict[str, str | int | float]:
    source = HARM / f"{trait}.standardized.tsv.gz"
    output = OUT / f"{trait}.supergnova.sumstats.gz"
    metrics = {key: 0 for key in [
        "raw_rows", "valid_rows", "duplicate_rows_removed", "missing_required_rows",
        "invalid_numeric_rows", "invalid_allele_rows", "p_outside_0_1",
        "z_positive", "z_negative", "z_zero", "reference_overlap", "allele_compatible",
        "allele_incompatible",
    ]}
    metrics["n_min"] = math.inf
    metrics["n_max"] = -math.inf
    seen: set[str] = set()
    with gzip.open(source, "rt", newline="") as src, gzip.open(output, "wt", compresslevel=1, newline="") as dst:
        header = src.readline().rstrip("\n").split("\t")
        index = {name: header.index(name) for name in ["SNP", "A1", "A2", "BETA", "SE", "P", "N"]}
        dst.write("SNP\tA1\tA2\tZ\tN\n")
        for line in src:
            metrics["raw_rows"] += 1
            fields = line.rstrip("\n").split("\t")
            vals = {name: fields[index[name]].strip() if index[name] < len(fields) else "" for name in index}
            if any(value == "" or value.upper() in {"NA", "NAN", "NULL"} for value in vals.values()):
                metrics["missing_required_rows"] += 1
                continue
            snp, a1, a2 = vals["SNP"], vals["A1"].upper(), vals["A2"].upper()
            if not all(is_number(vals[name]) for name in ["BETA", "SE", "P", "N"]) or float(vals["SE"]) <= 0:
                metrics["invalid_numeric_rows"] += 1
                continue
            if a1 not in COMP or a2 not in COMP or a1 == a2:
                metrics["invalid_allele_rows"] += 1
                continue
            if snp in seen:
                metrics["duplicate_rows_removed"] += 1
                continue
            seen.add(snp)
            z = float(vals["BETA"]) / float(vals["SE"])
            if not math.isfinite(z):
                metrics["invalid_numeric_rows"] += 1
                continue
            p = float(vals["P"])
            n = float(vals["N"])
            metrics["p_outside_0_1"] += int(p < 0 or p > 1)
            metrics["z_positive"] += int(z > 0)
            metrics["z_negative"] += int(z < 0)
            metrics["z_zero"] += int(z == 0)
            metrics["n_min"] = min(metrics["n_min"], n)
            metrics["n_max"] = max(metrics["n_max"], n)
            dst.write(f"{snp}\t{a1}\t{a2}\t{z:.12g}\t{vals['N']}\n")
            metrics["valid_rows"] += 1
            if snp in reference:
                metrics["reference_overlap"] += 1
                r1, r2 = reference[snp]
                if (a1, a2) in {(r1, r2), (r2, r1), (COMP[r1], COMP[r2]), (COMP[r2], COMP[r1])}:
                    metrics["allele_compatible"] += 1
                else:
                    metrics["allele_incompatible"] += 1
    metrics["reference_snps"] = len(reference)
    metrics["reference_overlap_fraction"] = metrics["reference_overlap"] / metrics["valid_rows"] if metrics["valid_rows"] else 0
    if metrics["n_min"] is math.inf:
        metrics["n_min"] = math.nan
        metrics["n_max"] = math.nan
    return {
        "trait": trait,
        "input_path": str(source.relative_to(ROOT)),
        "output_path": str(output.relative_to(ROOT)),
        "build": metadata["build"],
        "ancestry": metadata["ancestry"],
        "raw_N_manifest": metadata["N"],
        "effective_N_manifest": metadata["effective_N"],
        "scalar_N_used": str(round(float(metadata["effective_N"]))),
        **metrics,
        "created": now(),
    }


def main() -> None:
    freeze = load_freeze()
    reference = load_reference()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [prepare(trait, freeze[trait], reference) for trait in TRAITS]
    with QC.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
