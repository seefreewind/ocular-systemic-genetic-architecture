#!/usr/bin/env python3
"""Create the documented EUR LD-pruned SNP set for PLACO+ sensitivity."""

from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "data" / "phase1C" / "SUPERGNOVA_reference" / "bfiles"
OUT = ROOT / "results" / "phase2A" / "reference_freq"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    paths = []
    for chrom in range(1, 23):
        prefix = PANEL / f"eur_chr{chrom}_SNPmaf5"
        out_prefix = OUT / f"eur_chr{chrom}.ldprune2"
        prune = Path(str(out_prefix) + ".prune.in")
        if not prune.exists():
            subprocess.run(
                [
                    "plink2", "--bfile", str(prefix),
                    "--set-all-var-ids", "@:#:$r:$a",
                    "--rm-dup", "force-first",
                    "--indep-pairwise", "1000kb", "1", "0.1",
                    "--threads", "4", "--memory", "4096",
                    "--out", str(out_prefix),
                ], check=True, cwd=str(ROOT), stdout=subprocess.DEVNULL,
            )
        paths.append(prune)
    combined = OUT / "eur_ld_pruned_r2_0.1_1000kb.in"
    ids = []
    for path in paths:
        ids.extend(line.strip() for line in path.open() if line.strip())
    with combined.open("w") as handle:
        handle.write("\n".join(dict.fromkeys(ids)) + "\n")
    print(f"LD-pruned SNPs: {len(set(ids))}")


if __name__ == "__main__":
    main()
