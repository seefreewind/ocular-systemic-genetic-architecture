#!/usr/bin/env python3
"""Generate Figure 3A-3C using the Phase 3A frozen outputs."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures/phase3A"
OUT.mkdir(parents=True, exist_ok=True)

PALETTE = {
    "blue_main": "#0F4D92", "blue_secondary": "#3775BA", "green": "#8BCF8B",
    "red": "#B64342", "neutral": "#CFCECE", "teal": "#42949E", "violet": "#9A4D8E",
}


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def style():
    plt.rcParams.update({
        "font.family": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 9, "axes.spines.right": False, "axes.spines.top": False,
        "axes.linewidth": 1.2, "legend.frameon": False, "pdf.fonttype": 42,
    })


def figure3a():
    data = rows(ROOT / "results/phase3A/LOCUS_BIOLOGICAL_MATRIX.tsv")
    data.sort(key=lambda x: (0 if "NEGATIVE" in x["direction_class"] else 1, x["pair"], int(x["chr"]), x["lead_SNP"]))
    columns = ["negative", "positive", "Level1", "Level2", "Level3", "GTEx lookup", "multi-pair region"]
    matrix = np.zeros((len(data), len(columns)))
    recurrence = rows(ROOT / "results/phase3A/RECURRENT_CROSS_PAIR_LOCI.tsv")
    region_pairs = {r["genomic_region"]: int(r["n_pairs"]) for r in recurrence}
    for i, row in enumerate(data):
        matrix[i, 0] = "NEGATIVE" in row["direction_class"]
        matrix[i, 1] = "POSITIVE" in row["direction_class"]
        matrix[i, 2] = row["positional_priority"] == "LEVEL1"
        matrix[i, 3] = row["positional_priority"] == "LEVEL2"
        matrix[i, 4] = row["positional_priority"] == "LEVEL3"
        matrix[i, 5] = row["qtl_lookup_status"] == "LOOKUP_RETURNED"
        matrix[i, 6] = region_pairs.get(row["SUPERGNOVA_region"], 1) > 1
    fig, ax = plt.subplots(figsize=(7.0, 8.5))
    cmap = ListedColormap(["#F7F7F7", PALETTE["blue_main"]])
    ax.imshow(matrix, aspect="auto", interpolation="none", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(columns)), columns, rotation=35, ha="right")
    ax.set_yticks(range(len(data)), [x["lead_SNP"] for x in data], fontsize=5.8)
    ax.set_xlabel("Descriptive annotation feature")
    ax.set_ylabel("Frozen lead variant (58 loci)")
    ax.set_title("Figure 3A | Locus-level annotation matrix", loc="left", fontweight="bold")
    ax.axhline(sum("NEGATIVE" in x["direction_class"] for x in data) - 0.5, color=PALETTE["red"], lw=1.2)
    ax.text(-0.8, 15, "negative-rho", color=PALETTE["red"], rotation=90, va="center", ha="right", fontsize=8)
    ax.text(-0.8, 45, "positive-rho", color=PALETTE["blue_main"], rotation=90, va="center", ha="right", fontsize=8)
    ax.text(0, -3.5, "filled = present; blank = not observed in lookup", fontsize=7)
    fig.tight_layout(pad=1.2)
    fig.savefig(OUT / "Figure3A_locus_annotation_matrix.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def figure3b():
    sensitivity = rows(ROOT / "results/phase3A/MAJOR_LOCUS_ENRICHMENT_SENSITIVITY.tsv")
    primary = [r for r in sensitivity if r["analysis_set"] in {"NEGATIVE_ONE_LOCUS_ONE_GENE", "POSITIVE_ONE_LOCUS_ONE_GENE"}]
    labels = ["NEGATIVE\n(rho < 0)", "POSITIVE\n(rho > 0)"]
    values = []
    for name in ("NEGATIVE_ONE_LOCUS_ONE_GENE", "POSITIVE_ONE_LOCUS_ONE_GENE"):
        row = next((x for x in primary if x["analysis_set"] == name), None)
        values.append(int(row["n_terms_fdr_lt_0_05"]) if row else 0)
    fig, (ax, tx) = plt.subplots(1, 2, figsize=(8.5, 3.6), gridspec_kw={"width_ratios": [1.0, 1.6]})
    bars = ax.bar(labels, values, color=[PALETTE["red"], PALETTE["blue_main"]], edgecolor="black", linewidth=0.8)
    ax.set_ylabel("GO-BP terms with custom-background FDR < 0.05")
    ax.set_ylim(0, max(1, max(values) + 1))
    ax.set_title("Figure 3B | Direction-specific pathways", loc="left", fontweight="bold")
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.04, str(value), ha="center", va="bottom", fontsize=11)
    ax.text(0.5, 0.55, "No direction-specific\nGO-BP enrichment passed\nFDR < 0.05", transform=ax.transAxes,
            ha="center", va="center", fontsize=10, color="#4D4D4D")
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.6)
    tx.axis("off")
    tx.text(0.02, 0.90, "Interpretation", fontsize=11, fontweight="bold")
    tx.text(0.02, 0.78, "Primary gene sets were defined by one-locus-one-gene\npriority (Level 1 → Level 2 → Level 3).\n\n"
            "The custom reference was restricted to genes\nobserved among 2,096 valid independent loci.\n\n"
            "Reactome lookups are reported separately as\ndescriptive default-background support and were\nnot used to claim functional convergence.",
            va="top", fontsize=9, linespacing=1.55)
    fig.tight_layout(pad=1.2)
    fig.savefig(OUT / "Figure3B_direction_specific_pathways.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def figure3c():
    data = rows(ROOT / "results/phase3A/RECURRENT_CROSS_PAIR_LOCI.tsv")
    data = sorted(data, key=lambda x: (-int(x["n_pairs"]), -int(x["n_loci"]), x["genomic_region"]))[:12]
    data = list(reversed(data))
    labels = [x["genomic_region"].replace(":", "\n", 1) for x in data]
    neg = np.array([int(x["n_negative_loci"]) for x in data])
    pos = np.array([int(x["n_positive_loci"]) for x in data])
    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    y = np.arange(len(data))
    ax.barh(y, neg, color=PALETTE["red"], label="negative-rho loci", edgecolor="black", linewidth=0.5)
    ax.barh(y, pos, left=neg, color=PALETTE["blue_main"], label="positive-rho loci", edgecolor="black", linewidth=0.5)
    for i, row in enumerate(data):
        ax.text(int(row["n_loci"]) + 0.1, i, f"{row['n_pairs']} pairs", va="center", fontsize=8)
    ax.set_yticks(y, labels, fontsize=7)
    ax.set_xlabel("Independent frozen loci in region")
    ax.set_title("Figure 3C | Recurrent cross-pair regions", loc="left", fontweight="bold")
    ax.legend(loc="lower right")
    ax.set_xlim(0, max(int(x["n_loci"]) for x in data) + 2.3)
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.6)
    ax.set_axisbelow(True)
    fig.tight_layout(pad=1.2)
    fig.savefig(OUT / "Figure3C_recurrent_cross_pair_loci.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    style()
    figure3a(); figure3b(); figure3c()
    print("Phase 3A figures written to", OUT)


if __name__ == "__main__":
    main()
