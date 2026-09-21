#!/usr/bin/env python3
"""Create the three prespecified Phase 1D atlas figures."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "phase1D"
FIGURES = ROOT / "figures" / "phase1D"
FIGURES.mkdir(parents=True, exist_ok=True)

BLUE = "#0F4D92"
BLUE2 = "#3775BA"
GREEN = "#8BCF8B"
RED = "#B64342"
NEUTRAL = "#CFCECE"
GOLD = "#FFD700"
TEAL = "#42949E"
VIOLET = "#9A4D8E"
GRAY = "#5A5A5A"


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#444444")
    ax.spines["bottom"].set_color("#444444")
    ax.tick_params(colors="#333333", labelsize=8)


def save(fig, stem: str):
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def cumulative_coordinates(frame: pd.DataFrame):
    chrom_max = frame.groupby("chr")["end"].max().to_dict()
    offsets = {}
    running = 0
    for chrom in range(1, 23):
        offsets[chrom] = running
        running += float(chrom_max.get(chrom, 0)) + 5_000_000
    x = frame["start"].astype(float).to_numpy() + frame["chr"].map(offsets).astype(float).to_numpy()
    centers = {chrom: offsets[chrom] + float(chrom_max.get(chrom, 0)) / 2 for chrom in offsets if chrom in chrom_max}
    return x, centers


def figure_2a(atlas: pd.DataFrame):
    pairs = list(atlas["pair"].drop_duplicates())
    x, centers = cumulative_coordinates(atlas)
    p = atlas["P_rho"].astype(float).clip(lower=1e-300, upper=1.0)
    y_signal = np.sign(atlas["rho"].astype(float)) * (-np.log10(p))
    y_positions = {pair: len(pairs) - 1 - i for i, pair in enumerate(pairs)}
    fig, ax = plt.subplots(figsize=(15, 10.5))
    for pair in pairs:
        mask = atlas["pair"].eq(pair).to_numpy()
        ypos = y_positions[pair]
        pair_x = x[mask]
        pair_y = ypos + y_signal.to_numpy()[mask] / max(1.0, float(np.nanmax(np.abs(y_signal.to_numpy()))) * 1.25)
        pair_rho = atlas.loc[mask, "rho"].to_numpy()
        pair_q = atlas.loc[mask, "q_atlas"].to_numpy()
        positive = pair_rho > 0
        negative = pair_rho < 0
        for direction, direction_mask, color in [("positive", positive, BLUE), ("negative", negative, RED)]:
            if not direction_mask.any():
                continue
            sig = direction_mask & (pair_q < 0.05)
            ax.scatter(pair_x[direction_mask], pair_y[direction_mask], s=5, color=color, alpha=0.28, linewidths=0)
            if sig.any():
                ax.scatter(pair_x[sig], pair_y[sig], s=22, color=color, edgecolor=GOLD, linewidth=0.7, alpha=0.95, zorder=3)
        ax.axhline(ypos, color="#E7E7E7", linewidth=0.6, zorder=0)
    ax.axhline((len(pairs) - 1) / 2, color="#BBBBBB", linewidth=0.8)
    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels(pairs, fontsize=8)
    ax.set_xlabel("Genomic position across chromosomes (GRCh37)")
    ax.set_ylabel("Pair-level row; vertical displacement is signed −log10(Pρ)")
    ax.set_title("Figure 2A. Signed local covariance atlas", loc="left", fontsize=12, fontweight="bold")
    ax.set_xticks(list(centers.values()))
    ax.set_xticklabels([str(chrom) for chrom in centers], fontsize=8)
    ax.set_ylim(-1, len(pairs))
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor=BLUE, markersize=5, label="Positive ρ"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=RED, markersize=5, label="Negative ρ"),
        Line2D([0], [0], marker="o", color=GOLD, markerfacecolor=TEAL, markersize=6, label="Atlas FDR q < 0.05"),
    ], frameon=False, loc="upper right", fontsize=8)
    style(ax)
    save(fig, "Figure2A_signed_local_covariance_atlas")


def figure_2b(classification: pd.DataFrame, burden: pd.DataFrame):
    frame = classification.merge(burden, on=["pair", "ocular_trait", "systemic_trait"], how="left")
    frame["balance_index"] = frame["signed_net"] / frame["absolute_total"].replace(0, np.nan)
    color_map = {
        "MIXED_DIRECTION": TEAL,
        "POSITIVE_ONLY": BLUE,
        "NEGATIVE_ONLY": RED,
        "NO_DETECTABLE_LOCAL_SIGNAL": NEUTRAL,
    }
    colors = [color_map.get(value, GRAY) for value in frame["local_class"]]
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    sizes = 35 + 8 * frame["local_total_atlas_fdr_count"].fillna(0).to_numpy()
    ax.scatter(frame["global_rg"], frame["balance_index"], s=sizes, c=colors, alpha=0.9, edgecolor="white", linewidth=0.8)
    for _, row in frame.iterrows():
        if row["hidden_mixed_local_sharing"] == "HIDDEN_MIXED_LOCAL_SHARING":
            ax.annotate(row["pair"], (row["global_rg"], row["balance_index"]), xytext=(5, 5), textcoords="offset points", fontsize=8, color=VIOLET)
    ax.axvline(0, color="#AAAAAA", linewidth=0.8)
    ax.axhline(0, color="#AAAAAA", linewidth=0.8)
    ax.set_xlabel("Global LDSC genetic correlation (rɡ)")
    ax.set_ylabel("Local signed balance = Σρ / Σ|ρ|")
    ax.set_title("Figure 2B. Global versus local balance", loc="left", fontsize=12, fontweight="bold")
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor=TEAL, markersize=7, label="Mixed local directions"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=BLUE, markersize=7, label="Positive-only local"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=RED, markersize=7, label="Negative-only local"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=NEUTRAL, markersize=7, label="No atlas local signal"),
    ], frameon=False, fontsize=8, loc="best")
    style(ax)
    save(fig, "Figure2B_global_vs_local_balance")


def figure_2c(direction: pd.DataFrame):
    frame = direction.copy()
    frame["label"] = frame["pair"]
    frame = frame.iloc[::-1]
    y = np.arange(len(frame))
    positive = frame["atlas_fdr_positive"].astype(float).to_numpy()
    negative = frame["atlas_fdr_negative"].astype(float).to_numpy()
    fig, ax = plt.subplots(figsize=(8.5, 8.2))
    ax.barh(y, -negative, color=RED, alpha=0.9, label="Negative local covariance")
    ax.barh(y, positive, color=BLUE, alpha=0.9, label="Positive local covariance")
    ax.axvline(0, color="#444444", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(frame["label"], fontsize=8)
    ax.set_xlabel("Number of atlas-FDR-significant regions (q < 0.05)")
    ax.set_title("Figure 2C. Directional local signal counts", loc="left", fontsize=12, fontweight="bold")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    style(ax)
    save(fig, "Figure2C_directional_region_counts")


def main():
    atlas = pd.read_csv(RESULTS / "SUPERGNOVA_OCULAR_SYSTEMIC_ATLAS.tsv", sep="\t")
    classification = pd.read_csv(RESULTS / "GLOBAL_LOCAL_ARCHITECTURE_CLASSIFICATION.tsv", sep="\t")
    burden = pd.read_csv(RESULTS / "LOCAL_COVARIANCE_BURDEN.tsv", sep="\t")
    direction = pd.read_csv(RESULTS / "DIRECTIONAL_LOCAL_SIGNAL_SUMMARY.tsv", sep="\t")
    figure_2a(atlas)
    figure_2b(classification, burden)
    figure_2c(direction)


if __name__ == "__main__":
    main()
