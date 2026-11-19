#!/usr/bin/env python3
"""Figure 2d (historical matplotlib renderer).

Canonical paper panel is plot_topmed_vs_aou_prevalence_loeuf.R (theme_paper).
Haplo + z < -3; genes with >=2 underliers in both cohorts. No headline.
"""
from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(
    0,
    "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts/scripts/helpers",
)
from _paths import FIG1C, FIG2C  # noqa: E402
from paper_mpl import PAPER_DPI, TIER_COLORS, apply_paper_style  # noqa: E402

MIN_N = 2
TIER_ORDER = ["Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF"]
TIER_LABELS = {
    "Constrained": "Constrained",
    "Mid-low": "Mid-low",
    "Mid-high": "Mid-high",
    "Tolerant": "Tolerant",
    "No LOEUF": "No LOEUF",
}


def sunflower(k: int, radius: float) -> tuple[np.ndarray, np.ndarray]:
    """Deterministic Vogel spiral packing of k points inside `radius`."""
    if k == 1:
        return np.zeros(1), np.zeros(1)
    i = np.arange(k, dtype=float)
    r = radius * np.sqrt((i + 0.5) / k)
    theta = i * np.pi * (3.0 - np.sqrt(5.0))
    return r * np.cos(theta), r * np.sin(theta)


def spread_lattice(d: pd.DataFrame) -> pd.DataFrame:
    """Fan genes that share an (n_topmed, n_new) cell out around that node.

    Prevalence is n/N, so genes with identical counts land on one pixel and
    the last tier drawn hides the rest. Offsets stay inside ~half the local
    lattice spacing, which shrinks as counts rise and nodes crowd together.
    """
    d = d.copy()
    d["x_jit"] = d.x
    d["y_jit"] = d.y
    for (n_tm, n_aou), idx in d.groupby(["n_topmed", "n_new"]).groups.items():
        k = len(idx)
        if k == 1:
            continue
        spacing = min(
            np.log10((n_tm + 1) / n_tm),
            np.log10((n_aou + 1) / n_aou),
        )
        radius = min(0.45 * spacing, 0.014 * np.sqrt(k))
        dx, dy = sunflower(k, radius)
        d.loc[idx, "x_jit"] = d.loc[idx, "x"].to_numpy() + dx
        d.loc[idx, "y_jit"] = d.loc[idx, "y"].to_numpy() + dy
    return d


def main() -> int:
    apply_paper_style()
    FIG2C.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(FIG1C / "topmed_vs_new_z_-3_prevalence_per_gene.tsv", sep="\t")
    d = d[(d.n_topmed >= MIN_N) & (d.n_new >= MIN_N)].copy()
    d["LOEUF_tier"] = d["LOEUF_tier"].fillna("No LOEUF")
    d.loc[~d.LOEUF_tier.isin(TIER_ORDER), "LOEUF_tier"] = "No LOEUF"
    d["x"] = -np.log10(d.prev_topmed)
    d["y"] = -np.log10(d.prev_new)
    sp = float(d.prev_topmed.corr(d.prev_new, method="spearman"))
    pe = float(np.corrcoef(d.x, d.y)[0, 1])
    d = spread_lattice(d)

    lo = float(min(d.x_jit.min(), d.y_jit.min()))
    hi = float(max(d.x_jit.max(), d.y_jit.max()))
    pad = 0.06 * (hi - lo)
    lim = (lo - pad, hi + pad)

    fig, ax = plt.subplots(figsize=(6.3, 6.3))
    ax.plot(lim, lim, ls="--", color="grey", lw=1.0, zorder=0)
    # One shuffled pass so no tier is systematically painted over the others.
    shuffled = d.sample(frac=1.0, random_state=1)
    ax.scatter(
        shuffled.x_jit, shuffled.y_jit, s=9,
        c=[TIER_COLORS[t] for t in shuffled.LOEUF_tier],
        alpha=0.8, edgecolor="white", linewidth=0.15, zorder=2,
    )
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$-\log_{10}$ TOPMed RNA-outlier prevalence")
    ax.set_ylabel(r"$-\log_{10}$ AoU RNA-outlier prevalence")
    print(
        f"caption stats: n_genes={len(d):,}; Spearman={sp:.2f}; "
        f"Pearson(-log10)={pe:.2f}"
    )
    handles = [
        Line2D(
            [0], [0], marker="o", color="none",
            markerfacecolor=TIER_COLORS[t], markeredgecolor="none",
            markersize=7, label=TIER_LABELS[t],
        )
        for t in TIER_ORDER
    ]
    ax.legend(
        handles=handles, loc="upper left", frameon=False, fontsize=8,
        title="LOEUF tier", title_fontsize=9, borderaxespad=0.2,
        labelspacing=0.45, handletextpad=0.4, markerscale=1.1,
    )
    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.12, top=0.98)
    out = FIG2C / "fig2_topmed_vs_aou_prevalence_loeuf"
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=PAPER_DPI)
        print("wrote", f"{out}.{ext}")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
