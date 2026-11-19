#!/usr/bin/env python3
"""Two-panel apples-to-apples de Goede vs TOPMed lncRNA comparison.

a  Effect size: de Goede median of outlier median_Z vs TOPMed median z
   among people with |z| > 2 (same 'outlier-only' selection).
b  Prevalence: both sides count people with |Z| > 2; -log10(n/N), >=2 each.

Overlap: Table S5/S6 genes in all-gene HVG that are constraint lncRNA.
Does not overwrite Source_Data_Fig3a_ptv_vs_underlier_gene.tsv or the
z < -2-only prevalence figure.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(
    0,
    str(
        Path("/net/topmed11/working/porchard/for-others/julie-and-esther/"
             "2025-11-17/sheyang/pre_Watershed_scripts/scripts/helpers")
    ),
)
from paper_mpl import NATURE_DOUBLE_IN, OI, PAPER_DPI, apply_paper_style  # noqa: E402

ROOT = Path("/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts")
DRAFT = Path("/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026")
SRC = DRAFT / "source_panels"
SUPP = DRAFT / "supplemental_data"
OUT_PIPE = ROOT / "allgenes_pipeline_out"
PUB_FIG3A = SUPP / "Source_Data_Fig3a_ptv_vs_underlier_gene.tsv"
OVERLAP = OUT_PIPE / "degoede_overlap_topmed_median_z.tsv"
STEM = "Fig3a_degoede_vs_topmed_lncrna_zmatch"

N_TOPMED = 5489
N_DEGOEDE = 838 - 39
MIN_N = 2
DIR_COLORS = {
    "Under (Z < 0)": OI["blue"],
    "Over (Z > 0)": OI["orange"],
    "Mixed / 0": "#8A8A8A",
}


def sunflower(k: int, radius: float) -> tuple[np.ndarray, np.ndarray]:
    if k == 1:
        return np.zeros(1), np.zeros(1)
    i = np.arange(k, dtype=float)
    r = radius * np.sqrt((i + 0.5) / k)
    theta = i * np.pi * (3.0 - np.sqrt(5.0))
    return r * np.cos(theta), r * np.sin(theta)


def spread_lattice(d: pd.DataFrame, n_x: str, n_y: str) -> pd.DataFrame:
    d = d.copy()
    d["x_jit"] = d.x
    d["y_jit"] = d.y
    for (_, _), idx in d.groupby([n_x, n_y]).groups.items():
        k = len(idx)
        if k == 1:
            continue
        nx = int(d.loc[idx[0], n_x])
        ny = int(d.loc[idx[0], n_y])
        spacing = min(np.log10((nx + 1) / nx), np.log10((ny + 1) / ny))
        radius = min(0.45 * spacing, 0.014 * np.sqrt(k))
        dx, dy = sunflower(k, radius)
        d.loc[idx, "x_jit"] = d.loc[idx, "x"].to_numpy() + dx
        d.loc[idx, "y_jit"] = d.loc[idx, "y"].to_numpy() + dy
    return d


def direction(z: float) -> str:
    if z < 0:
        return "Under (Z < 0)"
    if z > 0:
        return "Over (Z > 0)"
    return "Mixed / 0"


def main() -> int:
    apply_paper_style()
    g = pd.read_csv(OVERLAP, sep="\t")
    g["direction"] = g.median_Z_degoede.map(direction)
    g["prev_topmed"] = g.n_abs_z_gt2 / N_TOPMED
    g["prev_degoede"] = g.n_degoede / N_DEGOEDE

    a = g[np.isfinite(g.median_Z_degoede) & np.isfinite(g.median_z_among_abs2)].copy()
    b = g[(g.n_degoede >= MIN_N) & (g.n_abs_z_gt2 >= MIN_N)].copy()
    b["x"] = -np.log10(b.prev_topmed)
    b["y"] = -np.log10(b.prev_degoede)

    spa = float(a.median_z_among_abs2.corr(a.median_Z_degoede, method="spearman"))
    pea = float(np.corrcoef(a.median_z_among_abs2, a.median_Z_degoede)[0, 1])
    spb = float(b.prev_topmed.corr(b.prev_degoede, method="spearman"))
    peb = float(np.corrcoef(b.x, b.y)[0, 1])
    print(f"panel a effect size n={len(a)}  Spearman={spa:.2f}  Pearson={pea:.2f}")
    print(f"panel b |Z|>2 prevalence n={len(b)}  Spearman={spb:.2f}  Pearson(-log10)={peb:.2f}")
    print("a direction", a.direction.value_counts().to_dict())
    print("b direction", b.direction.value_counts().to_dict())

    b = spread_lattice(b, "n_abs_z_gt2", "n_degoede")

    fig, axes = plt.subplots(
        1, 2, figsize=(NATURE_DOUBLE_IN, NATURE_DOUBLE_IN * 0.52)
    )

    # --- a: median Z vs median Z ---
    ax = axes[0]
    lo = float(min(a.median_z_among_abs2.min(), a.median_Z_degoede.min()))
    hi = float(max(a.median_z_among_abs2.max(), a.median_Z_degoede.max()))
    pad = 0.08 * (hi - lo)
    lim = (lo - pad, hi + pad)
    ax.plot(lim, lim, ls="--", color="grey", lw=1.0, zorder=0)
    ax.axhline(0, color="0.85", lw=0.6, zorder=0)
    ax.axvline(0, color="0.85", lw=0.6, zorder=0)
    sh = a.sample(frac=1.0, random_state=1)
    ax.scatter(
        sh.median_z_among_abs2,
        sh.median_Z_degoede,
        s=18,
        c=[DIR_COLORS[t] for t in sh.direction],
        alpha=0.85,
        edgecolor="white",
        linewidth=0.2,
        zorder=2,
    )
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("TOPMed median $z$ among $|z| > 2$")
    ax.set_ylabel("de Goede median of outlier median $Z$")
    ax.set_title(
        f"n={len(a)}  |  Spearman={spa:.2f}",
        loc="left", fontsize=9, fontweight="bold", pad=4,
    )
    ax.text(-0.12, 1.08, "a", transform=ax.transAxes,
            fontweight="bold", fontsize=14, va="bottom")

    # --- b: |Z|>2 prevalence ---
    ax = axes[1]
    blo = float(min(b.x_jit.min(), b.y_jit.min()))
    bhi = float(max(b.x_jit.max(), b.y_jit.max()))
    bpad = 0.06 * (bhi - blo)
    blim = (blo - bpad, bhi + bpad)
    ax.plot(blim, blim, ls="--", color="grey", lw=1.0, zorder=0)
    sh = b.sample(frac=1.0, random_state=1)
    ax.scatter(
        sh.x_jit,
        sh.y_jit,
        s=18,
        c=[DIR_COLORS[t] for t in sh.direction],
        alpha=0.85,
        edgecolor="white",
        linewidth=0.2,
        zorder=2,
    )
    ax.set_xlim(*blim)
    ax.set_ylim(*blim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$-\log_{10}$ TOPMed prevalence ($|z| > 2$)")
    ax.set_ylabel(r"$-\log_{10}$ de Goede prevalence ($|Z| > 2$)")
    ax.set_title(
        f"n={len(b)}  |  Spearman={spb:.2f}, Pearson(-log10)={peb:.2f}",
        loc="left", fontsize=9, fontweight="bold", pad=4,
    )
    ax.text(-0.12, 1.08, "b", transform=ax.transAxes,
            fontweight="bold", fontsize=14, va="bottom")

    handles = [
        Line2D([0], [0], linestyle="--", color="grey", lw=1.0, label="y = x"),
    ] + [
        Line2D(
            [0], [0], marker="o", color="none",
            markerfacecolor=DIR_COLORS[t], markeredgecolor="none",
            markersize=7, label=t,
        )
        for t in DIR_COLORS
        if (a.direction == t).any()
    ]
    axes[1].legend(
        handles=handles, loc="upper left", frameon=False, fontsize=8,
        borderaxespad=0.2, labelspacing=0.4, handletextpad=0.4,
    )
    fig.tight_layout()

    SRC.mkdir(parents=True, exist_ok=True)
    SUPP.mkdir(parents=True, exist_ok=True)
    src_tsv = SUPP / f"Source_Data_{STEM}.tsv"
    if src_tsv.resolve() == PUB_FIG3A.resolve():
        raise SystemExit("Refusing to overwrite published Source_Data_Fig3a")

    for folder in (DRAFT, SRC, OUT_PIPE):
        folder.mkdir(parents=True, exist_ok=True)
        for ext in ("png", "pdf"):
            fig.savefig(folder / f"{STEM}.{ext}", dpi=PAPER_DPI, bbox_inches="tight")
    plt.close(fig)

    g.assign(
        N_topmed=N_TOPMED,
        N_degoede=N_DEGOEDE,
        in_panel_a=np.isfinite(g.median_Z_degoede) & np.isfinite(g.median_z_among_abs2),
        in_panel_b=(g.n_degoede >= MIN_N) & (g.n_abs_z_gt2 >= MIN_N),
    ).to_csv(src_tsv, sep="\t", index=False)
    print("wrote", DRAFT / f"{STEM}.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
