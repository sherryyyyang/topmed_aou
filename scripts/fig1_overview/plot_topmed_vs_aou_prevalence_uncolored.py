#!/usr/bin/env python3
"""Figure 1d: TOPMed vs AoU RNA-outlier prevalence (haplo + z < -3 only).

No LOEUF coloring. No on-figure headline. One marker per distinct
(n_TOPMed, n_AoU) lattice point; marker area stays modest for 1–10
stacked genes and lifts 25 / 49 so the legend jump is readable.
"""
from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

sys.path.insert(
    0,
    "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts/scripts/helpers",
)
from _paths import FIG1C  # noqa: E402
from paper_mpl import COL_PTV, NATURE_DOUBLE_IN, PAPER_DPI, apply_paper_style  # noqa: E402

BASE_S = 10.0
# Knots: n^0.45 through 10 genes; designed areas at 25 and 49.
_SIZE_N = np.array([1.0, 2.0, 5.0, 10.0, 25.0, 49.0])
_SIZE_A = np.array([
    BASE_S * (1.0 ** 0.45),
    BASE_S * (2.0 ** 0.45),
    BASE_S * (5.0 ** 0.45),
    BASE_S * (10.0 ** 0.45),
    62.0,
    125.0,
])
CUT = 3.0
MIN_N = 2
POINT_COLOR = "#C8C8C8"
CORNER_EDGE = COL_PTV


def load_panel(path) -> tuple[pd.DataFrame, pd.DataFrame, int, int, float, float]:
    d = pd.read_csv(path, sep="\t")
    n_top = int(d.N_topmed.max())
    n_aou = int(d.N_new.max())
    d = d[(d.n_topmed >= MIN_N) & (d.n_new >= MIN_N)].copy()
    d["x"] = -np.log10(d.prev_topmed)
    d["y"] = -np.log10(d.prev_new)
    pos = (
        d.groupby(["n_topmed", "n_new", "x", "y"], as_index=False)
        .agg(n_genes=("gene_nv", "size"))
    )
    sp = float(d.prev_topmed.corr(d.prev_new, method="spearman"))
    pe = float(np.corrcoef(d.x, d.y)[0, 1])
    return d, pos, n_top, n_aou, sp, pe


def size_area(n: np.ndarray | pd.Series) -> np.ndarray:
    # 1–10 stay on the compressed n^0.45 scale. Above 10, interpolate to
    # larger designed areas so 25 and 49 read as a jump (still << old 18*n).
    return np.interp(np.asarray(n, dtype=float), _SIZE_N, _SIZE_A)


def draw_panel(ax, genes, pos, sp: float, pe: float) -> None:
    lo = float(min(pos.x.min(), pos.y.min()))
    hi = float(max(pos.x.max(), pos.y.max()))
    pad = 0.06 * (hi - lo)
    lim = (lo - pad, hi + pad)

    ax.plot(lim, lim, ls="--", color="grey", lw=1.0, zorder=0)
    ax.add_patch(
        Rectangle(
            (CUT, CUT), lim[1] - CUT, lim[1] - CUT,
            facecolor=COL_PTV, alpha=0.07, edgecolor=CORNER_EDGE,
            linewidth=1.1, zorder=0,
        )
    )
    ax.axvline(CUT, ls=":", lw=0.9, color=CORNER_EDGE, zorder=0)
    ax.axhline(CUT, ls=":", lw=0.9, color=CORNER_EDGE, zorder=0)
    rest = pos[(pos.x < CUT) | (pos.y < CUT)]
    corner = pos[(pos.x >= CUT) & (pos.y >= CUT)]
    ax.scatter(
        rest.x, rest.y, s=size_area(rest.n_genes),
        c=POINT_COLOR, alpha=0.42, edgecolor="white", linewidth=0.25, zorder=2,
    )
    ax.scatter(
        corner.x, corner.y, s=size_area(corner.n_genes),
        c=POINT_COLOR, alpha=0.55, edgecolor=CORNER_EDGE, linewidth=0.7, zorder=3,
    )
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$-\log_{10}$ TOPMed RNA-outlier prevalence")
    ax.set_ylabel(r"$-\log_{10}$ AoU RNA-outlier prevalence")


def main() -> int:
    apply_paper_style()
    z3_genes, z3_pos, n_top, n_aou, z3_sp, z3_pe = load_panel(
        FIG1C / "topmed_vs_new_z_-3_prevalence_per_gene.tsv"
    )

    fig, ax = plt.subplots(figsize=(NATURE_DOUBLE_IN * 0.72, 5.6))
    draw_panel(ax, z3_genes, z3_pos, z3_sp, z3_pe)

    max_stack = int(z3_pos.n_genes.max())
    counts = [c for c in (1, 2, 5, 10, 25) if c <= max_stack]
    if max_stack not in counts:
        counts.append(max_stack)
    size_handles = [
        Line2D(
            [0], [0], marker="o", color="none",
            markerfacecolor=POINT_COLOR, markeredgecolor="none",
            alpha=0.55, markersize=np.sqrt(size_area([c])[0]),
            label="1 gene" if c == 1 else f"{c} genes",
        )
        for c in counts
    ]
    fig.legend(
        handles=size_handles, loc="center left", bbox_to_anchor=(0.99, 0.48),
        fontsize=8, frameon=False, labelspacing=2.35, handletextpad=0.8,
        handlelength=2.0, title="genes stacked\non one point", title_fontsize=8,
        prop={"family": "Arial"},
    )
    fig.tight_layout(rect=(0, 0, 0.82, 1.0))

    out_stem = FIG1C / "fig1_topmed_vs_aou_prevalence_uncolored"
    for ext in ("png", "pdf"):
        fig.savefig(f"{out_stem}.{ext}", dpi=PAPER_DPI, bbox_inches="tight")
        print("wrote", f"{out_stem}.{ext}")
    plt.close(fig)
    print(
        f"caption stats: TOPMed n={n_top:,}, AoU n={n_aou:,}; "
        f"n_genes={len(z3_genes):,}; Spearman={z3_sp:.2f}; "
        f"Pearson(-log10)={z3_pe:.2f}; n_corner={int(z3_pos[(z3_pos.x>=CUT)&(z3_pos.y>=CUT)].n_genes.sum())}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
