#!/usr/bin/env python3
"""Study-design schematic for Figure 1a. Short banner, not a tall hero panel.
Headlines belong in the Word caption, not on the graphic.
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(
    0,
    "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts/scripts/helpers",
)
from _paths import FIG1A  # noqa: E402
from paper_mpl import (  # noqa: E402
    COL_AOU,
    COL_PTV,
    COL_TOPMED,
    NATURE_DOUBLE_IN,
    PAPER_DPI,
    apply_paper_style,
    tint,
)

OUT = str(FIG1A / "fig1_study_design_topmed_aou")


def box(ax, x, y, w, h, fc, ec, lw=1.2):
    p = FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.02",
        facecolor=fc, edgecolor=ec, linewidth=lw, transform=ax.transAxes,
        clip_on=False,
    )
    ax.add_patch(p)
    return p


def arrow(ax, x1, y1, x2, y2):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1), (x2, y2), transform=ax.transAxes,
            arrowstyle="-|>", mutation_scale=10, linewidth=1.15,
            color="#666666", shrinkA=0, shrinkB=0, clip_on=False,
        )
    )


def main() -> int:
    apply_paper_style()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    # Two-line cards: keep type readable without a tall banner.
    fig, ax = plt.subplots(figsize=(NATURE_DOUBLE_IN, 0.72))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    box(ax, 0.008, 0.10, 0.292, 0.80, tint(COL_TOPMED), COL_TOPMED, 1.2)
    ax.text(0.154, 0.64, "TOPMed  ·  n = 5,489", ha="center", va="center",
            fontsize=11, fontweight="bold", color=COL_TOPMED,
            transform=ax.transAxes, fontfamily="Arial")
    ax.text(0.154, 0.32, "discovery  ·  WGS + blood RNA", ha="center", va="center",
            fontsize=9, color="#333333", transform=ax.transAxes, fontfamily="Arial")

    box(ax, 0.312, 0.10, 0.292, 0.80, tint(COL_AOU), COL_AOU, 1.2)
    ax.text(0.458, 0.64, "All of Us  ·  n = 8,895", ha="center", va="center",
            fontsize=11, fontweight="bold", color=COL_AOU,
            transform=ax.transAxes, fontfamily="Arial")
    ax.text(0.458, 0.32, "replication  ·  WGS + blood RNA", ha="center", va="center",
            fontsize=9, color="#333333", transform=ax.transAxes, fontfamily="Arial")

    arrow(ax, 0.610, 0.50, 0.658, 0.50)

    box(ax, 0.665, 0.10, 0.327, 0.80, "#F7F7F7", "#333333", 1.15)
    ax.text(0.828, 0.64, r"$Z < -3$  and  log$_2$CPM < mean $-$ 1",
            ha="center", va="center", fontsize=10, fontweight="bold", color=COL_PTV,
            transform=ax.transAxes, fontfamily="Arial")
    ax.text(0.828, 0.32, "same underlier call  ·  \u2248 het. LoF dosage",
            ha="center", va="center", fontsize=9, color="#333333",
            transform=ax.transAxes, fontfamily="Arial")

    fig.subplots_adjust(left=0.004, right=0.996, top=0.96, bottom=0.04)
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}.{ext}", dpi=PAPER_DPI, bbox_inches="tight", pad_inches=0.02)
        print("wrote", f"{OUT}.{ext}")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
