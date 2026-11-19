"""Matplotlib counterpart of paper_theme.R (Arial, Word-readable, Okabe–Ito)."""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.colors as mcolors

PAPER_FONT = "Arial"
NATURE_DOUBLE_IN = 183 / 25.4
PAPER_DPI = 300

# Okabe–Ito; sequential colour uses viridis.
OI = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
}
OI_CYCLE = list(OI.values())
COL_TOPMED = OI["blue"]
COL_AOU = OI["purple"]
COL_PTV = OI["vermillion"]
COL_NON = "#8A8A8A"
TIER_COLORS = {
    "Constrained": OI["vermillion"],
    "Mid-low": OI["orange"],
    "Mid-high": OI["sky"],
    "Tolerant": OI["blue"],
    "No LOEUF": OI["black"],
}


def tint(hex_color: str, white: float = 0.82) -> str:
    rgb = mcolors.to_rgb(hex_color)
    mixed = tuple(c + (1.0 - c) * white for c in rgb)
    return mcolors.to_hex(mixed)


def apply_paper_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [PAPER_FONT, "Helvetica", "DejaVu Sans"],
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "legend.title_fontsize": 9,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "figure.dpi": PAPER_DPI,
            "savefig.dpi": PAPER_DPI,
            "savefig.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.prop_cycle": mpl.cycler(color=OI_CYCLE),
            "image.cmap": "viridis",
        }
    )
