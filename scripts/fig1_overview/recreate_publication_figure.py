"""Create an editable, publication-style SVG/PDF from the supplied figure.

Run: python3 recreate_publication_figure.py

Counts and labels live in the CONFIG section below. The TOPMed cohort labels
currently sum to 6,108 while the stated matched total is 5,489. The figure
preserves both supplied values and marks that discrepancy for confirmation.
"""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, PathPatch, FancyBboxPatch, Rectangle
from matplotlib.path import Path as MplPath


# ------------------------------- CONFIG ------------------------------------
OUT = Path(__file__).resolve().parent / "assets"
TOPMED_MATCHED = 5489
TOPMED_COHORTS = {
    "FHS": 792,
    "GALAII": 1827,
    "SAGE": 675,
    "SPIROMICS": 1538,
    "WHI": 1276,
}
MISSING_EXCLUSIONS = "~500 excluded for missing data"
AOU_MATCHED = 9043
# ---------------------------------------------------------------------------

NAVY = "#24316F"
TEAL = "#209A7A"
INK = "#252A34"
MUTED = "#667085"
LIGHT = "#E9EDF2"
PALETTE = ["#5AAE9A", "#238F78", "#75B7A6", "#2F806F", "#9ACBBC"]


def ribbon(ax, x0, y0, x1, y1, thickness, color, alpha=0.9):
    """Curved flow ribbon between two vertical positions."""
    c = (x1 - x0) * 0.48
    verts = [
        (x0, y0 + thickness / 2), (x0 + c, y0 + thickness / 2),
        (x1 - c, y1 + thickness / 2), (x1, y1 + thickness / 2),
        (x1, y1 - thickness / 2), (x1 - c, y1 - thickness / 2),
        (x0 + c, y0 - thickness / 2), (x0, y0 - thickness / 2),
        (x0, y0 + thickness / 2),
    ]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color,
                           edgecolor="none", alpha=alpha, zorder=1))


def draw_person(ax, x, y, color, scale=1.0):
    """Small flat participant pictogram centered at x,y."""
    r = 0.105 * scale
    ax.add_patch(Circle((x, y + 0.46 * scale), r, color=color, ec="none", zorder=3))
    body = FancyBboxPatch((x - 0.13 * scale, y + 0.02 * scale),
                          0.26 * scale, 0.34 * scale,
                          boxstyle=f"round,pad=0.02,rounding_size={0.08*scale}",
                          color=color, ec="none", zorder=2)
    ax.add_patch(body)
    ax.plot([x - 0.07 * scale, x - 0.07 * scale], [y + 0.03 * scale, y - 0.31 * scale],
            color=color, lw=5 * scale, solid_capstyle="round", zorder=2)
    ax.plot([x + 0.07 * scale, x + 0.07 * scale], [y + 0.03 * scale, y - 0.31 * scale],
            color=color, lw=5 * scale, solid_capstyle="round", zorder=2)
    ax.plot([x - 0.11 * scale, x - 0.24 * scale], [y + 0.27 * scale, y + 0.02 * scale],
            color=color, lw=5 * scale, solid_capstyle="round", zorder=2)
    ax.plot([x + 0.11 * scale, x + 0.24 * scale], [y + 0.27 * scale, y + 0.02 * scale],
            color=color, lw=5 * scale, solid_capstyle="round", zorder=2)


def draw_dna(ax, x, y, color=NAVY):
    import numpy as np
    t = np.linspace(0, 2.5 * np.pi, 120)
    xx = x + 0.43 * t / (2.5 * np.pi)
    yy1 = y + 0.11 * np.sin(t)
    yy2 = y - 0.11 * np.sin(t)
    ax.plot(xx, yy1, color=color, lw=2.2, solid_capstyle="round")
    ax.plot(xx, yy2, color="#36A6D8", lw=2.2, solid_capstyle="round")
    for i in range(0, len(t), 9):
        ax.plot([xx[i], xx[i]], [yy1[i], yy2[i]], color="#7A91BD", lw=1.1)


def main():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9,
        "svg.fonttype": "none",  # keep text editable in Illustrator/Inkscape
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(12.4, 5.4), facecolor="white")
    ax = fig.add_axes([0.025, 0.04, 0.95, 0.92])
    ax.set_xlim(0, 12.4)
    ax.set_ylim(0, 5.4)
    ax.axis("off")

    # Panel titles and subtle divider
    ax.text(0.25, 5.08, "TOPMed", fontsize=15, weight="bold", color=INK, va="center")
    ax.text(0.25, 4.78, "Matched WGS + bulk RNA-seq samples", fontsize=9.5, color=MUTED, va="center")
    ax.text(6.45, 5.08, "All of Us", fontsize=15, weight="bold", color=INK, va="center")
    ax.text(6.45, 4.78, "Matched WGS + bulk RNA-seq samples", fontsize=9.5, color=MUTED, va="center")
    ax.plot([6.15, 6.15], [0.35, 4.95], color=LIGHT, lw=1)

    # TOPMed panel: source cohorts flow into a shared bulk RNA-seq set.
    ax.text(0.28, 4.30, "WGS cohorts", color=MUTED, fontsize=8, weight="bold")
    ax.text(4.58, 4.30, "Whole-blood RNA-seq", color=MUTED, fontsize=8,
            weight="bold", ha="center")
    total_labels = sum(TOPMED_COHORTS.values())
    max_val = max(TOPMED_COHORTS.values())
    ys = [3.75, 3.10, 2.45, 1.80, 1.15]
    target_y = 2.47
    scale = 1.18 / max_val
    for (name, count), y, color in zip(TOPMED_COHORTS.items(), ys, PALETTE):
        th = max(0.045, count * scale)
        ax.text(0.28, y, f"{name}  {count:,}", ha="left", va="center", color=INK, fontsize=9)
        ax.add_patch(Rectangle((2.28, y - 0.075), 0.07, 0.15, color="#A8AFB8", ec="none", zorder=3))
        ribbon(ax, 2.35, y, 4.30, target_y + (y - target_y) * 0.15, th, color)
    ax.add_patch(Rectangle((4.30, 1.58), 0.12, 1.78, color="#A8AFB8", ec="none", zorder=3))
    ax.text(4.52, 2.48, f"Listed cohorts\n{total_labels:,}", color=INK, fontsize=9,
            va="center", ha="left", weight="bold")
    ax.text(0.28, 0.59, f"Stated matched total  N = {TOPMED_MATCHED:,}",
            fontsize=12, color=NAVY, weight="bold")
    ax.text(0.28, 0.29, MISSING_EXCLUSIONS, fontsize=8.5, color=MUTED)
    if total_labels != TOPMED_MATCHED:
        ax.text(5.88, 0.29, f"Check counts: listed cohorts exceed matched N by {total_labels-TOPMED_MATCHED:,}",
                fontsize=7.2, color="#A45B00", ha="right")

    # All of Us panel: concise acquisition-to-match sequence.
    draw_person(ax, 7.18, 2.55, "#65A982", 0.78)
    draw_person(ax, 7.55, 2.55, "#E16A83", 0.78)
    draw_person(ax, 7.92, 2.55, "#43A0B6", 0.78)
    ax.text(7.55, 1.95, "Participants", ha="center", color=INK, fontsize=9, weight="bold")
    ax.annotate("", xy=(8.62, 2.75), xytext=(8.30, 2.75),
                arrowprops=dict(arrowstyle="-|>", color="#9AA3AF", lw=1.5))

    # Two omics inputs
    ax.add_patch(FancyBboxPatch((8.72, 3.00), 1.23, 0.78,
                                boxstyle="round,pad=0.08,rounding_size=0.08",
                                facecolor="#F4F6FB", edgecolor="#DCE2EF", lw=1))
    draw_dna(ax, 8.84, 3.48)
    ax.text(9.42, 3.39, "WGS", fontsize=9, weight="bold", color=INK, va="center")
    ax.text(9.42, 3.17, "Genome", fontsize=7.5, color=MUTED, va="center")
    ax.add_patch(FancyBboxPatch((8.72, 1.93), 1.23, 0.78,
                                boxstyle="round,pad=0.08,rounding_size=0.08",
                                facecolor="#F2FAF8", edgecolor="#D7ECE5", lw=1))
    # RNA glyph
    import numpy as np
    t = np.linspace(0, 2 * np.pi, 50)
    ax.plot(8.83 + 0.37 * t/(2*np.pi), 2.36 + 0.075*np.sin(t), color=TEAL, lw=2)
    ax.text(9.42, 2.35, "Bulk RNA-seq", fontsize=8, weight="bold", color=INK, va="center")
    ax.text(9.42, 2.13, "Transcriptome", fontsize=7.2, color=MUTED, va="center")
    ax.annotate("", xy=(10.42, 2.75), xytext=(10.05, 2.75),
                arrowprops=dict(arrowstyle="-|>", color="#9AA3AF", lw=1.5))
    ax.add_patch(FancyBboxPatch((10.55, 2.32), 1.42, 0.88,
                                boxstyle="round,pad=0.10,rounding_size=0.10",
                                facecolor=NAVY, edgecolor="none"))
    ax.text(11.26, 2.86, f"N = {AOU_MATCHED:,}", ha="center", va="center",
            fontsize=14, weight="bold", color="white")
    ax.text(11.26, 2.57, "matched samples", ha="center", va="center",
            fontsize=7.5, color="white")
    ax.text(6.55, 0.59, "Whole-genome sequencing + bulk RNA-seq", fontsize=9, color=INK)
    ax.text(6.55, 0.29, "Participant illustration and logos omitted for a clean, scalable layout",
            fontsize=7.5, color=MUTED)

    fig.savefig(OUT / "publication_figure.svg", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / "publication_figure.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / "publication_figure.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Created {OUT / 'publication_figure.svg'}")
    print(f"Created {OUT / 'publication_figure.pdf'}")
    print(f"Created {OUT / 'publication_figure.png'}")


if __name__ == "__main__":
    main()
