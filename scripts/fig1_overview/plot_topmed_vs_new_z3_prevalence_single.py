#!/usr/bin/env python3
"""Standalone version of panel B (haplo + z<-3) from the two-panel figure
topmed_vs_new_haplo_prevalence_scatter_refined.png.

Reuses the merge/QC logic in compare_topmed_vs_new_haplo_prevalence.py so the
numbers are identical to the panel it is taken from. Adds secondary axes labeled
as "1 in N" carriers, since -log10(prevalence) increases toward rarer genes.
"""
from __future__ import annotations

import os

import numpy as np

from compare_topmed_vs_new_haplo_prevalence import (
    JOBS,
    PLOT_ORDER,
    TIER_COLORS,
    TIER_LEGEND_LABELS,
    TIER_LEVELS,
    load_loeuf_tiers,
    merge_job,
)

import sys
sys.path.insert(0, "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts/scripts/helpers")
from _paths import ROOT, FIG1A, FIG1B, FIG1C, FIG4D, SUP_GN, DATA_CORE, DATA_UNDERLIER, AOU_RNA
OUT_BASE = str(FIG1C / "topmed_vs_new_z_-3_prevalence_scatter_single")


def main() -> int:
    loeuf = load_loeuf_tiers()
    pan = merge_job(JOBS[1], loeuf)          # JOBS[1] is haplo + z < -3
    df = pan["df"]

    ratio = df["ratio_new_over_topmed"].replace([np.inf, -np.inf], np.nan).dropna()
    print(f"\npanel B: n={pan['n']} genes  Spearman={pan['sp']:.3f}  "
          f"Pearson(-log10)={pan['pe']:.3f}")
    print(f"AoU/TOPMed prevalence ratio: median={ratio.median():.3g}  "
          f"IQR=[{ratio.quantile(.25):.3g}, {ratio.quantile(.75):.3g}]")
    print(f"prevalence floors from the >=2 underlier filter: "
          f"TOPMed 2/{pan['n_topmed_subj']:,} = -log10 {-np.log10(2 / pan['n_topmed_subj']):.2f}; "
          f"AoU 2/{pan['n_new_subj']:,} = -log10 {-np.log10(2 / pan['n_new_subj']):.2f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    x = -np.log10(df["prev_topmed"].to_numpy())
    y = -np.log10(df["prev_new"].to_numpy())
    lo, hi = float(min(x.min(), y.min())), float(max(x.max(), y.max()))
    pad = 0.05 * (hi - lo)

    fig, ax = plt.subplots(figsize=(7.6, 6.6))
    ax.plot([lo, hi], [lo, hi], ls="--", color="grey", lw=0.8, zorder=0)
    counts = df["LOEUF_tier"].value_counts()
    for tier in PLOT_ORDER:
        sub = df[df["LOEUF_tier"] == tier]
        if sub.empty:
            continue
        ax.scatter(
            -np.log10(sub["prev_topmed"]), -np.log10(sub["prev_new"]),
            s=18, alpha=0.55, color=TIER_COLORS[tier], edgecolor="none",
            zorder=1 + PLOT_ORDER.index(tier),
        )
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, which="both", ls=":", lw=0.3, alpha=0.5)
    ax.set_xlabel(r"$-\log_{10}$ TOPMed RNA-outlier prevalence")
    ax.set_ylabel(r"$-\log_{10}$ AoU RNA-outlier prevalence")
    ax.set_title(
        f"RNA-outlier prevalence per gene, haplo + z < \u22123\n"
        f"TOPMed n={pan['n_topmed_subj']:,} vs AoU n={pan['n_new_subj']:,}; "
        f"{pan['n']:,} genes with \u22652 underliers in both  |  "
        f"Spearman={pan['sp']:.2f}, Pearson={pan['pe']:.2f}",
        fontsize=10,
    )

    # secondary axes restating the log scale as "1 in N" carriers
    ident = (lambda v: v, lambda v: v)
    ticks = [t for t in (1.5, 2.0, 2.5, 3.0, 3.5) if lo - pad <= t <= hi + pad]
    labels = [("1 in %s" % f"{10 ** t:,.0f}") if float(t).is_integer()
              else f"1 in {10 ** t:,.0f}" for t in ticks]
    for loc, setter in (("top", "secondary_xaxis"), ("right", "secondary_yaxis")):
        sec = getattr(ax, setter)(loc, functions=ident)
        if setter == "secondary_xaxis":
            sec.set_xticks(ticks); sec.set_xticklabels(labels, fontsize=7.5)
        else:
            sec.set_yticks(ticks); sec.set_yticklabels(labels, fontsize=7.5)
        sec.tick_params(length=2)

    handles = [Line2D([0], [0], ls="--", color="grey", lw=0.8, label="y = x")]
    for tier in TIER_LEVELS:
        n = int(counts.get(tier, 0))
        handles.append(Line2D([0], [0], marker="o", color="none",
                              markerfacecolor=TIER_COLORS[tier], markersize=8,
                              label=f"{TIER_LEGEND_LABELS[tier]}  n={n}"))
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.16, 0.5),
              fontsize=8, frameon=False, borderaxespad=0.0,
              labelspacing=0.55, handletextpad=0.7, markerscale=1.2, handlelength=1.8)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = f"{OUT_BASE}.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print("wrote", out)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
