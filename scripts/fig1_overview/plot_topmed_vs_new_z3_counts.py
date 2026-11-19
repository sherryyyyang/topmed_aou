#!/usr/bin/env python
"""Count-aware redraw of the TOPMed-vs-AoU haplo+z<-3 prevalence panel.

Prevalence is n_underliers / N_subjects with integer n, so every gene lands on a
lattice point and the plain scatter hides how many genes share a point. Here one
marker is drawn per distinct (n_topmed, n_new) pair, sized and coloured by the
number of genes on it, so the top-right corner reads as 28 positions holding 624
genes rather than as ~28 genes.

Same gene set as panel B of topmed_vs_new_haplo_prevalence_scatter_refined.png:
>= 2 underliers in both cohorts.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

import sys
sys.path.insert(0, "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts/scripts/helpers")
from _paths import ROOT, FIG1A, FIG1B, FIG1C, FIG4D, SUP_GN, DATA_CORE, DATA_UNDERLIER, AOU_RNA
HERE = str(FIG1C)
SRC = str(FIG1C / "topmed_vs_new_z_-3_prevalence_per_gene.tsv")
CUT = 3.0          # -log10 prevalence threshold defining the rare-outlier corner
MIN_N = 2          # underliers required in both cohorts

d = pd.read_csv(SRC, sep="\t")
N_TOP, N_AOU = int(d.N_topmed.max()), int(d.N_new.max())
d = d[(d.n_topmed >= MIN_N) & (d.n_new >= MIN_N)].copy()
d["x"] = -np.log10(d.n_topmed / N_TOP)
d["y"] = -np.log10(d.n_new / N_AOU)

pos = (d.groupby(["n_topmed", "n_new", "x", "y"])
         .agg(n_genes=("gene_nv", "size"))
         .reset_index())
corner = pos[(pos.x >= CUT) & (pos.y >= CUT)]
n_corner_genes = int(corner.n_genes.sum())

print(f"genes plotted (>= {MIN_N} underliers in both): {len(d):,} at {len(pos)} lattice positions")
print(f"corner (-log10 >= {CUT} on both axes): {n_corner_genes} genes at {len(corner)} positions")
print(f"most crowded position overall: {int(pos.n_genes.max())} genes; in the corner: "
      f"{int(corner.n_genes.max())} genes")
print(f"positions holding a single gene: {(pos.n_genes == 1).sum()} of {len(pos)}")

# marker area scaled to gene count; sqrt keeps area, not radius, proportional
size = lambda n: 12 + 34 * np.sqrt(n)
norm = LogNorm(vmin=1, vmax=pos.n_genes.max())
CMAP = "viridis"

fig, axes = plt.subplots(1, 2, figsize=(15.5, 7.0),
                         gridspec_kw=dict(width_ratios=[1.25, 1]))

ax = axes[0]
sc = ax.scatter(pos.x, pos.y, s=size(pos.n_genes), c=pos.n_genes, cmap=CMAP,
                norm=norm, alpha=0.85, edgecolor="white", linewidth=0.3)
lim = [min(pos.x.min(), pos.y.min()) - 0.08, max(pos.x.max(), pos.y.max()) + 0.08]
ax.plot(lim, lim, ls="--", lw=0.9, color="grey", zorder=0, label="y = x")
ax.axvline(CUT, ls=":", lw=0.9, color="firebrick")
ax.axhline(CUT, ls=":", lw=0.9, color="firebrick")
ax.add_patch(plt.Rectangle((CUT, CUT), lim[1] - CUT, lim[1] - CUT, facecolor="firebrick",
                           alpha=0.05, edgecolor="none", zorder=0))
ax.text(lim[1] - 0.03, lim[0] + 0.03,
        f"shaded box = rare-outlier corner:\n{n_corner_genes} genes on {len(corner)} points\n"
        f"({(pos.n_genes == 1).sum()} of {len(pos)} positions panel-wide hold 1 gene)",
        ha="right", va="bottom", fontsize=8.5, color="firebrick")
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_title("A. One marker per distinct outlier-count pair\n"
             f"{len(d):,} genes collapse onto {len(pos)} positions", fontsize=11, loc="left")
fig.colorbar(sc, ax=ax, label="genes at this position", fraction=0.046, pad=0.02)

ax = axes[1]
ax.scatter(corner.x, corner.y, s=size(corner.n_genes) * 2.6, c=corner.n_genes, cmap=CMAP,
           norm=norm, alpha=0.9, edgecolor="white", linewidth=0.4)
for _, r in corner.iterrows():
    ax.annotate(int(r.n_genes), (r.x, r.y), ha="center", va="center", fontsize=7.5,
                color="white" if r.n_genes > 20 else "black", fontweight="bold")
xt = sorted(corner.x.unique()); yt = sorted(corner.y.unique())
ax.set_xticks(xt); ax.set_yticks(yt)
ax.set_xticklabels([f"{v:.2f}\n{int(round(10**-v * N_TOP))}/{N_TOP:,}" for v in xt], fontsize=7.5)
ax.set_yticklabels([f"{v:.2f}  {int(round(10**-v * N_AOU))}/{N_AOU:,}" for v in yt], fontsize=7.5)
ax.plot([CUT, max(xt) + 0.05], [CUT, max(xt) + 0.05], ls="--", lw=0.9, color="grey", zorder=0)
ax.set_title("B. The corner, labelled with gene counts\n"
             f"{len(corner)} positions, {n_corner_genes} distinct genes "
             "(the 3-dot / 4-dot figure gene set)", fontsize=11, loc="left")
ax.set_xlim(min(xt) - 0.05, max(xt) + 0.05)
ax.set_ylim(min(yt) - 0.04, max(yt) + 0.06)

for ax in axes:
    ax.set_xlabel(f"-log10 TOPMed RNA-outlier prevalence  (n = {N_TOP:,})")
    ax.set_ylabel(f"-log10 AoU RNA-outlier prevalence  (n = {N_AOU:,})")
    ax.grid(True, ls=":", lw=0.3, alpha=0.5)
    ax.set_axisbelow(True)

fig.suptitle("Why the rare-outlier corner looks sparse: prevalence is quantised, so genes stack exactly "
             "(haplo + z<-3, >= 2 underliers in both cohorts)", fontsize=12, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.95))
for ext in ("png", "pdf"):
    out = str(FIG1C / f"topmed_vs_new_z3_corner_counts.{ext}")
    fig.savefig(out, dpi=200 if ext == "png" else None, bbox_inches="tight")
    print("wrote", out)

pos.sort_values("n_genes", ascending=False).to_csv(
    os.path.join(HERE, "topmed_vs_new_z3_lattice_counts.tsv"), sep="\t", index=False)
print("wrote topmed_vs_new_z3_lattice_counts.tsv")
