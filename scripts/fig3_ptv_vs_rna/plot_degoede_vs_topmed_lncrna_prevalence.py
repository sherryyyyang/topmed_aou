#!/usr/bin/env python3
"""Fig 2-style scatter: underlier prevalence in all subjects, overlapping lncRNAs.

Table S6 median_Z_score is treated as a per-person z.
Underlier = z < -2.
de Goede prev = n people with z < -2 / 799 (all GTEx subjects after
global-outlier removal).
TOPMed prev = n kept subjects with z < -2 / n subjects assayed for that gene
(5489). No log2CPM drop.

Universe: S6 genes that are TOPMed all-gene HVG constraint lncRNA.
Scatter: >=1 underlier person on each side.
"""
from __future__ import annotations

import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
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
SHEYANG = ROOT.parent
DRAFT = Path("/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026")
SRC = DRAFT / "source_panels"
SUPP = DRAFT / "supplemental_data"
OUT_PIPE = ROOT / "allgenes_pipeline_out"
PUB_FIG3A = SUPP / "Source_Data_Fig3a_ptv_vs_underlier_gene.tsv"
STEM = "Fig3a_degoede_vs_topmed_lncrna_prevalence"

ZIP_S6 = Path(
    "/net/wonderland/home/sheyang/.cursor/projects/"
    "net-topmed11-working-porchard-for-others-julie-and-esther/uploads/"
    "3018f0bc-ed93-4a42-8e00-bac9bc0c1ab9/"
    "ScienceDirect_files_16Oct2026_20-47-41.902.zip"
)
XLSX_S6 = "1-s2.0-S0092867421003810-mmc6.xlsx"

N_TOPMED = 5489
# GTEx v8 838 donors minus 39 global outliers (de Goede STAR Methods).
# Per-gene N (individuals with >=5 tissues and the gene tested) is not in
# Table S6; a shared denominator shifts -log10 by a constant and leaves
# Spearman / Pearson unchanged.
N_DEGOEDE = 838 - 39
Z_CUT = -2
MIN_N = 2
COL_UNDER = OI["purple"]
COL_OVER = OI["orange"]
TM_COUNT = OUT_PIPE / "underlier_gene_counts_z_-2_zonly_allgenes_pc2005_kept.tsv"


def nv(x: str) -> str:
    return str(x).split(".", 1)[0]


def sunflower(k: int, radius: float) -> tuple[np.ndarray, np.ndarray]:
    if k == 1:
        return np.zeros(1), np.zeros(1)
    i = np.arange(k, dtype=float)
    r = radius * np.sqrt((i + 0.5) / k)
    theta = i * np.pi * (3.0 - np.sqrt(5.0))
    return r * np.cos(theta), r * np.sin(theta)


def spread_lattice(d: pd.DataFrame) -> pd.DataFrame:
    d = d.copy()
    d["x_jit"] = d.x
    d["y_jit"] = d.y
    for (_, _), idx in d.groupby(["n_topmed", "n_degoede"]).groups.items():
        k = len(idx)
        if k == 1:
            continue
        n_tm = int(d.loc[idx[0], "n_topmed"])
        n_dg = int(d.loc[idx[0], "n_degoede"])
        spacing = min(
            np.log10((n_tm + 1) / n_tm),
            np.log10((n_dg + 1) / n_dg),
        )
        radius = min(0.45 * spacing, 0.014 * np.sqrt(k))
        dx, dy = sunflower(k, radius)
        d.loc[idx, "x_jit"] = d.loc[idx, "x"].to_numpy() + dx
        d.loc[idx, "y_jit"] = d.loc[idx, "y"].to_numpy() + dy
    return d


S6_CACHE = OUT_PIPE / "degoede_tableS6_lincRNA_outlier_gene_counts.tsv"


def load_table_s6() -> pd.DataFrame:
    cached = S6_CACHE.exists() and "n_degoede_z_lt2" in (
        S6_CACHE.read_text().splitlines()[0] if S6_CACHE.exists() else ""
    )
    if cached:
        return pd.read_csv(S6_CACHE, sep="\t")
    with zipfile.ZipFile(ZIP_S6) as z, z.open(XLSX_S6) as fh:
        wb = openpyxl.load_workbook(fh, read_only=True, data_only=True)
        ws = wb["Table_S6_lincRNA_outliers"]
        rows = ws.iter_rows(values_only=True)
        header = [str(c) for c in next(rows)]
        recs = [dict(zip(header, r)) for r in rows]
        wb.close()
    inds: dict[str, set[str]] = defaultdict(set)
    under: dict[str, set[str]] = defaultdict(set)
    over: dict[str, set[str]] = defaultdict(set)
    z_lt2: dict[str, set[str]] = defaultdict(set)
    z_lt3: dict[str, set[str]] = defaultdict(set)
    wide: dict[str, set[str]] = defaultdict(set)
    for r in recs:
        g = nv(r["ensembl_gene_id"])
        ind = str(r["individual"])
        inds[g].add(ind)
        z = r["median_Z_score"]
        if z is None:
            continue
        z = float(z)
        if z < 0:
            under[g].add(ind)
        elif z > 0:
            over[g].add(ind)
        if z < -2:
            z_lt2[g].add(ind)
        if z < -3:
            z_lt3[g].add(ind)
        if r["widely_detected"] == 1:
            wide[g].add(ind)
    out = pd.DataFrame(
        {
            "gene_nv": list(inds),
            "n_degoede_any": [len(inds[g]) for g in inds],
            "n_degoede_z_lt2": [len(z_lt2[g]) for g in inds],
            "n_degoede_z_lt3": [len(z_lt3[g]) for g in inds],
            "n_degoede_under": [len(under[g]) for g in inds],
            "n_degoede_over": [len(over[g]) for g in inds],
            "n_degoede_widely": [len(wide[g]) for g in inds],
        }
    )
    OUT_PIPE.mkdir(parents=True, exist_ok=True)
    out.to_csv(S6_CACHE, sep="\t", index=False)
    return out


def panel_df(g: pd.DataFrame, n_dg: str, n_tm: str) -> pd.DataFrame:
    d = g.copy()
    d["n_degoede"] = d[n_dg].astype(int)
    d["n_topmed"] = d[n_tm].astype(int)
    d["N_topmed_gene"] = d["n_topmed_subj"].astype(int)
    d["prev_topmed"] = d.n_topmed / d.N_topmed_gene
    d["prev_degoede"] = d.n_degoede / N_DEGOEDE
    d = d[(d.n_topmed >= MIN_N) & (d.n_degoede >= MIN_N)].copy()
    d["x"] = -np.log10(d.prev_topmed)
    d["y"] = -np.log10(d.prev_degoede)
    return d


def draw_panel(ax, d: pd.DataFrame, color: str, xlab: str, ylab: str, tag: str) -> None:
    d = spread_lattice(d)
    lo = float(min(d.x_jit.min(), d.y_jit.min()))
    hi = float(max(d.x_jit.max(), d.y_jit.max()))
    pad = 0.06 * (hi - lo)
    lim = (lo - pad, hi + pad)
    sp = float(d.prev_topmed.corr(d.prev_degoede, method="spearman"))
    pe = float(np.corrcoef(d.x, d.y)[0, 1])
    ax.plot(lim, lim, ls="--", color="grey", lw=1.0, zorder=0)
    sh = d.sample(frac=1.0, random_state=1)
    ax.scatter(
        sh.x_jit, sh.y_jit, s=22, c=color, alpha=0.85,
        edgecolor="white", linewidth=0.2, zorder=2,
    )
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(
        f"n={len(d)}  |  Spearman={sp:.2f}, Pearson(-log10)={pe:.2f}",
        loc="left", fontsize=9, fontweight="bold", pad=4,
    )
    ax.text(-0.12, 1.08, tag, transform=ax.transAxes,
            fontweight="bold", fontsize=14, va="bottom")
    ax.legend(
        handles=[Line2D([0], [0], linestyle="--", color="grey", lw=1.0, label="y = x")],
        loc="upper left", frameon=False, fontsize=8, borderaxespad=0.2,
    )
    print(f"panel {tag} n={len(d)} Spearman={sp:.2f} Pearson(-log10)={pe:.2f}")


def main() -> int:
    apply_paper_style()
    overlap = OUT_PIPE / "degoede_overlap_topmed_median_z.tsv"
    if not overlap.exists():
        raise SystemExit(f"missing {overlap}")
    g = pd.read_csv(overlap, sep="\t")
    print(f"S6 ∩ TOPMed HVG lncRNA={len(g)}")

    under = panel_df(g, "n_Z_lt2", "n_z_lt2")
    over = panel_df(g, "n_Z_gt2", "n_z_gt2")
    if under.empty or over.empty:
        raise SystemExit("empty under or over panel")

    fig, axes = plt.subplots(
        1, 2, figsize=(NATURE_DOUBLE_IN, NATURE_DOUBLE_IN * 0.52)
    )
    draw_panel(
        axes[0], under, COL_UNDER,
        r"$-\log_{10}$ TOPMed prevalence ($z < -2$, all subjects)",
        r"$-\log_{10}$ de Goede prevalence (person $z < -2$)",
        "a",
    )
    draw_panel(
        axes[1], over, COL_OVER,
        r"$-\log_{10}$ TOPMed prevalence ($z > 2$, all subjects)",
        r"$-\log_{10}$ de Goede prevalence (person $z > 2$)",
        "b",
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
        N_degoede=N_DEGOEDE,
        prev_topmed_z_lt2=g.n_z_lt2 / g.n_topmed_subj,
        prev_degoede_z_lt2=g.n_Z_lt2 / N_DEGOEDE,
        prev_topmed_z_gt2=g.n_z_gt2 / g.n_topmed_subj,
        prev_degoede_z_gt2=g.n_Z_gt2 / N_DEGOEDE,
        in_panel_a=(g.n_z_lt2 >= MIN_N) & (g.n_Z_lt2 >= MIN_N),
        in_panel_b=(g.n_z_gt2 >= MIN_N) & (g.n_Z_gt2 >= MIN_N),
    ).to_csv(src_tsv, sep="\t", index=False)
    print("wrote", DRAFT / f"{STEM}.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
