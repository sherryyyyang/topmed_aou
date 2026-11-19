#!/usr/bin/env python3
"""TOPMed vs AoU prevalence colored/shaped by clinical dosage class.

Panel a: same -log10 prevalence scatter as Fig 2 LOEUF panel, but marker
shape / outline by exclusive clinical class:
  ClinGen haploinsufficient > ACMG SF > Other disease (OMIM) > Remaining
Panel b: Spearman correlation of TOPMed vs AoU prevalence within each class
(one number per class; clearer than a stacked 'is correlated?' bar).

Does not overwrite the LOEUF-colored Fig 2 panel.
"""
from __future__ import annotations

import re
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
from _paths import FIG1C, FIG2C, ROOT, SHEYANG  # noqa: E402
from paper_mpl import NATURE_DOUBLE_IN, OI, PAPER_DPI, apply_paper_style  # noqa: E402

MIN_N = 2
STEM = "fig2_topmed_vs_aou_prevalence_clinical_class"
DRAFT = Path("/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026")
SRC = DRAFT / "source_panels"
SUPP = DRAFT / "supplemental_data"
WRITING = SHEYANG / "writing_and_documentations" / "figures"

CLINGEN_CSV = ROOT / "data" / "references" / "ClinGen_HI_Sufficient_evidence.csv"
ACMG_CSV = SHEYANG / "annotation" / "variant_annotations" / "acmg_sf_v3.2.csv"
OMIM_TXT = SHEYANG / "annotation" / "genemap2.txt"

CLASS_ORDER = [
    "Remaining genes",
    "Other disease genes",
    "ACMG secondary findings",
    "ClinGen haploinsufficient",
]
CLASS_COLORS = {
    "Remaining genes": "#BDBDBD",
    "Other disease genes": OI["sky"],
    "ACMG secondary findings": OI["orange"],
    "ClinGen haploinsufficient": OI["vermillion"],
}
CLASS_MARKERS = {
    "Remaining genes": "o",
    "Other disease genes": "s",
    "ACMG secondary findings": "D",
    "ClinGen haploinsufficient": "^",
}


def parse_clingen_symbol(x: str) -> str:
    x = re.sub(r"HGNC:\d+$", "", str(x))
    m = re.match(r"^([A-Za-z0-9\-@.]+)", x)
    return m.group(1).upper() if m else x.upper()


def load_clingen() -> set[str]:
    cg = pd.read_csv(CLINGEN_CSV)
    cg = cg[cg.Type.astype(str).str.lower().eq("gene")].copy()
    cg = cg[cg["HI Score"].astype(str).str.contains("Sufficient", case=False, na=False)]
    return set(cg["Gene Symbol /Region Name"].map(parse_clingen_symbol))


def load_acmg() -> set[str]:
    acmg = pd.read_csv(ACMG_CSV)
    return set(acmg["Gene"].astype(str).str.upper())


def load_omim() -> tuple[set[str], set[str]]:
    syms: set[str] = set()
    ens: set[str] = set()
    with open(OMIM_TXT) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 13:
                continue
            approved = parts[8].strip().upper()
            ensembl = parts[10].strip().split(".", 1)[0]
            phenotypes = parts[12].strip()
            if not phenotypes:
                continue
            if approved:
                syms.add(approved)
            if ensembl:
                ens.add(ensembl)
    return syms, ens


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
    for (_, _), idx in d.groupby(["n_topmed", "n_new"]).groups.items():
        k = len(idx)
        if k == 1:
            continue
        n_tm = int(d.loc[idx[0], "n_topmed"])
        n_aou = int(d.loc[idx[0], "n_new"])
        spacing = min(np.log10((n_tm + 1) / n_tm), np.log10((n_aou + 1) / n_aou))
        # Allow mild neighbor-bin bleed so coincident count pairs separate
        radius = min(1.25 * spacing, 0.06 * np.sqrt(k) + 0.01)
        dx, dy = sunflower(k, radius)
        d.loc[idx, "x_jit"] = d.loc[idx, "x"].to_numpy() + dx
        d.loc[idx, "y_jit"] = d.loc[idx, "y"].to_numpy() + dy
    return d


def assign_class(row: pd.Series, clingen: set[str], acmg: set[str],
                 omim_sym: set[str], omim_ens: set[str]) -> str:
    sym = row.symbol_u
    if sym in clingen:
        return "ClinGen haploinsufficient"
    if sym in acmg:
        return "ACMG secondary findings"
    if sym in omim_sym or row.gene_nv in omim_ens:
        return "Other disease genes"
    return "Remaining genes"


def main() -> int:
    apply_paper_style()
    clingen = load_clingen()
    acmg = load_acmg()
    omim_sym, omim_ens = load_omim()

    prev = pd.read_csv(FIG1C / "topmed_vs_new_z_-3_prevalence_per_gene.tsv", sep="\t")
    d = prev[(prev.n_topmed >= MIN_N) & (prev.n_new >= MIN_N)].copy()
    d["symbol_u"] = d.symbol.astype(str).str.upper()
    d["clin_class"] = d.apply(
        assign_class, axis=1,
        clingen=clingen, acmg=acmg, omim_sym=omim_sym, omim_ens=omim_ens,
    )
    d["clin_class"] = pd.Categorical(d.clin_class, categories=CLASS_ORDER, ordered=True)
    print(d.clin_class.value_counts().to_string())

    d["x"] = -np.log10(d.prev_topmed)
    d["y"] = -np.log10(d.prev_new)
    d = spread_lattice(d)

    corr_rows = []
    for c in CLASS_ORDER:
        sub = d[d.clin_class == c]
        if len(sub) < 3:
            sp = pe = np.nan
        else:
            sp = float(sub.prev_topmed.corr(sub.prev_new, method="spearman"))
            pe = float(np.corrcoef(
                -np.log10(sub.prev_topmed), -np.log10(sub.prev_new)
            )[0, 1])
        corr_rows.append({"clin_class": c, "n": len(sub), "spearman": sp, "pearson_log10": pe})
        print(f"{c}: n={len(sub)} Spearman={sp:.3f}")
    corr = pd.DataFrame(corr_rows)

    lo = float(min(d.x_jit.min(), d.y_jit.min()))
    hi = float(max(d.x_jit.max(), d.y_jit.max()))
    pad = 0.06 * (hi - lo)
    lim = (lo - pad, hi + pad)

    fig, axes = plt.subplots(
        1, 2, figsize=(NATURE_DOUBLE_IN, NATURE_DOUBLE_IN * 0.52),
        gridspec_kw={"width_ratios": [1.35, 1]},
    )

    ax = axes[0]
    ax.plot(lim, lim, ls="--", color="grey", lw=1.0, zorder=0)
    # draw Remaining first, clinical classes on top
    for c in CLASS_ORDER:
        sub = d[d.clin_class == c]
        if sub.empty:
            continue
        edge = "white" if c == "Remaining genes" else "black"
        lw = 0.08 if c == "Remaining genes" else 0.22
        s = 2.5 if c == "Remaining genes" else 7
        alpha = 0.22 if c == "Remaining genes" else 0.70
        ax.scatter(
            sub.x_jit, sub.y_jit,
            s=s, marker=CLASS_MARKERS[c],
            c=CLASS_COLORS[c],
            alpha=alpha, edgecolor=edge, linewidth=lw, zorder=2 if c == "Remaining genes" else 3,
            label=f"{c} (n={len(sub)})",
        )
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$-\log_{10}$ TOPMed RNA-outlier prevalence")
    ax.set_ylabel(r"$-\log_{10}$ AoU RNA-outlier prevalence")
    ax.set_title("a", loc="left", fontweight="bold", fontsize=14, pad=8)
    ax.legend(
        loc="upper left", frameon=False, fontsize=7,
        borderaxespad=0.15, labelspacing=0.35, handletextpad=0.3,
        markerscale=1.8,
    )

    ax = axes[1]
    y = np.arange(len(CLASS_ORDER))
    vals = [corr.loc[corr.clin_class == c, "spearman"].iloc[0] for c in CLASS_ORDER]
    ns = [int(corr.loc[corr.clin_class == c, "n"].iloc[0]) for c in CLASS_ORDER]
    colors = [CLASS_COLORS[c] for c in CLASS_ORDER]
    ax.barh(y, vals, color=colors, edgecolor="black", linewidth=0.4, height=0.65)
    ax.axvline(0, color="grey", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c}\n(n={n})" for c, n in zip(CLASS_ORDER, ns)], fontsize=7)
    ax.set_xlabel("Spearman correlation\n(TOPMed vs AoU prevalence)")
    ax.set_xlim(-0.05, max(0.75, max(v for v in vals if np.isfinite(v)) + 0.08))
    ax.set_title("b", loc="left", fontweight="bold", fontsize=14, pad=8)
    for yi, v in zip(y, vals):
        if np.isfinite(v):
            ax.text(v + 0.015, yi, f"{v:.2f}", va="center", fontsize=8)

    fig.tight_layout()

    out_dirs = [FIG2C, DRAFT, SRC, WRITING]
    for folder in out_dirs:
        folder.mkdir(parents=True, exist_ok=True)
        for ext in ("png", "pdf"):
            fig.savefig(folder / f"{STEM}.{ext}", dpi=PAPER_DPI, bbox_inches="tight")
    plt.close(fig)

    SUPP.mkdir(parents=True, exist_ok=True)
    d.sort_values(["clin_class", "prev_topmed"], ascending=[True, False]).to_csv(
        SUPP / f"Source_Data_{STEM}.tsv", sep="\t", index=False
    )
    corr.to_csv(SUPP / f"Source_Data_{STEM}_spearman.tsv", sep="\t", index=False)
    print("wrote", DRAFT / f"{STEM}.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
