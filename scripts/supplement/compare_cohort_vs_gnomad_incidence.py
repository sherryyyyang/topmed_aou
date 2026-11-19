import sys
sys.path.insert(0, "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts/scripts/helpers")
from _paths import ROOT, FIG1A, FIG1B, FIG1C, FIG4D, SUP_GN, DATA_CORE, DATA_UNDERLIER, AOU_RNA
#!/usr/bin/env python3
"""Compare cohort-derived recessive incidence (ClinVar P/LP, Hardy-Weinberg)
against an independent gnomAD-based AR incidence table.

Cohort incidence per gene (pan-cohort HWE):
    q_gene      = sum of cohort allele frequencies of ClinVar P/LP SNV/indel
                  alleles assigned to the gene
    incidence   = q_gene ** 2          (prob. of biallelic pathogenic genotype)
    carrier_fr  = 1 - (1 - q_gene)**2  (prob. of >=1 pathogenic allele)

Cohort allele frequencies come from the ClinVar n cohort intersect VCFs
(AC/AN), so only P/LP alleles actually observed in the cohort contribute
(cohort therefore under-estimates vs a large reference like gnomAD).

Outputs:
  cohort_vs_gnomad_incidence_per_gene.tsv
  cohort_vs_gnomad_incidence_scatter.png
"""
from __future__ import annotations

import glob
import math
import os
import subprocess
import sys

import numpy as np
import pandas as pd

CLINVAR = "/net/wonderland/home/sheyang/topmed11/annotation/variant_annotations/clinvar.vcf.gz"
INTERSECT_GLOB = "/net/wonderland/home/sheyang/topmed11/topmed_outlier_analysis/quick_clinvar_overlap_slurm/vcf/clinvar_x_cohort.chr*.vcf.gz"
XLSX = "Copy of IE_AR_USA_gnomAD_20250205.xlsx"
OUT_TSV = "cohort_vs_gnomad_incidence_per_gene.tsv"
OUT_PNG = "cohort_vs_gnomad_incidence_scatter.png"
OUT_PNG_REFINED = "cohort_vs_gnomad_incidence_scatter_refined.png"

ACGT = set("ACGT")


def norm_chrom(c: str) -> str:
    return c[3:] if c.startswith("chr") else c


def is_small(ref: str, alt: str) -> bool:
    # SNV or short indel; both alleles plain ACGT and <= 50 bp.
    return (
        1 <= len(ref) <= 50
        and 1 <= len(alt) <= 50
        and set(ref) <= ACGT
        and set(alt) <= ACGT
    )


def main() -> int:
    print("[1/5] Loading gnomAD AR incidence table:", XLSX, flush=True)
    xl = pd.read_excel(XLSX)
    xl = xl[["gene", "incidence_rate_low", "incidence_rate_base", "incidence_rate_high"]].copy()
    xl["gene_uc"] = xl["gene"].str.upper()
    excel_genes = set(xl["gene_uc"])
    print(f"    {len(excel_genes)} AR genes", flush=True)

    print("[2/5] Building ClinVar P/LP key->gene map (genome-wide) ...", flush=True)
    # Pathogenic / Likely_pathogenic, excluding Conflicting; keep all genes,
    # but only retain variants touching an Excel gene to bound memory.
    expr = 'CLNSIG ~ "[Pp]athogenic" && CLNSIG !~ "onflicting"'
    fmt = r"%CHROM\t%POS\t%REF\t%ALT\t%INFO/GENEINFO\n"
    p = subprocess.Popen(
        ["bcftools", "query", "-i", expr, "-f", fmt, CLINVAR],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    )
    key2genes: dict[str, set[str]] = {}
    n_plp = 0
    assert p.stdout is not None
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        chrom, pos, ref, alt, geneinfo = f[0], f[1], f[2].upper(), f[3].upper(), f[4]
        if not is_small(ref, alt):
            continue
        n_plp += 1
        if geneinfo in (".", ""):
            continue
        genes = set()
        for tok in geneinfo.split("|"):
            sym = tok.split(":")[0].upper()
            if sym in excel_genes:
                genes.add(sym)
        if not genes:
            continue
        key = f"{norm_chrom(chrom)}:{pos}:{ref}:{alt}"
        key2genes.setdefault(key, set()).update(genes)
    p.wait()
    print(f"    P/LP SNV/indel records scanned: {n_plp:,}", flush=True)
    print(f"    P/LP alleles assigned to Excel genes: {len(key2genes):,}", flush=True)

    print("[3/5] Reading cohort AC/AN from intersect VCFs ...", flush=True)
    key_af: dict[str, float] = {}
    files = sorted(glob.glob(INTERSECT_GLOB))
    for fp in files:
        q = subprocess.Popen(
            ["bcftools", "query", "-f", r"%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n", fp],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
        )
        assert q.stdout is not None
        for line in q.stdout:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            chrom, pos, ref = norm_chrom(f[0]), f[1], f[2].upper()
            alts = f[3].upper().split(",")
            acs = f[4].split(",")
            an = f[5]
            try:
                an_i = int(an)
            except ValueError:
                continue
            if an_i <= 0:
                continue
            for alt, ac in zip(alts, acs):
                key = f"{chrom}:{pos}:{ref}:{alt}"
                if key not in key2genes:
                    continue
                try:
                    af = int(ac) / an_i
                except ValueError:
                    continue
                # keep max if duplicated across rows
                if af > key_af.get(key, -1.0):
                    key_af[key] = af
        q.wait()
    print(f"    P/LP alleles observed in cohort: {len(key_af):,}", flush=True)

    print("[4/5] Aggregating per-gene cohort incidence (HWE q^2) ...", flush=True)
    # n P/LP alleles assigned per gene (in Excel set), and observed sum(AF).
    gene_nplp: dict[str, int] = {}
    gene_nobs: dict[str, int] = {}
    gene_q: dict[str, float] = {}
    for key, genes in key2genes.items():
        af = key_af.get(key, 0.0)
        for g in genes:
            gene_nplp[g] = gene_nplp.get(g, 0) + 1
            if af > 0:
                gene_nobs[g] = gene_nobs.get(g, 0) + 1
                gene_q[g] = gene_q.get(g, 0.0) + af

    rows = []
    for g in sorted(set(gene_nplp)):
        q = gene_q.get(g, 0.0)
        rows.append(
            dict(
                gene_uc=g,
                n_plp_alleles_clinvar=gene_nplp[g],
                n_plp_alleles_in_cohort=gene_nobs.get(g, 0),
                cohort_q=q,
                cohort_incidence=q * q,
                cohort_carrier_freq=1 - (1 - q) ** 2,
            )
        )
    cohort = pd.DataFrame(rows)

    print("[5/5] Merging with gnomAD table + stats + plot ...", flush=True)
    m = xl.merge(cohort, on="gene_uc", how="left")
    for c in ["n_plp_alleles_clinvar", "n_plp_alleles_in_cohort"]:
        m[c] = m[c].fillna(0).astype(int)
    for c in ["cohort_q", "cohort_incidence", "cohort_carrier_freq"]:
        m[c] = m[c].fillna(0.0)
    m["ratio_cohort_over_gnomad"] = np.where(
        m["incidence_rate_base"] > 0, m["cohort_incidence"] / m["incidence_rate_base"], np.nan
    )
    m = m.sort_values("incidence_rate_base", ascending=False)
    m.to_csv(OUT_TSV, sep="\t", index=False)
    print("    wrote", OUT_TSV, flush=True)

    obs = m[(m["cohort_incidence"] > 0) & (m["incidence_rate_base"] > 0)].copy()
    n_excel = len(m)
    n_with_plp = int((m["n_plp_alleles_clinvar"] > 0).sum())
    n_obs = len(obs)
    print(f"    Excel AR genes: {n_excel}", flush=True)
    print(f"    ... with >=1 ClinVar P/LP allele: {n_with_plp}", flush=True)
    print(f"    ... with >=1 P/LP allele observed in cohort (nonzero incidence): {n_obs}", flush=True)

    def analyze_and_plot(df, out_png, title):
        if len(df) < 3:
            print(f"    [{out_png}] too few genes ({len(df)}) to plot.", flush=True)
            return
        lx = np.log10(df["incidence_rate_base"].values)
        ly = np.log10(df["cohort_incidence"].values)
        sp = pd.Series(lx).corr(pd.Series(ly), method="spearman")
        pe = pd.Series(lx).corr(pd.Series(ly), method="pearson")
        ratio = df["ratio_cohort_over_gnomad"].replace([np.inf, -np.inf], np.nan).dropna()
        log10_ratio = np.log10(ratio[ratio > 0])
        print(f"    [{out_png}] n={len(df)}  Spearman={sp:.3f}  Pearson(log)={pe:.3f}", flush=True)
        print(
            f"        cohort/gnomAD ratio: median={ratio.median():.3g}  "
            f"IQR=[{ratio.quantile(.25):.3g}, {ratio.quantile(.75):.3g}]  "
            f"log10-SD={log10_ratio.std():.3f}",
            flush=True,
        )

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(7.2, 7.0))
        ax.scatter(
            df["incidence_rate_base"], df["cohort_incidence"],
            s=18, alpha=0.55, color="#2166AC", edgecolor="none",
        )
        lo = min(df["incidence_rate_base"].min(), df["cohort_incidence"].min())
        hi = max(df["incidence_rate_base"].max(), df["cohort_incidence"].max())
        ax.plot([lo, hi], [lo, hi], ls="--", color="grey", lw=0.8, label="y = x")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("gnomAD AR incidence (base estimate)")
        ax.set_ylabel("Cohort recessive incidence  (HWE q$^2$, ClinVar P/LP)")
        ax.set_title(
            f"{title}\n"
            f"n={len(df)} genes  |  Spearman={sp:.2f}, Pearson(log)={pe:.2f}",
            fontsize=10,
        )
        ax.legend(loc="upper left", fontsize=8)
        ax.grid(True, which="both", ls=":", lw=0.3, alpha=0.5)
        fig.tight_layout()
        fig.savefig(out_png, dpi=200)
        plt.close(fig)
        print("    wrote", out_png, flush=True)

    # Full set.
    analyze_and_plot(obs, OUT_PNG, "Cohort vs gnomAD recessive incidence per gene")

    # Refined: drop genes resting on a single observed P/LP allele (floor-pinned).
    refined = obs[obs["n_plp_alleles_in_cohort"] >= 2].copy()
    n_dropped = n_obs - len(refined)
    print(f"    Refined: dropped {n_dropped} single-observed-allele (floor-pinned) genes; "
          f"{len(refined)} remain.", flush=True)
    analyze_and_plot(
        refined, OUT_PNG_REFINED,
        "Cohort vs gnomAD recessive incidence per gene (>=2 observed P/LP alleles)",
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
