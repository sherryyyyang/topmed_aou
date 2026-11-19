#!/usr/bin/env python3
"""Max gnomAD AF_exomes among TOPMed-overlapping PTVs, per gene.

For each of the 624 shared rare RNA-outlier genes:
  1. take unique PTV sites from ptv_all.tsv (chrom, pos, ref, alt, gene)
  2. look those sites up in the local gnomAD joint v4.1 exomes VCF
  3. keep AF_exomes for exact chrom/pos/ref/alt matches
  4. collapse to max(AF_exomes) per gene

The gnomAD VCF has one variant per row and no gene annotation, so gene
identity comes from the TOPMed PTV table (VEP Ensembl id).
"""
import gzip
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GENES = os.path.join(HERE, "shared_rare_outlier_genes_4dot.tsv")
PTV = "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/ptv_all.tsv"
VCF = ("/net/topmed11/working/porchard/for-others/julie-and-esther/"
       "2025-11-17/sheyang/annotation/gnomad/gnomad.joint.v4.1.sites.exomes.vcf.bgz")
OUT = os.path.join(HERE, "max_af_exome_topmed_overlap.tsv")
SITES_OUT = os.path.join(HERE, "topmed_ptv_sites_624genes_with_af_exome.tsv")


def parse_af(info: str, n_alt: int):
    for field in info.split(";"):
        if field.startswith("AF_exomes="):
            vals = field.split("=", 1)[1].split(",")
            out = []
            for v in vals:
                try:
                    out.append(float(v))
                except ValueError:
                    out.append(None)
            if len(out) < n_alt:
                out.extend([None] * (n_alt - len(out)))
            return out[:n_alt]
    return [None] * n_alt


def main() -> None:
    genes = pd.read_csv(GENES, sep="\t", usecols=["gene_nv", "symbol"])
    gene_set = set(genes.gene_nv)

    ptv = pd.read_csv(PTV, sep="\t", usecols=["chrom", "pos", "ref", "alt", "symbol", "gene"])
    ptv["gene_nv"] = ptv.gene.str.replace(r"\..*$", "", regex=True)
    sites = (
        ptv[ptv.gene_nv.isin(gene_set)]
        .drop_duplicates(["chrom", "pos", "ref", "alt", "gene_nv"])
        .copy()
    )
    sites["pos"] = sites["pos"].astype(int)
    print(f"unique TOPMed PTV sites in {len(gene_set)} genes: {len(sites)} "
          f"({sites.gene_nv.nunique()} genes with >=1 site)", flush=True)

    key_to_rows = {}
    for i, r in sites.iterrows():
        key_to_rows.setdefault((str(r.chrom), int(r.pos), str(r.ref), str(r.alt)), []).append(i)

    wanted = set(key_to_rows)
    sites["AF_exomes"] = pd.NA
    n_hit = 0
    n_scan = 0
    last_chrom = None
    with gzip.open(VCF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            n_scan += 1
            chrom, pos_s, _id, ref, alt, _qual, _filt, info = line.rstrip("\n").split("\t")[:8]
            pos = int(pos_s)
            if chrom != last_chrom:
                print(f"  scanning {chrom} (scanned {n_scan:,} records, hits {n_hit})", flush=True)
                last_chrom = chrom
            alts = alt.split(",")
            afs = None
            for a, allele in enumerate(alts):
                key = (chrom, pos, ref, allele)
                if key not in wanted:
                    continue
                if afs is None:
                    afs = parse_af(info, len(alts))
                af = afs[a] if a < len(afs) else None
                for idx in key_to_rows[key]:
                    sites.at[idx, "AF_exomes"] = af
                n_hit += 1

    n_matched = int(sites["AF_exomes"].notna().sum())
    print(f"scanned {n_scan:,} gnomAD records")
    print(f"TOPMed sites matched in gnomAD: {n_matched} / {len(sites)}")
    print(f"genes with >=1 matched AF_exomes: "
          f"{sites.loc[sites.AF_exomes.notna(), 'gene_nv'].nunique()}")

    matched = sites.dropna(subset=["AF_exomes"])
    if len(matched):
        per_match = matched.groupby("gene_nv", as_index=False).agg(
            n_matched_af_exome=("AF_exomes", "size"),
            max_AF_exomes=("AF_exomes", "max"),
            median_AF_exomes=("AF_exomes", "median"),
        )
    else:
        per_match = pd.DataFrame(columns=["gene_nv", "n_matched_af_exome",
                                          "max_AF_exomes", "median_AF_exomes"])
    n_all = sites.groupby("gene_nv").size().rename("n_topmed_ptv_sites")
    per_gene = (
        genes.merge(n_all, on="gene_nv", how="left")
        .merge(per_match, on="gene_nv", how="left")
    )
    per_gene["n_topmed_ptv_sites"] = per_gene["n_topmed_ptv_sites"].fillna(0).astype(int)
    per_gene["n_matched_af_exome"] = per_gene["n_matched_af_exome"].fillna(0).astype(int)

    sites.sort_values(["chrom", "pos"]).to_csv(SITES_OUT, sep="\t", index=False)
    per_gene.to_csv(OUT, sep="\t", index=False)
    print(f"wrote {OUT}")
    print(f"wrote {SITES_OUT}")
    have = per_gene["max_AF_exomes"].notna()
    print(f"genes with max AF_exomes: {int(have.sum())} / {len(per_gene)}")
    if have.any():
        m = per_gene.loc[have, "max_AF_exomes"]
        print(f"max AF_exomes: median={m.median():.4g}  min={m.min():.4g}  max={m.max():.4g}")


if __name__ == "__main__":
    main()
