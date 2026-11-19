#!/usr/bin/env python
"""Per-gene pLoF carrier prevalence from the gnomAD v4.1 exome sites VCFs.

Same quantity as fetch_gnomad_plof_af.py but read straight out of the public
gnomAD VCFs with remote tabix instead of the browser API, which is far faster:
one tabix call per gene, streaming only the blocks that overlap it. The .tbi
indexes are downloaded once and passed back to htslib with the file##idx##index
syntax, otherwise every call would re-fetch the index over the network.

Both sides are GRCh38: gnomAD v4.1 is GRCh38 and GENCODE v30 is Ensembl 96 /
GRCh38, and both use chr-prefixed contig names.

Variant selection, matching the cohort PTV definition as closely as the VCF allows:
  FILTER = PASS, LOFTEE LoF = HC with an empty LoF_filter, attributed to the gene
  by the Ensembl gene id in the VEP annotation, allele frequency AF <= 0.5%.

  q            = 1 - prod(1 - AF_i)   cumulative pLoF allele frequency
  carrier_prev = 1 - (1 - q)^2        Hardy-Weinberg, ~2q for rare alleles

Frequencies are exomes-only (AN ~ 1.2M alleles). The API script uses the joint
exome+genome callset, so its numbers run slightly higher; do not mix the two
within one estimate.
"""
import argparse
import gzip
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "shared_rare_outlier_genes_3dot.tsv")
GFF = os.path.join(HERE, "..", "annotation", "gencode.v30.annotation.gff3.gz")
VCF = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/"
       "exomes/gnomad.exomes.v4.1.sites.chr{chrom}.vcf.bgz")
OUT = os.path.join(HERE, "gnomad_plof_prevalence_tabix.tsv")
AF_CAP = 0.005     # same ceiling as the cohort PTV call set
AF_RARE = 1e-4     # "rare only" sensitivity estimate
PAD = 1000         # bp added either side of the gene body
N_COHORT = 5489    # TOPMed subjects, for an expected-carrier count


def gene_coords(gene_ids: set[str]) -> pd.DataFrame:
    """Gene body coordinates for unversioned Ensembl ids, from GENCODE v30."""
    rows = []
    opener = gzip.open if GFF.endswith(".gz") else open
    with opener(GFF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "gene":
                continue
            gid = ""
            for kv in f[8].split(";"):
                if kv.startswith("gene_id="):
                    gid = kv[8:].split(".")[0]
                    break
            if gid in gene_ids:
                rows.append((gid, f[0], int(f[3]), int(f[4])))
    return pd.DataFrame(rows, columns=["gene_nv", "chrom", "start", "end"])


def csq_index(chrom: str) -> dict[str, int]:
    hdr = subprocess.run(["tabix", "-H", VCF.format(chrom=chrom)],
                         capture_output=True, text=True, check=True).stdout
    for line in hdr.splitlines():
        if "ID=vep" in line and "Format:" in line:
            fmt = line.split("Format:")[1].strip().rstrip('">').strip()
            return {name: i for i, name in enumerate(fmt.split("|"))}
    raise RuntimeError(f"no vep Format line in the chr{chrom} header")


def parse_chrom(chrom: str, genes: pd.DataFrame, tmpdir: str) -> list[dict]:
    idx = csq_index(chrom)
    i_gene, i_lof, i_filt = idx["Gene"], idx["LoF"], idx["LoF_filter"]
    i_canon, i_sym = idx["CANONICAL"], idx["SYMBOL"]
    wanted = set(genes.gene_nv)
    bed = os.path.join(tmpdir, f"regions.chr{chrom}.bed")
    with open(bed, "w") as fh:
        for _, r in genes.iterrows():
            fh.write(f"{r.chrom}\t{max(r.start - PAD - 1, 0)}\t{r.end + PAD}\n")

    acc: dict[str, dict] = {g: dict(afs_canon=[], afs_any=[], n_nonpass=0,
                                   n_over_cap=0, symbol="") for g in wanted}
    proc = subprocess.Popen(["tabix", "-R", bed, VCF.format(chrom=chrom)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    n_rec = 0
    for line in proc.stdout:
        n_rec += 1
        f = line.rstrip("\n").split("\t", 8)
        if len(f) < 8 or f[4] == "*" or "," in f[4]:
            continue
        info = f[7]
        # LOFTEE values sit in positional CSQ subfields, so there is no "LoF"
        # substring to pre-filter on; only the vep field itself can be required
        if "vep=" not in info:
            continue
        passing = f[6] in ("PASS", ".")
        ac = an = af = None
        vep = ""
        for kv in info.split(";"):
            if kv.startswith("AF="):
                af = float(kv[3:].split(",")[0])
            elif kv.startswith("AC="):
                ac = float(kv[3:].split(",")[0])
            elif kv.startswith("AN="):
                an = float(kv[3:])
            elif kv.startswith("vep="):
                vep = kv[4:]
        if af is None:
            af = (ac / an) if (ac and an) else 0.0
        if not vep or af <= 0:
            continue
        # a variant can be pLoF for more than one gene; attribute per annotation
        hits: dict[str, tuple[bool, str]] = {}
        for tx in vep.split(","):
            c = tx.split("|")
            if len(c) <= i_filt:
                continue
            g = c[i_gene].split(".")[0]
            if g not in wanted or c[i_lof] != "HC" or c[i_filt] not in ("", "."):
                continue
            canon = c[i_canon] == "YES"
            prev = hits.get(g)
            hits[g] = (canon or bool(prev and prev[0]), c[i_sym] or (prev[1] if prev else ""))
        for g, (canon, sym) in hits.items():
            a = acc[g]
            if sym and not a["symbol"]:
                a["symbol"] = sym
            if not passing:
                a["n_nonpass"] += 1
                continue
            if af > AF_CAP:
                a["n_over_cap"] += 1
                continue
            a["afs_any"].append(af)
            if canon:
                a["afs_canon"].append(af)
    proc.stdout.close()
    err = proc.stderr.read()
    if proc.wait() != 0:
        raise RuntimeError(f"tabix failed on chr{chrom}: {err[:300]}")
    print(f"  chr{chrom}: {len(genes)} genes, {n_rec:,} VCF records scanned", flush=True)

    cum = lambda afs: 1.0 - math.prod(1.0 - a for a in afs)
    rows = []
    for _, r in genes.iterrows():
        a = acc[r.gene_nv]
        q, q_any = cum(a["afs_canon"]), cum(a["afs_any"])
        q_rare = cum([x for x in a["afs_canon"] if x < AF_RARE])
        prev = 1.0 - (1.0 - q) ** 2
        rows.append(dict(
            gene_nv=r.gene_nv, symbol_gnomad=a["symbol"], chrom=r.chrom,
            start=r.start, end=r.end,
            n_plof=len(a["afs_canon"]), n_plof_any_tx=len(a["afs_any"]),
            n_plof_rare=sum(1 for x in a["afs_canon"] if x < AF_RARE),
            n_dropped_nonpass=a["n_nonpass"], n_dropped_af_over_cap=a["n_over_cap"],
            max_variant_af=max(a["afs_canon"], default=0.0),
            gnomad_plof_af=q, gnomad_plof_carrier_prev=prev,
            gnomad_plof_carrier_prev_any_tx=1.0 - (1.0 - q_any) ** 2,
            gnomad_plof_carrier_prev_rare=1.0 - (1.0 - q_rare) ** 2,
            expected_carriers_in_cohort=prev * N_COHORT))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=8, help="chromosomes in flight")
    ap.add_argument("--chrom", default="", help="restrict to one chromosome (testing)")
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()

    genes = pd.read_csv(SRC, sep="\t")[["gene_nv", "symbol"]].drop_duplicates()
    coords = gene_coords(set(genes.gene_nv))
    missing = sorted(set(genes.gene_nv) - set(coords.gene_nv))
    print(f"{len(genes)} genes, coordinates found for {len(coords)}"
          + (f", missing {len(missing)}: {missing[:5]}" if missing else ""), flush=True)
    if args.chrom:
        coords = coords[coords.chrom == args.chrom]
        print(f"restricted to {args.chrom}: {len(coords)} genes", flush=True)

    groups = [(c.replace("chr", ""), sub) for c, sub in coords.groupby("chrom")]
    groups.sort(key=lambda kv: -len(kv[1]))
    tmpdir = os.path.join(HERE, ".gnomad_tabix_tmp")
    os.makedirs(tmpdir, exist_ok=True)
    rows: list[dict] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(parse_chrom, c, sub, tmpdir) for c, sub in groups]
        for fut in futures:
            rows += fut.result()

    d = pd.DataFrame(rows).merge(genes, on="gene_nv", how="right")
    d.to_csv(args.out, sep="\t", index=False)
    ok = d[d.n_plof.notna()]
    med = ok.gnomad_plof_carrier_prev.median()
    print(f"\nwrote {args.out}: {len(d)} genes")
    print(f"  pLoF variants per gene (canonical tx): median {ok.n_plof.median():.0f}, "
          f"max {ok.n_plof.max():.0f}; genes with none: {(ok.n_plof == 0).sum()}")
    print(f"  carrier prevalence: median {med:.3g}"
          + (f" (1 in {round(1 / med):,})" if med > 0 else ""))
    print(f"  expected carriers in {N_COHORT:,} subjects: median "
          f"{ok.expected_carriers_in_cohort.median():.1f}")
    print(f"  dropped: {int(ok.n_dropped_nonpass.sum())} non-PASS, "
          f"{int(ok.n_dropped_af_over_cap.sum())} above AF {AF_CAP:.1%}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
