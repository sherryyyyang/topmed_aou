#!/usr/bin/env python
"""Per-gene pLoF carrier prevalence in gnomAD v4.1, for the 624 shared rare-outlier genes.

The cohort PTV dot in shared_rare_outlier_genes_3dot.png counts TOPMed subjects
carrying a LOFTEE high-confidence PTV with AF <= 0.5%. This script builds the
external counterpart from gnomAD v4.1: for each gene, sum the allele frequencies
of LOFTEE high-confidence pLoF variants and convert to a carrier prevalence.

gnomAD publishes no bulk pLoF table for v4, so the numbers come from the browser
GraphQL API one gene at a time (batched). Where gnomAD has manually curated a
pLoF variant, the verdict is used to drop variants judged not to be real LoF.

  q            = 1 - prod(1 - AF_i) over kept pLoF variants  (cumulative pLoF AF)
  carrier_prev = 1 - (1 - q)^2                               (HWE, ~2q when small)

Output is appended per batch, so the script can be re-run to resume.
"""
import argparse
import json
import math
import os
import time
import urllib.error
import urllib.request

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
API = "https://gnomad.broadinstitute.org/api"
SRC = os.path.join(HERE, "shared_rare_outlier_genes_3dot.tsv")
OUT = os.path.join(HERE, "gnomad_plof_af_per_gene.tsv")
AF_CAP = 0.005          # same AF ceiling as the cohort PTV call set
NOT_LOF = {"not_lof", "likely_not_lof"}

FIELDS = "lof lof_filter flags exome{ac an} genome{ac an} joint{ac an} lof_curation{verdict}"


def query(gene_ids: list[str], tries: int = 5) -> dict:
    parts = [f'g{i}: gene(gene_id:"{g}", reference_genome:GRCh38)'
             f'{{ gene_id symbol variants(dataset:gnomad_r4){{ {FIELDS} }} }}'
             for i, g in enumerate(gene_ids)]
    payload = json.dumps({"query": "{" + " ".join(parts) + "}"}).encode()
    delay = 5.0
    for attempt in range(tries):
        try:
            req = urllib.request.Request(
                API, data=payload, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=300) as fh:
                body = json.loads(fh.read())
            if body.get("errors") and not body.get("data"):
                raise RuntimeError(json.dumps(body["errors"])[:300])
            return body
        except urllib.error.HTTPError as exc:
            msg = exc.read()[:300].decode(errors="replace")
            # the API caps query cost at 25 and one gene's variant list costs 11,
            # so more than 2 genes per request is rejected outright: not retryable
            if "too expensive" in msg:
                raise RuntimeError(f"reduce --batch: {msg}") from None
            if attempt == tries - 1:
                raise RuntimeError(f"HTTP {exc.code}: {msg}") from None
            print(f"    retry {attempt + 1} after {delay:.0f}s (HTTP {exc.code})", flush=True)
            time.sleep(delay)
            delay *= 2
        except (urllib.error.URLError, TimeoutError,
                RuntimeError, json.JSONDecodeError) as exc:
            if attempt == tries - 1:
                raise
            print(f"    retry {attempt + 1} after {delay:.0f}s ({type(exc).__name__}: "
                  f"{str(exc)[:120]})", flush=True)
            time.sleep(delay)
            delay *= 2
    raise AssertionError("unreachable")


def af_of(v: dict) -> tuple[float, int, str]:
    """Allele frequency, preferring the joint v4.1 callset over exomes alone."""
    for key in ("joint", "exome", "genome"):
        f = v.get(key)
        if f and f.get("an"):
            return f["ac"] / f["an"], f["an"], key
    return 0.0, 0, "none"


def summarise(gene_id: str, node: dict) -> dict:
    variants = node["variants"]
    hc = [v for v in variants if v.get("lof") == "HC" and not v.get("lof_filter")]
    kept, dropped_af, dropped_curated, curated = [], 0, 0, 0
    for v in hc:
        af, an, src = af_of(v)
        if af <= 0:
            continue
        verdict = (v.get("lof_curation") or {}).get("verdict")
        if verdict:
            curated += 1
            if verdict in NOT_LOF:
                dropped_curated += 1
                continue
        if af > AF_CAP:
            dropped_af += 1
            continue
        kept.append((af, an, src))
    cum = lambda afs: 1.0 - math.prod(1.0 - a for a in afs)
    afs = [af for af, _, _ in kept]
    q = cum(afs)
    # a single pLoF variant sitting just under the 0.5% cap dominates the sum, so
    # also report the estimate restricted to genuinely rare alleles (AF < 0.01%),
    # which is the regime the cohort PTV counts actually sample
    rare = [a for a in afs if a < 1e-4]
    q_rare = cum(rare)
    ans = sorted(an for _, an, _ in kept)
    return dict(
        gene_nv=gene_id, symbol=node["symbol"],
        n_hc_plof=len(hc), n_used=len(kept), n_curated=curated,
        n_dropped_curated_not_lof=dropped_curated, n_dropped_af_over_cap=dropped_af,
        gnomad_plof_af=q, gnomad_plof_carrier_prev=1.0 - (1.0 - q) ** 2,
        max_variant_af=max(afs) if afs else 0.0,
        n_used_rare=len(rare), gnomad_plof_af_rare=q_rare,
        gnomad_plof_carrier_prev_rare=1.0 - (1.0 - q_rare) ** 2,
        median_an=ans[len(ans) // 2] if ans else 0,
        freq_source=kept[0][2] if kept else "none",
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    # 2 is the ceiling: the API rejects queries costing more than 25 and each
    # gene's variant list costs 11
    ap.add_argument("--batch", type=int, default=2, help="genes per API request")
    ap.add_argument("--sleep", type=float, default=2.5, help="seconds between requests")
    ap.add_argument("--limit", type=int, default=0, help="stop after N genes (testing)")
    ap.add_argument("--slice", default="1/1",
                    help="k/n: take every nth gene starting at k, for parallel workers")
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()

    genes = pd.read_csv(SRC, sep="\t")[["gene_nv", "symbol"]].drop_duplicates()
    k, n = (int(x) for x in args.slice.split("/"))
    if n > 1:
        genes = genes.iloc[k - 1::n]
        print(f"slice {k}/{n}: {len(genes)} genes", flush=True)
    done: set[str] = set()
    if os.path.exists(args.out):
        done = set(pd.read_csv(args.out, sep="\t").gene_nv)
        print(f"resuming: {len(done)} genes already in {os.path.basename(args.out)}", flush=True)
    todo = [g for g in genes.gene_nv if g not in done]
    if args.limit:
        todo = todo[:args.limit]
    print(f"{len(todo)} genes to fetch, {args.batch} per request", flush=True)

    header_needed = not os.path.exists(args.out) or os.path.getsize(args.out) == 0
    for i in range(0, len(todo), args.batch):
        chunk = todo[i:i + args.batch]
        body = query(chunk)
        rows, missing = [], []
        for j, gid in enumerate(chunk):
            node = (body.get("data") or {}).get(f"g{j}")
            if not node:
                missing.append(gid)
                continue
            rows.append(summarise(gid, node))
        if missing:
            # genes gnomAD does not carry under this id; recorded as absent so the
            # resume logic does not retry them forever
            nan = float("nan")
            rows += [dict(gene_nv=g, symbol="", n_hc_plof=0, n_used=0, n_curated=0,
                          n_dropped_curated_not_lof=0, n_dropped_af_over_cap=0,
                          gnomad_plof_af=nan, gnomad_plof_carrier_prev=nan,
                          max_variant_af=nan, n_used_rare=0, gnomad_plof_af_rare=nan,
                          gnomad_plof_carrier_prev_rare=nan,
                          median_an=0, freq_source="absent") for g in missing]
        df = pd.DataFrame(rows)
        df.to_csv(args.out, sep="\t", index=False, mode="a", header=header_needed)
        header_needed = False
        print(f"  [{i + len(chunk)}/{len(todo)}] " +
              ", ".join(f"{r['symbol'] or r['gene_nv']}:{r['n_used']}" for r in rows),
              flush=True)
        if i + args.batch < len(todo):
            time.sleep(args.sleep)

    d = pd.read_csv(args.out, sep="\t")
    print(f"\n{len(d)} genes in {args.out}")
    ok = d[d.freq_source != "absent"]
    print(f"  genes with gnomAD data: {len(ok)}  |  absent: {len(d) - len(ok)}")
    print(f"  pLoF variants used per gene: median {ok.n_used.median():.0f}, "
          f"max {ok.n_used.max()}")
    print(f"  carrier prevalence: median {ok.gnomad_plof_carrier_prev.median():.3g}")
    print(f"  variants dropped by curation verdict: {int(ok.n_dropped_curated_not_lof.sum())}"
          f"  |  by AF > {AF_CAP:.1%}: {int(ok.n_dropped_af_over_cap.sum())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
