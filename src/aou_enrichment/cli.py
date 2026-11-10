"""Command line entry point for participant-level phenotype enrichment."""
import argparse
from pathlib import Path
import pandas as pd
from .analysis import run_enrichment


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--eligible", required=True)
    parser.add_argument("--underliers", required=True)
    parser.add_argument("--phenotypes", required=True)
    parser.add_argument("--mapping", default="data/gene_phenotype_map.csv")
    parser.add_argument("--selected-genes", nargs="+", required=True)
    parser.add_argument("--include-proxies", action="store_true")
    parser.add_argument("--alternative", choices=["greater", "two-sided"], default="greater")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    mapping = pd.read_csv(args.mapping, dtype=str)
    mapping = mapping[mapping.gene.isin(args.selected_genes)]
    if not args.include_proxies:
        mapping = mapping[mapping.outcome_tier.eq("primary")]
    if mapping.empty:
        raise SystemExit("No primary mappings for selected genes. Review mapping_coverage.csv or pass --include-proxies.")
    results = run_enrichment(
        pd.read_csv(args.eligible, dtype=str),
        pd.read_csv(args.underliers, dtype=str),
        pd.read_csv(args.phenotypes, dtype={"subject_id": str, "category": str, "phenoname": str}),
        mapping[["gene", "category", "phenoname"]].drop_duplicates(),
        alternative=args.alternative,
    )
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(out, index=False)
    print(f"Wrote {len(results)} prespecified gene-phenotype tests to {out}")


if __name__ == "__main__":
    main()
