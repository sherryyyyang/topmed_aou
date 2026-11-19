#!/usr/bin/env python3
"""Mean single-cell expression per cell type, restricted to requested genes.

The reference TSV is genes x cells. Column names (except GeneSymbol) are cell
types; there is no separate donor column. Values are raw counts. Duplicate
GeneSymbol rows are averaged.
"""
from __future__ import annotations

import argparse
import gzip
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

GLOBINS = {
    "HBA1",
    "HBA2",
    "HBB",
    "HBD",
    "HBE1",
    "HBG1",
    "HBG2",
    "HBM",
    "HBQ1",
    "HBZ",
}


def read_symbols(path: Path) -> set[str]:
    df = pd.read_csv(path, sep="\t")
    col = "symbol" if "symbol" in df.columns else df.columns[0]
    return {str(x) for x in df[col].dropna() if str(x) and str(x) != "nan"}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--sc", required=True)
    p.add_argument("--symbols", required=True, help="TSV with a 'symbol' column")
    p.add_argument("--out-signature", required=True)
    p.add_argument("--out-overlap", required=True)
    p.add_argument("--out-cells", required=True)
    args = p.parse_args()

    want = read_symbols(Path(args.symbols))
    opener = gzip.open if args.sc.endswith(".gz") else open
    sums: dict[str, np.ndarray] = {}
    n_rows: dict[str, int] = {}
    n_sc_rows = 0
    sc_symbols: set[str] = set()

    with opener(args.sc, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header[0] != "GeneSymbol":
            raise SystemExit(f"unexpected first column: {header[0]!r}")
        types = np.array(header[1:], dtype=object)
        uniq, inv = np.unique(types, return_inverse=True)
        n_type = len(uniq)
        n_cells = np.bincount(inv, minlength=n_type).astype(np.int64)

        for line in fh:
            n_sc_rows += 1
            gene, _, rest = line.partition("\t")
            sc_symbols.add(gene)
            if gene not in want:
                continue
            vals = np.fromstring(rest, sep="\t", dtype=np.float64)
            if vals.size != types.size:
                raise SystemExit(
                    f"{gene}: expected {types.size} values, got {vals.size}"
                )
            acc = np.bincount(inv, weights=vals, minlength=n_type)
            if gene in sums:
                sums[gene] += acc
                n_rows[gene] += 1
            else:
                sums[gene] = acc
                n_rows[gene] = 1

    if not sums:
        raise SystemExit("no requested symbols found in the single-cell file")

    genes = sorted(sums)
    mat = np.vstack([sums[g] / n_cells / n_rows[g] for g in genes])
    sig = pd.DataFrame(mat, index=genes, columns=list(uniq))
    sig.index.name = "symbol"
    Path(args.out_signature).parent.mkdir(parents=True, exist_ok=True)
    sig.to_csv(args.out_signature, sep="\t")

    cells = pd.DataFrame({"cell_type": uniq, "n_cells": n_cells})
    cells.to_csv(args.out_cells, sep="\t", index=False)

    overlap = pd.DataFrame(
        {
            "item": [
                "n_requested_symbols",
                "n_sc_rows",
                "n_sc_unique_symbols",
                "n_overlap_symbols",
                "n_requested_missing_from_sc",
                "n_sc_duplicate_symbol_rows_in_overlap",
                "n_overlap_globins",
                "overlap_globins",
            ],
            "value": [
                len(want),
                n_sc_rows,
                len(sc_symbols),
                len(genes),
                len(want - sc_symbols),
                int(sum(n - 1 for n in n_rows.values())),
                sum(g in GLOBINS for g in genes),
                ",".join(g for g in genes if g in GLOBINS),
            ],
        }
    )
    overlap.to_csv(args.out_overlap, sep="\t", index=False)
    print(
        f"signature {sig.shape[0]} genes x {sig.shape[1]} cell types; "
        f"overlap globins={overlap.loc[overlap.item=='overlap_globins','value'].iloc[0]}"
    )


if __name__ == "__main__":
    main()
