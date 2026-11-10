"""Unadjusted, prespecified gene-phenotype enrichment with explicit eligibility."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio


def require_columns(frame, columns, label):
    missing = set(columns) - set(frame.columns)
    if missing:
        raise ValueError(f"{label}: missing columns {sorted(missing)}")


def clean_keys(frame, columns):
    frame = frame.copy()
    for column in columns:
        if frame[column].isna().any():
            raise ValueError(f"Missing {column}")
        frame[column] = frame[column].astype(str).str.strip()
        if frame[column].eq("").any():
            raise ValueError(f"Empty {column}")
    return frame


def binary(series, label, allow_missing=False):
    # Missing phenotype status means excluded/unknown, never control.
    numeric = pd.to_numeric(series, errors="coerce")
    if (series.notna() & numeric.isna()).any() or not numeric.dropna().isin([0, 1]).all():
        raise ValueError(f"{label} must contain only 0/1" + (" or missing" if allow_missing else ""))
    if not allow_missing and numeric.isna().any():
        raise ValueError(f"Missing {label}")
    return numeric


def bh_adjust(values):
    """BH across all planned rows; non-estimable hypotheses count as p=1."""
    p = np.asarray(values, dtype=float)
    ranked_p = np.where(np.isnan(p), 1.0, p)
    order = np.argsort(ranked_p)
    adjusted = np.minimum.accumulate(
        (ranked_p[order] * len(p) / np.arange(1, len(p) + 1))[::-1]
    )[::-1].clip(0, 1)
    result = np.empty(len(p))
    result[order] = adjusted
    result[np.isnan(p)] = np.nan
    return result


def run_enrichment(eligible, underliers, phenotypes, mapping, alternative="greater"):
    """One test per (gene, category, phenoname); no pooling across genes.

    eligible: subject_id,gene -- all callable, QC-passing gene-person pairs.
    underliers: subject_id,gene -- positive underlier calls, subset of eligible.
    phenotypes: subject_id,category,phenoname,status -- 1=case,0=control,NA=exclude.
    mapping: gene,category,phenoname -- clinically selected hypotheses.
    Repeated diagnosis rows must be reduced to statuses before calling this API.
    """
    if alternative not in ("greater", "two-sided"):
        raise ValueError("alternative must be greater or two-sided")
    require_columns(eligible, ["subject_id", "gene"], "eligible")
    require_columns(underliers, ["subject_id", "gene"], "underliers")
    require_columns(phenotypes, ["subject_id", "category", "phenoname", "status"], "phenotypes")
    require_columns(mapping, ["gene", "category", "phenoname"], "mapping")
    eligible = clean_keys(eligible, ["subject_id", "gene"]).drop_duplicates(["subject_id", "gene"])
    underliers = clean_keys(underliers, ["subject_id", "gene"]).drop_duplicates(["subject_id", "gene"])
    mapping = clean_keys(mapping, ["gene", "category", "phenoname"])
    phenotypes = clean_keys(phenotypes, ["subject_id", "category", "phenoname"])
    if eligible.empty or mapping.empty:
        raise ValueError("eligible and mapping must not be empty")
    outside = underliers.merge(eligible[["subject_id", "gene"]], how="left", indicator=True)
    if outside["_merge"].eq("left_only").any():
        raise ValueError("Underlier calls contain gene-person pairs absent from eligible")
    if set(underliers.gene) - set(mapping.gene):
        raise ValueError("Some underlier genes have no selected phenotype mapping")
    if set(mapping.gene) - set(eligible.gene):
        raise ValueError("Some mapped genes have no eligible participants")
    phenotypes["status"] = binary(phenotypes.status, "phenotype status", allow_missing=True)
    keys = ["subject_id", "category", "phenoname"]
    if phenotypes.duplicated(keys).any():
        raise ValueError("Phenotypes must have one status per subject/category/phenoname")
    rows = []
    for gene, category, pheno in mapping[["gene", "category", "phenoname"]].drop_duplicates().itertuples(index=False, name=None):
        population = eligible.loc[eligible.gene.eq(gene), ["subject_id"]].copy()
        positive = set(underliers.loc[underliers.gene.eq(gene), "subject_id"])
        population["underlier"] = population.subject_id.isin(positive)
        statuses = phenotypes.loc[phenotypes.category.eq(category) & phenotypes.phenoname.eq(pheno), ["subject_id", "status"]]
        joined = population.merge(statuses, how="left", on="subject_id", validate="one_to_one")
        missing = joined.status.isna()
        excluded_underliers = int((missing & joined.underlier).sum())
        excluded_others = int((missing & ~joined.underlier).sum())
        data = joined.loc[~missing]
        a = int((data.underlier & data.status.eq(1)).sum())
        b = int((data.underlier & data.status.eq(0)).sum())
        c = int((~data.underlier & data.status.eq(1)).sum())
        d = int((~data.underlier & data.status.eq(0)).sum())
        table = [[a, b], [c, d]]
        estimable = min(a+b, c+d, a+c, b+d) > 0
        p = float(fisher_exact(table, alternative=alternative).pvalue) if estimable else np.nan
        sample_or = (a*d / (b*c)) if b*c else (np.inf if a*d else np.nan)
        conditional_or = low = high = np.nan
        if estimable:
            estimate = odds_ratio(table, kind="conditional")
            conditional_or = float(estimate.statistic)
            ci = estimate.confidence_interval(confidence_level=0.95)
            low, high = float(ci.low), float(ci.high)
        rows.append(dict(gene=gene, category=category, phenoname=pheno,
                         underlier_cases=a, underlier_controls=b,
                         other_cases=c, other_controls=d,
                         underlier_case_fraction=a/(a+b) if a+b else np.nan,
                         other_case_fraction=c/(c+d) if c+d else np.nan,
                         excluded_underliers=excluded_underliers, excluded_others=excluded_others,
                         odds_ratio=sample_or, conditional_odds_ratio=conditional_or,
                         ci95_low=low, ci95_high=high, p_value=p,
                         alternative=alternative, test_status="tested" if estimable else "not_estimable"))
    result = pd.DataFrame(rows)
    result["q_value_bh"] = bh_adjust(result.p_value)
    return result.sort_values(["q_value_bh", "p_value", "gene"], na_position="last").reset_index(drop=True)
