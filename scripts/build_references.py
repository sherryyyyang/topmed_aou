"""Build public references from official workbook exports; no participant data.

Usage: python scripts/build_references.py --acmg-workbook public.xlsx
       --aou-workbook All_by_All_phenotypes_v7.xlsx
"""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import openpyxl
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ACMG_URL = "https://docs.google.com/spreadsheets/d/1ecSUe0bJ-2gAMXBsrf5U1wMEDOmiIjtR_A5388uj7HY/edit"
AOU_URL = "https://docs.google.com/spreadsheets/d/1sRxX4nzEPCKnzCgjJiwnt9h3HNrL_-zwpTjqZsikn34/edit?gid=1330947904"

# Explicit disease-label rules, not fuzzy matches. Primary means a disease-level
# outcome, not proof of inherited etiology. Secondary outcomes are proxies or
# manifestations; they are excluded unless explicitly requested in the notebook.
# Source disease spellings are preserved in the output, including source typos.
RULES = [
    ("adrenoleukodystrophy", "GE_962.392", "primary", "Disease-level phenotype"),
    ("thoracic aortic aneurysm", "CV_438.12", "secondary", "Thoracic aneurysm proxy; does not establish familial etiology"),
    ("hypertrophic cardiomyopathy", "CV_414.1", "primary", "Disease-level phenotype"),
    ("hemorrhagic telangiectasia", "GE_975.2", "primary", "Disease-level phenotype"),
    ("adenomatous polyposis", "CA_101.41", "secondary", "Colorectal cancer manifestation; not a polyposis diagnosis"),
    ("hypercholesterolemia", "GE_965.4", "primary", "Familial hypercholesterolemia phenotype"),
    ("wilson disease", "GE_967.1", "primary", "Disease-level phenotype"),
    ("dilated cardiomyopathy", "CV_414.2", "primary", "Disease-level phenotype"),
    ("dliated cardiomyopathy", "CV_414.2", "primary", "Source spelling corrected for mapping only"),
    ("myofibrillar myopathy", "NS_340", "secondary", "Broad myopathy proxy; many other causes"),
    ("juvenile polyposis", "CA_101.41", "secondary", "Colorectal cancer manifestation; not a polyposis diagnosis"),
    ("hereditary breast and ovarian cancer", "CA_105", "secondary", "Breast cancer manifestation; not hereditary syndrome diagnosis"),
    ("hereditary breast and ovarian cancer", "CA_106.3", "secondary", "Ovarian cancer manifestation; not hereditary syndrome diagnosis"),
    ("hereditary breast cancer", "CA_105", "secondary", "Breast cancer manifestation; not hereditary syndrome diagnosis"),
    ("biotinidase deficiency", "GE_967.6", "primary", "Disease-level phenotype"),
    ("malignant hyperthermia", "GE_981.8", "primary", "Recorded disease; susceptibility may remain clinically silent"),
    ("long-qt syndrome", "CV_416.6", "primary", "Long QT phenotype; not subtype specific"),
    ("long qt syndrome", "CV_416.6", "primary", "Long QT phenotype; not subtype specific"),
    ("catecholaminergic polymorphic ventricular tachycardia", "CV_416.12", "secondary", "Broad ventricular tachycardia proxy; not CPVT-specific"),
    ("ehlers-danlos syndrome, vascular", "GE_978.222", "primary", "Vascular Ehlers-Danlos phenotype"),
    ("arrhythmogenic right ventricular cardiomyopathy", "CV_414", "secondary", "Broad cardiomyopathy proxy; not ARVC-specific"),
    ("marfan syndrome", "GE_978.21", "primary", "Disease-level phenotype"),
    ("pompe disease", "GE_964.32", "primary", "Disease-level phenotype"),
    ("fabry disease", "GE_963.11", "primary", "Disease-level phenotype"),
    ("hemochromatosis", "EM_247.71", "secondary", "Hemochromatosis proxy; not specific to HFE genotype"),
    ("maturity-onset of diabetes", "EM_202", "secondary", "Broad diabetes proxy; not MODY-specific"),
    ("paraganglioma-pheochromocytoma", "CA_114.6", "secondary", "Pheochromocytoma manifestation; not hereditary syndrome diagnosis"),
    ("paraganglioma-pheochromocytoma", "CA_114.5", "secondary", "Paraganglioma manifestation; not hereditary syndrome diagnosis"),
    ("multiple endocrine neoplasia type 1", "GE_961.21", "primary", "Disease-level phenotype"),
    ("multiple endocrine neoplasia type 2a", "GE_961.22", "primary", "Disease-level phenotype"),
    ("multiple endocrine neoplasia type 2b", "GE_961.23", "primary", "Disease-level phenotype"),
    ("lynch syndrome", "CA_101.41", "secondary", "Colorectal cancer manifestation; not Lynch syndrome diagnosis"),
    ("mutyh-associated polyposis", "CA_101.41", "secondary", "Colorectal cancer manifestation; not MUTYH polyposis diagnosis"),
    ("nf2-related schwannomatosis", "GE_961.112", "primary", "Neurofibromatosis type 2 is the legacy diagnostic label"),
    ("ornithine transcarbamylase", "GE_962.54", "secondary", "Broad ornithine metabolism disorder proxy"),
    ("pten hamartoma tumor", "GE_961.14", "primary", "Disease-level phenotype"),
    ("rpe65-related retinopathy", "SO_374.4", "secondary", "Broad retinopathy proxy; not RPE65-specific"),
    ("loeys-dietz", "CV_438.12", "secondary", "Thoracic aneurysm manifestation; not Loeys-Dietz diagnosis"),
    ("tuberous sclerosis", "GE_961.12", "primary", "Disease-level phenotype"),
    ("hereditary transthyretin", "EM_249", "secondary", "Broad amyloidosis proxy; includes other etiologies"),
    ("von hippel-lindau", "GE_961.13", "primary", "Disease-level phenotype"),
]


def extract_sheet(path, sheet, width):
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    rows = [tuple(r[:width]) for r in workbook[sheet].values if r and r[0] is not None]
    return pd.DataFrame(rows[1:], columns=rows[0])


def build(acmg_workbook, aou_workbook):
    data = ROOT / "data"
    source = data / "sources"
    source.mkdir(parents=True, exist_ok=True)
    genes = extract_sheet(acmg_workbook, "ACMG SF v3.3 List", 8)
    genes.columns = ["gene", "gene_mim", "disease", "disease_mim", "disease_category", "inheritance", "introduced_in_sf_version", "variant_reporting_guidance"]
    for col in ["gene_mim", "disease_mim"]:
        genes[col] = genes[col].map(lambda v: str(int(v)) if isinstance(v, (int, float)) and pd.notna(v) else str(v or ""))
    if genes.gene.nunique() != 84:
        raise ValueError("Expected the pinned ACMG SF v3.3 list with 84 genes")
    genes["list_version"] = "ACMG SF v3.3"
    genes["source_url"] = ACMG_URL
    genes.to_csv(data / "gene_disease_associations.csv", index=False)
    genes[["gene"]].drop_duplicates().sort_values("gene").to_csv(data / "acmg_genes.csv", index=False)
    catalog = pd.concat([
        extract_sheet(aou_workbook, name, 9).assign(catalog_sheet=name)
        for name in ["Analyzed_phenotypes", "Excluded_phenotypes"]
    ], ignore_index=True)
    # Keep only catalog metadata, not published ancestry case counts.
    catalog = catalog.loc[catalog.category.eq("phecodeX"), ["phenoname", "category", "description", "catalog_sheet"]].drop_duplicates()
    catalog.to_csv(source / "aou_v7_phecodex_catalog.csv", index=False)
    info = pd.read_csv(source / "phecodeX_info.csv", encoding="latin1", dtype={"phecode": str})
    info = info.rename(columns={"phecode": "phenoname", "phecode_string": "description"})
    info.to_csv(source / "phecodex_v1_info.csv", index=False)
    valid = info.set_index("phenoname")
    mappings = []
    coverage = []
    for row in genes.to_dict("records"):
        matches = [rule for rule in RULES if rule[0] in row["disease"].lower()]
        for term, code, tier, note in matches:
            if code not in valid.index:
                raise ValueError(f"Unknown PhecodeX definition: {code}")
            mappings.append({**row, "category": "phecodeX", "phenoname": code,
                             "phenotype_description": valid.loc[code, "description"],
                             "outcome_tier": tier, "mapping_note": note,
                             "mapping_source": "Explicit research crosswalk v1; review before confirmatory analysis",
                             "in_aou_v7_catalog": code in set(catalog.phenoname)})
        coverage.append({"gene": row["gene"], "disease": row["disease"],
                         "coverage": "primary" if any(r[2] == "primary" for r in matches) else "secondary_only" if matches else "unmapped",
                         "note": "" if matches else "No sufficiently specific mapping selected; add a reviewed ICD/OMOP definition"})
    mapping = pd.DataFrame(mappings).drop_duplicates(["gene", "disease", "category", "phenoname"])
    mapping.to_csv(data / "gene_phenotype_map.csv", index=False)
    pd.DataFrame(coverage).to_csv(data / "mapping_coverage.csv", index=False)
    summary = pd.DataFrame(coverage).groupby("gene").coverage.agg(lambda s: "primary" if "primary" in set(s) else "secondary_only" if "secondary_only" in set(s) else "unmapped")
    manifest = {"built_utc": datetime.now(timezone.utc).isoformat(), "acmg_source": ACMG_URL,
                "aou_catalog_source": AOU_URL, "acmg_genes": int(genes.gene.nunique()),
                "gene_disease_rows": len(genes), "mapping_rows": len(mapping),
                "gene_coverage": summary.value_counts().to_dict(),
                "phecodex_version": "1.0",
                "phecodex_commit": "237f7c58daedd0188e4f5ca0201e4fc83083e53a",
                "sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(data.glob("*.csv"))}}
    (data / "reference_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: manifest[k] for k in ["acmg_genes", "gene_disease_rows", "mapping_rows", "gene_coverage"]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--acmg-workbook", required=True)
    parser.add_argument("--aou-workbook", required=True)
    args = parser.parse_args()
    build(args.acmg_workbook, args.aou_workbook)
