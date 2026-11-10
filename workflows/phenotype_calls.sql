-- Run inside the All of Us Researcher Workbench BigQuery environment.
-- Supply underlier and eligible subject IDs as query parameters/temporary tables.
-- Returns explicit 1/0/missing phenotype status; never label absent EHR as control.
-- Replace {{...}} placeholders with validated values in the Workbench notebook.
WITH eligible AS (
  SELECT DISTINCT subject_id, gene FROM {{eligible_table}}
), underliers AS (
  SELECT DISTINCT subject_id, gene FROM {{underlier_table}}
), phecode_icd AS (
  SELECT DISTINCT phenoname, ICD, vocabulary_id
  FROM {{phecodex_mapping_table}}
  WHERE phenoname IN UNNEST(@selected_phenonames)
), condition_rows AS (
  SELECT CAST(co.person_id AS STRING) AS subject_id,
         co.condition_start_date,
         co.condition_source_value,
         co.condition_source_concept_id,
         c.vocabulary_id
  FROM {{cdm_dataset}}.condition_occurrence co
  JOIN {{cdm_dataset}}.concept c
    ON c.concept_id = co.condition_source_concept_id
  WHERE co.person_id IN (SELECT DISTINCT SAFE_CAST(subject_id AS INT64) FROM eligible)
), phenotypes AS (
  SELECT e.subject_id, m.phenoname,
         COUNTIF(REGEXP_REPLACE(UPPER(cr.condition_source_value), r'\.', '') =
                 REGEXP_REPLACE(UPPER(m.ICD), r'\.', '')) AS matching_records
  FROM eligible e
  CROSS JOIN (SELECT DISTINCT phenoname FROM phecode_icd) p
  JOIN phecode_icd m USING (phenoname)
  LEFT JOIN condition_rows cr
    ON cr.subject_id=e.subject_id
   AND cr.vocabulary_id=m.vocabulary_id
   AND REGEXP_REPLACE(UPPER(cr.condition_source_value), r'\.', '') =
       REGEXP_REPLACE(UPPER(m.ICD), r'\.', '')
  GROUP BY e.subject_id, m.phenoname
)
SELECT e.gene, p.phenoname, e.subject_id,
       IF(p.matching_records > 0, 1, 0) AS status
FROM eligible e
JOIN phenotypes p USING (subject_id)
JOIN (SELECT DISTINCT phenoname FROM phecode_icd) m USING (phenoname)
;
