# DATA folder

The files in this folder are available upon request. The contents of the folder are

- **01-Jan-2022-ClinicalEvidenceSummaries_Annotated.xlsx**: Retrieved from the january 2022 release here
  <https://civicdb.org/releases>. Some manual curation of values in `variant` and `drugs` column was performed and
  a column `category` was added which classifies variants into any combination of `m` (mutations), `f` (fusions) and `c`
  (copy-number alterations).

  NOTE: manually curated values are highlighted in green.

- **01-Jan-2022-GeneSummaries.tsv**: Retrieved from the january 2022 release here
  <https://civicdb.org/releases>

- **CIViC_Curation_And_Rules_Mutation.xlsx**: The sheets in this excel workbook define values/rows to be excluded/included
  during the annotation of mutations as well as rules for matching between values of `Variant_Classification` in the
  input MAF table and values of `variant` in the civic table of clinical evidence summaries.

  NOTE: This sheet is a joint work between PhD S. Nikolaev and PhDc Y. Pradat.

- **Table_Correspondence_Tumor_Type.xlsx**: A table for converting between different classifications of tumors. This
  table contains TCGA classification, MSKCC classification and CIViC classification of tumors.
