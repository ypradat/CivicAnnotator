#!/bin/bash

maf=examples/data/example_maf.tsv
civ=data/01-Jan-2022-ClinicalEvidenceSummaries_Annotated.xlsx
rul=data/CIViC_Curation_And_Rules_Mutation.xlsx
category=mut
tumor_types="Lung Cancer|Lung Carcinoma|Lung Non-small Cell Carcinoma|Lung Adenocarcinoma|Solid Tumor|Cancer"
out=examples/data/example_maf_annotated.tsv
log=examples/logs/example_maf_annotated.log

python civic_annotator.py \
    --input ${maf} \
    --civic ${civ} \
    --rules ${rul} \
    --category ${category} \
    --tumor_types "${tumor_types}" \
    --output ${out} &> ${log}
