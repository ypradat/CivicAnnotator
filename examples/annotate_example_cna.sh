#!/bin/bash

cna=examples/data/example_cna.tsv
civ=data/01-Jan-2022-ClinicalEvidenceSummaries_Annotated.xlsx
rul=data/CIViC_Curation_And_Rules_Mutation.xlsx
category=cna
tumor_types="Lung Cancer|Lung Carcinoma|Lung Non-small Cell Carcinoma|Lung Adenocarcinoma|Solid Tumor|Cancer"
out=examples/data/example_cna_annotated.tsv
log=examples/logs/example_cna_annotated.log

python civic_annotator.py \
    --input ${cna} \
    --civic ${civ} \
    --rules ${rul} \
    --category ${category} \
    --tumor_types "${tumor_types}" \
    --output ${out} &> ${log}