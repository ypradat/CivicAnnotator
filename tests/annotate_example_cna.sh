#!/bin/bash

cna=tests/data/example_cna.tsv
civ=data/01-Jan-2022-ClinicalEvidenceSummaries_Annotated.xlsx
rul=data/CIViC_Curation_And_Rules_Mutation.xlsx
category=cna
out=tests/data/example_cna_annotated.tsv
log=tests/logs/example_cna_annotated.log

python civic.py \
    --input ${cna} \
    --civic ${civ} \
    --rules ${rul} \
    --category ${category} \
    --output ${out} &> ${log}
