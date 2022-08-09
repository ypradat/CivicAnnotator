#!/bin/bash

alt=tests/data/example_alt.tsv
civ=data/01-Jan-2022-ClinicalEvidenceSummaries_Annotated.xlsx
rul=data/CIViC_Curation_And_Rules_Mutation.xlsx
category=alt
out=tests/data/example_alt_annotated.tsv
log=tests/logs/example_alt_annotated.log

python civic.py \
    --input ${alt} \
    --civic ${civ} \
    --rules ${rul} \
    --category ${category} \
    --output ${out} &> ${log}
