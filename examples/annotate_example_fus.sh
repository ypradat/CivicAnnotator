#!/bin/bash

fus=examples/data/example_fus.tsv
civ=data/01-Jan-2022-ClinicalEvidenceSummaries_Annotated.xlsx
rul=data/CIViC_Curation_And_Rules_Mutation.xlsx
category=fus
out=examples/data/example_fus_annotated.tsv
log=examples/logs/example_fus_annotated.log

python civic_annotator.py \
    --input ${fus} \
    --civic ${civ} \
    --rules ${rul} \
    --category ${category} \
    --output ${out} &> ${log}
