#!/bin/bash

maf=tests/data/example_maf_egfr.tsv
civ=data/01-Jan-2022-ClinicalEvidenceSummaries_Annotated.xlsx
rul=data/CIViC_Curation_And_Rules_Mutation.xlsx
category=mut
out=tests/data/example_maf_eggr_annotated.tsv
log=tests/logs/example_maf_egfr_annotated.log

python civic.py \
    --input ${maf} \
    --civic ${civ} \
    --rules ${rul} \
    --category ${category} \
    --output ${out} &> ${log}
