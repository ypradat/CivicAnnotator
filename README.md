# CivicAnnotator
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![codecov](https://codecov.io/gh/ypradat/CivicAnnotator/branch/master/graph/badge.svg?token=YVVOSTUAWS)](https://codecov.io/gh/ypradat/CivicAnnotator)

This repository hosts some scripts for performing non-exhaustive annotations of mutations, fusions or copy-number events
using the [CIViC database](https://civicdb.org/home). The script will leave the number of rows untouched and will add
extra columns in the following format

- `[Type]:[Direction]:[Level]`

with

- **Type**: one of *Predictive* or *Prognostic*.
- **Direction**: either *N* (Negative) or *P* (Positive).
- **Level**: one of the 5 levels *A, B, C, D* and *E* used by CIViC.

**NOTE**: not all possible combinations of these three tags are automatically added to the input table in order to avoid
completely empty columns.

The script will also add the following 5 columns that detail how the alteration was matched in CIViC in order to
allow for a manual review of the annotations.

- `CIViC_Matching_Disease`
- `CIViC_Matching_Type`
- `CIViC_Matching_Gene_Variant`
- `CIViC_Matching_Evidence_Id`
- `CIViC_Matching_Citation`

# Requirements

In order to run the analysis, ensure the `conda` command is available and ensure you have all files in the `data` folder
available.

## Conda

Please read the [official conda
documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) for instructions on how to
install the command. 

NOTE: Miniconda is a more lightweight version of Anaconda, you may opt for Miniconda if you are not interested in all the
extra tools and libraries shipped with Anaconda.

Any basic python environment should be enough to run the script `civic_annotator.py`. For the sake of completeness, an
example python environment is provided in `envs/civic-annotator.yaml`. You may install it locally through the command

```
conda env create -f envs/civic-annotator.yaml
```

## Data

You need to have the following files in your local `data` folder in order for the script to have the resources it needs

- **01-Jan-2022-ClinicalEvidenceSummaries_Annotated.xlsx**
- **01-Jan-2022-GeneSummaries.tsv**
- **CIViC_Curation_And_Rules_Mutation.xlsx**

Please read the README `data/README.md` for more details about what these files are.

**NOTE**: it is not mandatory to use the CIViC release of 01-Jan-2022, you may use any other release.

# Examples

Examples of commands to annotate mutations, fusions and copy-number events are provided in the `tests` folder. The
commands are described in the bash scripts while the input files to these commands are provided in `tests/data`.

# References

1.	Griffith, M. et al. CIViC is a community knowledgebase for expert crowdsourcing the clinical interpretation of
variants in cancer. Nat. Genet. 49, 170–174 (2017).
