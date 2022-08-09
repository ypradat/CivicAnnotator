from dataclasses import dataclass
from civic import main

@dataclass
class Args:
    input: str
    output: str
    category: str
    civic: str = "data/01-Jan-2022-ClinicalEvidenceSummaries_Annotated.xlsx"
    rules: str = "data/CIViC_Curation_And_Rules_Mutation.xlsx"
    tumor_types: str = ""


def test_annotate_maf():
    args = Args(input="tests/data/example_maf.tsv", output="tests/data/example_maf_annotated.tsv", category="mut")
    main(args)


def test_annotate_cna():
    args = Args(input="tests/data/example_cna.tsv", output="tests/data/example_cna_annotated.tsv", category="cna")
    main(args)


def test_annotate_fus():
    args = Args(input="tests/data/example_fus.tsv", output="tests/data/example_fus_annotated.tsv", category="fus")
    main(args)


def test_annotate_alt():
    args = Args(input="tests/data/example_alt.tsv", output="tests/data/example_alt_annotated.tsv", category="alt")
    main(args)


def test_annotate_maf_egfr():
    args = Args(input="tests/data/example_maf_egfr.tsv", output="tests/data/example_maf_egfr_annotated.tsv", category="mut")
    main(args)
