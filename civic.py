# -*- coding: utf-8 -*-
"""
@created: Jan 14 2022
@modified: May 13 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Use CIViC database to annotate variants.
"""

import argparse
from functools import reduce
import numpy as np
import pandas as pd
import re

from utils import explode_df, build_alt_identifier, clean_alt_identifier
from civic_preprocessor import CivicPreprocessor
from civic_annotator import CivicAnnotator

# functions ============================================================================================================

def annotate_tumor_types(df_alt, civic, tumor_types, category, rules):
    # load CIViC database
    df_civ = pd.read_excel(civic)

    # split tumor types
    tumor_types_split = tumor_types.split("|")

    # apply a series of filters on CIViC table to select only relevant lines
    preprocessor = CivicPreprocessor(tumor_types=tumor_types_split, category=category, rules=rules)
    df_civ = preprocessor.run(df_civ)

    # perform the annotation
    annotator = CivicAnnotator(category=category, rules=rules)
    df_ann = annotator.run(df_civ, df_alt)
    df_ann = clean_alt_identifier(df_ann, category)

    return df_ann


def get_columns_civic():
    return ["CIViC_Matching_Disease", "CIViC_Matching_Type", "CIViC_Matching_Gene_Variant",
            "CIViC_Matching_Evidence_Id", "CIViC_Matching_Citation",
            "Predictive:N:A","Predictive:N:B","Predictive:N:C","Predictive:N:D","Predictive:N:E",
            "Predictive:P:A","Predictive:P:B","Predictive:P:C","Predictive:P:D","Predictive:P:E",
            "Diagnostic:N:A","Diagnostic:N:B","Diagnostic:N:C","Diagnostic:N:D","Diagnostic:N:E",
            "Diagnostic:P:A","Diagnostic:P:B","Diagnostic:P:C","Diagnostic:P:D","Diagnostic:P:E",
            "Prognostic:N:A","Prognostic:N:B","Prognostic:N:C","Prognostic:N:D","Prognostic:N:E",
            "Prognostic:P:A","Prognostic:P:B","Prognostic:P:C","Prognostic:P:D","Prognostic:P:E"]

def main(args):
    # load input alterations table
    df_alt = pd.read_table(args.input)

    if df_alt.shape[0]==0:
        df_ann = df_alt.copy()
    else:
        cols_old = df_alt.columns.tolist()

        # build identifier
        df_alt = build_alt_identifier(df_alt, args.category)

        # checks
        if len(args.tumor_types) > 0:
            print("-INFO: using tumor types from option --tumor_types")
            print("-INFO: processing %d lines from tumor type %s" % (len(df_alt), args.tumor_types))
            df_ann = annotate_tumor_types(df_alt=df_alt, civic=args.civic, tumor_types=args.tumor_types,
                                          category=args.category, rules=args.rules)
        else:
            if "Civic_Disease" not in df_alt:
                raise ValueError("-ERROR: you have not specified a value for --tumor_types and the column" +
                                 " Civic_Disease is absent from the input table %s." % args.input)
            elif df_alt["Civic_Disease"].isnull().sum() > 0:
                raise ValueError("-ERROR: you have not specified a value for --tumor_types and the column" +
                                 " Civic_Disease from the input table %s contains NaN." % args.input)

            tumor_types_all = df_alt["Civic_Disease"].unique().tolist()

            # annotate group of samples sharing same tumor types
            # NOTE: the value of --tumor_types has priority over the presence the values of Civic_Disease.
            print("-INFO: the civic annotator will process %d group(s) of samples from %d different tumor type(s)..." %
                  (len(tumor_types_all), len(tumor_types_all)))

            dfs_ann = []
            for tumor_types in tumor_types_all:
                df_alt_sub = df_alt.loc[df_alt["Civic_Disease"]==tumor_types].copy()
                print("="*40)
                print("-INFO: processing %d lines from tumor type %s" % (len(df_alt_sub), tumor_types))
                df_ann_sub = annotate_tumor_types(df_alt=df_alt_sub, civic=args.civic, tumor_types=tumor_types,
                                                  category=args.category, rules=args.rules)
                dfs_ann.append(df_ann_sub)

            # concatenate annotated samples per group of tumor type
            df_ann = pd.concat(dfs_ann)

        if args.category != "alt":
            if df_ann.shape[0]>0 and "CIViC_Matching_Disease" in df_ann:
                # order columns
                cols_civ = get_columns_civic()
                cols_new = list(set(df_ann.columns.tolist()).difference(set(cols_old)))
                cols_new_ord = [x for x in cols_civ if x in cols_new]
                df_ann = df_ann[cols_old + cols_new_ord]
            else:
                # empty dataframe
                df_ann = df_ann.iloc[:0,:]

    # save
    print("="*40)
    df_ann.to_csv(args.output, sep="\t", index=False)
    print("-INFO: annotated table saved at %s" % args.output)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform annotation of alterations using CIViC database.")
    parser.add_argument('--input', type=str, help='Path to input table of variants.',
                        default="tests/data/example_maf.tsv")
    parser.add_argument('--civic', type=str, help='Path to CIViC database of clinical evidence summaries.',
                        default="data/01-Jan-2022-ClinicalEvidenceSummaries_Annotated.xlsx")
    parser.add_argument('--rules', type=str, help='Path to table of rules for cleaning the database and matching.',
                        default="data/CIViC_Curation_And_Rules_Mutation.xlsx")
    parser.add_argument('--category', type=str, default="mut",
                        help="Choose one of alt, cna, mut or fus. alt considers only annotations spanning 2 types" + \
                         " of data or more.")
    parser.add_argument('--tumor_types', type=str, default='',
                        help="Tumor type designations in CIViC, separated with |. In case the input table" +
                        " contains Civic_Disease column, the value of this option has priority over the values of"
                        " Civic_Disease column.")
    parser.add_argument('--output', type=str, help='Path to output table of variants.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
