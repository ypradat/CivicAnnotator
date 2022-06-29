# -*- coding: utf-8 -*-
"""
@created: Jan 14 2022
@modified: May 11 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Useful functions.
"""

import argparse
import numpy as np
import pandas as pd
import re

AMINO_ACIDS = {"Ala": "A",
               "Arg": "R",
               "Asn": "N",
               "Asp": "D",
               "Cys": "C",
               "Gln": "Q",
               "Glu": "E",
               "Gly": "G",
               "His": "H",
               "Ile": "I",
               "Leu": "L",
               "Lys": "K",
               "Met": "M",
               "Phe": "F",
               "Pro": "P",
               "Pyl": "O",
               "Ser": "S",
               "Sec": "U",
               "Thr": "T",
               "Trp": "W",
               "Tyr": "Y",
               "Val": "V",
               "Asx": "B",
               "Glx": "Z",
               "Xaa": "X",
               "Xle": "J",
               "Ter": "*"}

AMINO_ACIDS_REV = {v: k for k,v in AMINO_ACIDS.items()}


# functions ============================================================================================================

def explode_df(df, cols, sep=',', fill_value='', preserve_index=False):
    # transform comma-separated to list
    df = df.assign(**{col:df[col].str.split(sep) for col in cols}).copy()
    if (cols is not None and len(cols) > 0 and not isinstance(cols, (list, tuple, np.ndarray, pd.Series))):
        cols = [cols]
    # calculate lengths of lists
    lens = df[cols[0]].str.len()
    # format NaN to [NaN] and strip unwanted characters
    for col in cols:
        df.loc[df[col].isnull(), col] = df.loc[df[col].isnull(), col].apply(lambda x: [np.nan])
        df.loc[lens > 1, col] = df.loc[lens > 1, col].apply(lambda x: [y.strip() for y in x])
    # all columns except `cols`
    idx_cols = df.columns.difference(cols)
    # preserve original index values    
    idx = np.repeat(df.index.values, lens)
    # create "exploded" DF
    df_expanded = (pd.DataFrame({col:np.repeat(df[col].values, lens) for col in idx_cols},
                index=idx).assign(**{col:np.concatenate(df.loc[lens>0, col].values) for col in cols}))
    # append those rows that have empty lists
    if (lens == 0).any():
        # at least one list in cells is empty
        df_expanded = (df_expanded.append(df.loc[lens==0, idx_cols], sort=False).fillna(fill_value))
    # revert the original index order
    df_expanded = df_expanded.sort_index()
    # reset index if requested
    if not preserve_index:
        df_expanded = df_expanded.reset_index(drop=True)
    return df_expanded


def shorten_aa(x):
    for aa_l, aa_s in AMINO_ACIDS.items():
        x = x.replace(aa_l, aa_s)
    return x


def expand_aa(x):
    x_split = list(x)
    x_expan = []
    for s in x_split:
        if s in AMINO_ACIDS_REV:
            x_expan.append(AMINO_ACIDS_REV[s])
        else:
            x_expan.append(s)
    return "".join(x_expan)



def build_row_identifier(df, mode=1, verbose=False):
    if mode==1:
        cols_id = ["evidence_id", "variant_id", "gene_id"]
    else:
        cols_id = ["gene", "variant", "disease", "drugs", "evidence_type", "evidence_direction", "evidence_level",
                   "clinical_significance"]
    df["Row_Identifier"] = df[cols_id].fillna("-").astype(str).apply("_".join, axis=1)
    n_unq = df["Row_Identifier"].nunique()
    n_row = df.shape[0]
    if verbose:
        if n_unq < n_row:
            print("-warning! %d unique values of Row_Identifier < %d rows" % (n_unq, n_row))
    return df


def build_mut_identifier(df, col_gene="Hugo_Symbol", col_start="Start_Position", col_end="End_Position",
                         col_ref="Reference_Allele", col_alt="Tumor_Seq_Allele2"):
    dfc = df.copy()
    cols = [col_gene, col_start, col_end, col_ref, col_alt]
    dfc[col_start] =  dfc[col_start].apply(lambda x: "%d" % x if not np.isnan(x) else "")
    dfc[col_end] =  dfc[col_end].apply(lambda x: "%d" % x if not np.isnan(x) else "")
    dfc["Mut_Identifier"] = dfc[cols].fillna("-").astype(str).apply("_".join, axis=1)
    df["Mut_Identifier"] = dfc["Mut_Identifier"]
    return df


def build_fus_identifier(df, col_gene_1="Gene_1", col_gene_2="Gene_2"):
    dfc = df.copy()
    cols = [col_gene_1, col_gene_2]
    dfc["Fus_Identifier"] = dfc[cols].fillna("N/A").astype(str).apply("--".join, axis=1)
    df["Fus_Identifier"] = dfc["Fus_Identifier"]
    return df


def build_cna_identifier(df, col_gene="Hugo_Symbol", col_alt="Alteration"):
    dfc = df.copy()
    cols = [col_gene, col_alt]
    dfc[col_alt] = dfc[col_alt].str.upper()
    dfc["Cna_Identifier"] = dfc[cols].fillna("N/A").astype(str).apply("--".join, axis=1)
    df["Cna_Identifier"] = dfc["Cna_Identifier"]
    return df


def build_alt_identifier(df, category, **kwargs):
    if category=="mut":
        return build_mut_identifier(df, **kwargs)
    elif category=="fus":
        return build_fus_identifier(df, **kwargs)
    elif category=="cna":
        return build_cna_identifier(df, **kwargs)
    elif category=="alt":
        df_c = df.loc[df["Alteration_Category"]=="c"].copy()
        df_f = df.loc[df["Alteration_Category"]=="f"].copy()
        df_m = df.loc[df["Alteration_Category"]=="m"].copy()
        df_c = build_cna_identifier(df_c, **kwargs)
        df_f = build_fus_identifier(df_f, **kwargs)
        df_m = build_mut_identifier(df_m, **kwargs)
        return pd.concat((df_c, df_f, df_m), axis=0)
    else:
        raise ValueError("-ERROR! Unsupported value %s for category." %  category)


def clean_alt_identifier(df, category):
    if category=="mut":
        del df["Mut_Identifier"]
    elif category=="fus":
        del df["Fus_Identifier"]
    elif category=="cna":
        del df["Cna_Identifier"]
    elif category=="alt":
        for col in ["Mut_Identifier", "Fus_Identifier", "Cna_Identifier"]:
            if col in df:
                del df[col]

    return df
