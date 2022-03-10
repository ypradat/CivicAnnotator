# -*- coding: utf-8 -*-
"""
@created: Jan 14 2022
@modified: Mar 08 2022
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
    else:
        raise ValueError("-ERROR! Unsupported value %s for category." %  category)


def clean_alt_identifier(df, category):
    if category=="mut":
        del df["Mut_Identifier"]
    elif category=="fus":
        del df["Fus_Identifier"]
    elif category=="cna":
        del df["Cna_Identifier"]

    return df


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


class CivicPreprocessor(object):
    def __init__(self, tumor_types, category, rules):
        self.tumor_types = tumor_types
        self.category = category
        self.rules = rules


    def _check_args(self, df_civ):
        # check tumor types
        if self.tumor_types is not None and len(self.tumor_types)>0:
            civ_tumor_types = df_civ["disease"].unique()
            mis_tumor_types = set(self.tumor_types).difference(set(civ_tumor_types))
            if len(mis_tumor_types)>0:
                print("-WARNING! the following tumor types are absent from CIViC database!")
                print("\t" + "\n\t".join(list(mis_tumor_types)))


        # check category
        allowed_categories = ["mut", "cna", "fus"]
        if self.category not in allowed_categories:
            raise ValueError("-the value of category %s is not recognized! Choose one of %s" % \
                             (category, ",".join(allowed_categories)))


    def _filter_tumor_types(self, df_civ):
        if args.tumor_types is not None and len(args.tumor_types)>0:
            mask = df_civ["disease"].isin(args.tumor_types)
            print("-INFO: selected %s/%s evidence lines corresponding to the following tumor types: \n\t%s" \
                  % (sum(mask), len(mask), "\n\t".join(args.tumor_types)))
            return df_civ.loc[mask].copy()
        else:
            return df_civ


    def _filter_category(self, df_civ):
        long2short= {"mut": "m", "cna": "c", "fus": "f"}
        mask = df_civ["category"]==long2short[self.category]
        print("-INFO: selected %s/%s evidence lines corresponding to the %s category" % \
              (sum(mask), len(mask), self.category))
        return df_civ.loc[mask].copy()


    def _filter_using_rules(self, df_civ):
        df_sub = df_civ.copy()

        # exclude rows for which the column `col` contains any of the values in `vals`
        df_rul = pd.read_excel(self.rules, sheet_name="CIViC_Exclude_Values")
        df_rul = df_rul.fillna("N/A")
        for col in df_rul:
            vals = df_rul[col].unique()
            mask = df_sub[col].fillna("N/A").isin(vals)
            print("-INFO: excluded %s/%s evidence lines having the following values of %s: %s" % \
                  (sum(mask), len(mask), col, ",".join(vals)))
            df_sub = df_sub.loc[~mask].copy()

        # exclude rows manually excluded
        df_rul = pd.read_excel(self.rules, sheet_name="CIViC_Exclude_Lines")
        df_rul = build_row_identifier(df_rul, mode=1)
        df_sub = build_row_identifier(df_sub, mode=1)

        col = "Row_Identifier"
        vals = df_rul["Row_Identifier"].unique()
        mask = df_sub[col].isin(vals)
        print("-INFO: excluded %s/%s evidence lines manually reviewed" % (sum(mask), len(mask)))
        df_sub = df_sub.loc[~mask].copy()

        return df_sub


    def _match_specific_variant(self, x, variants):
        if type(x)!=str:
            return False
        else:
            if x.lower() in [v.lower() for v in variants]:
                return True
            else:
                return False


    def _split_hgvs_pc(self, x):
        if "(c." in x:
            hgvsp = x.split("(c.")[0].strip()
            hgvsc_match = re.search("(?<=\()[\.\>\w\+\-\*]+(?=\))", x)
            if hgvsc_match is not None:
                hgvsc = hgvsc_match.group(0)
            else:
                hgvsc = np.nan
        elif " c." in x:
            hgvsp = x.split(" c.")[0].strip()
            hgvsc = "c." + x.split(" c.")[1].strip()

        if hgvsp in ["Splice Site", "Splice Region", "3'UTR alteration", "Intronic deletion"]:
            hgvsp = np.nan
        else:
            hgvsp = "p." + expand_aa(hgvsp)

        return [hgvsp, hgvsc]


    def _reformat_mut_variant_aggregated(self, x):
        regex = "^[A-Z]+[0-9]+[A-Z]+/[A-Z]+$"
        match = re.search(regex, x)

        if x in ["G12/G13", "V600/K601E"]:
            return x
        elif x=="V600E/K and Amplification":
            return "V600E and Amplification/V600K and Amplification"
        elif match is not None:
            # the variant has the format [AA][XX][AA]/[AA] which should be transformed into [AA][XX][AA]/[AA][XX][AA]
            x_split = x.split("/")
            x_base = re.search("^[A-Z]+[0-9]+", x_split[0]).group(0)
            return "/".join([x_split[0], x_base+x_split[1]])
        else:
            print("-warning! the variant %s has unrecognized specific format." % x)
            return x


    def _reformat_mut_variant(self, df_civ):
        # where possible, build correct HGVSp and HGVSc entries. HGVS_Type will serve for specific matching rules
        df_civ["variant_reformatted"] = False
        df_civ["HGVSp"] = np.nan
        df_civ["HGVS_Type"] = np.nan
        df_civ["HGVSc"] = np.nan

        # split variant aggregated by the "/" separator
        mask = df_civ["variant"].apply(lambda x: "/" in x)
        df_civ_a = df_civ.loc[mask].copy()
        df_civ_b = df_civ.loc[~mask].copy()
        if df_civ_a.shape[0]>0:
            df_civ_a["variant"] = df_civ_a["variant"].apply(self._reformat_mut_variant_aggregated)
            df_civ_a = explode_df(df_civ_a, cols=["variant"], sep="/")
        df_civ = pd.concat((df_civ_a, df_civ_b), axis=0)

        # don't reformat specific variants
        df_rul = pd.read_excel(self.rules, sheet_name="CIViC_Variant_Matching")
        mask = df_civ["variant"].apply(self._match_specific_variant, variants=df_rul["Variant"].unique())
        df_civ.loc[mask, "variant_reformatted"] = True

        mask = df_civ["variant"].apply(lambda x: ("(c." in x or " c." in x) if type(x)==str else False)
        if sum(mask)>0:
            df_split = df_civ.loc[mask, "variant"].apply(self._split_hgvs_pc).apply(pd.Series)
            df_civ.loc[mask, "HGVSp"] = df_split.iloc[:,0]
            df_civ.loc[mask, "HGVSc"] = df_split.iloc[:,1]
            df_civ.loc[mask, "variant"] = df_civ.loc[mask, "variant"].apply(lambda x: x.split("(c.")[0].strip())
            df_civ.loc[mask, "variant_reformatted"] = True

        # reformat partial HGVS protein format
        regex_par = "^[A-Z]+[0-9]+$"
        mask = df_civ["variant"].apply(lambda x: re.search(regex_par, x) is not None)
        df_civ.loc[mask, "HGVSp"] = "p." + df_civ.loc[mask, "variant"].apply(expand_aa)
        df_civ.loc[mask, "HGVS_Type"] = "partial"
        df_civ.loc[mask, "variant_reformatted"] = True

        # reformat HGVS protein format sub
        regex_sub = "^[A-Z\*]+[0-9]+[A-Z\*]+$"
        mask = df_civ["variant"].apply(lambda x: re.search(regex_sub, x) is not None)
        df_civ.loc[mask, "HGVSp"] = "p." + df_civ.loc[mask, "variant"].apply(expand_aa)
        df_civ.loc[mask, "HGVS_Type"] = "complete"
        df_civ.loc[mask, "variant_reformatted"] = True

        # reformat HGVS protein del format
        mask = ~df_civ["variant_reformatted"]
        df_civ.loc[mask, "variant"] = df_civ.loc[mask, "variant"].str.replace("DEL", "del")

        regex_del = "(^[A-Z]+[0-9]+del$)|(^[A-Z]+[0-9]+_[A-Z]+[0-9]+del$)"
        mask = df_civ["variant"].apply(lambda x: re.search(regex_del, x) is not None)
        df_civ.loc[mask, "HGVSp"] = "p." + df_civ.loc[mask, "variant"].apply(expand_aa)
        df_civ.loc[mask, "HGVS_Type"] = "complete"
        df_civ.loc[mask, "variant_reformatted"] = True

        # reformat HGVS protein dup format
        mask = ~df_civ["variant_reformatted"]
        df_civ.loc[mask, "variant"] = df_civ.loc[mask, "variant"].str.replace("DUP", "dup")
        regex_dup = "(^[A-Z]+[0-9]+dup$)|(^[A-Z]+[0-9]+_[A-Z]+[0-9]+dup$)"
        mask = df_civ["variant"].apply(lambda x: re.search(regex_dup, x) is not None)
        df_civ.loc[mask, "HGVSp"] = "p." + df_civ.loc[mask, "variant"].apply(expand_aa)
        df_civ.loc[mask, "HGVS_Type"] = "complete"
        df_civ.loc[mask, "variant_reformatted"] = True

        # reformat HGVS protein ins format
        mask = ~df_civ["variant_reformatted"]
        df_civ.loc[mask, "variant"] = df_civ.loc[mask, "variant"].str.replace("DUP", "ins")
        regex_ins = "(^[A-Z]+[0-9]+ins[A-Z]+$)|(^[A-Z]+[0-9]+_[A-Z]+[0-9]+ins[A-Z]+$)"
        mask = df_civ["variant"].apply(lambda x: re.search(regex_ins, x) is not None)
        df_civ.loc[mask, "HGVSp"] = "p." + df_civ.loc[mask, "variant"].apply(expand_aa)
        df_civ.loc[mask, "HGVS_Type"] = "complete"
        df_civ.loc[mask, "variant_reformatted"] = True

        # reformat HGVS protein delins format
        regex_delins = "(^[A-Z]+[0-9]+delins[A-Z]+$)|(^[A-Z]+[0-9]+_[A-Z]+[0-9]+delins[A-Z]+$)"
        mask = df_civ["variant"].apply(lambda x: re.search(regex_delins, x) is not None)
        df_civ.loc[mask, "HGVSp"] = "p." + df_civ.loc[mask, "variant"].apply(expand_aa)
        df_civ.loc[mask, "HGVS_Type"] = "complete"
        df_civ.loc[mask, "variant_reformatted"] = True

        # reformat HGVS protein frameshift format
        mask = ~df_civ["variant_reformatted"]
        df_civ.loc[mask, "variant"] = df_civ.loc[mask, "variant"].str.replace("FS", "fs")
        regex_fs = "(^[A-Z]+[0-9]+[A-Z]*fs[\*]*[A-Z]*[0-9]*$)"
        mask = df_civ["variant"].apply(lambda x: re.search(regex_fs, x) is not None)
        df_civ.loc[mask, "HGVSp"] = "p." + df_civ.loc[mask, "variant"].apply(expand_aa)
        df_civ.loc[mask, "HGVS_Type"] = "complete"
        df_civ.loc[mask, "variant_reformatted"] = True

        # reformat HGVSc
        mask = df_civ["variant"].apply(lambda x: x.startswith("c."))
        df_civ.loc[mask, "HGVSc"] = df_civ.loc[mask, "variant"]
        df_civ.loc[mask, "HGVS_Type"] = "complete"
        df_civ.loc[mask, "variant_reformatted"] = True

        return df_civ.reset_index(drop=True)


    def _reformat_fus_variant(self, df_civ):
        # where possible, build gene_1 and gene_2 columns. If variant is "Fusions", gene_1 is set to gene
        # and gene_2 is set to "Any"
        df_civ["variant_reformatted"] = False
        df_civ["gene_1"] = np.nan
        df_civ["gene_2"] = np.nan
        df_civ["gene_1_exon"] = np.nan
        df_civ["gene_2_exon"] = np.nan

        # trim white spaces
        df_civ["variant"] = df_civ["variant"].str.strip()

        # process "Fusion" variant
        regex_fus = "^Fusion$"
        mask = df_civ["variant"].apply(lambda x: re.search(regex_fus, x, re.IGNORECASE) is not None)
        if sum(mask) > 0:
            df_civ.loc[mask, "gene_1"] = df_civ.loc[mask, "gene"]
            df_civ.loc[mask, "gene_2"] = "Any"
            df_civ.loc[mask, "variant_reformatted"] = True

        # process variants in format gene_1-gene_2
        regex_ggs = "^[A-Z0-9]+-[A-Z0-9]+$"
        mask = df_civ["variant"].apply(lambda x: re.search(regex_ggs, x) is not None)
        if sum(mask) > 0:
            df_ggs = df_civ.loc[mask,"variant"].apply(lambda x: x.split("-")).apply(pd.Series)
            df_civ.loc[mask,"gene_1"] = df_ggs[0]
            df_civ.loc[mask,"gene_2"] = df_ggs[1]
            df_civ.loc[mask, "variant_reformatted"] = True

        # process variants in format gene_1-gene_2 exon_1-exon_2
        regex_ggee = "^[A-Z0-9]+-[A-Z0-9]+\s+[eE]{1}\d+-[eE]{1}\d+$"
        mask = df_civ["variant"].apply(lambda x: re.search(regex_ggee, x) is not None)
        if sum(mask)>0:
            ggee_split = lambda x: x.split()[0].split("-") + x.split()[1].split("-")
            df_ggee = df_civ.loc[mask,"variant"].apply(ggee_split).apply(pd.Series)
            df_ggee[2] = df_ggee[2].apply(lambda x: re.sub("e", "", x, re.IGNORECASE))
            df_ggee[3] = df_ggee[3].apply(lambda x: re.sub("e", "", x, re.IGNORECASE))
            df_civ.loc[mask,"gene_1"] = df_ggee[0]
            df_civ.loc[mask,"gene_2"] = df_ggee[1]
            df_civ.loc[mask,"gene_1_exon"] = df_ggee[2]
            df_civ.loc[mask,"gene_2_exon"] = df_ggee[3]
            df_civ.loc[mask, "variant_reformatted"] = True

        return df_civ.reset_index(drop=True)


    def _reformat_cna_variant(self, df_civ):
        df_civ["variant_reformatted"] = False
        cna_values = ["DELETION", "LOSS", "LOH", "CN LOH", "AMPLIFICATION"]

        # trim white spaces
        df_civ["variant"] = df_civ["variant"].str.strip()

        # harmonized between capital and small letters
        df_civ["variant"] = df_civ["variant"].str.upper()

        # rename some
        old2new = {"COPY NUMBER VARIATION": "CNV",
                   "COPY-NEUTRAL LOSS OF HETEROZYGOSITY": "CN LOH"}
        df_civ["variant"] = df_civ["variant"].replace(old2new)

        # split CNV into DELETION, LOSS, CN LOH, AMPLIFICATION
        if df_civ.shape[0] > 0:
            df_civ["variant"] = df_civ["variant"].replace({"CNV": "/".join(cna_values)})
            df_civ = explode_df(df_civ, cols=["variant"], sep="/")

        # expected values
        mask = df_civ["variant"].isin(cna_values)
        df_civ.loc[mask, "variant_reformatted"] = True

        return df_civ.reset_index(drop=True)


    def _reformat_variant(self, df_civ):
        if self.category=="mut":
            return self._reformat_mut_variant(df_civ)
        elif self.category=="fus":
            return self._reformat_fus_variant(df_civ)
        elif self.category=="cna":
            return self._reformat_cna_variant(df_civ)


    def _add_variant_identifier(self, df_civ):
        if self.category=="mut":
            return build_alt_identifier(df_civ, self.category, col_gene="gene", col_start="start", col_end="stop",
                                        col_ref="reference_bases", col_alt="variant_bases")
        elif self.category=="fus":
            return build_alt_identifier(df_civ, self.category, col_gene_1="gene_1", col_gene_2="gene_2")
        elif self.category=="cna":
            return build_alt_identifier(df_civ, self.category, col_gene="gene", col_alt="variant")


    def run(self, df_civ):
        self._check_args(df_civ)
        df_civ = self._filter_tumor_types(df_civ)
        df_civ = self._filter_category(df_civ)
        df_civ = self._filter_using_rules(df_civ)
        df_civ = self._reformat_variant(df_civ)
        df_civ = self._add_variant_identifier(df_civ)
        return df_civ



class CivicAnnotator(object):
    def __init__(self, category, rules):
        self.category = category
        self.rules = rules

    def _is_match_mut(self, x, y):
        df_rul = pd.read_excel(self.rules, sheet_name="CIViC_Variant_Matching")
        x_is_null = x.isnull()

        start_ok = not x_is_null["start"]
        stop_ok = not x_is_null["stop"]
        alleles_ok = not (x_is_null["reference_bases"] & x_is_null["variant_bases"])
        hgvsp_ok = not x_is_null["HGVSp"]
        hgvsc_ok = not x_is_null["HGVSc"]
        hgvsp_partial = x["HGVS_Type"]=="partial"
        hgvsp_complete = x["HGVS_Type"]=="complete"
        variant_specific = x["variant"].lower() in [x.lower() for x in df_rul["Variant"].values]

        # extract values required to assess the match
        x_var = x["variant"]
        x_start = x["start"]
        y_start = y["Start_Position"]
        x_stop = x["stop"]
        y_stop = y["End_Position"]
        x_ref = x["reference_bases"]
        y_ref = y["Reference_Allele"]
        x_alt = x["reference_bases"]
        y_alt = y["Reference_Allele"]
        x_mid = x["Mut_Identifier"]
        y_mid = y["Mut_Identifier"]
        x_hgvsc = x["HGVSc"]
        y_hgvsc = y["HGVSc"]
        x_hgvsp = x["HGVSp"]
        y_hgvsp = y["HGVSp"]
        y_vc = y["Variant_Classification"]
        y_hgvsp_all = y["all_effects"]

        # possible matches are 
        #  - A -> if all of start, stop, reference_bases and tumor_bases are present. Build an identifier from these 4
        #    columns and match on this identifier. Be aware of INDEL where one of start or stop will be empty.
        #  - B -> start and stop are present but none of reference_bases and tumor_bases. Any match is conditional
        #    to the candidate having start AND stop within the start and stop of the database.
        #    * B.1 the variant is in the list of specific variants. Use the corresponding values of 
        #       Variant_Classification to match.
        #    * B.2 the variant has partial HGVS protein format p.[AA][XX]. Then any alteration where at least one of all 
        #       the transcripts has a protein consequence starting with this pattern will be a match.
        #    * B.3 the variant has complete HGVS protein format
        #      -- substitution p.[AA][XX][AA]
        #      -- deletion p.[AA][XX](_[AA][XX])del
        #      -- duplication p.[AA][XX](_[AA][XX])dup
        #      -- insertion p.[AA][XX](_[AA][XX])ins[AA]
        #      -- deletion-insertion p.[AA][XX](_[AA][XX])delins[AA]
        #      -- frameshift p.[AA][XX][AA]fs[*][XX]
        #       Then any alteration where at least one of all the transcripts has a protein consequence containing the
        #       pattern will be a match.
        #    * B.4 the variant has complete HGVS dna format
        #      -- substitution p.[XX][N]>[N]
        #      -- deletion p.[XX](_[XX])del
        #      -- duplication p.[XX](_[XX])dup
        #      -- insertion p.[XX](_[XX])ins[N]
        #      -- invertion p.[XX](_[XX])inv[N]
        #      -- deletion-insertion p.[XX](_[XX])delins[N]
        #  - C -> start and stop are not present
        #    * C.2 = B.2
        #    * C.3 = B.3
        #    * C.4 = B.4
        if start_ok and stop_ok and alleles_ok:
            return (x_mid==y_mid, "A")
        elif start_ok and stop_ok:
            # TODO: is and good? would or be better or lead to false-matches?
            if y_start >= x_start and y_stop <= x_stop:
                if variant_specific:
                    vcs = df_rul.loc[df_rul["Variant"].str.lower()==x_var.lower(),"Variant_Classifications"].iloc[0]
                    vcs = vcs.split(",")
                    return (y_vc in vcs, "B.1")
                elif hgvsp_partial:
                    return (x_hgvsp in y_hgvsp_all, "B.2")
                elif hgvsp_complete:
                    return (x_hgvsp in y_hgvsp_all, "B.3")
                elif hgvsc_ok:
                    return (x_hgvsc==y_hgpvc, "B.4")
        else:
            if hgvsp_partial:
                return (x_hgvsp in y_hgvsp_all, "C.2")
            elif hgvsp_complete:
                return (x_hgvsp in y_hgvsp_all, "C.3")
            elif hgvsc_ok:
                return (x_hgvsc==y_hgpvc, "C.4")

        return (False, "F")


    def _is_match_fus(self, x, y):
        # extract values required to assess the match
        x_gene_1 = x["gene_1"]
        x_gene_2 = x["gene_2"]
        y_gene_1 = y["Gene_1"]
        y_gene_2 = y["Gene_2"]

        # possible matches are 
        #  - A -> gene 1 matches exactly and gene 2 matches exactly
        #  - B -> gene_1 matches exactly and gene 2 matches by Any or the other way around
        if x_gene_1 == y_gene_1 and x_gene_2 == y_gene_2:
            return (True, "A")
        else:
            if x_gene_1 == "Any":
                return (y_gene_2==x_gene_2 or y_gene_1==x_gene_2, "B")
            elif x_gene_2 == "Any":
                return (y_gene_2==x_gene_1 or y_gene_1==x_gene_1, "B")

        return (False, "F")


    def _is_match_cna(self, x, y):
        # extract values required to assess the match
        x_gene = x["gene"]
        y_gene = y["Hugo_Symbol"]
        x_alt = x["variant"]
        y_alt = y["Alteration"].upper()

        # possible matches are 
        #  - A -> gene matches exactly and variant matches exactly
        #  - B -> gene matches exactly and variant matches by transitivity
        #     * DELETION matches to LOSS and DELETION
        #     * AMPLIFICATION matches to AMPLIFICATION
        #     * LOH matches to LOH and CN LOH
        if x_gene == y_gene and x_alt == y_alt:
            return (True, "A")
        else:
            if x_gene == y_gene:
                if y_alt == "DELETION" and x_alt in ["DELETION", "LOSS"]:
                    return (True, "B")
                elif y_alt == "LOH" and x_alt in ["LOH", "CN LOH"]:
                    return (True, "B")

        return (False, "F")


    def _add_annotations(self, x, y, match_type):
        # add matching type 
        col = "CIViC_Matching_Type"
        if col not in y.index:
            y[col] = match_type
        else:
            y[col] = y[col] + ";" + match_type

        # add identity of matching disease
        col = "CIViC_Matching_Disease"
        if col not in y.index:
            y[col] = x["disease"]
        else:
            y[col] = y[col] + ";" + x["disease"]

        # add identity of matching variant
        col = "CIViC_Matching_Variant"

        if col not in y.index:
            y[col] = x["variant"]
        else:
            y[col] = y[col] + ";" + x["variant"]

        # add type of annotation (Predictive, Diagnostic, Prognostic, Oncogenic)
        # and direction of annotation (Positive or Negative)
        if x["evidence_type"]!="Oncogenic":
            is_positive = x["clinical_significance"] in ["Sensitivity/Response", "Positive", "Better Outcome"]
            is_negative = x["clinical_significance"] in ["Resistance", "Poor Outcome", "Reduced Sensitivity",
                                                         "Adverse Response", "Negative", "Dominant Negative",
                                                         "Pathogenic", "Likely Pathogenic"]
            if is_positive:
                col = "%s:P:%s" % (x["evidence_type"], x["evidence_level"])
            elif is_negative:
                col = "%s:N:%s" % (x["evidence_type"], x["evidence_level"])
            else:
                raise ValueError("-warning! could not determine evidence direction for %s with matching variant %s" %
                      (y["Mut_Identifier"], x["variant"]))
                # print("-warning! could not determine evidence direction for %s with matching variant %s" %
                #       (y["Mut_Identifier"], x["variant"]))
                return y
        else:
            col = "%s:%s" % (x["evidence_type"], x["evidence_level"])

        if x["evidence_type"]=="Predictive":
            if x["drug_interaction_type"]=="Combination":
                val = x["drugs"].replace(",","+")
            elif x["drug_interaction_type"]=="Substitutes":
                val = x["drugs"].replace(",","|")
            elif x["drug_interaction_type"]=="Sequential":
                val = x["drugs"].replace(",","|")
            else:
                val = x["drugs"]
        else:
            val = "Y"

        if col not in y.index:
            y[col] = val
        else:
            y[col] = y[col] + ";" + val

        return y


    def _annotate_mut(self, df_civ, df_alt):
        x_genes = set(df_civ["gene"].dropna())
        y_genes = set(df_alt["Hugo_Symbol"].dropna())

        if len(x_genes.intersection(y_genes))==0:
            print("-CIViC database and input mutations table have no gene in common, no annotation.")
            return df_alt.copy()
        else:
            ys = []
            for _, y in df_alt.iterrows():
                y_gene = y["Hugo_Symbol"]
                if y_gene in x_genes:
                    for _, x in df_civ.loc[df_civ["gene"]==y_gene].iterrows():
                        is_match, match_type = self._is_match_mut(x, y)
                        if is_match:
                            y = self._add_annotations(x, y, match_type)
                ys.append(y.to_frame().T)

            return pd.concat(ys, axis=0)


    def _annotate_fus(self, df_civ, df_alt):
        x_genes_1 = set(df_civ["gene_1"].dropna())
        x_genes_2 = set(df_civ["gene_2"].dropna())
        y_genes_1 = set(df_alt["Gene_1"].dropna())
        y_genes_2 = set(df_alt["Gene_2"].dropna())

        if len(x_genes_1.intersection(y_genes_1))==0 & len(x_genes_2.intersection(y_genes_2))==0:
            print("-CIViC database and input fusions table have no gene in common, no annotation.")
            return df_alt.copy()
        else:
            ys = []
            for _, y in df_alt.iterrows():
                y_gene_1 = y["Gene_1"]
                y_gene_2 = y["Gene_2"]
                if y_gene_1 in x_genes_1 or y_gene_2 in x_genes_2:
                    for _, x in df_civ.loc[(df_civ["gene_1"]==y_gene_1) | (df_civ["gene_2"]==y_gene_2)].iterrows():
                        is_match, match_type = self._is_match_fus(x, y)
                        if is_match:
                            y = self._add_annotations(x, y, match_type)
                ys.append(y.to_frame().T)

            return pd.concat(ys, axis=0)


    def _annotate_cna(self, df_civ, df_alt):
        x_genes = set(df_civ["gene"].dropna())
        y_genes = set(df_alt["Hugo_Symbol"].dropna())

        if len(x_genes.intersection(y_genes))==0:
            print("-CIViC database and input cnas table have no gene in common, no annotation.")
            return df_alt.copy()
        else:
            ys = []
            for _, y in df_alt.iterrows():
                y_gene = y["Hugo_Symbol"]
                if y_gene in x_genes:
                    for _, x in df_civ.loc[df_civ["gene"]==y_gene].iterrows():
                        is_match, match_type = self._is_match_cna(x, y)
                        if is_match:
                            y = self._add_annotations(x, y, match_type)
                ys.append(y.to_frame().T)

            return pd.concat(ys, axis=0)


    def run(self, df_civ, df_alt):
        if self.category=="mut":
            df_ann = self._annotate_mut(df_civ, df_alt)
        elif self.category=="fus":
            df_ann = self._annotate_fus(df_civ, df_alt)
        elif self.category=="cna":
            df_ann = self._annotate_cna(df_civ, df_alt)

        # print info about numbers of variants annotated
        if "CIViC_Matching_Variant" in df_ann.columns:
            mask = ~df_ann["CIViC_Matching_Variant"].isnull()
        else:
            mask = pd.Series(False, index=df_ann.index)
        print("-INFO: %s/%s lines from the input table of alterations were matched in CIViC" % (sum(mask), len(mask)))

        cols_alt = list(df_alt.columns)
        cols_ann = list(df_ann.columns)
        cols_new = list(set(cols_ann).difference(set(cols_alt)))
        return df_ann[cols_alt + sorted(cols_new)]


def main(args):
    # load input alterations table
    df_alt = pd.read_table(args.input)

    if df_alt.shape[0]==0:
        df_ann = df_alt.copy()
    else:
        # build identifier
        df_alt = build_alt_identifier(df_alt, args.category)

        # transform input tumor types to list
        args.tumor_types = args.tumor_types.split("|")

        # load and process CIViC table
        df_civ = pd.read_excel(args.civic)

        # apply a series of filters on CIViC table to select only relevant lines
        preprocessor = CivicPreprocessor(tumor_types=args.tumor_types, category=args.category, rules=args.rules)
        df_civ = preprocessor.run(df_civ)

        # perform the annotation
        annotator = CivicAnnotator(category=args.category, rules=args.rules)
        df_ann = annotator.run(df_civ, df_alt)
        df_ann = clean_alt_identifier(df_ann, args.category)

    # save
    df_ann.to_csv(args.output, sep="\t", index=False)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add Tumor_Sample and Normal_Sample fields.")
    parser.add_argument('--input', type=str, help='Path to input table of variants.',
                        default="examples/data/example_cna.tsv")
    parser.add_argument('--civic', type=str, help='Path to CIViC database of clinical evidence summaries.',
                        default="data/01-Jan-2022-ClinicalEvidenceSummaries_Annotated.xlsx")
    parser.add_argument('--rules', type=str, help='Path to table of rules for cleaning the database and matching.',
                        default="data/CIViC_Curation_And_Rules_Mutation.xlsx")
    parser.add_argument('--category', type=str, help='Choose one of cna, mut or fus.', default='cna')
    parser.add_argument('--tumor_types', type=str, help='Tumor type designations in CIViC, separated with |.',
            default='Lung Cancer|Lung Carcinoma|Lung Non-small Cell Carcinoma|Lung Adenocarcinoma|Solid Tumor|Cancer')
    parser.add_argument('--output', type=str, help='Path to output table of variants.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
