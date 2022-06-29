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

Class for preprocessing the CIViC database (variant filtering and reformatting).
"""

import numpy as np
import pandas as pd
import re

from utils import explode_df, build_alt_identifier, build_row_identifier, expand_aa, shorten_aa

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
        allowed_categories = ["mut", "cna", "fus", "alt"]
        if self.category not in allowed_categories:
            raise ValueError("-the value of category %s is not recognized! Choose one of %s" % \
                             (category, ",".join(allowed_categories)))


    def _filter_tumor_types(self, df_civ):
        if self.tumor_types is not None and len(self.tumor_types)>0:
            mask = df_civ["disease"].isin(self.tumor_types)
            print("-INFO: selected %s/%s evidence lines corresponding to the following tumor types: \n\t%s" \
                  % (sum(mask), len(mask), "\n\t".join(self.tumor_types)))
            return df_civ.loc[mask].copy()
        else:
            return df_civ


    def _filter_category(self, df_civ):
        long2short= {"mut": ["m"], "cna": ["c"], "fus": ["f"], "alt": ["m", "f", "c", "o"]}
        if self.category != "alt":
            mask = df_civ["category"].isin(long2short[self.category])
        else:
            mask = ~df_civ["category"].isin(long2short[self.category])
        print("-INFO: selected %s/%s evidence lines corresponding to the %s category" % \
              (sum(mask), len(mask), self.category))
        return df_civ.loc[mask].copy()


    def _filter_using_rules(self, df_civ):
        df_sub = df_civ.copy()

        # exclude rows for which the column `col` contains any of the values in `vals`
        df_rul = pd.read_excel(self.rules, sheet_name="CIViC_Exclude_Values")
        for col in df_rul:
            vals = df_rul[col].dropna().unique()
            mask = df_sub[col].fillna("Is_N/A").isin(vals)
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

        # trim white spaces
        df_civ["variant"] = df_civ["variant"].str.strip()

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


    def _aggregate_values(self, x, sep=" and "):
        x_nna = x.fillna("")
        return sep.join(x_nna)


    def _reformat_cfm_variant(self, df_civ, cols_exp, cols_old, sep):
        df_civ_c = df_civ.loc[df_civ["category"]=="c"].copy()
        df_civ_f = df_civ.loc[df_civ["category"]=="f"].copy()
        df_civ_m = df_civ.loc[df_civ["category"]=="m"].copy()
        df_civ_m["index_original"] = df_civ_m.index
        df_civ_c = self._reformat_cna_variant(df_civ_c)
        df_civ_f = self._reformat_fus_variant(df_civ_f)
        df_civ_m = self._reformat_mut_variant(df_civ_m)

        # for mutations, add dedoubling number in case a variant like V600E/K was split into 2 lines
        # for V600E and V600K
        col_dn = "dedoubling_number"
        counts_seen = {i: 0 for i in df_civ_m["index_original"].unique()}
        df_civ_m[col_dn] = np.nan
        for i in df_civ_m.index:
            index_original_i = df_civ_m.loc[i, "index_original"]
            count_seen_i = counts_seen[index_original_i]
            df_civ_m.loc[i, col_dn] = count_seen_i
            counts_seen[index_original_i] += 1
        del df_civ_m["index_original"]

        # left merge in order to double any row in df_civ_c for df_civ_c from an evidence id with dedoubled rows
        df_civ_c = df_civ_c.merge(df_civ_m[["evidence_id", col_dn]].drop_duplicates(), how="left", on="evidence_id")
        df_civ_f = df_civ_f.merge(df_civ_m[["evidence_id", col_dn]].drop_duplicates(), how="left", on="evidence_id")

        df_civ = pd.concat((df_civ_c, df_civ_f, df_civ_m), axis=0)
        df_civ[col_dn] = df_civ[col_dn].fillna(0)

        cols_new = [x for x in df_civ.columns.tolist() if x not in cols_old + [col_dn]]
        dt_agg_old = {x: "first" for x in cols_old if x not in ["evidence_id"] + [cols_exp]}
        dt_agg_new = {x: lambda y: self._aggregate_values(y, sep=sep) for x in cols_new + cols_exp}
        dt_agg = {**dt_agg_old, **dt_agg_new}
        df_civ = df_civ.groupby(["evidence_id", col_dn]).agg(dt_agg).reset_index(drop=False)
        del df_civ[col_dn]

        return df_civ


    def _reformat_alt_variant(self, df_civ):
        df_civ["variant_reformatted"] = False

        # trim white spaces
        df_civ["variant"] = df_civ["variant"].str.strip()

        # remove multiple consecutive spaces
        df_civ["variant"] = df_civ["variant"].apply(lambda x: " ".join(x.split()))

        # if variant is prefixed by the gene symbol, remove it
        mask = df_civ["variant"].apply(lambda x: x.split(" ")[0] if type(x)==str else x)==df_civ["gene"]
        df_civ.loc[mask, "variant"] = df_civ.loc[mask, "variant"].apply(lambda x: " ".join(x.split(" ")[1:]))

        # reformat variant from f and m category to add "and" tags
        mask_f_and_m = df_civ["category"].isin(["f and m", "f and m and m"])
        df_civ.loc[mask_f_and_m, "variant"] = df_civ.loc[mask_f_and_m, "variant"].apply(lambda x: " and ".join(x.split()))

        # process combinations
        mask = df_civ["variant"].apply(lambda x: "and" in x)
        df_civ = df_civ.loc[mask]
        sep = " and "
        cols_exp = ["category", "variant"]
        cols_old = df_civ.columns.tolist()
        df_civ = explode_df(df_civ, cols=cols_exp, sep=sep)
        df_civ = self._reformat_cfm_variant(df_civ=df_civ, cols_exp=cols_exp, cols_old=cols_old, sep=sep)
        cols_new = [x for x in df_civ if x not in cols_old]

        # split on the "and" tag
        df_civ = explode_df(df_civ, cols=cols_exp + cols_new, sep=sep)
        df_civ = df_civ.replace("", np.nan)

        return df_civ.reset_index(drop=True)


    def _reformat_variant(self, df_civ):
        if self.category=="mut":
            return self._reformat_mut_variant(df_civ)
        elif self.category=="fus":
            return self._reformat_fus_variant(df_civ)
        elif self.category=="cna":
            return self._reformat_cna_variant(df_civ)
        elif self.category=="alt":
            return self._reformat_alt_variant(df_civ)


    def _add_variant_identifier(self, df_civ):
        if self.category=="mut":
            return build_alt_identifier(df_civ, self.category, col_gene="gene", col_start="start", col_end="stop",
                                        col_ref="reference_bases", col_alt="variant_bases")
        elif self.category=="fus":
            return build_alt_identifier(df_civ, self.category, col_gene_1="gene_1", col_gene_2="gene_2")
        elif self.category=="cna":
            return build_alt_identifier(df_civ, self.category, col_gene="gene", col_alt="variant")
        elif self.category=="alt":
            df_civ_c = df_civ.loc[df_civ["category"]=="c"].copy()
            df_civ_f = df_civ.loc[df_civ["category"]=="f"].copy()
            df_civ_m = df_civ.loc[df_civ["category"]=="m"].copy()
            df_civ_c =  build_alt_identifier(df_civ_c, "cna", col_gene="gene", col_alt="variant")
            df_civ_f =  build_alt_identifier(df_civ_f, "fus", col_gene_1="gene_1", col_gene_2="gene_2")
            df_civ_m =  build_alt_identifier(df_civ_m, "mut", col_gene="gene", col_start="start", col_end="stop",
                                             col_ref="reference_bases", col_alt="variant_bases")
            return pd.concat((df_civ_c, df_civ_f, df_civ_m), axis=0)


    def run(self, df_civ):
        self._check_args(df_civ)
        df_civ = self._filter_tumor_types(df_civ)
        df_civ = self._filter_category(df_civ)
        df_civ = self._filter_using_rules(df_civ)
        df_civ = self._reformat_variant(df_civ)
        df_civ = self._add_variant_identifier(df_civ)
        return df_civ
