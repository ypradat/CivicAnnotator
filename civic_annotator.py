# -*- coding: utf-8 -*-
"""
@created: Jan 14 2022
@modified: Aug 09 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Class for annotating a set of alterations using the preprocessed CIViC database (variant filtered and reformatted).
"""

import numpy as np
import pandas as pd

def convert_num_to_str(x):
    try:
        y = "%d" % int(x)
    except:
        try:
            y = "%f" % float(x)
            if y=="nan":
                y = x
        except:
            y = x
    return y


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
        if start_ok and stop_ok and alleles_ok and x_mid==y_mid:
            return (True, "A")

        if start_ok and stop_ok:
            # TODO: is "and" good? would "or" be better or lead to false-matches?
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

        if hgvsp_partial and x_hgvsp is str:
            return (x_hgvsp in y_hgvsp_all, "C.2")
        elif hgvsp_complete and x_hgvsp is str:
            return (x_hgvsp in y_hgvsp_all, "C.3")
        elif hgvsc_ok:
            return (x_hgvsc==y_hgvsc, "C.4")

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

        # add identity of matching gene variant
        col = "CIViC_Matching_Gene_Variant"
        if col not in y.index:
            y[col] = x["gene"] + " " + x["variant"]
        else:
            y[col] = y[col] + ";" + x["gene"] + " " + x["variant"]

        # add identity of evidence
        col = "CIViC_Matching_Evidence_Id"
        if col not in y.index:
            y[col] = str(x["evidence_id"])
        else:
            y[col] = y[col] + ";" + str(x["evidence_id"])

        # add citation
        col = "CIViC_Matching_Citation"
        if col not in y.index:
            y[col] = x["citation"]
        else:
            y[col] = y[col] + ";" + x["citation"]

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


    def _annotate_alt(self, df_civ, df_alt):
        x_genes = set(df_civ["gene"].dropna())
        x_genes_1 = set(df_civ["gene_1"].dropna())
        x_genes_2 = set(df_civ["gene_2"].dropna())
        y_genes = set(df_alt["Hugo_Symbol"].dropna())
        y_genes_1 = set(df_alt["Gene_1"].dropna())
        y_genes_2 = set(df_alt["Gene_2"].dropna())

        old2news = {}
        cols_keeps =  {}
        for cat, nam in zip(["m", "c", "f"], ["Mut", "Cna", "Fus"]):
            old2new = {"Tumor_Sample_Barcode": "Tumor_Sample_Barcode_%s" % nam,
                       "Sample_Id": "Sample_Id_%s" % nam,
                       "Matched_Norm_Sample_Barcode": "Matched_Norm_Sample_Barcode_%s" % nam,
                       "CIViC_Matching_Type": "CIViC_Matching_Type_%s" % nam,
                       "CIViC_Matching_Gene_Variant": "CIViC_Matching_Gene_Variant_%s" % nam}
            old2news[cat] = old2new
            cols_keeps[cat] = list(old2new.values()) + ["%s_Identifier" % nam]

        cols_com = ["Subject_Id", "CIViC_Matching_Disease", "CIViC_Matching_Evidence_Id", "CIViC_Matching_Citation"]
        tags_civ = ["Predictive:", "Prognostic:", "Diagnostic:"]

        # check that there is a potential match by checking that there is 1 least gene in common
        if len(x_genes.intersection(y_genes))==0 & len(x_genes_1.intersection(y_genes_1))==0 & \
           len(x_genes_2.intersection(y_genes_2))==0:
            print("-CIViC database and input alts table have no gene in common, no annotation.")
            return pd.DataFrame({c: [] for c in cols_com})
        else:
            ys = [pd.DataFrame({c: [] for c in cols_com})]
            for evidence_id in df_civ["evidence_id"].unique():
                # select rows in civic db for evidence id
                df_civ_ei = df_civ.loc[df_civ["evidence_id"]==evidence_id].copy()
                x_genes = set(df_civ_ei["gene"].dropna())
                x_genes_1 = set(df_civ_ei["gene_1"].dropna())
                x_genes_2 = set(df_civ_ei["gene_2"].dropna())

                for subject_id in df_alt["Subject_Id"].unique():
                    # select rows in alterations table for subject id
                    df_alt_si = df_alt.loc[df_alt["Subject_Id"]==subject_id].copy()
                    y_genes = set(df_alt_si["Hugo_Symbol"].dropna())
                    y_genes_1 = set(df_alt_si["Gene_1"].dropna())
                    y_genes_2 = set(df_alt_si["Gene_2"].dropna())

                    # check that there is a potential match by checking that there is 1 least gene in common
                    if len(x_genes.intersection(y_genes))==0 & len(x_genes_1.intersection(y_genes_1))==0 & \
                       len(x_genes_2.intersection(y_genes_2))==0:
                        pass
                    else:
                        matches = {var: False for var in df_civ_ei["variant"].unique()}
                        match_cats = {var: None for var in df_civ_ei["variant"].unique()}
                        match_types = {var: None for var in df_civ_ei["variant"].unique()}
                        match_xs = {var: None for var in df_civ_ei["variant"].unique()}
                        match_ys = {var: None for var in df_civ_ei["variant"].unique()}

                        for _, x, in df_civ_ei.iterrows():
                            if x["category"] in ["m", "c"]:
                                mask_cat = df_alt_si["Alteration_Category"]==x["category"]
                                mask_gen = df_alt_si["Hugo_Symbol"]==x["gene"]
                            elif x["category"] == "f":
                                mask_cat = df_alt_si["Alteration_Category"]==x["category"]
                                if x["gene_1"]=="Any" or x["gene_2"]=="Any":
                                    f_genes = [x["gene_1"], x["gene_2"]]
                                    mask_gen = (df_alt_si["Gene_1"].isin(f_genes)) | (df_alt_si["Gene_2"].isin(f_genes))
                                else:
                                    mask_gen = (df_alt_si["Gene_1"]==x["gene_1"]) | (df_alt_si["Gene_2"]==x["gene_2"])
                            df_alt_si_gen = df_alt_si.loc[mask_cat & mask_gen]

                            for _, y in df_alt_si_gen.iterrows():
                                if x["category"]=="m":
                                    is_match, match_type = self._is_match_mut(x, y)
                                elif x["category"]=="c":
                                    is_match, match_type = self._is_match_cna(x, y)
                                elif x["category"]=="f":
                                    is_match, match_type = self._is_match_fus(x, y)
                                if is_match:
                                    matches[x["variant"]] = is_match
                                    match_cats[x["variant"]] = x["category"]
                                    match_types[x["variant"]] = match_type
                                    match_xs[x["variant"]] = x
                                    match_ys[x["variant"]] = self._add_annotations(x, y, match_type)

                        if all(matches.values()):
                            y_vars = []
                            for var, y in match_ys.items():
                                cat = match_cats[var]
                                old2new = old2news[cat]
                                cols_keep = cols_keeps[cat]
                                y = y.to_frame().T
                                y = y.rename(columns=old2new)

                                cols_civ = [c for c in y if any(c.startswith(t) for t in tags_civ)]
                                cols_keep = cols_com + cols_keep + cols_civ
                                cols_keep = [x for x in cols_keep if x in y]
                                y_vars.append(y[cols_keep].set_index(cols_com+cols_civ))

                            ys.append(pd.concat(y_vars, axis=1).reset_index(drop=False))

            # concatenate and collapse annotations from the same combination of (subject, samples, alterations)
            df_ann = pd.concat(ys, axis=0)

            if df_ann.shape[0] > 0:
                cols_agg = ["CIViC_Matching_Disease", "CIViC_Matching_Evidence_Id", "CIViC_Matching_Citation"]
                cols_agg += [x for x in df_ann if x.startswith("CIViC_Matching_Type")]
                cols_agg += [x for x in df_ann if x.startswith("CIViC_Matching_Gene_Variant")]
                cols_agg += sorted([x for x in df_ann if any(x.startswith(t) for t in tags_civ)])
                cols_gby = [x for x in df_ann if x not in cols_agg]
                df_ann[cols_gby] = df_ann[cols_gby].fillna("N/A")
                df_ann = df_ann.groupby(cols_gby).agg({c: lambda x: self._aggregate_str(x) for c in cols_agg})
                df_ann = df_ann.reset_index().replace("N/A", np.nan)

            return df_ann


    def _aggregate_str(self, x, sep=";"):
        x_nna = x.dropna().tolist()
        if len(x_nna)==0:
            return np.nan
        else:
            x_nna = [convert_num_to_str(e) for e in x_nna]
            return sep.join(x_nna)


    def run(self, df_civ, df_alt):
        if self.category=="mut":
            df_ann = self._annotate_mut(df_civ, df_alt)
        elif self.category=="fus":
            df_ann = self._annotate_fus(df_civ, df_alt)
        elif self.category=="cna":
            df_ann = self._annotate_cna(df_civ, df_alt)
        elif self.category=="alt":
            df_ann = self._annotate_alt(df_civ, df_alt)

        if self.category != "alt":
            # print info about numbers of variants annotated
            if "CIViC_Matching_Gene_Variant" in df_ann.columns:
                mask = ~df_ann["CIViC_Matching_Gene_Variant"].isnull()
            else:
                mask = pd.Series(False, index=df_ann.index)
            print("-INFO: %s/%s lines from the input table of alterations were matched in CIViC" % (sum(mask), len(mask)))

            cols_alt = list(df_alt.columns)
            cols_ann = list(df_ann.columns)
            cols_new = list(set(cols_ann).difference(set(cols_alt)))
            return df_ann[cols_alt + sorted(cols_new)]
        else:
            print("-INFO: %d combinations from the input table of alterations were matched in CIViC" % len(df_ann))
            return df_ann
