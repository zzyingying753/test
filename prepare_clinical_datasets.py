import pandas as pd
import os
import numpy as np
import pysam
import requests
import json
from collections import defaultdict
from multiprocessing import Pool
import sys
sys.path.append("/home/yz2296/bin")
import all_function_py3
import logging
import copy
from os import path
from glob import glob
import subprocess
import matplotlib.pyplot as plt
dirname = os.path.dirname(os.path.abspath(__file__))
sys.path.append(dirname + "/NetFlow3D")
import funcs
from Bio import SeqIO
import gzip
from residuecontact import get_alphafold_pdbresidue_to_uniprot_map, build_PDB_residues_connection_graph
from joblib import Parallel, delayed
import xml.etree.ElementTree as ET
import random
import networkx as nx
from lifelines import CoxPHFitter

def convert_vital_status(vital_status):
    if vital_status == "dead":
        return 1
    elif vital_status == "alive":
        return 0
    else:
        return np.nan
    
def convert_cursor_length(cursor_length):
    try:
        return float(cursor_length)
    except:
        return np.nan
    
def subnetwork_filter(uniprots, subnetworks):
    for x in subnetworks:
        if set(uniprots.split(",")) - set(x.split(",")) == set():
            return True
    return False

def get_tumor_status(tumor_status):
    if tumor_status == "WITH TUMOR":
        return 1
    elif tumor_status == "TUMOR FREE":
        return 0
    else:
        return np.nan


# df_subnetwork = pd.read_csv("/NFS4/storage/yz2296/cancer_hotspot_new/website/website/revision/NetFlow3D/output/TCGA_subnetworks.txt", sep = "\t")
# df_subnetwork = df_subnetwork[df_subnetwork["Subnetwork_size"] > 5]
# subnetwork_uniprots = set()
# for x in df_subnetwork["Subnetwork_UniProts"]:
#     subnetwork_uniprots = subnetwork_uniprots.union(set(x.split(",")))

file = "/local/storage/yz2296/cancer_hotspot_new/website/paper/survival/clinic_info.txt"
# file = "/local/storage/yz2296/cancer_hotspot_new/reference_data/clinic_info.txt"
df_clinic = pd.read_csv(file, sep = "\t", dtype = str).rename(columns = {"submitter_id": "bcr_patient_barcode"})
# df_clinic["OS"] = df_clinic["vital_status"].apply(convert_vital_status)
# df_clinic["OS.time"] = df_clinic["cursor_length"].apply(convert_cursor_length)

lof = {"Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site"}
inframe = {"In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"}
df_mut = pd.read_csv("/NFS4/storage/yz2296/cancer_hotspot_new/website/website/revision/mutations/TCGA_preprocessed_single.maf", sep = "\t", dtype = str)
df_lof = df_mut[df_mut["Variant_Classification"].isin(lof)]
df_inframe = df_mut[df_mut["Variant_Classification"].isin(inframe)]

pval_cutoff = 1e-3

df_intra1 = pd.read_csv("/NFS4/storage/yz2296/cancer_hotspot_new/website/website/revision/NetFlow3D/output/TCGA/PDB_intra_pvalue.txt", sep = "\t", dtype = str)
df_intra2 = pd.read_csv("/NFS4/storage/yz2296/cancer_hotspot_new/website/website/revision/NetFlow3D/output/TCGA/AlphaFold2_intra_pvalue_pLDDT0.txt", sep = "\t", dtype = str)
df_intra1 = df_intra1[df_intra1["Raw_pvalue"].astype(float) < pval_cutoff / df_intra2.shape[0]]
df_intra2 = df_intra2[df_intra2["Raw_pvalue"].astype(float) < pval_cutoff / df_intra2.shape[0]]
df_intra = pd.concat([df_intra1, df_intra2])

df_inter1 = pd.read_csv("/NFS4/storage/yz2296/cancer_hotspot_new/website/website/revision/NetFlow3D/output/TCGA/PDB_inter_pvalue.txt", sep = "\t", dtype = str)
df_inter2 = pd.read_csv("/NFS4/storage/yz2296/cancer_hotspot_new/website/website/revision/NetFlow3D/output/TCGA/PIONEER_inter_pvalue.txt", sep = "\t", dtype = str)
df_inter1 = df_inter1[df_inter1["Raw_pvalue"].astype(float) < pval_cutoff / df_inter1.shape[0]]
df_inter2 = df_inter2[df_inter2["Raw_pvalue"].astype(float) < pval_cutoff / df_inter2.shape[0]]
df_inter = pd.concat([df_inter1, df_inter2])

f = open("/NFS4/storage/yz2296/cancer_hotspot_new/website/website/revision/survival/cox_regression.txt", "w")
f.write("\t".join(["cancer", "HR", "coef", "se(coef)", "95CI_lower", "95CI_higher", "pval"]) + "\n")
for cancer in df_inframe["TUMORTYPE"].unique():
    patients = set(df_inframe[df_inframe["TUMORTYPE"] == cancer]["Tumor_Sample_Barcode"])
    specific_genes = {"AR", "CDH4", "EGFR", "EPHA3", "ERBB4", "FGFR2", "FLT3", "FOXA1", "FOXA2", "MECOM", "MIR142", "MSH4", "PDGFRA", "SOX1", "SOX9", "SOX17", "TBX3", "WT1"}
    expr_genes = pd.read_csv("/NFS4/storage/yz2296/cancer_hotspot_new/website/website/revision/transcriptome_profiles/parsed_files/" + cancer + "_tumor.txt", sep = "\t", header = None)[0].tolist()
    prioritized_uniprots = set(df_inframe[df_inframe["Hugo_Symbol"].isin(specific_genes)]["UniProt"].tolist() + df_inframe[df_inframe["Gene"].isin(expr_genes)]["UniProt"].tolist())

    # expression filter
    print("before expr filter", df_intra.shape[0], df_inter.shape[0])
    df_intra_tmp = df_intra[df_intra["Uniprots"].isin(prioritized_uniprots)]
    df_inter_tmp = df_inter[df_inter["Uniprots"].apply(lambda x: set(x.split(",")) - prioritized_uniprots == set())]
    print("after expr filter", df_intra_tmp.shape[0], df_inter_tmp.shape[0])
    
#     # subnetwork filter
#     df_intra_tmp = df_intra_tmp[df_intra_tmp["Uniprots"].isin(subnetwork_uniprots)]
#     df_inter_tmp = df_inter_tmp[df_inter_tmp["Uniprots"].apply(lambda x: subnetwork_filter(x, df_subnetwork["Subnetwork_UniProts"].tolist()))]
      
    uniprot2res = defaultdict(set)
    for residues in df_intra_tmp["Residues"].tolist() + df_inter_tmp["Residues"].tolist():
        for res in residues.split(","):
            uniprot, pos = res.split("_")
            uniprot2res[uniprot] = uniprot2res[uniprot].union({int(pos)})
    
    df_tmp = df_inframe[df_inframe["TUMORTYPE"] == cancer]
    case_patient = set(df_tmp[df_tmp.apply(lambda x: set(all_function_py3.split_consecutive_pos(x["Protein_position"])).intersection(uniprot2res[x["UniProt"]]) != set(), axis = 1)]["Tumor_Sample_Barcode"])
    print(cancer + "##################")
    print(len(case_patient), len(patients - case_patient))
    survival_type = "OS"
    df_clinic_tmp = pd.DataFrame({"bcr_patient_barcode": sorted(case_patient) + sorted(patients - case_patient), "label": ([1] * len(case_patient)) + ([0] * len(patients - case_patient))})
    df_clinic_tmp["bcr_patient_barcode"] = df_clinic_tmp["bcr_patient_barcode"].apply(lambda x: "-".join(x.split("-")[0:3]))
    df_clinic_tmp = pd.merge(df_clinic.dropna(subset = [survival_type, survival_type + ".time"]), df_clinic_tmp, on = "bcr_patient_barcode")
    if df_clinic_tmp.shape[0] > 0:
        df_clinic_tmp[survival_type] = df_clinic_tmp[survival_type].astype(int)
        df_clinic_tmp[survival_type + ".time"] = df_clinic_tmp[survival_type + ".time"].astype(float)     
        df_survival = df_clinic_tmp.copy()
        df_clinic_tmp["tumor_status"] = df_clinic_tmp["tumor_status"].apply(get_tumor_status)
        df_clinic_tmp = df_clinic_tmp.dropna(subset = ["tumor_status"])
        if df_clinic_tmp.shape[0] > 0:
            df_clinic_tmp["tumor_status"] = df_clinic_tmp["tumor_status"].astype(int)
            for col in ["tumor_status", "label", survival_type]:
                df_clinic_tmp[col] = df_clinic_tmp[col].astype("category")
            to_fit = pd.get_dummies(df_clinic_tmp[[survival_type, survival_type + ".time", "tumor_status", "label"]], drop_first = True)
            cph = CoxPHFitter()
            formula = " + ".join([x for x in to_fit.columns if x not in {survival_type + ".time", survival_type + "_1"}])
            try:
                cph.fit(to_fit, duration_col = survival_type + '.time', event_col = survival_type + '_1', formula = formula)
                lower_95ci = np.exp(cph.confidence_intervals_.loc["label_1", "95% lower-bound"])
                upper_95ci = np.exp(cph.confidence_intervals_.loc["label_1", "95% upper-bound"])
                print(cph.summary["exp(coef)"]["label_1"], np.exp((cph.confidence_intervals_.loc["label_1", "95% lower-bound"] + cph.confidence_intervals_.loc["label_1", "95% upper-bound"]) / 2))
                HR = cph.summary["exp(coef)"]["label_1"]
                f.write("\t".join([cancer, str(HR), str(cph.summary["coef"]["label_1"]), str(cph.summary["se(coef)"]["label_1"]), str(lower_95ci), str(upper_95ci), str(cph.summary["p"]["label_1"])]) + "\n")
                print(cph.summary["p"])
            except:
                logging.error(cancer)
                HR = None
        else:
            logging.warning(cancer)
            HR = None
        all_function_py3.simple_survival_curve(df_survival, "/NFS4/storage/yz2296/cancer_hotspot_new/website/website/revision/figures/" + cancer + "_" + survival_type + ".png", column_name="label", case_label=1, control_label=0, vital_label=survival_type, duration_label=survival_type + ".time", title=";".join([cancer, "case:" + str(df_survival[df_survival["label"] == 1].shape[0]), "control:" + str(df_survival[df_survival["label"] == 0].shape[0])]), \
                                               HR=HR, tick_bins=4, xlabel="Time (days)", ylabel=survival_type + " Survival", figsize = (5, 5), fontsize = 16, case_color = "red", control_color = "gray")
f.close()

# survival_type = "PFI"
# to_fit = pd.get_dummies(df_lm[[survival_type, survival_type + ".time", "tumor_status", "ajcc_pathologic_tumor_stage", "log_intra", "log_inter"]], drop_first = True)
# df_lm["histological_type"] = df_lm["histological_type"].astype("category")
# stages = sorted(set(df_lm["ajcc_pathologic_tumor_stage"].tolist()))
# df_lm["tumor_stage"] = df_lm["ajcc_pathologic_tumor_stage"].apply(lambda x: stages.index(x)).astype(float)
# df_lm["gender"] = df_lm["gender"].astype("category")
# df_lm["race"] = df_lm["race"].astype("category")
# df_lm["histological_grade"] = df_lm["histological_grade"].astype("category")
# df_lm["age_at_initial_pathologic_diagnosis"] = df_lm["age_at_initial_pathologic_diagnosis"].astype(float)


# from lifelines import CoxPHFitter
#     df_lm = pd.merge(df_final, df_clinic, on = ["bcr_patient_barcode"])
#     df_lm["tumor_status"] = df_lm["tumor_status"].astype("category")
#     df_lm["ajcc_pathologic_tumor_stage"] = df_lm["ajcc_pathologic_tumor_stage"].astype("category")
#     df_lm["PFI"] = df_lm['PFI'].astype("category")
#     df_lm["OS"] = df_lm["OS"].astype("category")
#     df_lm["DFI"] = df_lm["DFI"].astype("category")
#     df_lm["DSS"] = df_lm["DSS"].astype("category")

#     survival_type = "PFI"
#     to_fit = pd.get_dummies(df_lm[[survival_type, survival_type + ".time", "tumor_status", "ajcc_pathologic_tumor_stage", "log_intra", "log_inter"]], drop_first = True)


#     # df_lm["histological_type"] = df_lm["histological_type"].astype("category")
#     # stages = sorted(set(df_lm["ajcc_pathologic_tumor_stage"].tolist()))
#     # df_lm["tumor_stage"] = df_lm["ajcc_pathologic_tumor_stage"].apply(lambda x: stages.index(x)).astype(float)
#     # df_lm["gender"] = df_lm["gender"].astype("category")
#     # df_lm["race"] = df_lm["race"].astype("category")



#     # # df_lm["age_at_initial_pathologic_diagnosis"] = df_lm["age_at_initial_pathologic_diagnosis"].astype(float)
#     # df_lm["histological_grade"] = df_lm["histological_grade"].astype("category")

#     cph = CoxPHFitter()
#     formula = " + ".join([item for item in to_fit.columns if item not in {survival_type + ".time", survival_type + "_1.0"}]) + " + log_intra * log_inter"
#     cph.fit(to_fit.dropna(), duration_col = survival_type + '.time', event_col = survival_type + '_1.0', formula = formula)
#     print(cph.summary)

# cph = CoxPHFitter()
# formula = " + ".join([item for item in to_fit.columns if item not in {survival_type + ".time", survival_type + "_1.0"}])
# cph.fit(to_fit, duration_col = survival_type + '.time', event_col = survival_type + '_1.0', formula = formula)
# pval = cph.summary["p"]["NO_of_mutations"]
# HR = cph.summary["exp(coef)"]["NO_of_mutations"]
# coef = cph.summary["coef"]["NO_of_mutations"]
# se = cph.summary["se(coef)"]["NO_of_mutations"]


# #---------------------------------------------------------------------------------------------------------------------
# # Clinical data
# command = "/fs/cbsuhyfs1/storage1/yz2296/gdc-client download -d /fs/cbsuhyfs1/storage1/yz2296/TCGA_raw/clinical -m " + dirname + "/clinical_data/gdc_manifest.2023-08-19.txt"
# os.system(command)

# def extract_information(element, depth=0, path=""):
#	 """Recursively extract information from the XML element and its children."""
	
#	 # Extract the current element's tag, text, and attributes
#	 tag = element.tag
#	 text = element.text.strip() if element.text else ""
#	 attributes = element.attrib

#	 # Print or store the extracted details. 
#	 # The depth and path parameters help in showing the depth and the full path of the element in the tree.
#	 print(f"{'  ' * depth}Path: {path}/{tag}")
#	 print(f"{'  ' * depth}Tag: {tag}")
#	 print(f"{'  ' * depth}Text: {text}")
#	 print(f"{'  ' * depth}Attributes: {attributes}")
#	 print()

#	 # Recursively apply for all child elements
#	 for child in element:
#		 extract_information(child, depth + 1, f"{path}/{tag}")

# # Parse the XML document
# tree = ET.parse('/fs/cbsuhyfs1/storage1/yz2296/TCGA_raw/clinical/f88ddd98-a2f9-469c-a7d1-2a99ba9510e5/nationwidechildrens.org_clinical.TCGA-25-1630.xml')
# root = tree.getroot()

# # Extract information from the root and its children
# extract_information(root)
