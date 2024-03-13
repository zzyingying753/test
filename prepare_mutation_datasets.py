import requests
from bs4 import BeautifulSoup
import logging
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

def retrieve_disease_name(row):
    if type(row["NCI_CODE"]) == str:
        response = requests.get("https://ncit.nci.nih.gov/ncitbrowser/ConceptReport.jsp?dictionary=NCI%20Thesaurus&ns=null&code=" + row["NCI_CODE"])
        # for EFO terms, use http://www.ebi.ac.uk/efo/EFO_0000222
        response.raise_for_status()  # Raise an error for bad responses
        soup = BeautifulSoup(response.text, 'html.parser')
        disease_name_tag = soup.find('td', class_='texttitle-blue')
#         disease_name_tag = soup.find('title', id="pageTitle") # EFO terms
        if disease_name_tag:
            return disease_name_tag.get_text().split(" (")[0]
        else:
            return "Unknown"
    else:
        return " ".join([row["PRIMARY_SITE"], row["HISTOLOGY_SUBTYPE_1"], row["PRIMARY_HISTOLOGY"]])
    
    
def get_tcga(disease_name, tcga2keyword = {"LAML": {"leukemia", "myelo"}, "ACC": "adrenocortical", "BLCA": "bladder", "LGG": "brain", "GBM": {"glio", "astrocyt"}, "BRCA": "breast", "CESC": "cervical", "CHOL": "cholangio", "COAD": "colon", "ESCA": "esophag", "KIRC": "kidney", "KIRP": "kidney", "KICH": "kidney", "LIHC": "liver",  "LUSC": ["lung", "squam"], "LUAD": "lung", "DLBC": "lymphoma", "MESO": "mesothelioma", "OV": "ovar", "PAAD": "pancrea", "PCPG": "pheoch",  "PRAD": "prostate", "READ": "rect", "UVM": ["uveal", "melanoma"], "SKCM": {"melanoma", "skin"}, "STAD": "stomach", "TGCT": "testi", "THYM": "thym", "THCA": "thyro", "HNSC": "head and neck", "UCS": "uteri", "UCEC": "uteri",  "SARC": "sarcoma"}):
    disease_name = disease_name.lower().replace("renal", "kidney").replace("bile", "cholangio").replace("gastric", "stomach").replace("hepato", "liver")
    for tcga, keyword in tcga2keyword.items():
        if type(keyword) == str:
            if disease_name.find(keyword) != -1:
                return "TCGA-" + tcga
        elif type(keyword) == list:
            count = 0
            for term in keyword:
                if disease_name.find(term) != -1:
                    count = count + 1
            if count == len(keyword):
                return "TCGA-" + tcga
        else:
            for term in keyword:
                if disease_name.find(term) != -1:
                    return "TCGA-" + tcga
    return "Unknown"

def fix_tsb(tsb):
    if len(tsb) != len("TCGA-02-0003-01"):
        return tsb + "-01"
    else:
        return tsb

def get_tumor_type(other, other2tcga, tcga_cancers):
    if other in other2tcga:
        return "TCGA-" + other2tcga[other].upper()
    elif other in tcga_cancers:
        return "TCGA-" + other.upper()
    else:
        return other

def get_sample_types(tsb, types = ["01", "03", "09"]):
    if tsb.split("-")[3][0:2] in types:
        return True
    else:
        return False

#---------------------------------------------------------------------------------------------------------------------
# COSMIC

# Get primary tumors
df_sample = pd.read_csv(dirname + "/mutations/COSMIC/Cosmic_Sample_v98_GRCh37.tsv.gz", sep = "\t")
df_sample = df_sample[(df_sample["TUMOUR_SOURCE"] == "primary") & ~df_sample["SAMPLE_TYPE"].isin({"cell-line", "xenograft", "organoid culture", "short-term culture"})]
df_sample = df_sample[["COSMIC_SAMPLE_ID", "TUMOUR_ID"]].drop_duplicates()

# Get mutations
df1 = pd.read_csv(dirname + "/mutations/COSMIC/Cosmic_GenomeScreensMutant_v98_GRCh37.tsv.gz", sep = "\t", dtype = str)
df1 = df1[df1["SAMPLE_NAME"].apply(lambda x: x[0:4] != "TCGA")]
df1 = pd.merge(df1, df_sample, on = "COSMIC_SAMPLE_ID")
df1 = df1[["TUMOUR_ID", "COSMIC_PHENOTYPE_ID", "CHROMOSOME", "GENOME_START", "COSMIC_STUDY_ID", "GENOMIC_WT_ALLELE", "GENOMIC_MUT_ALLELE"]].rename(columns = {"TUMOUR_ID": "Tumor_Sample_Barcode", "CHROMOSOME": "Chromosome", "GENOME_START": "Start_Position", "GENOMIC_WT_ALLELE": "Reference_Allele", "GENOMIC_MUT_ALLELE": "Tumor_Seq_Allele2"})
df1 = df1.dropna(subset = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode"]).drop_duplicates()
df1["vep"] = df1.apply(lambda x: " ".join([x["Chromosome"], x["Start_Position"], str(int(x["Start_Position"]) + len(x["Reference_Allele"]) - 1), x["Reference_Allele"] + "/" + x["Tumor_Seq_Allele2"]]), axis = 1)
df1["Start_Position"] = df1["Start_Position"].astype(int)
df1 = df1.sort_values(by = ["Chromosome", "Start_Position"])
df1[["vep"]].drop_duplicates().to_csv(dirname + "/mutations/COSMIC/VEP_input.txt", sep = "\t", header = None, index = None)

# Annotate tumor type
df2 = pd.read_csv(dirname + "/mutations/COSMIC/Cosmic_Classification_v98_GRCh37.tsv.gz", sep = "\t", dtype = str)
df2 = df2[df2["COSMIC_PHENOTYPE_ID"].isin(set(df1["COSMIC_PHENOTYPE_ID"]))]
df2["disease_name"] = df2.apply(retrieve_disease_name, axis = 1)
df2["TUMORTYPE"] = df2["disease_name"].apply(get_tcga)
df_merge = pd.merge(df1.drop(columns = ["vep"]), df2[["COSMIC_PHENOTYPE_ID", "TUMORTYPE"]].drop_duplicates(), on = "COSMIC_PHENOTYPE_ID")
df_merge.to_csv(dirname + "/mutations/COSMIC/GenomeScreensMutant.v98.maf", sep = "\t", header = True, index = None)

VEP_output = {
  "all_aa": "http://grch37.ensembl.org/Homo_sapiens/Download/Tools/VEP\?format\=vcf\;tl\=4xvI4afG5a3Lct7I\-9543392",
  "all_ab": "http://grch37.ensembl.org/Homo_sapiens/Download/Tools/VEP\?format\=vcf\;tl\=O2lMA0rtubxAFBqe\-9543399",
  "all_ac": "http://grch37.ensembl.org/Homo_sapiens/Download/Tools/VEP\?format\=vcf\;tl\=VvY18MEaYtHMK8A3\-9543401",
  "all_ad": "http://grch37.ensembl.org/Homo_sapiens/Download/Tools/VEP\?format\=vcf\;tl\=etLSM2EbMPQmzAjL\-9543404",
  "single_aa": "http://grch37.ensembl.org/Homo_sapiens/Download/Tools/VEP\?format\=vcf\;tl\=Ex4XjFzwheQvya4t\-9543390",
  "single_ab": "http://grch37.ensembl.org/Homo_sapiens/Download/Tools/VEP\?format\=vcf\;tl\=apj1KfgEZ89mZqB1\-9543393",
  "singla_ac": "http://grch37.ensembl.org/Homo_sapiens/Download/Tools/VEP\?format\=vcf\;tl\=9ZKK4AzQlroLoVBI\-9543400",
  "single_ad": "http://grch37.ensembl.org/Homo_sapiens/Download/Tools/VEP\?format\=vcf\;tl\=oAiIpzDUzwv6z1iW\-9543403"}

for x in VEP_output:
  os.system("wget -c -O " + dirname + "/mutations/COSMIC/VEP_output/" + x + ".vcf " + VEP_output[x])

for mapping in ["single", "all"]:
    count = 0
    for file in os.listdir(dirname + "/mutations/COSMIC/VEP_output"):
        if file.split("_")[0] == mapping:
            count = count + 1
            if count == 1:
                mode = "w"
            else:
                mode = "a+"
            all_function_py3.parse_VEP_output(dirname + "/mutations/COSMIC/VEP_output/" + file, \
                dirname + "/mutations/COSMIC/parsed_" + mapping + ".txt", mode = mode, cols = ["SWISSPROT", "Consequence", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "HGVSc", "HGVSp", "Protein_position", "Amino_acids", "Codons", "ENSP", "SIFT", "PolyPhen", "gnomADe_AF", "CADD_PHRED", "CADD_RAW", "MPC"])




#---------------------------------------------------------------------------------------------------------------------
# TCGA and non-TCGA datasets from Chang et al. (https://www.nature.com/articles/nbt.3391)

# load old and new datasets
df1 = pd.read_csv(dirname + "/mutations/pancan_unfiltered.maf", sep = "\t", dtype = str)
df1["Tumor_Sample_Barcode"] = df1["Tumor_Sample_Barcode"].apply(fix_tsb)
df2 = pd.read_csv(dirname + "/mutations/mc3.v0.2.8.PUBLIC.maf.gz", sep = "\t", dtype = str).dropna(subset = ["Tumor_Seq_Allele2"])
df2["vep"] = df2.apply(lambda x: " ".join([x["Chromosome"], x["Start_Position"], x["End_Position"], x["Reference_Allele"] + "/" + x["Tumor_Seq_Allele2"]]), axis = 1)
df2["Start_Position"] = df2["Start_Position"].astype(int)
df2.sort_values(by = ["Chromosome", "Start_Position"])[["vep"]].drop_duplicates().to_csv(dirname + \
  "/mutations/VEP_input_MC3.txt", sep = "\t", header = None, index = None)

# load all samples with somatic mutation data in TCGA
with open(dirname + "/mutations/biospecimen.cases_selection.2023-10-03.json") as f:
    case = json.load(f)
tcga_samples = set()
for x in case:
    for y in x["samples"]:
        if y["submitter_id"][0:4] == "TCGA":
            tcga_samples = tcga_samples.union({y["submitter_id"][0:len("TCGA-02-0003-01")]})

# load case information
with open(dirname + "/mutations/cases.2023-07-30.json") as f:
    data = json.load(f)
    case2project = defaultdict(lambda: np.nan)
    for x in data:
        case2project[x["submitter_id"]] = x["project"]["project_id"]

# split TCGA and non-TCGA samples
new_samples = {x[0:len("TCGA-02-0003-01")] for x in df2["Tumor_Sample_Barcode"]}
old_samples = set()
non_tcga_samples = set()
for x in df1["Tumor_Sample_Barcode"].unique():
    if x[0:4] == "TCGA":
        old_samples = old_samples.union({x})
    else:
        non_tcga_samples = non_tcga_samples.union({x})
added_samples = old_samples.intersection(tcga_samples) - new_samples
df_vep = df1[df1["Tumor_Sample_Barcode"].isin(added_samples.union(non_tcga_samples))]
df_vep["vep_input"] = df_vep.apply(lambda x: x["Chromosome"] + " " + x["Start_Position"] + " " + x["End_Position"] + " " + x["Reference_Allele"] + "/" + x["Tumor_Seq_Allele2"], axis = 1)
df_vep["Start_Position"] = df_vep["Start_Position"].astype(int)
df_vep = df_vep.sort_values(by = ["Chromosome", "Start_Position"])
df_vep[["vep_input"]].drop_duplicates().to_csv(dirname + "/mutations/VEP_input_Chang.txt", sep = "\t", header = None, index = None)
cols = ["NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Strand", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode"]

# retain primary tumors in TCGA
df_tcga = pd.concat([df1[df1["Tumor_Sample_Barcode"].isin(added_samples)][cols], df2[cols]]).drop_duplicates()
df_tcga = df_tcga[df_tcga["Tumor_Sample_Barcode"].apply(get_sample_types)]

# annotate tumor types
df_tcga["TUMORTYPE"] = df_tcga["Tumor_Sample_Barcode"].apply(lambda x: case2project["-".join(x.split("-")[0:3])])
df_tcga.dropna(subset = ["TUMORTYPE"]).to_csv(dirname + "/mutations/TCGA.maf", sep = "\t", header = True, index = None)
df_non_tcga = df1[df1["Tumor_Sample_Barcode"].isin(non_tcga_samples)][cols]
other2tcga = {'hgg': 'gbm', 'pias': 'lgg', 'lusm': 'lusc', 'mcl': 'dlbc', 'all': 'dlbc', 'lymbc': 'dlbc', 'cll': 'dlbc', 'cscc': 'skcm', 'coadread': 'coad', 'mmyl': 'laml', 'mds': 'laml', 'acyc': 'hnsc', 'npc': 'hnsc', 'gbc': 'chol'}
tcga_cancers = {x.split(".")[0].split("-")[1].lower().split("_")[0] for x in os.listdir(dirname + "/transcriptome_profiles/parsed_files")}
df_non_tcga["TUMORTYPE"] = df_non_tcga["TUMORTYPE"].apply(lambda x: get_tumor_type(x, other2tcga, tcga_cancers))
df_non_tcga.to_csv(dirname + "/mutations/non_TCGA.maf", sep = "\t", header = True, index = None)

VEP_output = {
"all_chang": "http://grch37.ensembl.org/Homo_sapiens/Download/Tools/VEP\?format\=vcf\;tl\=CqH9IoEytqzMOuP8\-9542574",
"single_chang": "http://grch37.ensembl.org/Homo_sapiens/Download/Tools/VEP\?format\=vcf\;tl\=7As86uz9147Zhz1J\-9542572",
"all_aa": "http://grch37.ensembl.org/Homo_sapiens/Download/Tools/VEP\?format\=vcf\;tl\=0sxrT52oxIvzKtu0\-9375321",
"all_ab": "http://grch37.ensembl.org/Homo_sapiens/Download/Tools/VEP\?format\=vcf\;tl\=kVKjTKdU55TaXrnh\-9375323",
"single_aa": "http://grch37.ensembl.org/Homo_sapiens/Download/Tools/VEP\?format\=vcf\;tl\=NGAWJneRZkoN9ywu\-9375325",
"single_ab": "http://grch37.ensembl.org/Homo_sapiens/Download/Tools/VEP\?format\=vcf\;tl\=kaJIT7Zj06Q34XL3\-9375330"
}

for x in VEP_output:
  os.system("wget -c -O " + dirname + "/mutations/" + x + ".vcf " + VEP_output[x])

for mapping in ["single", "all"]:
  count = 0
  for file in os.listdir(dirname + "/mutations/VEP_output"):
      if file.split("_")[0] == mapping:
          count = count + 1
          if count == 1:
              mode = "w"
          else:
              mode = "a+"
          all_function_py3.parse_VEP_output(dirname + "/mutations/VEP_output/" + file, \
              dirname + "/mutations/parsed_" + mapping + ".txt", mode = mode, cols = ["SWISSPROT", "Consequence", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "HGVSc", "HGVSp", "Protein_position", "Amino_acids", "Codons", "ENSP", "SIFT", "PolyPhen", "gnomADe_AF", "CADD_PHRED", "CADD_RAW", "MPC"])
  df = pd.read_csv(dirname + "/mutations/parsed_" + mapping + ".txt", sep = "\t", dtype = str).drop_duplicates()
  df.to_csv(dirname + "/mutations/parsed_" + mapping + ".txt", sep = "\t", header = True, index = None)
