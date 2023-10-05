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

def get_hugo_symbol(row, ensg2symbol):
	if row["Gene"] in ensg2symbol and ensg2symbol[row["Gene"]] != set():
		return list(ensg2symbol[row["Gene"]])[0]
	else:
		return row["SYMBOL"]

def mutation_preprocessing(maf_file, vep_file, mapping_flag, af_cutoff, output):
	'''
	mapping_flag: 1. single; 2. all.
	'''

	# definition
	lof = {"Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site"}
	inframe = {"In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"}
	specific_genes = {"AR", "CDH4", "EGFR", "EPHA3", "ERBB4", "FGFR2", "FLT3", "FOXA1", "FOXA2", "MECOM", "MIR142", "MSH4", "PDGFRA", "SOX1", "SOX9", "SOX17", "TBX3", "WT1"}

	# load raw mutation data
	df = pd.read_csv(maf_file, sep = "\t", dtype = str)
	df = df.dropna(subset = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"])
	df["#Uploaded_variation"] = df.apply(lambda x: x["Chromosome"] + "_" + str(x["Start_Position"]) + "_" + x["Reference_Allele"] + "/" + x["Tumor_Seq_Allele2"], axis = 1)

	# load VEP output
	df_vep = pd.read_csv(vep_file, sep = "\t", dtype = str).rename(columns = {"Feature": "Transcript_ID"})
	df_vep["Variant_Classification"] = df_vep["Consequence"].apply(all_function_py3.consequence2varclass_coding)
	ensg2symbol = all_function_py3.gene_ID_conversion(df_vep["Gene"].dropna().unique().tolist(), convert_from = "ensembl.gene", convert_to = "symbol", species = "human")
	df_vep["Hugo_Symbol"] = df_vep.apply(lambda x: get_hugo_symbol(x, ensg2symbol), axis = 1)
	df_map = pd.merge(df, df_vep, on = "#Uploaded_variation")

	# exclude unexpressed mutations
	tumor2gene = {}
	gene2cancers = defaultdict(set)
	for x in os.listdir(dirname + "/transcriptome_profiles/parsed_files/"):
		if x[-len("_tumor.txt"):] == "_tumor.txt":
			cancer = x.split("_")[0]
			tumor2gene[cancer] = pd.read_csv(dirname + "/transcriptome_profiles/parsed_files/" + x, sep = "\t", header = None)[0].tolist()
			for gene in tumor2gene[cancer]:
				gene2cancers[gene] = gene2cancers[gene].union({cancer})
	pancan_genes = {gene for gene in gene2cancers if len(gene2cancers[gene]) / len(tumor2gene) >= 0.8}
	df1 = df_map[df_map["Hugo_Symbol"].isin(specific_genes)]
	df2 = df_map[~df_map["Hugo_Symbol"].isin(specific_genes)]
	df2_known = df2[df2["TUMORTYPE"].isin(tumor2gene)]
	df2_known = df2_known[df2_known.apply(lambda x: x["Gene"] in tumor2gene[x["TUMORTYPE"]], axis = 1)]
	df2_unknown = df2[~df2["TUMORTYPE"].isin(tumor2gene)]
	df2_unknown = df2_unknown[df2_unknown["Gene"].isin(pancan_genes)]
	df_map = pd.concat([df1, df2_known, df2_unknown])

	# exclude variants in gnomAD
	df_map = df_map[df_map["gnomADe_AF"].astype(float) <= af_cutoff]
	df_map.to_csv(output, sep = "\t", header = True, index = None)

	# expressed genes
	expr_genes = []
	for tumor in df2_known["TUMORTYPE"].unique():
		expr_genes.extend(tumor2gene[tumor])
	expr_genes = set(expr_genes).union(pancan_genes)
	pd.DataFrame({"Gene": sorted(expr_genes)}).to_csv(output.replace(".txt", "expr.txt"), sep = "\t", header = None, index = None)


for mapping_flag in ["single", "all"]:
	mutation_preprocessing(dirname + "/mutations/TCGA.maf", dirname + "/mutations/parsed_" + mapping_flag + ".txt", mapping_flag, 0, dirname + "/mutations/TCGA_preprocessed_" + mapping_flag + ".maf")
	mutation_preprocessing(dirname + "/mutations/non_TCGA.maf", dirname + "/mutations/parsed_" + mapping_flag + ".txt", mapping_flag, 0, dirname + "/mutations/non_TCGA_preprocessed_" + mapping_flag + ".maf")
	mutation_preprocessing(dirname + "/mutations/COSMIC/GenomeScreensMutant.v98.maf", dirname + "/mutations/COSMIC/parsed_" + mapping_flag + ".txt", mapping_flag, 0, dirname + "/mutations/COSMIC/GenomeScreensMutant_preprocessed_" + mapping_flag + ".maf")
