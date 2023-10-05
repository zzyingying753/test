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

def get_single_uniprot(row, ensp2uniprot, enst2uniprot, ensg2uniprot):
	if type(row["ENSP"]) == str and row["ENSP"].split(".")[0] in ensp2uniprot:
		return ensp2uniprot[row["ENSP"].split(".")[0]]
	elif type(row["Transcript_ID"]) == str and row["Transcript_ID"].split(".")[0] in enst2uniprot:
		return enst2uniprot[row["Transcript_ID"].split(".")[0]][0]
	elif type(row["Gene"]) == str and row["Gene"].split(".")[0] in ensg2uniprot:
		return ensg2uniprot[row["Gene"].split(".")[0]][0]
	else:
		return np.nan

def get_multiple_uniprots(row, ensp2uniprot, enst2uniprot, ensg2uniprot):
	if type(row["ENSP"]) == str and len(ensp2uniprot[row["ENSP"].split(".")[0]]) > 0:
		uniprot_output = ensp2uniprot[row["ENSP"].split(".")[0]]
	elif type(row["Transcript_ID"]) == str and len(enst2uniprot[row["Transcript_ID"].split(".")[0]]) > 0:
		uniprot_output = enst2uniprot[row["Transcript_ID"].split(".")[0]]
	elif type(row["Gene"]) == str and len(ensg2uniprot[row["Gene"].split(".")[0]]) > 0:
		uniprot_output = ensg2uniprot[row["Gene"].split(".")[0]]
	else:
		uniprot_output = set()

	if len(uniprot_output) > 0:
		new_rows = []
		for uniprot in uniprot_output:
			new_row = row.copy()
			new_row["UniProt"] = uniprot
			new_rows.append(new_row)
	else:
		new_row = row.copy()
		new_row["UniProt"] = np.nan
		new_rows = [new_row]

	return pd.DataFrame(new_rows)


def uniprot_sanity_check(row, uniprot2seq):
	pos = int(row["Protein_position"].split("-")[0])
	if row["UniProt"] in uniprot2seq and len(uniprot2seq[row["UniProt"]]) >= pos and uniprot2seq[row["UniProt"]][pos - 1] == row["Amino_acids"].split("/")[0][0]:
		return True
	else:
		return False

def eliminate_redundancy(key, df_one):
	df_processed = df_one[["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "TUMORTYPE", "Variant_Classification", "Gene", "Hugo_Symbol", "UniProt", "Protein_position", "Amino_acids", "Codons"]].drop_duplicates()
	if df_processed.shape[0] > 1:
		if df_processed.dropna(subset = ["Protein_position"]).shape[0] > 0:
			return df_processed.dropna(subset = ["Protein_position"]).sort_values(by = "Protein_position")[0:1]
		else:
			return df_processed.sort_values(by = ["Gene"])[0:1]
	else:
		return df_processed

def eliminate_redundancy_multirun(args):
	return eliminate_redundancy(*args)

def mutation_preprocessing(maf_file, vep_file, mapping_flag, af_cutoff, output):
	'''
	mapping_flag: 1. single; 2. all.
	'''

	# definition
	lof = {"Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site"}
	inframe = {"In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"}
	specific_genes = {"AR", "CDH4", "EGFR", "EPHA3", "ERBB4", "FGFR2", "FLT3", "FOXA1", "FOXA2", "MECOM", "MIR142", "MSH4", "PDGFRA", "SOX1", "SOX9", "SOX17", "TBX3", "WT1"}

	# # load raw mutation data
	# df = pd.read_csv(maf_file, sep = "\t", dtype = str)
	# df = df.dropna(subset = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"])
	# df["#Uploaded_variation"] = df.apply(lambda x: x["Chromosome"] + "_" + str(x["Start_Position"]) + "_" + x["Reference_Allele"] + "/" + x["Tumor_Seq_Allele2"], axis = 1)

	# # load VEP output
	# df_vep = pd.read_csv(vep_file, sep = "\t", dtype = str).rename(columns = {"Feature": "Transcript_ID"})
	# df_vep["Variant_Classification"] = df_vep["Consequence"].apply(all_function_py3.consequence2varclass_coding)
	# ensg2symbol = all_function_py3.gene_ID_conversion(df_vep["Gene"].dropna().unique().tolist(), convert_from = "ensembl.gene", convert_to = "symbol", species = "human")
	# df_vep["Hugo_Symbol"] = df_vep.apply(lambda x: get_hugo_symbol(x, ensg2symbol), axis = 1)
	# df_map = pd.merge(df, df_vep, on = "#Uploaded_variation")

	# # exclude unexpressed mutations
	# tumor2gene = {}
	# gene2cancers = defaultdict(set)
	# for x in os.listdir(dirname + "/transcriptome_profiles/parsed_files/"):
	# 	if x[-len("_tumor.txt"):] == "_tumor.txt":
	# 		cancer = x.split("_")[0]
	# 		tumor2gene[cancer] = pd.read_csv(dirname + "/transcriptome_profiles/parsed_files/" + x, sep = "\t", header = None)[0].tolist()
	# 		for gene in tumor2gene[cancer]:
	# 			gene2cancers[gene] = gene2cancers[gene].union({cancer})
	# pancan_genes = {gene for gene in gene2cancers if len(gene2cancers[gene]) / len(tumor2gene) >= 0.8}
	# df1 = df_map[df_map["Hugo_Symbol"].isin(specific_genes)]
	# df2 = df_map[~df_map["Hugo_Symbol"].isin(specific_genes)]
	# df2_known = df2[df2["TUMORTYPE"].isin(tumor2gene)]
	# df2_known = df2_known[df2_known.apply(lambda x: x["Gene"] in tumor2gene[x["TUMORTYPE"]], axis = 1)]
	# df2_unknown = df2[~df2["TUMORTYPE"].isin(tumor2gene)]
	# df2_unknown = df2_unknown[df2_unknown["Gene"].isin(pancan_genes)]
	# df_map = pd.concat([df1, df2_known, df2_unknown])

	# # expressed genes
	# expr_genes = []
	# for tumor in df2_known["TUMORTYPE"].unique():
	# 	expr_genes.extend(tumor2gene[tumor])
	# expr_genes = set(expr_genes).union(pancan_genes)
	# pd.DataFrame({"Gene": sorted(expr_genes)}).to_csv(output.replace(".maf", ".expr"), sep = "\t", header = None, index = None)

	# # exclude variants in gnomAD
	# df_map = df_map[df_map["gnomADe_AF"].astype(float) <= af_cutoff]

	# # map to UniProt
	# if mapping_flag == "single":
	# 	with open(dirname + "/NetFlow3D/metadata/ensp2uniprot.json") as f:
	# 		ensp2uniprot = json.load(f)
	# 	with open(dirname + "/NetFlow3D/metadata/ensg2uniprot.json") as f:
	# 		ensg2uniprot = json.load(f)
	# 	with open(dirname + "/NetFlow3D/metadata/enst2uniprot.json") as f:
	# 		enst2uniprot = json.load(f)
	# 	df_map["UniProt"] = df_map.apply(lambda x: get_single_uniprot(x, ensp2uniprot, enst2uniprot, ensg2uniprot), axis = 1)
	# else:
	# 	ensp2uniprot = all_function_py3.gene_ID_conversion({x.split(".")[0] for x in df_map["ENSP"].dropna().unique()}, convert_from = "ensembl.protein", convert_to = "uniprot", species = "human")
	# 	enst2uniprot = all_function_py3.gene_ID_conversion({x.split(".")[0] for x in df_map["Transcript_ID"].dropna().unique()}, convert_from = "ensembl.transcript", convert_to = "uniprot", species = "human")
	# 	ensg2uniprot = all_function_py3.gene_ID_conversion({x.split(".")[0] for x in df_map["Gene"].dropna().unique()}, convert_from = "ensembl.gene", convert_to = "uniprot", species = "human")
	# 	df_map = pd.concat(df_map.apply(lambda x: get_multiple_uniprots(x, ensp2uniprot, enst2uniprot, ensg2uniprot), axis = 1).tolist(), ignore_index = True)
	# df_map = df_map.dropna(subset = ["UniProt"])
	
	# sanity check
	df_map = pd.read_csv(output, sep = "\t", dtype = str)
	with open(dirname + "/NetFlow3D/metadata/uniprot2proseq.json") as f:
		uniprot2seq = json.load(f)
	df_inframe = df_map[df_map["Variant_Classification"].isin(inframe)].dropna(subset = ["Protein_position"])
	df_inframe = df_inframe[df_inframe.apply(lambda x: uniprot_sanity_check(x, uniprot2seq), axis = 1)]
	df_lof = df_map[df_map["Variant_Classification"].isin(lof)]
	function_input = []
	for x in ["df_inframe", "df_lof"]:
		for key, df_one in eval(x).groupby(["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "UniProt"]):
			function_input.append([key, df_one])
	p = Pool(50)
	df_map = p.map(eliminate_redundancy_multirun, function_input)
	p.close()
	df_map = pd.concat(df_map)
	df_map.to_csv(output, sep = "\t", header = True, index = None)
	return


for mapping_flag in ["all"]:#["single", "all"]:
	mutation_preprocessing(dirname + "/mutations/TCGA.maf", dirname + "/mutations/parsed_" + mapping_flag + ".txt", mapping_flag, 0, dirname + "/mutations/TCGA_preprocessed_" + mapping_flag + ".maf")
	mutation_preprocessing(dirname + "/mutations/non_TCGA.maf", dirname + "/mutations/parsed_" + mapping_flag + ".txt", mapping_flag, 0, dirname + "/mutations/non_TCGA_preprocessed_" + mapping_flag + ".maf")
	mutation_preprocessing(dirname + "/mutations/COSMIC/GenomeScreensMutant.v98.maf", dirname + "/mutations/COSMIC/parsed_" + mapping_flag + ".txt", mapping_flag, 0, dirname + "/mutations/COSMIC/GenomeScreensMutant_preprocessed_" + mapping_flag + ".maf")
