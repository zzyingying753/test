import pandas as pd
import os
import numpy as np
import pysam
import requests
import json
from collections import defaultdict
from multiprocessing import Pool
import logging
import copy
from os import path
from glob import glob
import subprocess
import matplotlib.pyplot as plt
dirname = os.path.dirname(os.path.abspath(__file__))
from Bio import SeqIO
import gzip
from joblib import Parallel, delayed
import xml.etree.ElementTree as ET
import random
import networkx as nx

def gene_tissue_expr(biospecimen_file, manifest_file, sample_type, output_path, data_path = "/fs/cbsuhyfs1/storage1/yz2296/TCGA_raw/transcriptome/"):
	'''
	sample_type: 1. normal; 2. tumor;
	bispecimen_file: "/fs/cbsuhy01/storage/yz2296/cancer_hotspot_new/website/website/revision/clinical_data/biospecimen.cases_selection.2023-07-29/sample.tsv"
	'''
	def parse_one_file(row):
		df_tmp = pd.read_csv(data_path + row["id"] + "/" + row["filename"], sep = "\t", comment = "#")
		for gene in df_tmp[df_tmp["tpm_unstranded"] >= 1]["gene_id"]:
			gene2sample[gene.split(".")[0]].append(row["Tumor_Sample_Barcode"][0:-1])
		return True

	df_manifest = pd.read_csv(manifest_file, sep = "\t")
	df_manifest["Tumor_Sample_Barcode"] = df_manifest["Tumor_Sample_Barcode"].apply(lambda x: x[0:16])
	df = pd.read_csv(biospecimen_file, sep = "\t")
	df = df[df["project_id"].apply(lambda x: x[0:4] == "TCGA")]
	for key,df_one in df.groupby("project_id"):
		if key in {"TCGA-LAML", "TCGA-DLBC"}:
			if sample_type == "normal":
				df_one = df_one[df_one["sample_type"].isin({"Solid Tissue Normal", "Blood Derived Normal", "Bone Marrow Normal"})]
			else:
				df_one = df_one[df_one["sample_type"].isin({"Primary Tumor", "Primary Blood Derived Cancer - Peripheral Blood"})]
		else:
			if sample_type == "normal":
				df_one = df_one[df_one["sample_type"] == "Solid Tissue Normal"]
			else:
				df_one = df_one[df_one["sample_type"] == "Primary Tumor"]
		df_manifest_one = df_manifest[df_manifest["Tumor_Sample_Barcode"].isin(df_one["sample_submitter_id"].tolist())]
		if df_manifest_one.shape[0] > 0:
			gene2sample = defaultdict(list)
			df_manifest_one["status"] = df_manifest_one.apply(parse_one_file, axis = 1)
			all_samples = set(df_manifest_one["Tumor_Sample_Barcode"].apply(lambda x: x[0:-1]))
			expr_genes = []
			for gene in gene2sample:
				if len(set(gene2sample[gene])) / len(all_samples) >= 0.8:
					expr_genes.append(gene)
			pd.DataFrame({"Gene": expr_genes}).to_csv(output_path + key + "_" + sample_type + ".txt", sep = "\t", header = None, index = None)
	
	return

# Transcriptome profile
url = 'https://api.gdc.cancer.gov/files/'
manifest_df = pd.read_csv(dirname + "/transcriptome_profiles/gdc_manifest.2023-07-29.txt", sep = "\t")

tumor_sample_barcodes = []
for file_id in manifest_df['id']:
	 response = requests.get(url + file_id, params={'fields': 'associated_entities.entity_submitter_id'})
	 data = json.loads(response.text)
	 tumor_sample_barcodes.append(data['data']['associated_entities'][0]['entity_submitter_id'])
	 if len(data['data']['associated_entities']) > 1:
	 	print("yoyoyoyoyoyo" + file_id)

manifest_df['Tumor_Sample_Barcode'] = tumor_sample_barcodes
manifest_df.to_csv(dirname + '/transcriptome_profiles/gdc_manifest.2023-07-29.annotated.txt', sep = '\t', header = True, index = None)

for sample_type in ["normal", "tumor"]:
	gene_tissue_expr(dirname + "/transcriptome_profiles/biospecimen.cases_selection.2023-07-29/sample.tsv", \
		dirname + "/transcriptome_profiles/gdc_manifest.2023-07-29.annotated.txt", sample_type, \
		dirname + "/transcriptome_profiles/parsed_files/")
