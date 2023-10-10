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

command = "python " + dirname + "/NetFlow3D/NetFlow3D.py -m " + dirname + "/mutations/TCGA_preprocessed_single.maf \
-R low -I TCGA -t 50 -X " + dirname + "/mutations/TCGA_preprocessed_single.expr"

os.system(command)

df_mut = pd.read_csv(dirname + "/mutations/TCGA_preprocessed_single.maf", sep = "\t")
all_genes = set(df_mut["Hugo_Symbol"])
cancer_genes = all_function_py3.get_cgc_genes()

pval_cutoff = [float("1e-" + str(x)) for x in range(0, 31)]
for file in ["PDB_intra_pvalue.txt", "PDB_inter_pvalue.txt", "AlphaFold2_intra_pvalue_pLDDT0.txt", "PIONEER_inter_pvalue.txt"]:#, "All_intra_LoF_pvalue.txt"]:
	df = pd.read_csv(dirname + "/NetFlow3D/output/TCGA/" + file, sep = "\t")
	enr_list = []
	tick_list = []
	df_all = []
	for x in pval_cutoff:
		df_tmp = df[df["Raw_pvalue"] < float(x)]
		uniprot2res = defaultdict(set)
		for y in df_tmp["Residues"]:
			for res in y.split(","):
				uniprot, pos = res.split("_")
				uniprot2res[uniprot] = uniprot2res[uniprot].union({int(pos)})
		df_output = {"UniProt": [], "Protein_position": []}
		for uniprot in uniprot2res:
			df_output["UniProt"].extend([uniprot] * len(uniprot2res[uniprot]))
			df_output["Protein_position"].extend(sorted(uniprot2res[uniprot]))
		df_output = pd.DataFrame(df_output)
		output_file = dirname + "/pval_sensitivity/" + "_".join(file.split("_")[0:2]) + "_" + str(x) + ".txt"
		df_output.to_csv(output_file, sep = "\t", header = True, index = None)
		
		# enrichment of known cancer genes
		candidate_genes = set(df_mut[df_mut["UniProt"].isin(uniprot2res)]["Hugo_Symbol"])
		c1 = len(candidate_genes.intersection(cancer_genes))
		n1 = len(candidate_genes)
		c2 = len(all_genes.intersection(cancer_genes))
		n2 = len(all_genes)
		enr, pval = all_function_py3.calc_enr(c1, n1, c2, n2)
		enr_list.append(enr)
		tick_list.append(-np.log10(x))
		print(-np.log10(x), len(candidate_genes), enr, pval)

		# enrichment on functional sites
		command = "python /local/storage/resources/tools/functional_sites.py -m " + output_file + " -o " + output_file
		os.system(command)
		df_output = pd.read_csv(output_file, sep = "\t")
		df_output["-log10(pval_cutoff)"] = -np.log10(x)
		df_all.append(df_output)

	# enrichment on functional sites
	df_all = pd.concat(df_all)
	for site_type, df_one in df_all.groupby("Type"):
		if site_type != "ENPS_dna_binding_sites":
			df_one = df_one[df_one["Proteins_of_interest"] == "input"]
			fig, ax = all_function_py3.scatter_plot(df_one["-log10(pval_cutoff)"].tolist(), df_one["Enrichment"].tolist(), yerr = 0, yerr_capsize = 2, yerr_width = 1, point_size = 20, color = "#525252", markers = 'o', \
													point_alpha = 1, point_edge = 0, edge_color = "#737373", xlim = None, ylim = None, title = file, xlabel="-log10(p value)", ylabel="Enrichment of\n" + site_type, fontsize = 14, \
													xtick_rotation = 0, xtick_align = "center", figsize = (6,3), labelpad = 10)
			ax.axvline(x = -np.log10(0.05 / df.shape[0]), color = 'r', linestyle = '--')
			fig.savefig(dirname + "/figures/" + "_".join(file.split("_")[0:2]) + "_" + site_type + ".png")
	
	# enrichment of known cancer genes 
	fig, ax = all_function_py3.scatter_plot(tick_list, enr_list, yerr = 0, output = None, yerr_capsize = 2, yerr_width = 1, point_size = 20, color = "#525252", markers = 'o', \
		point_alpha = 1, point_edge = 0, edge_color = "#737373", xlim = None, ylim = None, title = file, xlabel="-log10(p value)", ylabel="Enrichment of\nknown cancer genes", fontsize = 14, \
		xtick_rotation = 0, xtick_align = "center", figsize = (6,3), labelpad = 10)
	# ax.plot(tick_list, enr_list)
	ax.axvline(x = -np.log10(0.05 / df.shape[0]), color = 'r', linestyle = '--')
	fig.savefig(dirname + "/figures/pval_sensitivity_" + "_".join(file.split("_")[0:2]) + ".png")

