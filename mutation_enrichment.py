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

def get_pathway_enr(module_proteins, mut_all, mut_cluster, total_prolen, prolen_dict):
    module_len = 0
    for x in set(module_proteins):
        module_len = module_len + prolen_dict[x]
    
    # enrichment of clustered muts
    cluster_one = mut_cluster[mut_cluster["UniProt"].isin(module_proteins)]
    c1 = cluster_one.shape[0]
    n1 = mut_cluster.shape[0]
    enr1, pval1 = all_function_py3.calc_enr(c1, n1, module_len, total_prolen)
    
    # enrichment of all muts
    all_one = mut_all[mut_all["UniProt"].isin(module_proteins)]
    c2 = all_one.shape[0]
    n2 = mut_all.shape[0]
    enr2, pval2 = all_function_py3.calc_enr(c2, n2, module_len, total_prolen)
    
    return [enr1, pval1, c1, n1, enr2, pval2, c2, n2]


cgc_genes = all_function_py3.get_cgc_genes()
df_mut = pd.read_csv(dirname + "/independent_datasets/TCGA_preprocessed_single.maf", sep = "\t")
TP_uniprots = set(df_mut[df_mut["Hugo_Symbol"].isin(cgc_genes)]["UniProt"])
df_mut = df_mut[df_mut["Variant_Classification"].isin({"Missense_Mutation", "In_Frame_Del", "In_Frame_Ins"})]

with open(dirname + "/NetFlow3D/metadata/uniprot2prolen.json") as f:
    prolen = json.load(f)
with open(dirname + "/NetFlow3D/metadata/unused/gene2uniprot.json") as f:
    gene2uniprot = json.load(f)

total_len = 0
all_proteins = set(df_mut["UniProt"])
for x in all_proteins:
    total_len = total_len + prolen[x]

df = pd.read_csv(dirname + "/NetFlow3D/output/TCGA_drivers.txt", sep = "\t")
df = df[df["Variant_Classification"].isin({"Missense_Mutation", "In_Frame_Del", "In_Frame_Ins"})]

# initialize output data frame
df_box = {"Group": [], "Type": [], "Pathway": [], "Enrichment": [], "Pvalue": [], "Module_mutations": [], \
          "All_mutations": []}

# known cancer pathways
df_pathway = pd.read_excel(dirname + "/table_s1_to_s57/table_s8.xlsx")
cancer_pathways = df_pathway["Cancer_pathway"].tolist()
for i,x in enumerate(df_pathway["Genes"]):
    module_proteins = {gene2uniprot[y.strip(" ")] for y in x.split(",") if y in gene2uniprot}
    if len(module_proteins) > 0:
        enr1, pval1, c1, n1, enr2, pval2, c2, n2 = get_pathway_enr(module_proteins, df_mut, df, total_len, prolen)
        df_box["Group"].extend(["Known cancer pathways", "Known cancer pathways"])
        df_box["Type"].extend(["Spatially clustered mutations", "All mutations"])
        df_box["Pathway"].extend([cancer_pathways[i], cancer_pathways[i]])
        df_box["Enrichment"].extend([enr1, enr2])
        df_box["Pvalue"].extend([pval1, pval2])
        df_box["Module_mutations"].extend([c1, c2])
        df_box["All_mutations"].extend([n1, n2])
        if enr1 < 1:
            print(x, enr1)

# Our subnetworks
tag2file = {"NetFlow3D": dirname + "/NetFlow3D/output/TCGA_subnetworks_intercept1.0_lowres_edgeweightTrue.txt", \
# "Standard": dirname + "/NetFlow3D/output/TCGA_standard_network_subnetworks.txt", \
# "Mix": ""
}
# for x in os.listdir(dirname + "/NetFlow3D/output/"):
#     if x[0:len("TCGA_subnetworks")] == "TCGA_subnetworks":
#         tag2file[x.replace("TCGA_subnetworks_", "").replace(".txt", "")] = dirname + "/NetFlow3D/output/" + x
for tag,file in tag2file.items():
    df_pathway = pd.read_csv(file, sep = "\t")
    df_pathway = df_pathway[df_pathway["Subnetwork_size"] > 5]
    for key in df_pathway["Subnetwork_UniProts"]:
        enr1, pval1, c1, n1, enr2, pval2, c2, n2 = get_pathway_enr(key.split(","), df_mut, df, total_len, prolen)
        df_box["Type"].extend(["Spatially clustered mutations", "All mutations"])
        df_box["Pathway"].extend([key, key])
        df_box["Enrichment"].extend([enr1, enr2])
        df_box["Pvalue"].extend([pval1, pval2])
        df_box["Module_mutations"].extend([c1, c2])
        df_box["All_mutations"].extend([n1, n2])
        if tag == "NetFlow3D": 
            if set(key.split(",")).intersection(TP_uniprots) != set():
                df_box["Group"].extend([tag + "_KCG", tag + "_KCG"])
            else:
                df_box["Group"].extend([tag + "_NKCG", tag + "_NKCG"])
        else:
            df_box["Group"].extend([tag, tag])
    
# GO biological processes
df_pathway = pd.read_excel(dirname + "/table_s1_to_s57/table_s9.xlsx", sep = "\t")
df_pathway = df_pathway[df_pathway["Genes"].apply(lambda x: set(x.split(",")).intersection(cgc_genes) == set())]
all_pathways = df_pathway["Biological_process"].tolist()
for i,x in enumerate(df_pathway["Genes"]):
    module_proteins = {gene2uniprot[y] for y in x.split(",") if y in gene2uniprot}
    if len(module_proteins) > 0:
        enr1, pval1, c1, n1, enr2, pval2, c2, n2 = get_pathway_enr(module_proteins, df_mut, df, total_len, prolen)
        df_box["Group"].extend(["GO biological processes", "GO biological processes"])
        df_box["Type"].extend(["Spatially clustered mutations", "All mutations"])
        df_box["Pathway"].extend([all_pathways[i], all_pathways[i]])
        df_box["Enrichment"].extend([enr1, enr2])
        df_box["Pvalue"].extend([pval1, pval2])
        df_box["Module_mutations"].extend([c1, c2])
        df_box["All_mutations"].extend([n1, n2])
df_box = pd.DataFrame(df_box)
df_box["log(Enrichment)"] = df_box["Enrichment"].apply(lambda x: np.log(x))
df_box.to_csv(dirname + "/table_s1_to_s57/pathway_enrichment.txt", sep = "\t", header = True, index = None)

df_fc = {"Group": [], "Fold_change": [], "Pathway": []}#, "Pvalue": []}
for group, df_one in df_box.groupby("Group"):
    for pathway, df_path in df_one.groupby("Pathway"):
        df_path = df_path.sort_values(by = "Type")
        c2, c1 = df_path["Module_mutations"].tolist()
        n2, n1 = df_path["All_mutations"].tolist()
        enr, pval = all_function_py3.calc_enr(c1, n1, c2, n2)
        df_fc["Fold_change"].append(enr)
        df_fc["Pathway"].append(pathway)
        df_fc["Group"].append(group)
        # df_fc["Pvalue"].append(pval)
df_fc = pd.DataFrame(df_fc)
df_fc.to_csv(dirname + "/revision/table_s1_to_s57/pathway_FC.txt", sep = "\t", header = True, index = None)


# ----------------------------------------------------------
# draw the figure

# df_box = pd.read_csv(dirname + "/table_s1_to_s57/pathway_enrichment.txt", sep = "\t")
df_box = df_box[df_box["Group"].isin({"Known cancer pathways", "NetFlow3D_KCG", "NetFlow3D_NKCG", "GO biological processes"})]
color = {"Spatially clustered mutations": "red", "All mutations": "#bdbdbd"}
for x,df_one in df_box.groupby("Group"):
    df1 = df_one[df_one["Type"] == "Spatially clustered mutations"]
    df2 = df_one[df_one["Type"] == "All mutations"]
    _, pval = scipy.stats.mannwhitneyu(df1["Enrichment"].tolist(), df2["Enrichment"].tolist())
    print(x, pval, df1.shape[0], df1["log(Enrichment)"].median(), df2["log(Enrichment)"].median())
fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw = {'width_ratios': [3, 1]}, figsize = (16, 7))
_, _ = all_function_py3.box_plot(df_box, "Group", "log(Enrichment)", color, contrast_col = "Type", barWidth = 0.7, \
    fontsize = 20, xlim = (-4,4), tick_rotation = 0, tick_align = "right", legend = True, tick_bins = 6,\
    orientation = "horizontal", contrast_widthfrac = 0.85, notch = False, outlier = False, ax = ax1, subplot = True)

# df_fc = pd.read_csv(dirname + "/table_s1_to_s57/pathway_FC.txt", sep = "\t")
df_fc = df_fc[df_fc["Group"].isin({"Known cancer pathways", "NetFlow3D_KCG", "NetFlow3D_NKCG", "GO biological processes"})]
for x, df_one in df_fc.groupby("Group"):
    print(x, df_one["Fold_change"].median())
for x in ["NetFlow3D_KCG", "NetFlow3D_NKCG"]:
    for y in ["Known cancer pathways", "GO biological processes"]:
        _, pval = scipy.stats.mannwhitneyu(df_fc[df_fc["Group"] == x]["Fold_change"].dropna().tolist(), \
                                           df_fc[df_fc["Group"] == y]["Fold_change"].dropna().tolist())
        print(x, y, pval, df_fc[df_fc["Group"] == x].shape[0], df_fc[df_fc["Group"] == y].shape[0])
ax2 = plt.subplot(1,2,2, sharey = ax1)
_, ax = all_function_py3.violin_plot(df_fc, "Group", "Fold_change", {x: "#73c2fb" for x in df_fc["Group"].unique()},
                                       x_col_order = ["Known cancer pathways", "NetFlow3D_KCG", "NetFlow3D_NKCG", "GO biological processes"], \
                                       # x_col_order = df_box["Group"].unique().tolist(), 
                                       xlim = (-1, 9), linewidth = 3, tick_bins = 5, fontsize = 20,
                                       tick_rotation = 0, tick_align = "right", orientation = "horizontal", \
                                       ax = ax2, subplot = True)
ax2.tick_params('y', labelleft=False)
ymin, ymax = ax.get_ylim()
ax.plot([1,1], [0,len(df_fc["Group"].unique())], color = "red", linestyle = "dashed")
ax.set_ylim(ymin, ymax)
fig.savefig(dirname + "/figures/mutation_enrichment.png")
