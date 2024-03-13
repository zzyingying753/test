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


def get_map(row):
    pdb_pos = [x + ":" + row["chain"] for x in all_function_py3.unzip_interface_residue(row["MappableResInPDBChainOnPDBBasis"])]
    uniprot_pos = all_function_py3.unzip_interface_residue(row["MappableResInPDBChainOnUniprotBasis"])
    df = pd.DataFrame({"structure_id": [row["PDBId"]] * len(pdb_pos), "residues": pdb_pos, "Uniprot": [row["hugo"]] * len(pdb_pos), "Uniprot_Position": uniprot_pos})
    return df

def get_HOTMAPS_input(shortversion_mutation_data, tissue, output_prefix):
    def get_path(pdb):
        if path.exists("/fs/cbsuhy01/storage/NetFlow3D/Data/PDB/" + pdb[1:-1] + "/pdb" + pdb + ".ent.gz"):
            return "/fs/cbsuhy01/storage/NetFlow3D/Data/PDB/" + pdb[1:-1] + "/pdb" + pdb + ".ent.gz"
        elif path.exists("/fs/cbsuhy01/storage/NetFlow3D/Data/PDBlike/" + pdb[1:-1] + "/pdb" + pdb + ".ent.gz"):
            return "/fs/cbsuhy01/storage/NetFlow3D/Data/PDBlike/" + pdb[1:-1] + "/pdb" + pdb + ".ent.gz"
        else:
            return np.nan

    df_sifts = pd.read_csv("/local/storage/resources/sifts/parsed_files/pdbresiduemapping.txt", sep = "\t", dtype = str).dropna()
    df = pd.read_csv(shortversion_mutation_data, sep = "\t", dtype = str)[["Uniprot", "Uniprot_Position", "Mutations"]].rename(columns = {"Mutations": "occurrence"})
    df_sifts = df_sifts[df_sifts["UniProt"].isin(set(df["Uniprot"]))].rename(columns = {"PDB": "PDBId", "Chain": "chain", "UniProt": "hugo"})
    df_sifts["PDBId"] = df_sifts["PDBId"].apply(lambda x: x.lower())
    df_sifts["path"] = df_sifts["PDBId"].apply(get_path)
    logging.info(shortversion_mutation_data.split("/")[-2] + " before filtering: " + str(df_sifts.shape[0]))
    df_sifts = df_sifts.dropna(subset = ["path"])
    logging.info(shortversion_mutation_data.split("/")[-2] + " after filtering: " + str(df_sifts.shape[0]))
    df_sifts[["PDBId", "chain", "hugo", "path"]].to_csv(output_prefix + "_1.txt", sep = "\t", header = True, index = None)
    df_map = []
    for uniprot,df_tmp in df_sifts.groupby("hugo"):
        df_map.extend(df_tmp.apply(get_map, axis = 1).tolist())
    df_map = pd.concat(df_map)
    df = pd.merge(df, df_map, on = ["Uniprot", "Uniprot_Position"])
    df["tissue"] = tissue
    df[["structure_id", "tissue", "residues", "occurrence"]].to_csv(output_prefix + "_2.txt", sep = "\t", header = True, index = None)
    return

def get_3dhotspot_input(maf_file, output):
    df = pd.read_csv(maf_file, sep = "\t", dtype = str, usecols = ["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", "HGVSp", "Transcript_ID", "Gene", "Protein_position", "ENSP", "SWISSPROT", "TUMORTYPE"])
    logging.info(maf_file.split("/")[-1] + ";" + str(df.shape[0]) + ";raw")
    df = df[df["Variant_Classification"] == "Missense_Mutation"].dropna(subset = ["Protein_position"])
    logging.info(maf_file.split("/")[-1] + ";" + str(df.shape[0]) + ";missense")
    df["CODE"] = df["TUMORTYPE"].apply(lambda x: x.replace("TCGA-", ""))
    df["HGVSp_Short"] = df["HGVSp"].apply(all_function_py3.get_HGVSp_short)
    df[["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", "HGVSp_Short", "Transcript_ID", "Gene", "Protein_position", "ENSP", "SWISSPROT", "CODE"]].to_csv(output, sep = "\t", header = True, index = None)

def get_HOTMAPS_input_multirun(args):
    return get_HOTMAPS_input(*args)

def get_3dhotspot_input_multirun(args):
    return get_3dhotspot_input(*args)


function_input_hotmaps = []
function_input_3dhotspot = []
for maf_file in [dirname + "/mutations/COSMIC/GenomeScreensMutant_preprocessed_single.maf", dirname + "/mutations/TCGA_preprocessed_single.maf", dirname + "/mutations/non_TCGA_preprocessed_single.maf"]:
    command = "python /local/storage/yz2296/cancer_hotspot_new/website/website/revision/NetFlow3D/NetFlow3D.py -m " + \
    maf_file + " -R low -I " + maf_file.split("/")[-1].replace("_preprocessed_single.maf", "") + " -t 50 -X " + \
    maf_file.replace(".maf", ".expr") + " -N"
    os.system(command)
    function_input_hotmaps.append([dirname + "/NetFlow3D/output/" + maf_file.split("/")[-1].replace("_preprocessed_single.maf", "") + "/ShortVersion_mutation_data.txt", "all", maf_file.replace(".maf", "_HOTMAPS")])
    function_input_3dhotspot.append([maf_file, maf_file.replace(".maf", "_3dhotspot.txt")])

p = Pool(3)
output = p.map(get_HOTMAPS_input_multirun, function_input_hotmaps)
p.close()

p = Pool(3)
output = p.map(get_3dhotspot_input_multirun, function_input_3dhotspot)
p.close()
