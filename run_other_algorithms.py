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
    pdb_pos = all_function_py3.unzip_interface_residue(row["MappableResInPDBChainOnPDBBasis"])
    uniprot_pos = all_function_py3.unzip_interface_residue(row["MappableResInPDBChainOnUniprotBasis"])
    df = pd.DataFrame({"pdb_id": [row["PDB"]] * len(pdb_pos), "chain": [row["Chain"]] * len(pdb_pos), "residue": pdb_pos, "UniProt": [row["UniProt"]] * len(pdb_pos), "Reference Codon Position": uniprot_pos})
    return df

def get_HOTMAPS_input(maf_file, output):
    df_sifts = pd.read_csv("/local/storage/resources/sifts/parsed_files/pdbresiduemapping.txt", sep = "\t", dtype = str).dropna()
    df = pd.read_csv(maf_file, sep = "\t", dtype = str)
    df = df[df["Variant_Classification"] == "Missense_Mutation"].dropna(subset = ["UniProt", "Protein_position"])
    df = df[df["Protein_position"].apply(lambda x: x.find("-") == -1)]
    df["Transcript_ID"] = df["Transcript_ID"].apply(lambda x: x.split(".")[0])
    df["Amino_acids"] = df["Amino_acids"].apply(lambda x: x.split("/")[0])
    df = df[["UniProt","Hugo_Symbol", "Transcript_ID", "Protein_position", "Amino_acids", "Start_Position", "Chromosome"]].rename(columns = {"Hugo_Symbol": "HUGO symbol", "Transcript_ID": "Reference Transcript", "Protein_position": "Reference Codon Position", "Amino_acids": "Reference AA", "Start_Position": "Reference Genomic Position"})

    df_map = []
    for uniprot in set(df["UniProt"]).intersection(set(df_sifts["UniProt"])):
        df_tmp = df_sifts[df_sifts["UniProt"] == uniprot]
        df_map.extend(df_tmp.apply(get_map, axis = 1).tolist())
    df_map = pd.concat(df_map)
    df = pd.merge(df, df_map, on = ["UniProt", "Reference Codon Position"])[['pdb_id', 'chain', 'residue', 'HUGO symbol', 'Reference Transcript', 'Reference Codon Position', 'Reference AA', 'Reference Genomic Position', 'Chromosome']]
    df.to_csv(output, sep = "\t", header = True, index = None)
    return

def get_3dhotspot_input(maf, output):
    df = pd.read_csv(maf_file, sep = "\t", dtype = str, usecols = ["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", "HGVSp", "Transcript_ID", "Gene", "Protein_position", "ENSP", "SWISSPROT", "TUMORTYPE"])
    df = df[df["Variant_Classification"] == "Missense_Mutation"].dropna(subset = ["Protein_position"])
    df["CODE"] = df["TUMORTYPE"].apply(lambda x: x.replace("TCGA-", ""))
    df["HGVSp_Short"] = df["HGVSp"].apply(all_function_py3.get_HGVSp_short)
    df[["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", "HGVSp_Short", "Transcript_ID", "Gene", "Protein_position", "ENSP", "SWISSPROT", "CODE"]].to_csv(output, sep = "\t", header = True, index = None)

def get_HOTMAPS_input_multirun(args):
    return get_HOTMAPS_input(*args)

def get_3dhotspot_input_multirun(args):
    return get_3dhotspot_input(*args)


function_input_hotmaps = []
function_input_3dhotspot = []
for maf_file in [dirname + "/mutations/TCGA_preprocessed_single.maf", dirname + "/mutations/non_TCGA_preprocessed_single.maf", dirname + "/mutations/COSMIC/GenomeScreensMutant_preprocessed_single.maf"]:
    function_input_hotmaps.append([maf_file, maf_file.replace(".maf", "_HOTMAPS.txt")])
    function_input_3dhotspot.append([maf_file, maf_file.replace(".maf", "_3dhotspot.txt")])

# p = Pool(3)
# output = p.map(get_HOTMAPS_input_multirun, function_input_hotmaps)
# p.close()

p = Pool(3)
output = p.map(get_3dhotspot_input_multirun, function_input_3dhotspot)
p.close()
