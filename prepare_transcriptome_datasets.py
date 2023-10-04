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
