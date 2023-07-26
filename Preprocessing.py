import sys
import pandas as pd
import numpy as np
import os
sys.path.append("/home/yz2296/bin")
import all_function_py3

def format_var(upload):
    chrom, pos, mut = upload.split("_")
    ref, alt = mut.split("/")
    return [chrom, pos, ref, alt]

def format_id_mapping(file):
    df_map = pd.read_csv(file, sep = "\t", dtype = str)
    df_map_final = []
    for key,df_one in df_map.groupby("From"):
        if df_one.shape[0] > 1:
            if df_one[df_one["Reviewed"] == "reviewed"].shape[0] > 0:
                df_one = df_one[df_one["Reviewed"] == "reviewed"].sort_values(by = "Entry")[0:1]
                df_map_final.append(df_one)
            else:
                df_map_final.append(df_one.sort_values(by = "Entry")[0:1])
        else:
            df_map_final.append(df_one)
    df_map = pd.concat(df_map_final)
    return df_map

def convert_hgvs(hgvs):
    if type(hgvs) == str:
        return hgvs.split(":")[1]
    else:
        return hgvs

def get_expr_genes(expr, expr_level = 1, sample_frac = 0.8):
    i = 0
    for x in expr.split(","):
        if float(x) >= expr_level:
            i = i + 1
    if i / len(expr.split(",")) >= sample_frac:
        return True
    else:
        return False

def get_expressed_PPInetwork(cancers, output, specific_uniprot_file = "/local/storage/yz2296/cancer_hotspot_new/website/website/input/specific_uniprots.txt"):
    """
    [NOTE] Expression cutoff: >=1 FPKM. Sample fraction: >= 80%.
    """
    expr_prot = set(pd.read_csv(specific_uniprot_file, sep = "\t", header = None)[0].tolist())
    expr_nodes = []
    for x in cancers:
        df = pd.read_csv("/local/storage/yz2296/TCGA_expression/parsed_files/" + x + "_normal.txt", sep = "\t", dtype = str)
        expr_nodes.extend(df[df["Expression"].apply(get_expr_genes)]["ENSG"].tolist())
    expr_nodes = set(expr_nodes)
    ensg2uniprot = all_function_py3.gene_ID_conversion(expr_nodes, convert_from = "ensembl.gene", convert_to = "uniprot", species = "human")
    for x in expr_nodes:
        expr_prot = expr_prot.union(ensg2uniprot[x])
    f = open(output, "w")
    f.write("\n".join(sorted(expr_prot)) + "\n")
    f.close()
    return

file_path = "/local/storage/yz2296/cancer_hotspot_new/website/website/input/"
    
# load raw mutation data
df = pd.read_csv(file_path + "pancan_unfiltered.maf", sep = "\t", dtype = {"Tumor_Sample_Barcode": str, "Chromosome": str}).sort_values(by = ["Chromosome", "Start_Position"])
df = df[["Center", "NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Strand", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2", "TUMORTYPE"]]
df["#Uploaded_variation"] = df.apply(lambda x: x["Chromosome"] + "_" + str(x["Start_Position"]) + "_" + x["Reference_Allele"] + "/" + x["Tumor_Seq_Allele2"], axis = 1)

# relabel tumors
index = (df["TUMORTYPE"] == "coadread") & df["Tumor_Sample_Barcode"].apply(lambda x: x[0:4] == "TCGA")
df_coadread = df[index].drop(columns = ["TUMORTYPE"])
df_other = df[~index]
df_coadread["case_submitter_id"] = df_coadread["Tumor_Sample_Barcode"].apply(lambda x: "-".join(x.split("-")[0:3]))
df_map = pd.read_csv(file_path + "coadread_samples.tsv", sep = "\t")
df_coadread = pd.merge(df_coadread, df_map, on = "case_submitter_id").drop(columns = ["case_submitter_id"])
df = pd.concat([df_coadread, df_other])

# # # prepare VEP input
# # df["vep_input"] = df.apply(lambda x: str(x["Chromosome"]) + " " + str(x["Start_Position"]) + " " + str(x["End_Position"]) + " " + x["Reference_Allele"] + "/" + x["Tumor_Seq_Allele2"], axis = 1)
# # df[["vep_input"]].drop_duplicates().to_csv(file_path + "VEP_input.txt", sep = "\t", header = None, index = None)

# load VEP output (one consequence per variant)
df_vep = pd.read_csv(file_path + "VEP_output.txt", sep = "\t", na_values = {"", "-"})
df_vep = df_vep[["#Uploaded_variation", "Consequence", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "SYMBOL_SOURCE", "HGNC_ID", "CCDS", "ENSP", "SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM", "SIFT", "PolyPhen", "AF", "AA_AF", "EA_AF", "gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF", "gnomAD_FIN_AF", "gnomAD_NFE_AF", "gnomAD_OTH_AF", "gnomAD_SAS_AF", "MPC", "CADD_PHRED", "CADD_RAW"]]
df_vep["Variant_Classification"] = df_vep["Consequence"].apply(all_function_py3.consequence2varclass_coding)

# remove germline variants
af_cols = ["AA_AF", "EA_AF", "gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF", "gnomAD_FIN_AF", "gnomAD_NFE_AF", "gnomAD_OTH_AF", "gnomAD_SAS_AF"]
df_vep = df_vep[~(df_vep[af_cols].max(axis = 1, skipna = True) > 0)]

# get protein altering mutations
df_lof = df_vep[df_vep["Variant_Classification"].apply(lambda x: all_function_py3.whether_muttype(x, muttype = "LOF"))]
df_nt = df_vep[df_vep["Variant_Classification"].apply(lambda x: all_function_py3.whether_muttype(x, "NT"))]
df_vep = pd.concat([df_lof, df_nt])

# filter out mutations in unexpressed genes
df = pd.merge(df, df_vep, on = "#Uploaded_variation")
specific_genes = ["AR", "CDH4", "EGFR", "EPHA3", "ERBB4", "FGFR2", "FLT3", "FOXA1", "FOXA2", "MECOM", "MIR142", "MSH4", "PDGFRA", "SOX1", "SOX9", "SOX17", "TBX3", "WT1"]
df1 = df[df["SYMBOL"].isin(specific_genes)]
df2 = df[~df["SYMBOL"].isin(specific_genes)]
df_final = []
for key,df_one in df2.groupby("TUMORTYPE"):
    file = "/local/storage/yz2296/TCGA_expression/parsed_files/" + key + "_normal.txt"
    if os.path.exists(file):
        print("####################", key)
        df_normal = pd.read_csv(file, sep = "\t", dtype = str)
        expr_genes = set(df_normal[df_normal["Expression"].apply(get_expr_genes)]["ENSG"].tolist())
        df_final.append(df_one[df_one["Gene"].isin(expr_genes)])
    else:
        df_final.append(df_one)
        print(key)
df_final = pd.concat(df_final + [df1])

# prepare pipeline input
cancers = {'blca', 'brca', 'cesc', 'coad', 'read', 'gbm', 'hnsc', 'kich', 'kirc', 'kirp', 'lihc', 'luad', 'lusc', 'paad', 'prad', 'skcm', 'stad', 'thca', 'ucec'}
df_pancan = df_final[df_final["TUMORTYPE"].isin(cancers) & df_final["Tumor_Sample_Barcode"].apply(lambda x: x[0:4] == "TCGA")].drop(columns = ["#Uploaded_variation", "UNIPROT_ISOFORM", "AF", "AA_AF","EA_AF", "gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF", "gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF", "CADD_PHRED", "CADD_RAW", "MPC", "SWISSPROT", "TREMBL"])
df_pancan["HGVSc"] = df_pancan["HGVSc"].apply(convert_hgvs)
df_pancan["HGVSp"] = df_pancan["HGVSp"].apply(convert_hgvs)
df_pancan["Hugo_Symbol"] = df_pancan["SYMBOL"]
df_pancan.to_csv(file_path + "TCGA_" + str(len(cancers)) + ".maf", sep = "\t", header = True, index = None)
get_expressed_PPInetwork(cancers, file_path + "ExprGen_" + str(len(cancers)) + ".txt")
for x in cancers:
    df_pancan[df_pancan["TUMORTYPE"] == x].to_csv(file_path + "TCGA_" + x + ".maf", sep = "\t", header = True, index = None)
    get_expressed_PPInetwork({x}, file_path + "ExprGen_" + x + ".txt")
