# -*- coding: utf-8 -*-
"""
Created on Mon May 31 12:50:32 2021

@author: milil
"""

#Proteomic data.

from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
import argparse

# Arguments.
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--proteomic_dir', help = "/path/to/ proteomic files directory")
parser.add_argument('-o', '--output', help = "/path/to/ output directory")

args = parser.parse_args()

# Paths.
proteomic_files_path = args.proteomic_dir # Path to proteomic files.
proteomic_result_path = args.output # Path of output directory.

#%% Classes and Functions:

def peptides(pept, type_protein, values, door, dict_proteins_result, protein):
    if type_protein in ["omega-gliadin", "gamma-gliadin", "ATIs", "LTP"]:
       if all(value in protein["Description"].lower() for value in values):
           dict_proteins_result[type_protein] = dict_proteins_result.get(type_protein, 0) + 1
           door = "close"
    if type_protein == "alpha-gliadin":
        if all(value in protein["Description"].lower() for value in values[0]) or all(value in protein["Description"] for value in values[1]):
            dict_proteins_result[type_protein] = dict_proteins_result.get(type_protein, 0) + 1
            door = "close"
    if type_protein in ["HMW", "LMW"]:
        if any(value in protein["Description"].lower() for value in values):
            dict_proteins_result[type_protein] = dict_proteins_result.get(type_protein, 0) + 1
            door = "close"
    if type_protein in ["Globulin", "Triticin", "Serpin", "Avenin", "Hordein", "Secalin"]:
        if values in protein["Description"].lower():
            dict_proteins_result[type_protein] = dict_proteins_result.get(type_protein, 0) + 1
            door = "close"
    return(door)

#%%Read proteomic files.
proteomic_files = [f for f in listdir(proteomic_files_path) if isfile(join(proteomic_files_path, f))] # List with all files in directory.

#%%
# Word keys to search protein types in proteomic files.
# A previous work has to be done to identify all the keys for each type of protein.
# In other protein type will be classified all the proteins that did  not match with the keys.
dict_proteins = {"omega-gliadin": ("omega", "gliadin"),
                 "alpha-gliadin": (("alpha", "gliadin"), ("alfa", "gliadin")),
                 "gamma-gliadin": ("gamma", "gliadin"),
                 "HMW": ("high molecular weight", "high-molecular-weight", "hmw", "glu-a1x"),
                 "LMW": ("low molecular weight", "low-molecular-weight", "lmw"),
                 "ATIs": ("alpha", "amylase", "inhibitor"),
                 "LTP": ("lipid", "transfer", "protein"),
                 "Globulin": ("globulin"),
                 "Triticin": ("triticin"),
                 "Serpin": ("serpin"),
                 "Avenin": ("avenin"),
                 "Hordein": ("hordein"),
                 "Secalin": ("secalin"),
                 "Others": ()}


# Search in Triticum aestivum, Triticum durum and Hordeum vulgare.
# In this study, only proteins from these species are of interest.
species = ["triticum aestivum", "triticum turgidum", "hordeum"]

# Read proteomic files in xlsx format and count the number of unique peptides for each type of protein.
proteomic_dict = {} # Dict for generate the dataframe with peptide counts for each type of protein.
proteomic_other_dict = {} # Dict with the names of proteins considered as "Others".

for excel_proteomic_data in proteomic_files:
    NCBI = pd.read_excel(proteomic_files_path + excel_proteomic_data, sheet_name = "Protein_List")
    proteomic = pd.read_excel(proteomic_files_path + excel_proteomic_data, sheet_name = "Peptide_List")
    dict_proteins_result = {} # Dict with results for each proteomic file.
    dict_proteins_others = {} # Dict with results for other proteins that are not included in the search process in each proteomic file.
    NCBI_Acc = [] # List of NCBI Acc supported by only 1 peptide.
    for index, protein in NCBI.iterrows():
        if protein["Num_Pept"] != 1:
            NCBI_Acc.append(protein["NCBI_Acc"])
    list_peptides = []
    proteomic["NCBI_Acc"] = proteomic["NCBI_Acc"].fillna(method='ffill')
    proteomic["Description"] = proteomic["Description"].fillna(method='ffill')
    for index, protein in proteomic.iterrows():
        pept = (protein["Peptide_Sequence"], protein["Modifications"], protein["Pep_Score"], protein["m_z"], protein["Mr_Da"], protein["z_charge"], protein["Pep_expect"])
        if pept not in list_peptides:
            door = "open" # Close when the peptide belongs to one of the protein types.
            if protein["Pep_expect"] < 0.05 and any(specie in protein["Description"].lower() for specie in species) and any(AccNum.lower() in protein["NCBI_Acc"].lower() for AccNum in NCBI_Acc):
                for type_protein, values in  dict_proteins.items():
                    door = peptides(pept, type_protein, values, door, dict_proteins_result, protein)
                if door == "open":
                    dict_proteins_result["Others"] = dict_proteins_result.get("Others", 0) + 1
                    dict_proteins_others[protein["Description"]] = dict_proteins_others.get(protein["Description"], 0) + 1
                list_peptides.append(pept)
    for protein in list(dict_proteins.keys()):
        if protein in dict_proteins_result:            
            proteomic_dict.setdefault(excel_proteomic_data.split(".")[0], []).append(dict_proteins_result[protein])
        else:
            proteomic_dict.setdefault(excel_proteomic_data.split(".")[0], []).append(0)
    proteomic_other_dict[excel_proteomic_data] = dict_proteins_others
    
proteomic_table = pd.DataFrame(proteomic_dict, index = [list(dict_proteins.keys())])
proteomic_table.to_csv(proteomic_result_path + "unique_peptides.txt", sep = "\t")