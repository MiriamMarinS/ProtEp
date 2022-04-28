# -*- coding: utf-8 -*-
"""
Created on Mon May 31 12:50:32 2021

@author: Miriam Marin Sanz
"""

#Proteomic data.
from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
import argparse
from source.peptidesv2 import peptidesv2
from progressbar import *

# Arguments.
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--proteomic_dir', help = "/path/to/ proteomic files directory")
parser.add_argument('-o', '--output', help = "/path/to/ output directory")
args = parser.parse_args()

widgets = ['Test: ', Percentage(), ' ', Bar(marker='0',left='[',right=']'),
           ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options

# Paths.
proteomicFilesPath = args.proteomic_dir # Path to proteomic files.
proteomicResultPath = args.output # Path of output directory.

def main():
    #Read proteomic files.
    proteomicFiles = [f for f in listdir(proteomicFilesPath) if isfile(join(proteomicFilesPath, f))] # List with all files in directory.

    # Word keys to search protein types in proteomic files.
    # A previous work has to be done to identify all the keys for each type of protein.
    # In other protein type will be classified all the proteins that did  not match with the keys.
    dictProteins = {"omega-gliadin": ("omega", "gliadin"),
                    "alpha-gliadin": (("alpha", "gliadin"), ("alfa", "gliadin")),
                    "gamma-gliadin": ("gamma", "gliadin"),
                    "HMW": ("high molecular weight", "high-molecular-weight", "hmw", "glu-a1x"),
                    "LMW": ("low molecular weight", "low-molecular-weight", "lmw"),
                    "ATIs": ("alpha", "amylase", "trypsin", "inhibitor"),
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
    proteomicDict = {} # Dict for generate the dataframe with peptide counts for each type of protein.
    proteomicOther = {} # Dict with the names of proteins considered as "Others".

    count = 0
    pbar = ProgressBar(widgets=widgets, maxval=len(proteomicFiles))
    pbar.start()

    for excel_proteomic_data in proteomicFiles:
        NCBI = pd.read_excel(proteomicFilesPath + excel_proteomic_data, sheet_name = "Protein-List")
        proteomic = pd.read_excel(proteomicFilesPath + excel_proteomic_data, sheet_name = "Peptide-List")
        proteinsResult = {} # Dict with results for each proteomic file.
        proteinsOthers = {} # Dict with results for other proteins that are not included in the search process in each proteomic file.
        NCBI_Acc = [] # List of NCBI Acc supported by only 1 peptide.
        NCBI_Acc_annotation = {} # Dict with description and NCBI_Acc.
        for index, protein in NCBI.iterrows():
            if protein["NCBI_Acc"] not in list(NCBI_Acc_annotation.keys()):
                NCBI_Acc_annotation[protein["NCBI_Acc"]] = protein["Description"]
            if protein["Num_Pept"] != 1:
                NCBI_Acc.append(protein["NCBI_Acc"])
        listPeptides = {}
        for type_protein, values in dictProteins.items():
            listPeptides[type_protein] = []
        proteomic["NCBI_Acc"] = proteomic["NCBI_Acc"].fillna(method='ffill')
        proteomic["Description"] = proteomic["Description"].fillna(method='ffill')
        for index, protein in proteomic.iterrows():
            pept = (protein["Peptide_Sequence"], protein["Modifications"], protein["Pep_Score"], protein["m_z"],
            protein["Mr_Da"], protein["z_Charge"], protein["Pep_expect"])
            door = "open" # Close when the peptide belongs to one of the protein types.
            if protein["Pep_expect"] < 0.05 and any(specie in protein["Description"].lower() for specie in species) and any(AccNum.lower() in protein["NCBI_Acc"].lower() for AccNum in NCBI_Acc):
                for type_protein, values in  dictProteins.items():
                    door = peptidesv2(pept, type_protein, values, door, proteinsResult, protein, NCBI_Acc_annotation, listPeptides)
                if door == "open":
                    if pept not in listPeptides["Others"]:
                        proteinsResult["Others"] = proteinsResult.get("Others", 0) + 1
                        proteinsOthers[NCBI_Acc_annotation[protein["NCBI_Acc"]]] = proteinsOthers.get(NCBI_Acc_annotation[protein["NCBI_Acc"]], 0) + 1
                        listPeptides["Others"].append(pept)
        for protein in list(dictProteins.keys()):
            if protein in proteinsResult:            
                proteomicDict.setdefault(excel_proteomic_data.split(".")[0], []).append(proteinsResult[protein])
            else:
                proteomicDict.setdefault(excel_proteomic_data.split(".")[0], []).append(0)
        proteomicOther[excel_proteomic_data] = proteinsOthers
        count += 1
        pbar.update(count)
    pbar.finish()
    proteomic_table = pd.DataFrame(proteomicDict, index = [list(dictProteins.keys())])
    proteomic_table.to_csv(proteomicResultPath + "unique_peptidesv2.txt", sep = "\t")

if __name__ == "__main__":
    main()
