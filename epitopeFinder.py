# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 09:16:45 2021

@author: Miriam Marin Sanz
"""

# For analysis of proteomic data of wheat samples, find canonical epitopes and monoclonal antibody recognition sites.

from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
import argparse
from source.epitope import Epitope
from source.proteomic import Proteomic
from source.writeTable import writeTable

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--path_proteomic_files", help = "/path/to/ proteomic files directory.")
parser.add_argument("-o", "--path_output", help = "/path/to/ output files directory.", default="./")
parser.add_argument("-e", "--epitopes", help = "/path/to/epitope fasta file.")
parser.add_argument("-m", "--method", help = "Method is the type of epitopes. Example: epitopes with original sequences, method = original. Options = original|deaminated1|deaminated2|p31|wage|moAb|ias|all.")
parser.add_argument("-d", "--remove_epitope_duplication", help = "Remove epitope with duplicated sequences. Optios = yes | no.", default="yes")
parser.add_argument("-mm", "--mismatch", help = "Number of allowed mismatches to find epitopes. Optios = 0 | 1 or more.", type=int, default=0)
parser.add_argument("-t", "--type_table", help = "Table by epitope name: g12_1, a1_1, ..., or by epitope type: g12, a1.... Optios = name | type.")
args = parser.parse_args()

def main():
    proteomic_files_path = args.path_proteomic_files
    proteomic_files_path_results = args.path_output
    epitope_file = args.epitopes
    method = args.method
    remove_epitope_duplication = args.remove_epitope_duplication
    mismatch = args.mismatch
    type_table = args.type_table

    # Processing file of epitopes.
    print("You are working with " + method + " epitope sequences.")
    epitope = Epitope(epitope_file, method)
    epitope_sequences = epitope.library

    #Modify library for working without duplications.
    #If epitope of the same type (alpha, gamma, omega, ...) have the same sequence, the duplication will be removed.
    if remove_epitope_duplication == "yes":
        epitope.epitope_library_remove_duplication()
        epitope_sequences = epitope.epitope_sequences2

    # Find epitopes in peptide sequences.
    proteomic_files = [f for f in listdir(proteomic_files_path) if isfile(join(proteomic_files_path, f)) and f.endswith("xlsx")] # List with all files in directory.
    list_of_samples = [] #List of all posible genotype_enzyme combinations.
    new_dict_annotation_name_epitope_table = {} #Dict wiht all genotypes and enzymes for epitope name representation.
    NGP_list = [] #List of NGP with epitopes.

    for proteomic_file in proteomic_files:
        print("Processing file: ", proteomic_file)
        proteomic = Proteomic(proteomic_file, proteomic_files_path, epitope_sequences, mismatch)
        proteomic.epitope() # Find epitopes in peptide sequences.
        proteomic.annotation() # Annotate peptides.
        new_dict_annotation_name_epitope_table, NGP_list_with_epitopes = proteomic.type_representation_(new_dict_annotation_name_epitope_table)
        for NGP in NGP_list_with_epitopes:
            if NGP not in NGP_list:
                NGP_list.append(NGP) # Print for study the NGP with epitopes.
        if (proteomic.genotype, proteomic.enzyme, proteomic.type) not in list_of_samples:
            list_of_samples.append((proteomic.genotype, proteomic.enzyme, proteomic.type))

    writeTable(new_dict_annotation_name_epitope_table, type_table, list_of_samples, mismatch, proteomic_files_path_results, epitope_sequences, method)

if __name__ == "__main__":
    main()
