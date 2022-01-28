# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 09:16:45 2021

@author: milil
"""

# For analysis of proteomic data of wheat samples, find canonical epitopes and monoclonal antibody recognition sites.

from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
import argparse
import collections
from Bio import SeqIO
from itertools import product

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--path_proteomic_files", help = "/path/to/ proteomic files directory.")
parser.add_argument("-o", "--path_output", help = "/path/to/ output files directory.")
parser.add_argument("-e", "--epitopes", help = "/path/to/epitope fasta file.")
parser.add_argument("-m", "--method", help = "Method is the type of epitopes. Example: epitopes with original sequences, method = original. Options = original|deaminated1|deaminated2|p31|wage|moAb|ias|all.")
parser.add_argument("-d", "--remove_epitope_duplication", help = "Remove epitope with duplicated sequences. Optios = yes | no.")
parser.add_argument("-mm", "--mismatch", help = "Number of allowed mismatches to find epitopes. Optios = 0 | 1 or more.")
parser.add_argument("-t", "--type_table", help = "Table by epitope name: g12_1, a1_1, ..., or by epitope type: g12, a1.... Optios = name | type.")
args = parser.parse_args()

proteomic_files_path = args.path_proteomic_files
proteomic_files_path_results = args.path_output
epitope_file = args.epitopes
method = args.method
remove_epitope_duplication = args.remove_epitope_duplication
mismatch = args.mismatch
type_table = args.type_table

#%%Functions and classes:

#For epitopes with original sequence. Method is the type of epitope: with original sequence, deaminated1 or deaminated2 sequence, p31 sequence, wage sequence, moAb sequence or ias sequence.
class Epitope():
    def __init__(self, epitope_file, method):
        self.epitope_file = epitope_file
        self.method = method
        self.library = self.epitope_library()

    def epitope_library(self):
        self.epitope_sequences = {}
        type_epitopes = [("a", "alpha"), ("g", "gamma"), ("w", "omega"), ("L", "LMW"), ("H", "HMW"), ("hor", "hordein"), ("sec", "secalin"), ("ave", "avenin"), ("p31", "p31"), ("a1", "a1"), ("g12", "g12"), ("r5", "r5")]
        for record in SeqIO.parse(self.epitope_file, "fasta"):
            if str(record.id).split(";")[1] == self.method:
                if self.method == "original" or self.method == "deaminated1" or self.method == "deaminated2":
                    if "glia" in str(record.id) or "glut" in str(record.id):
                        print("glia or glut")
                        print(str(record.id))
                        epitope = [type_epitope[1] for type_epitope in type_epitopes if str(record.id).split(";")[0].split("_")[2].startswith(type_epitope[0])][0]
                        print("Done")
                    elif any(i in str(record.id) for i in ["hor", "sec", "ave"]):
                        print("hor u otro")
                        print(str(record.id))
                        epitope = [type_epitope[1] for type_epitope in type_epitopes if str(record.id).split(";")[0].split("_")[1].startswith(type_epitope[0])][0]
                    else:
                        print(str(record.id))
                        epitope = [type_epitope[1] for type_epitope in type_epitopes if str(record.id).split(";")[0].split("_")[0].startswith(type_epitope[0])][0]

                elif self.method == "moAb":
                    epitope = [type_epitope[1] for type_epitope in type_epitopes if str(record.id).split(";")[0].split("_")[0].strip() == type_epitope[0]][0]
                else:
                    epitope = "non_classified"
                self.epitope_sequences[str(record.id).split(";")[0]] = (str(record.seq), epitope)

        return(self.epitope_sequences)

    def epitope_library_remove_duplication(self):
        #In the case of original, deaminated1, deaminated2 and moAb epitopes methods, the epitopes are grouped in classes. If epitopes in the same class have the same sequences, they will be considered as one epitope. In other methods, the duplication will be removed in the group of all epitopes, because they are not grouped in classes.

        #Find duplicated sequences in epitope library.
        duplicated_sequences = {}
        for epitope_id, value in self.library.items():
            epitope_seq = value[0].strip()
            epitope_type = value[1]
            duplicated_sequences.setdefault((epitope_seq, epitope_type), []).append(epitope_id)

        #New dictionary with epitopes grouped by their sequences separated in original, deaminated1 and deaminated2 groups, or unclassified group for the other methods.
        self.new_epitope_library = {}
        for k, v in duplicated_sequences.items():
            self.new_epitope_library[';'.join(v)] = k


        #The name of epitopes with duplicated sequences are concatenated by ;. But the name is too long. So we kept the first epitope name, and it can be checked in epitope library printed.
        self.epitope_sequences2 = {}
        for k, v in self.new_epitope_library.items():
            if len(k.split(";")) > 1:
                self.epitope_sequences2[k.split(";")[0]] = v
            else:
                self.epitope_sequences2[k] = v


        return(self.epitope_sequences2)

def protein_type(epitope_name, epitope, epitope_annotation, new_dict_annotation, alpha_list, alfa_list, gamma_list, omega_list, HMW_list, LMW_list, NGP_list, AccNum_list, genotype, enzyme, flour):
    alpha = 0
    gamma = 0
    omega = 0
    HMW = 0
    LMW = 0
    NGP = 0

    for k, v in epitope_annotation.items():
        for epitope_from_dict in v:
            if epitope_from_dict.name == epitope_name:
                if any(specie in k.lower() for specie in ["triticum aestivum", "triticum turgidum", "hordeum"]):
                    if not any(AccNum in k.lower() for AccNum in AccNum_list):
                        door = "open"
                        if all(word in k.lower() for word in alpha_list) or all(word in k.lower() for word in alfa_list):
                            alpha = alpha + epitope_from_dict.number_epitope_in_peptides
                            door = "close"
                        elif all(word in k.lower() for word in gamma_list):
                            gamma = gamma + epitope_from_dict.number_epitope_in_peptides
                            door = "close"
                        elif all(word in k.lower() for word in omega_list):
                            omega = omega + epitope_from_dict.number_epitope_in_peptides
                            door = "close"
                        elif any(word in k.lower() for word in HMW_list):
                            HMW = HMW + epitope_from_dict.number_epitope_in_peptides
                            door = "close"
                        elif any(word in k.lower() for word in LMW_list):
                            LMW = LMW + epitope_from_dict.number_epitope_in_peptides
                            door = "close"
                        if door == "open":
                            #There are unclassifed gliadin and glutenin proteins, named only as gliadin and glutenin.
                            NGP = NGP + epitope_from_dict.number_epitope_in_peptides
                            if epitope_from_dict.number_epitope_in_peptides != 0:
                                NGP_list.append(k)

    new_dict_annotation.setdefault(epitope_name, []).append([genotype, enzyme, flour, alpha, gamma, omega, HMW, LMW, NGP])

    return(new_dict_annotation)

class Proteomic():
    def __init__(self, proteomic_file_name, proteomic_files_path, epitope_library, mismatch):
        self.proteomic_file_name = proteomic_file_name
        self.proteomic_files_path = proteomic_files_path
        self.epitope_library = epitope_library
        self.read_proteomic_file()
        self.fill_NCBI_Acc()
        self.significant_peptides()
        self.sample()
        self.peptides_per_NCBI_Acc()
        self.mismatch = mismatch

    def read_proteomic_file(self):
        self.proteomic_file = self.proteomic_files_path + self.proteomic_file_name
        self.proteomic = pd.read_excel(self.proteomic_file, sheet_name = "Peptide_List")

    #Fill cells of NCBI Acc without NaN values. Keep value of upper row with NACBI Acc value.
    def fill_NCBI_Acc(self):
        self.proteomic["NCBI_Acc"] = self.proteomic["NCBI_Acc"].fillna(method='ffill')
        self.proteomic["Description"] = self.proteomic["Description"].fillna(method='ffill')

     #Keep significant peptides: P-value < 0.05.
    def significant_peptides(self):
        self.proteomic_significant_peptides = self.proteomic[self.proteomic["Pep_expect"] < 0.05]

    def sample(self):
        self.genotype = self.proteomic["Organism"][0]
        self.enzyme = self.proteomic["Enzyme"][0]
        self.type = self.proteomic["Type"][0]

    def peptides_per_NCBI_Acc(self):
        NCBI = ""
        self.dict_peptides = {}
        list_peptides = []
        Peptide = collections.namedtuple("Peptide", ("seq", "Modifications", "PeptScore", "mz", "MrDaexp", "zCharge", "evalue"))
        for index, peptide in self.proteomic_significant_peptides.iterrows():
            NCBI = peptide["NCBI_Acc"]
            pept = (peptide["Peptide_Sequence"], peptide["Modifications"], peptide["Pep_Score"], peptide["m_z"], peptide["Mr_Da"], peptide["z_charge"], peptide["Pep_expect"])
            if pept not in list_peptides:
                if any(specie in peptide["Description"].lower() for specie in ["triticum aestivum", "triticum turgidum", "hordeum"]):
                    self.dict_peptides.setdefault(NCBI, []).append(Peptide(seq=peptide["Peptide_Sequence"], Modifications=peptide["Modifications"], PeptScore=peptide["Pep_Score"], mz=peptide["m_z"], MrDaexp=peptide["Mr_Da"], zCharge=peptide["z_charge"], evalue=peptide["Pep_expect"]))
                    list_peptides.append(pept)

    #Match epitope sequences in peptides, match perfect.
    def epitope(self):
        self.epitopes = {}
        Epitope = collections.namedtuple("Epitope", ("name", "protein", "number_epitope_in_peptides"))
        for k, v in self.dict_peptides.items():
            for epitope_name, epitope_sequence_protein in self.epitope_library.items():
                number_epitope_in_peptides = 0
                for peptide in v:
                    if mismatch == 0:
                        if epitope_sequence_protein[0] in peptide.seq:
                            number_epitope_in_seq = occurrences(peptide.seq.lower().strip(), epitope_sequence_protein[0].lower().strip())
                            number_epitope_in_peptides = number_epitope_in_peptides + number_epitope_in_seq
                    else:
                        number_epitope_in_seq = occurrences_mismatch(peptide.seq.lower().strip(), epitope_sequence_protein[0].lower().strip(), mismatch, self.epitope_library)
                        number_epitope_in_peptides = number_epitope_in_peptides + number_epitope_in_seq
                self.epitopes.setdefault(k, []).append(Epitope(name=epitope_name, protein=epitope_sequence_protein[1],number_epitope_in_peptides=number_epitope_in_peptides))

    def annotation(self):
        self.proteomic_annotation = pd.read_excel(self.proteomic_file, sheet_name = "Protein_List")
        self.NCBI_Acc_annotation = {}
        self.AccNum_list = []
        for index, protein in self.proteomic_annotation.iterrows():
            self.NCBI_Acc_annotation[protein["NCBI_Acc"]] = protein["Description"]
            if protein["Num_Pept"] == 1:
                if protein["NCBI_Acc"].lower().strip() not in self.AccNum_list:
                    self.AccNum_list.append(protein["NCBI_Acc"].lower().strip())
        self.epitope_annotation = {}
        for k, v in self.epitopes.items():
            self.epitope_annotation[k + ";" + self.NCBI_Acc_annotation[k]] = v

    def type_representation_(self, new_dict_annotation_name_epitope_table):
        #Generate list of proteins (alpha, gamma, ...) or species for representation.
        self.alpha_list = ["alpha", "gliadin"]
        self.alfa_list = ["alfa", "gliadin"]
        self.gamma_list = ["gamma", "gliadin"]
        self.omega_list = ["omega", "gliadin"]
        self.HMW_list = ["high molecular weight", "hmw", "high-molecular-weight"]
        self.LMW_list = ["low molecular weight", "lmw", "low-molecular-weight"]

        self.list_epitopes = []
        #List of type of protein in which epitopes are: belongs to alpha, gamma, omega, .... Only for original, deminated1 and deaminated2 methods.
        self.list_proteins = []
        for k, v in self.epitope_annotation.items():
            for epitope in v:
                if (epitope.name, epitope.protein) not in self.list_epitopes:
                    self.list_epitopes.append((epitope.name, epitope.protein))
                if epitope.protein not in self.list_proteins:
                    self.list_proteins.append(epitope.protein)

        #For representatino of epitopes for genotype and enzyme for each protein or species.
        self.NGP_list = []
        self.new_dict_annotation_name_epitope_table = new_dict_annotation_name_epitope_table
        for epitope in self.list_epitopes:
            epitope_name = epitope[0]
            self.new_dict_annotation_name_epitope_table = protein_type(epitope_name, epitope, self.epitope_annotation, self.new_dict_annotation_name_epitope_table, self.alpha_list, self.alfa_list, self.gamma_list, self.omega_list, self.HMW_list, self.LMW_list, self.NGP_list, self.AccNum_list, self.genotype, self.enzyme, self.type)

        return(self.new_dict_annotation_name_epitope_table, self.NGP_list)

def writeTable(dict_table, type_table, list_of_samples, mismatch, proteomic_files_path_results):
    proteins = {'alpha': 3, 'gamma': 4, 'omega': 5, 'HMW': 6, 'LMW': 7, 'NGP': 8}

    #Table to export.
    dict_to_table = {}
    for protein, value in proteins.items():
        for epitope, samples_types in dict_table.items():
            name = epitope
            included_samples = []
            dict_to_table.setdefault("protein", []).append(protein)
            dict_to_table.setdefault("epitope", []).append(name)
            for sample_type in samples_types:
                genotype = sample_type[0]
                enzyme = sample_type[1]
                flour = sample_type[2]
                dict_to_table.setdefault(genotype + " " + enzyme + " " + flour, []).append(sample_type[proteins[protein]])
                if (genotype, enzyme, flour) not in included_samples:
                    included_samples.append((genotype, enzyme, flour))

            for i in included_samples:
                if i not in list_of_samples:
                    genotype = i[0]
                    enzyme = i[1]
                    flour = i[2]
                    dict_to_table.setdefault(i[0] + " " + i[1] + " " + i[2], []).append(0)

    table_df = pd.DataFrame(dict_to_table)
    if type_table == "name":
        table_df_write = table_df
    else:
        proteins_epitopes = list(product(list(proteins.keys()), set(v[1] for k, v in epitope_sequences.items())))
        table_df_write = pd.DataFrame(index = range(0, len(proteins_epitopes)), columns = list(table_df.columns))

        for index, row in table_df_write.iterrows():
            table_df_write.loc[index]["protein"] = proteins_epitopes[index][0]
            table_df_write.loc[index]["epitope"] = proteins_epitopes[index][1]
            epitope_names = list(k for k, v in epitope_sequences.items() if v[1] == proteins_epitopes[index][1])
            for column in list(table_df.columns)[2:]:
                table_df_write.loc[index][column] = table_df.loc[(table_df['protein'] == proteins_epitopes[index][0]) & (table_df['epitope'].isin(epitope_names)), column].sum()

    table_df_write.to_csv(proteomic_files_path_results + 'table_epitopes_' + method + '_' + type_table + '_' + str(mismatch) + 'mismatch.txt', header=True, index=False, sep='\t', mode='a')


def occurrences(string, sub):
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count


def occurrences_mismatch(string, sub, mismatch_allowed, epitope_sequences):

    epitopes = []
    for k, v in epitope_sequences.items():
        if v[0] not in epitopes:
            epitopes.append(v[0])

    count = 0
    for i in range(0, len(string) - len(sub) + 1):
        mismatch = 0
        amplicon = string[i:i+len(sub)]

        #if any(epitope == amplicon for epitope in epitopes):
        #    continue
        for ii in range(0, len(sub)):
            if amplicon[ii] != sub[ii]:
                mismatch += 1
        if mismatch == mismatch_allowed:
            count += 1
    return(count)


#%% Processing file of epitopes.
print("You are working with " + method + " epitope sequences.")
epitope = Epitope(epitope_file, method)
epitope_sequences = epitope.library

#Modify library for working without duplications.
#If epitope of the same type (alpha, gamma, omega, ...) have the same sequence, the duplication will be removed.
if remove_epitope_duplication == "yes":
    epitope.epitope_library_remove_duplication()
    epitope_sequences = epitope.epitope_sequences2

#%%Find epitopes in peptide sequences.
proteomic_files = [f for f in listdir(proteomic_files_path) if isfile(join(proteomic_files_path, f))] # List with all files in directory.
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

#%%Table of epitopes.
writeTable(new_dict_annotation_name_epitope_table, type_table, list_of_samples, mismatch, proteomic_files_path_results)
