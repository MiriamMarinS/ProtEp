import pandas as pd
import numpy as np
import collections
from source.utils import occurrences, occurrences_mismatch

# For analyse peptides from proteomic files.

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
                    if self.mismatch == 0:
                        if epitope_sequence_protein[0] in peptide.seq:
                            number_epitope_in_seq = occurrences(peptide.seq.lower().strip(), epitope_sequence_protein[0].lower().strip())
                            number_epitope_in_peptides = number_epitope_in_peptides + number_epitope_in_seq
                    else:
                        number_epitope_in_seq = occurrences_mismatch(peptide.seq.lower().strip(), epitope_sequence_protein[0].lower().strip(), self.mismatch, self.epitope_library)
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
            self.new_dict_annotation_name_epitope_table = self.protein_type(epitope_name, epitope)

        return(self.new_dict_annotation_name_epitope_table, self.NGP_list)
    
    def protein_type(self, epitope_name, epitope):
        alpha, gamma, omega, HMW, LMW, NGP = 0, 0, 0, 0, 0, 0

        for k, v in self.epitope_annotation.items():
            for epitope_from_dict in v:
                if epitope_from_dict.name == epitope_name:
                    if any(specie in k.lower() for specie in ["triticum aestivum", "triticum turgidum", "hordeum"]):
                        if not any(AccNum in k.lower() for AccNum in self.AccNum_list):
                            door = "open"
                            if all(word in k.lower() for word in self.alpha_list) or all(word in k.lower() for word in self.alfa_list):
                                alpha = alpha + epitope_from_dict.number_epitope_in_peptides
                                door = "close"
                            elif all(word in k.lower() for word in self.gamma_list):
                                gamma = gamma + epitope_from_dict.number_epitope_in_peptides
                                door = "close"
                            elif all(word in k.lower() for word in self.omega_list):
                                omega = omega + epitope_from_dict.number_epitope_in_peptides
                                door = "close"
                            elif any(word in k.lower() for word in self.HMW_list):
                                HMW = HMW + epitope_from_dict.number_epitope_in_peptides
                                door = "close"
                            elif any(word in k.lower() for word in self.LMW_list):
                                LMW = LMW + epitope_from_dict.number_epitope_in_peptides
                                door = "close"
                            if door == "open":
                                #There are unclassifed gliadin and glutenin proteins, named only as gliadin and glutenin.
                                NGP = NGP + epitope_from_dict.number_epitope_in_peptides
                                if epitope_from_dict.number_epitope_in_peptides != 0:
                                    self.NGP_list.append(k)

        self.new_dict_annotation_name_epitope_table.setdefault(epitope_name, []).append([self.genotype, self.enzyme, self.type, alpha, gamma, omega, HMW, LMW, NGP])

        return(self.new_dict_annotation_name_epitope_table)