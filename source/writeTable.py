import pandas as pd
import numpy as np
from itertools import product

def writeTable(dict_table, type_table, list_of_samples, mismatch, proteomic_files_path_results, epitope_sequences, method):
    proteins = {'alpha': 3, 'gamma': 4, 'omega': 5, 'HMW': 6, 'LMW': 7, 'NGP': 8}

    # Table to export.
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