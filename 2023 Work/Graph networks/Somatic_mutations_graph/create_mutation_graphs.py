import pandas as pd
import numpy as np
import json
# Needed data
wes_nivolumab_file = pd.read_csv('2023 Work/Data/Preprocessed_data/clinic_and_RNA_and_WES_data_raw_NIVOLUMAB.csv')
wes_nivolumab_file = wes_nivolumab_file[['MAF_Tumor_ID','Hugo_Symbol']]
rna_clinic_nivolumab_file = pd.read_csv('2023 Work/Data/Preprocessed_data/clinic_and_RNA_data_raw_NIVOLUMAB.csv')

# Group patients by their id to get their mutated genes
uniques = wes_nivolumab_file['MAF_Tumor_ID'].value_counts()
filter = wes_nivolumab_file.groupby(['MAF_Tumor_ID'])

# Save each patient and their mutated genes in a dictionary to have better access
dict_genes={}
for index, value in uniques.items():
    dict_genes[index] = (filter.get_group(index)['Hugo_Symbol']).to_list()

with open('2023 Work/Graph networks/Somatic_mutations_graph/Data/genes_per_patient.json', 'w') as outfile:
    json.dump(dict_genes, outfile)

# Load dicctionary 
with open('2023 Work/Graph networks/Somatic_mutations_graph/Data/genes_per_patient.json', 'r') as json_file:
    dict_genes = json.load(json_file)

# Iterate patients to create a unique dataframe with all the mutation gene information
final_df = pd.DataFrame()
for patient in dict_genes:
    pat_genes = [] # Mutated genes of the current patient wich are inclded in RNA file (we know the sobreexpression value)
    not_include = 0 # Mutated genes not included in RNA file or are duplicated
    line_patient = rna_clinic_nivolumab_file[rna_clinic_nivolumab_file['MAF_Tumor_ID'] == patient] # Curent patient clinic and rna information
    for g in dict_genes[patient]:
        if g in rna_clinic_nivolumab_file.columns and g not in pat_genes: # Sobreexpression is known and it is not duplicated
            pat_genes.append(g)
        else:
            not_include = not_include + 1
    print('Genes not included for patient ' + str(patient) + ' : '+str(not_include)+' out of ' + str(len(dict_genes[patient])) + '.')
    # Select only mutated genes columns of rna and clinic dataframe
    line_patient = line_patient[pat_genes]
    # Save rna and maf id in order to identificate current patient
    line_patient.insert(0, "RNA_ID", rna_clinic_nivolumab_file[rna_clinic_nivolumab_file['MAF_Tumor_ID'] == patient]['RNA_ID'], True)
    line_patient.insert(0, "MAF_Tumor_ID", rna_clinic_nivolumab_file[rna_clinic_nivolumab_file['MAF_Tumor_ID'] == patient]['MAF_Tumor_ID'], True)
    # Add patient to the final dataframe
    final_df = pd.concat([final_df, line_patient], axis = 0)
# Save final dataframe
final_df.to_csv('2023 Work/Graph networks/Somatic_mutations_graph/Data/mutated_genes_per_patient.csv')