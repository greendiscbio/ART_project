############### GIVEN INITIAL EXCEL WITH PLENTY SHEETS ###############
import pandas as pd
import os
os.chdir("Sunitinib_Avelumab+Axitinib/Data")

# Unpack sheets from excel file
complete_file = pd.ExcelFile('41591_2020_1044_MOESM3_ESM.xlsx')
RNA_sheet = pd.read_excel(complete_file, 'S13_Gene_expression_TPM')
clinic_sheet = pd.read_excel(complete_file, 'S11_Clinical_data', header=0)

###### RNA FILE ######
RNA_sheet.columns = RNA_sheet.iloc[0]
RNA_sheet = RNA_sheet.drop([0], axis=0)
print(RNA_sheet.head()) 
print('Total of patients with RNA data: '+str(len(RNA_sheet.columns)))

###### CLINICAL FILE ######
clinic_sheet.columns = clinic_sheet.iloc[0]
clinic_sheet = clinic_sheet.drop([0], axis=0)
print('Total of patients with clinic data: '+str(len(clinic_sheet)))
clinic_sheet.to_csv("Preprocessed_data/clinic_data.csv")

###### Compare patients in RNA and CLINICAL FILE ######
not_included = 0
ax_patients = 0
for patient in range(len(clinic_sheet)):
    if clinic_sheet['ID'].iloc[patient] not in RNA_sheet.columns:
        not_included = not_included+1
        if clinic_sheet['TRT01P'].iloc[patient] == 'Avelumab+Axitinib':
            ax_patients = ax_patients+1

print('Patients who did not undergo the RNA study: '+str(not_included))
print('Patients who did not undergo the RNA study from Avelumab+Axitinib arm: '+str(ax_patients))
print('Patients who did not undergo the RNA study from Sunitinib arm: '+str(not_included-ax_patients))

###### Compare genes from current and old RNA FILES ######
old_rna_data = pd.read_csv('../../2023 Work/Data/Preprocessed_data/clinic_and_RNA_data_raw_NIVOLUMAB.csv', low_memory=False)
not_included = 0
not_included_genes = []
for gene in range(len(RNA_sheet)):
    if RNA_sheet['HUGO'].iloc[gene] not in old_rna_data.columns:
        not_included = not_included+1
        not_included_genes.append(RNA_sheet['HUGO'].iloc[gene])

print('New genes which were not include in last RNA file: '+str(not_included))
print(not_included_genes)


#### FILTER CLINICAL FILE (only Avelumab+Axitinib patients) ####
clinic_sheet_Ave_Axi = clinic_sheet[(clinic_sheet['TRT01P']=='Avelumab+Axitinib')]

###### Create combined files ######
RNA_sheet = RNA_sheet.set_index("HUGO")
RNA_sheet = RNA_sheet.T
print(RNA_sheet.head())
RNA_sheet.insert(0, "ID", RNA_sheet.index, True)
print(RNA_sheet.head())
RNA_sheet.to_csv("Preprocessed_data/RNA_data.csv")

merged_df_RNA_clinic_Ave_Axi = pd.merge(clinic_sheet_Ave_Axi, RNA_sheet, on='ID')
merged_df_RNA_clinic_Ave_Axi.to_csv("Preprocessed_data/clinic_and_RNA_data_raw_Ave_Axi.csv")
