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
print(f'Total of patients with RNA data: {len(RNA_sheet.columns)-1}')


###### CLINICAL FILE ######

clinic_sheet.columns = clinic_sheet.iloc[0]
clinic_sheet = clinic_sheet.drop([0], axis=0)
print(f'Total of patients with clinic data: {len(clinic_sheet)}  \n')
# clinic_sheet.to_csv("Preprocessed_data/clinic_data.csv")


###### Compare patients in RNA and CLINICAL FILE ######

rna_patients = clinic_sheet[~clinic_sheet['ID'].isin(RNA_sheet.columns.values)]
ax_patients = len(rna_patients[rna_patients['TRT01P']=='Avelumab+Axitinib'])
suni_patients = len(rna_patients[rna_patients['TRT01P']=='Sunitinib'])
print(f'Patients who did not undergo the RNA study: {len(rna_patients)}')
print(f'Patients who did not undergo the RNA study from Avelumab+Axitinib arm: {ax_patients}')
print(f'Patients who did not undergo the RNA study from Sunitinib arm: {suni_patients}  \n')

rna_patients = clinic_sheet[clinic_sheet['ID'].isin(RNA_sheet.columns.values)]
ax_patients = len(rna_patients[rna_patients['TRT01P']=='Avelumab+Axitinib'])
suni_patients = len(rna_patients[rna_patients['TRT01P']=='Sunitinib'])
print(f'Patients who underwent the RNA study: {len(rna_patients)}')
print(f'Patients who underwent the RNA study from Avelumab+Axitinib arm: {ax_patients}')
print(f'Patients who underwent the RNA study from Sunitinib arm: {suni_patients} \n')


###### Compare genes from current and old RNA FILES ######

old_rna_data = pd.read_csv('../../2023 Work/Data/Preprocessed_data/clinic_and_RNA_data_raw_NIVOLUMAB.csv', low_memory=False)
not_included_genes = RNA_sheet[~RNA_sheet['HUGO'].isin(old_rna_data.columns.values)]['HUGO']

print(f'New genes which were not include in last RNA file: {len(not_included_genes)}')
# print(not_included_genes)


#### FILTER CLINICAL FILE (only Avelumab+Axitinib patients) ####

clinic_sheet_Ave_Axi = clinic_sheet[(clinic_sheet['TRT01P']=='Avelumab+Axitinib')]
clinic_sheet_Sunitinib = clinic_sheet[(clinic_sheet['TRT01P']=='Sunitinib')]


###### Create combined files ######

RNA_sheet = RNA_sheet.set_index("HUGO")
RNA_sheet = RNA_sheet.T
print(RNA_sheet.head())
RNA_sheet.insert(0, "ID", RNA_sheet.index, True)
print(RNA_sheet.head())
# RNA_sheet.to_csv("Preprocessed_data/RNA_data.csv")

merged_df_RNA_clinic_Ave_Axi = pd.merge(clinic_sheet_Ave_Axi, RNA_sheet, on='ID')
# merged_df_RNA_clinic_Ave_Axi.to_csv("Preprocessed_data/clinic_and_RNA_data_raw_Ave_Axi.csv")

merged_df_RNA_clinic_Sunitinib = pd.merge(clinic_sheet_Sunitinib, RNA_sheet, on='ID')
merged_df_RNA_clinic_Sunitinib.to_csv("Preprocessed_data/clinic_and_RNA_data_raw_Sunitinib.csv")
