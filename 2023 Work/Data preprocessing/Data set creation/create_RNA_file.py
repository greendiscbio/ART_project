############### GIVEN INITIAL EXCEL WITH PLENTY SHEETS ############################################
import pandas as pd

# Unpack sheets from excel file
complete_file = pd.ExcelFile('2023 Work/Data/41591_2020_839_MOESM2_ESM.xlsx')
RNA_sheet = pd.read_excel(complete_file, 'S4A_RNA_Expression')
clinic_sheet = pd.read_excel(complete_file, 'S1_Clinical_and_Immune_Data', header=0)
print(RNA_sheet.head()) 

###### RNA FILE ######
RNA_sheet = RNA_sheet.set_index('Supplementary Table 4. RNA expression (normalized expression matrix).')
RNA_sheet = RNA_sheet.T
RNA_sheet.insert(0, "RNA_ID", RNA_sheet['gene_name'], True)
RNA_sheet=RNA_sheet.set_index('gene_name')
RNA_sheet.to_csv("2023 Work/Data/Preprocessed_data/RNA_data.csv")

##### CLINICAL FILE #####
clinic_sheet.columns = clinic_sheet.iloc[0]
clinic_sheet = clinic_sheet.drop([0], axis=0)
merged_df_RNA_clinic = pd.merge(clinic_sheet, RNA_sheet, on='RNA_ID')
merged_df_RNA_clinic.to_csv("2023 Work/Data/Preprocessed_data/clinic_and_RNA_data_raw.csv")

#### FILTER CLINICAL FILE (only NIVOLUMAB patients) ####
clinic_sheet_NIVOLUMAB = clinic_sheet[(clinic_sheet['Arm']=='NIVOLUMAB')]
merged_df_RNA_clinic_NIVOLUMAB = pd.merge(clinic_sheet_NIVOLUMAB, RNA_sheet, on='RNA_ID')
merged_df_RNA_clinic_NIVOLUMAB.to_csv("2023 Work/Data/Preprocessed_data/clinic_and_RNA_data_raw_NIVOLUMAB.csv")

#### FILTER CLINICAL FILE (only EVEROLIMUS patients) ####
clinic_sheet_EVEROLIMUS = clinic_sheet[(clinic_sheet['Arm']=='EVEROLIMUS')]
merged_df_RNA_clinic_EVEROLIMUS = pd.merge(clinic_sheet_EVEROLIMUS, RNA_sheet, on='RNA_ID')
merged_df_RNA_clinic_EVEROLIMUS.to_csv("2023 Work/Data/Preprocessed_data/clinic_and_RNA_data_raw_EVEROLIMUS.csv")