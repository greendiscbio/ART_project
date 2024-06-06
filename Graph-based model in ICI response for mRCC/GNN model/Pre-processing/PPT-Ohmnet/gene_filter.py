from unittest import case
import pandas as pd
import numpy as np

path ='C:/Users/sandr/Documents/ART_project/GNN model/Data/PPT-Ohmnet/mRCC_big_pool/mrcc genes second big pool.xlsx'
data = pd.read_excel(path)
print('Total genes: '+str(len(data)))

data = data.sort_values('Score_gda')
data = data.drop_duplicates(data.columns[~data.columns.isin(['Gene'])], keep='first')

print('Filtered genes: '+str(len(data)))

print(len(data[data['Score_gda']>0.3]))
