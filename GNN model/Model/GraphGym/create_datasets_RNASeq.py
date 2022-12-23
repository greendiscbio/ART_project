import os
import pickle
import numpy as np
import pandas as pd
import networkx as nx
import torch
from sklearn.model_selection import StratifiedKFold

def create_samples_graphs(data, G):

   '''
   Create a graph for each sample in the dataset, using nodes and edges 
   attributes obtained previously from genetic variants information.
   '''

   samples = list(data.index)

   print('Creating samples graphs...')

   graphs_list = []

   for sample in samples:

      sample_graph = G.copy()
      # Add graph features
      label = data.loc[sample]['Y']
      sample_graph.graph['graph_label'] = torch.tensor([label]) # no cambiar lo de graph_label

      sample_graph.graph['sampleID'] = sample # IMPORTANTE PARA EL SPLIT, no cambiar 'sampleID'
      
      # Add node and edge features
      for n in G.nodes:
         expression = data.loc[sample][n]
         sample_graph.nodes[n]['node_feature'] = torch.tensor([expression]) # no cambiar lo de node_feature
      
      # Delete components with less than 10 nodes - cambiar para lo que necesites!!!
      for component in list(nx.connected_components(sample_graph)):
         if len(component)<10:
               for node in component:
                  sample_graph.remove_node(node)

      graphs_list.append(sample_graph)
      print(f'{sample} graph: nodes = {nx.number_of_nodes(sample_graph)}; edges = {nx.number_of_edges(sample_graph)}; id = {label}')
   
   print()

   return graphs_list

def create_folds_stratified_cv(n):


   with open(f'output_graph_dataset.pkl', 'rb') as f: # it doesn't matter the network we use because samples are the same
      list_of_dicts = pickle.load(f)
   
   graphs = []
   for i in list_of_dicts:
      graph = nx.Graph(i)
      sample = graph.graph['sampleID']
      label  = graph.graph['graph_label'].item()
      graphs.append([sample, label])

   df = pd.DataFrame(graphs, columns=['sample', 'y'])
   print(df)
   df['trick'] = df['y']
   df.set_index('sample', inplace=True)
   y = df['y']
   print(df)
   skf = StratifiedKFold(n_splits=n)
   k = 1
   for train_idx, val_idx in skf.split(df, y):
      print(f'Fold -  {k}  |  train -  {np.bincount(y[train_idx])}  |  test -  {np.bincount(y[val_idx])}')

      split_dict = {}
      split_dict['train'] = list(train_idx)
      split_dict['valid'] = list(val_idx)

      # print(split_dict['valid'])

      f = open(f'GNN model/Model/GraphGym/{n}fold_CV_split/k{k}.pkl', 'wb')
      pickle.dump(split_dict, f)
      f.close()

      k += 1
    
   print()

if __name__ == "__main__":

   input_network = nx.read_edgelist('GNN model/Data/PPT-Ohmnet/Graph_PPT_Ohmn/Second big pool/AD_SNAP_PPI_kidney_RNA_30_nodes.edgelist') #  [*]
   input_data = pd.read_csv('GNN model/Data/PPT-Ohmnet/mRCC_big_pool/Second big pool/mrcc_protein_matrix_RNA_30_nodes_with_patient_id.csv', index_col='id') # 
   input_data = input_data.drop(['Unnamed: 0'], axis=1)
   print(input_data)
   '''

   input_file_network.edgelist: edgelist que conecta los genes (nodos)
   asegúrate de que la red de entrada no tiene:
   - nodos o componentes aisladas (quédate con la componente más grande, en la func de arriba tiene un paso extra para hacer eso)
   - self-loops: ejes entre un propio nodo; a veces pasa porque las proteínas interaccionan consigo mismas

   input_file_data.csv: CSV con datos RNA: pacientes (filas) x genes (columnas)
   asegúrate de que hay una columna (que en mi caso es "clusters") que tiene la etiqueta de los pacientes en el dataset
   
   '''

   print('NETWORK\t', input_network)
   biogrid_gd = create_samples_graphs(input_data, input_network)
   with open(f'GNN model/Model/GraphGym/output_graph_dataset.pkl', 'wb') as f:
      pickle.dump(biogrid_gd, f)


   '''
   En este punto se hace el split, saca N conjuntos train-val distintos
   No hay test porque en GraphGym funciona igual que el validation, no pasa nada solo queremos ver cómo entrena
   '''

   create_folds_stratified_cv(5) # numero de folds, haz menos que tu muestra es más chica