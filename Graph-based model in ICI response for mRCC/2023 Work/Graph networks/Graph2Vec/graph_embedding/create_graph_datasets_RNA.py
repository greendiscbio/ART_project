import pickle
import pandas as pd
import networkx as nx
from sklearn.preprocessing import StandardScaler
import os
os.chdir('2023 Work/Graph networks/Graph2Vec/graph_embedding')

def create_rna_file(ppi, level):

   rna_data = pd.read_csv('../../../Data/Preprocessed_data/clinic_and_RNA_data_raw_Nivolumab.csv')
   nodes=[];no_nodes=[]
   rna_data.columns=rna_data.columns.str.upper()
   for node in ppi.nodes():
         if node in rna_data.columns:
            nodes.append(node)
         else:
            no_nodes.append(node)
   for node in no_nodes:
      ppi.remove_node(node)
   print(ppi)
   print(nx.number_connected_components(ppi))
   exit()
   pfs = rna_data['PFS'].values
   booleanMask = rna_data.columns.isin(nodes)
   # saving the selected columns
   rna_data = rna_data[rna_data.columns[booleanMask]]
   rna_data.insert(0, 'PFS', pfs)
   rna_data.to_csv(f'Data/RNA_Disgenet_matrix_{level}_Nivolumab.csv')
   return rna_data

def create_samples_graphs(data, G):

   '''
   Create a graph for each sample in the dataset, using nodes and edges
   attributes obtained previously from genetic variants information.
   '''

   samples = list(data.index)

   print('Creating samples graphs...')

   graphs_list = []
   for sample in samples:

      label = data.loc[sample]['PFS']
      sample_graph = G.copy()
      sample_data = data.loc[sample]

      # Add graph features
      sample_graph.graph['graph_label'] = label
      sample_graph.graph['sample_id']   = sample

      # Add node and edge features
      for node in G.nodes:
            node_expr = sample_data[node]
            sample_graph.nodes[node]['node_attr'] = node_expr

      graphs_list.append(sample_graph)
      # print(sample_graph.nodes(data=True))
      # print(f'{sample} graph: nodes = {nx.number_of_nodes(sample_graph)}; edges = {nx.number_of_edges(sample_graph)}; label = {label}')

   return graphs_list


if __name__ == "__main__":
   print('Level of the edgelist: ')
   n = input()
   ppi = nx.read_edgelist('../ppi_creation/Data/biogrid_found_genes_level_'+str(n)+'_Nivolumab.edgelist')
   print(ppi)
   expression_data = create_rna_file(ppi, n)
   exit()
   # expression_data = pd.read_csv('Data/RNA_Disgenet_matrix.csv', index_col=0)

   # scaler = StandardScaler()
   # scaler.fit(expression_data.drop(columns='PFS'))



   # Graph datasets with PPI networks and expression value per node
   variants_gd = create_samples_graphs(expression_data, ppi)

   with open('Data/RNA_dataset_'+str(n)+'_Nivolumab.pkl', 'wb') as f:
      pickle.dump(variants_gd, f)