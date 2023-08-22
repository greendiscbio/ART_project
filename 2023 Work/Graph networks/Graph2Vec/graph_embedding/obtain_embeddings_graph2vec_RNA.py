import pickle
import pandas as pd
import networkx as nx
import seaborn as sns
from karateclub import Graph2Vec
import os
os.chdir('2023 Work/Graph networks/Graph2Vec/graph_embedding')

def load_graphs(graph_dataset):

   graphs_int = []
   graph_ids = []

   for g in graph_dataset:

      sample_id = g.graph['sample_id']
      graph_ids.append(sample_id)

      for n, data in g.nodes(data=True):
         data['feature'] = data.pop('node_attr')

      h = nx.relabel.convert_node_labels_to_integers(g, label_attribute='gene')
      # print(h.nodes(data=True))
      graphs_int.append(h)

   return graphs_int, graph_ids


if __name__ == "__main__":


   print()

   epochs = [900,1000,1200,1300,1400,150,1600,1700,1800,1900,2000,2500];
   wl = [2, 3, 4] # hyperparams Graph2Vec
   dim = [256, 512, 1028] # default is 128
   for e in epochs:
      for d in dim:
         for w in wl:
            # Load networkx graph datsets
            with open('Data/RNA_dataset_2_Nivolumab.pkl', 'rb') as f:
               original_graphs = pickle.load(f)

            processed_graphs, graphs_ids = load_graphs(original_graphs)

            # Graph embeddings
            # Graph2Vec
            graph2vec = Graph2Vec(dimensions=d, attributed=True,
                                    epochs=e, wl_iterations=w)
            graph2vec.fit(processed_graphs)
            g_emb = graph2vec.get_embedding()

            embeddings = pd.DataFrame(g_emb)
            embeddings['sample_id'] = graphs_ids
            embeddings.set_index('sample_id', inplace=True)
            embeddings.to_csv(f'Results_embedding/Nivolumab/RNA_graph2vec_{e}e_{w}WL_{d}dim.csv')
            print(f'RNA_graph2vec_{e}e_{w}WL_{d}dim')