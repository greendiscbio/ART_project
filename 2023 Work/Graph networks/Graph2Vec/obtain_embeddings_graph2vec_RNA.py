import pickle
import pandas as pd
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
from karateclub import Graph2Vec, SF, GraphWave
import os
os.chdir('2023 Work\Graph networks\Graph2Vec')
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

   epochs = [100,200,300,400,500,600,700,900,1200,1400]; 
   wl = [2,3] # hyperparams Graph2Vec
   dim = [128,256,512] # default is 128
   for e in epochs:
      for d in dim:
         for w in wl:
            # Load networkx graph datsets
            with open('Data/RNA_dataset_3.pkl', 'rb') as f:
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
            embeddings.to_csv(f'Results_3/RNA_graph2vec_{e}e_{w}WL_{d}dim.csv')
            print(f'RNA_graph2vec_{e}e_{w}WL_{d}dim')

   # result = pd.concat([embeddings, table_labels], axis=1)
   # print(embeddings)

   # sns.heatmap(embeddings)
   # plt.show()

   # sns.clustermap(embeddings, col_cluster=True, standard_scale=1)
   # plt.show()

   # for i in range(1, dimensions):

   #    plt.figure()
   #    sns.scatterplot(data=result, x=0, y=i, hue='DX')
   #    plt.xlabel("Dimension 0")
   #    plt.ylabel(f"Dimension {i}")
   #    plt.show()
