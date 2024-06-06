import torch
from torch_geometric.data import Data
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
from torch_geometric.loader import DataLoader

from sklearn import preprocessing

class GNN():

    def create_edge_index(self):
        path ='GNN model/Data/Kinase_proteins/network_edges_kinase_proteins.tsv'
        data = pd.read_csv(path, delimiter='\t')
        edge_index1=data['#node1'].to_numpy()
        edge_index2=data['node2'].to_numpy()
        le = preprocessing.LabelEncoder()
        le.fit(edge_index1)
        # print(len(list(le.classes_)))
        edge_index1 = le.transform(edge_index1)
        edge_index2 = le.transform(edge_index2)
        edge_index = [edge_index1]+[edge_index2]
        edge_index = np.array(edge_index)

        edge_index = torch.tensor(edge_index, dtype=torch.int64)
        return edge_index

    def create_dataset(self, edge_idex):
        genes = pd.read_csv('GNN model/Data/Kinase_proteins/Kinase_gene_matrix.csv')
        Y = genes.Y
        genes = genes.iloc[:,1:118] 

        list_data=[]

        for g in range(len(genes)):
            b=[]
            for i in genes.iloc[g].to_numpy():
                a=[]
                a.append(i)
                b.append(a)
            x = torch.tensor([b], dtype=torch.long).reshape([-1])
            edge_index = edge_index
            y = torch.tensor([Y.iloc[g]], dtype=torch.float).reshape([-1, 1])
            data = Data(x=x, edge_index=edge_index, y=y)
            list_data.append(data)
        return list_data
    
    def dataset_info(dataset):
        data = dataset[0]
        print(f'Number of nodes: {data.num_nodes}')
        print(f'Number of charcateristics per node: {data.num_features}')
        print(f'Number of edges: {data.num_edges}')
        print(f'Average node degree: {data.num_edges / data.num_nodes:.2f}')
        print(f'Has isolated nodes: {data.has_isolated_nodes()}')
        print(f'Has self-loops: {data.has_self_loops()}')
        print(f'Is undirected: {data.is_undirected()}')

    def train_test_split(dataset):
        torch.manual_seed(0)
        random.shuffle(dataset)
        train_dataset = dataset[0:154]
        test_dataset = dataset[154:182]
        print(f'Number of training graphs: {len(train_dataset)}')
        print(f'Number of test graphs: {len(test_dataset)}')
        train_loader = DataLoader(train_dataset, batch_size=64, shuffle=True)
        test_loader = DataLoader(test_dataset, batch_size=64, shuffle=False)
        for step, data in enumerate(train_loader):
            print(f'Step {step + 1}:')
            print('=======')
            print(f'Number of graphs in the current batch: {data.num_graphs}')
            print(data)
            print()
        return train_loader,test_loader

    

    # from torch import nn
    # from sklearn.metrics import roc_auc_score

    # def train(epoch):
    #     model.train()
    #     criterion = nn.BCELoss()
    #     loss_all = 0
    #     for data in train_loader:
    #         output = model(data.x, data.edge_index, data.batch)
    #         # print("OUTPUT")
    #         # print(output)
    #         loss = criterion(output, data.y.squeeze(1))
    #         optimizer.zero_grad()
    #         loss.backward()
    #         # print(loss.item())
    #         # print(data.num_graphs)
    #         # print(loss_all)
    #         optimizer.step()
    #         loss_all += loss.item() * data.num_graphs

    #     return loss_all / len(train_dataset)


    # def test(loader):
    #     model.eval()

    #     correct = 0
    #     for data in loader:
    #         data = data
    #         output = model(data.x, data.edge_index, data.batch)
    #         for i in range(len(output)):
    #             if output[i]>0.5:
    #                 output[i]=1
    #             else:
    #                 output[i]=0
    #             if output[i]==data.y[i]:
    #                 correct=correct+1
    #     # print("Correct: "+str(correct) +" of "+str(len(loader.dataset)))
    #     return correct / len(loader.dataset)

    # import matplotlib.pyplot as plt
    # model = Net(dim=117)
    # optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    # train_epoch=[]
    # test_epoch=[]
    # for epoch in range(1,100):
    #     loss = train(epoch)
    #     train_acc = test(train_loader)
    #     test_acc = test(test_loader)
    #     train_epoch.append(train_acc)
    #     test_epoch.append(test_acc)
    #     print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, '
    #         f'Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')

    # plt.plot(train_epoch, color="red", label='Train')
    # plt.plot(test_epoch, color="blue", label = 'Test')
    # plt.xlabel("Epoch")
    # plt.ylabel("Accuracy")
    # plt.legend()