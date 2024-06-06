import torch
from torch.nn import Linear
import torch.nn.functional as F
from torch_geometric.nn import GraphConv, global_add_pool
embed_dim = 117
from torch_geometric.nn import TopKPooling
from torch_geometric.nn import global_mean_pool as gap, global_max_pool as gmp


class Net(torch.nn.Module):
        def __init__(self, dim):
            super(Net, self).__init__()
            self.dim = dim
            super(Net, self).__init__()
            self.conv1 = GraphConv(embed_dim, dim)
            self.pool1 = TopKPooling(dim, ratio=0.8)
            self.conv2 = GraphConv(dim, dim)
            self.pool2 = TopKPooling(dim, ratio=0.8)
            self.item_embedding = torch.nn.Embedding(num_embeddings=390, embedding_dim=embed_dim)
            self.lin1 = torch.nn.Linear(234, 51)
            self.lin3 = torch.nn.Linear(51, 1)
            self.act1 = torch.nn.ReLU()
            print(self)

        def forward(self, x, edge_index, batch):
            x = torch.tensor(x).to(torch.int)
            x = self.item_embedding(x)
            x = x.squeeze(1)

            x = F.relu(self.conv1(x, edge_index))

            x, edge_index, _, batch, _, _ = self.pool1(x, edge_index, None, batch)
            x1 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

            x = F.relu(self.conv2(x, edge_index))

            x, edge_index, _, batch, _, _ = self.pool2(x, edge_index, None, batch)
            x2 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

            x = x1 + x2

            x = self.lin1(x)
            x = self.act1(x)

            x = F.dropout(x, p=0.5, training=self.training)
            x = torch.sigmoid(self.lin3(x)).squeeze(1)
            return x
