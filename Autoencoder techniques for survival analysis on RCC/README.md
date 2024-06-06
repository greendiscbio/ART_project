#### Autoencoder techniques for survival analysis on RCC

##### Implementation
Houses the mainframe written in Python that preprocesses the data, trains and evaluate the model
Within it we can find:
* Logic
  * Autoencoders: folder containing tabular and graph autoencoders, for both default and variational (implementation changes a bit) penalties
  * Losses: folder containing the classes that implement the loss applied to the loss function, as well as the different losses one can choose from
  * Results: folder where the results are placed in
  * CustomDataset: holds the transcriptomic, clinical and predicted data to make the handling easier
  * CustomKFoldScikit: class that inherits from KFold to add certain invariants to the folding procedure involving censoring
  * FoldObject: holds the data for each fold as well as the results they each obtain
  * GraphDataLoader: custom dataloader for the graph data
  * IterationObject: used by FoldObject to house the data for each fold
  * Main: start point of the program. Many hyperparameters are defined here
  * TabularDataLoader: custom dataloader for the tabular data
  * Trainer: class that trains and evaluates the models
  * TrainingModel: class holding all the information to train the model
* Preprocessing
  * bio_networks: file with auxiliary functions for the creation of the PPI network (with neighbor gathering)
  * create_graph_datasets_RNA: file that creates the graph network (and the tabular dataset is derived from it)

##### Notebooks
Contains .ipynb files that are external to the mainframe
Within it we can find:
* Data analytics: gives analytics about the data such as distributions for different parameters, genes, etc. Used as exploratory analysis at the beginning
* Data preprocessing: preprocesses the direct dataset from the trial to make it a little easier to work with. Whereas the preprocessing in Implementation filters and normalizes the data, here we simply make sure the formatting is correct and values are actually represented in i.e float.
* Graph_analysis: analysis for the graph network we use in the graph autoencoder. Bridge analysis, betweeness, centrality, etc.
* MATH and Histology: _(not used)_ not really useful because we didn't use MATH scores in this work in the end. However, it preprocesses the data related to these variables.
