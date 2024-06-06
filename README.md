#### Autoencoder techniques for survival analysis in RCC
A research evaluating the quality of the compression performed by different autoencoder models.
On the one hand, the tabular autoencoder takes the tabular representation of transcriptomic data and can apply a series of penalties to add to the loss function.
On the other hand, the graph autoencoder considers protein-protein interaction sources between genes to use as input graphs.

The penalties that these autoencoders can use are:
* denoising: improves generalization of the data
* sparse: finds meaningful latent representations
* variational: fits a particular distribution to make the model generative
* penalty combinations: we combined denoising and sparse to benefit from the advantages they both bring

These latent representations obtained by the autoencoders are then used as input in a COX PH model combined with Breslow's estimator to predict the PFS, as well as the risk (Area under ROC) of the patients.

Moreover, we added interpretability by finding the Mutual Information between the original transcriptomic data and the latent representations. The genes that were most expressed
in the representations were LRP2 and ACE2, among a few others.

#### Graph-based model in ICI response for mRCC
A graph-based prediction model that includes information about different protein-protein interaction sources between genes. 
This work considers two different immunotherapies (nivolumab and avelumab+axitinib) and maps the data into two different PPI networks (BioGrid and kidney-specific PPT-Ohmnet) resulting in four
different graph datasets, which are embedded with Graph2Vec. The embeddings are then evaluated as input for a Random Forest classification model using the binarized Progression Free Survival (PFS)
of each treatment dataset. Pathway Enrichment Analysis was used to identify differences and similarities in the overrepresented biological functions that define each of the networks and treatment sets.
