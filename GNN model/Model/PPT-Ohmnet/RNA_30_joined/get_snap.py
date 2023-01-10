import networkx as nx
import pandas as pd
import numpy as np
import mygene

def get_snap(genes_file):
    G = nx.read_edgelist('GNN model\Data\PPT-Ohmnet\Graph_PPT_Ohmn/PPT-Ohmnet_tissues-combined.edgelist', nodetype=int, data=(('tissue', str),))

    tissues_edgelist = pd.read_csv('GNN model/Data/PPT-Ohmnet/Graph_PPT_Ohmn/PPT-Ohmnet_tissues-combined.edgelist', sep='\t')
    kidney_specific = tissues_edgelist[tissues_edgelist['tissue'] == 'kidney']
    tissues = ['retina','renal_glomerulus','pons','duodenum','cerebellum','chondrocyte','tear_gland','glia','hepatocyte','t_lymphocyte','culture_condition_cd8_cell','embryo','ovary','vascular_endothelium','parietal_lobe','spleen','natural_killer_cell','telencephalon','adrenal_gland','forebrain','neuron','thymocyte','temporal_lobe','amygdala','cornea','trophoblast','bronchus','epidermis','dendritic_cell','skin_fibroblast','skeletal_muscle','basophil','eosinophil','dentate_gyrus','tonsil','corpus_striatum','occipital_pole','megakaryocyte','lymph_node','hair_follicle','eye','blood_platelet','blood_vessel','tooth','large_intestine','lymphocyte','cochlea','caudate_nucleus','gastrointestinal_tract','serum','keratinocyte','ileum','diencephalon','caudate_putamen','spermatid','macrophage','pancreas','trachea','stomach','small_intestine','prostate_gland','adipose_tissue','astrocyte','pancreatic_islet','mammary_epithelium','hypothalamus','ovarian_follicle','podocyte','locus_ceruleus','corpus_luteum','salivary_gland','bronchial_epithelial_cell','choroid','artery','occipital_lobe','uterine_endometrium','uterus','liver','fetus','medulla_oblongata','monocyte','b_lymphocyte','thyroid_gland','blood','thalamus','mast_cell','kidney','uroepithelium','granulocyte','midbrain','ear','frontal_lobe','cerebellar_cortex','osteoblast','cecum','corpus_callosum','esophagus','basal_ganglion','cardiac_muscle','brain','substantia_nigra','spermatocyte','blood_plasma','skin','colon','lung','vascular_endothelial_cell','hypophysis','spinal_cord','adrenal_cortex','umbilical_cord','heart','cerebral_cortex','vermiform_appendix','hippocampus','myometrium','renal_tubule','smooth_muscle','urinary_bladder','muscle','leukocyte','spermatogonium','uterine_cervix','placenta','cartilage','bone','subthalamic_nucleus','peripheral_nervous_system','testis','lens','nephron','hematopoietic_stem_cell','bone_marrow','jejunum','nervous_system','neutrophil','nucleus_accumbens','mononuclear_phagocyte','central_nervous_system','intestine','mammary_gland','umbilical_vein_endothelial_cell','oviduct','aorta']
    kidney_specific = tissues_edgelist[tissues_edgelist['tissue'].isin(tissues)]

    kidney_specific.to_csv('GNN model/Data/PPT-Ohmnet/Graph_PPT_Ohmn/RNA_30_joined/PPT-Ohmnet_tissues-kidney.edgelist', sep='\t', index=False)
    G_kidney = nx.read_edgelist('GNN model/Data/PPT-Ohmnet/Graph_PPT_Ohmn/RNA_30_joined/PPT-Ohmnet_tissues-kidney.edgelist', nodetype=int, data=(('tissue', str),))
    rna_genes = np.loadtxt('GNN model/Data/all_genes_db.txt', dtype=str,  delimiter=',', skiprows=0)
    print("Nodes: 6")
    print(G_kidney.number_of_nodes())
    # List of genes to search for
    infile = open(genes_file, "r")
    genes = infile.read().split("\n")
    infile.close()
    print('Quering '+ str(len(genes))+' genes.')

    # Genes in PPT-Ohmnet are Entrez IDs, it is necessary to convert them to gene Symbols.
    mg = mygene.MyGeneInfo()
    out = mg.querymany(genes, scopes='symbol', fields='entrezgene', species='human')
    print(len(out))
    entrezgenes = []
    no_entrezgenes = []
    mapping = {}
    for o in out:
        print(o)
        if 'entrezgene' in o and o['query'] in rna_genes: #['ACACA', 'ADRA2B', 'ATP5F1A', 'CCN1', 'CCN2', 'DACH1', 'H2AX', 'H4C1', 'H4C11', 'H4C12', 'H4C13', 'H4C14', 'H4C15', 'H4C2', 'H4C3', 'H4C4', 'H4C5', 'H4C6', 'H4C8', 'H4C9', 'HNF1B', 'MARCHF8', 'NSD2', 'ORAI1', 'PECAM1', 'PRKN', 'TWNK', 'TXNIP']:
            entrezgenes.append(int(o['entrezgene']))
            mapping[int(o['entrezgene'])] = o['query']
        else:
            no_entrezgenes.append(o['query'])
            print(o)
    print('Genes without entrezgene or not in rna list: ' + str(len(no_entrezgenes)))
    print('Number of nodes found in MyGene: '+ str(len(entrezgenes)))
    print(entrezgenes)

    A_kidney_frozen = G_kidney.subgraph(entrezgenes)
    A_kidney = nx.Graph(A_kidney_frozen)
    original = A_kidney.number_of_nodes()
    print('Number of nodes included in PPT-Ohmnet: ' + str(original))
    
get_snap('GNN model/Data/rna_30.txt')