import pandas as pd
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_classif
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from bioinfokit.analys import stat
import scipy.stats as stats
import scanpy as sc
import numpy as np
def import_data (path):
    # load data file
    df = pd.read_csv(path)
    y = df['Y']
    df = df.iloc[:,1:17826]
    print(df)
    return df, y

def anova_study(df, y, n_features):
    X = df
    if n_features != "all":
        fvalue_Best = SelectKBest(f_classif, k=int(n_features))
        X_kbest = fvalue_Best.fit_transform(X, y)
        cols = fvalue_Best.get_feature_names_out(df.columns.values)
        # print(X_kbest)
        print('Original number of features:', X.shape)
        print('Reduced number of features:', X_kbest.shape)

        fvalue, pvalue = stats.f_oneway(
            *X_kbest
        )
        print("F-vlaue: " +str(fvalue))
        print("P-value: " +str(pvalue))

        X_kbest = pd.DataFrame(X_kbest, columns=cols)

        # # reshape the d dataframe suitable for statsmodels package 
        df_melt = pd.melt(X_kbest.reset_index(), id_vars=['index'], value_vars=X_kbest.columns.values)
        # replace column names
        df_melt.columns = ['index', 'genes', 'value']

        if int(n_features)<30:
            ax = sns.boxplot(x='genes', y='value', data=df_melt, color='#99c2a2')
            ax = sns.swarmplot(x="genes", y="value", data=df_melt, color='#7d0013', size=3)
            plt.show()
        else:
            print("The number of feature selected is to big to be shown in a box plot")

        if int(n_features)<1000:
            # ANOVA table using bioinfokit v1.0.3 or later (it uses wrapper script for anova_lm)
            res = stat()
            res.anova_stat(df=df_melt, res_var='value', anova_model='value ~ C(genes)')
            print(res.anova_summary)
        else:
            print("The number of feature selected is to big to be shown in the sumary")

    else:
        n = df.shape[0]   # 181 patients
        k = df.shape[1]   # 17826 genes

        # 1. Calculate the mean for allgenes:
        df.loc['Group Means'] = df.mean()
        print(df)

        overall_mean = df.iloc[-1].mean()
        print("Overall mean: " + str(overall_mean))

        # 2. Calculate the sum of squares
        '''The sum of squares of all observation is calculated by deducting 
        each observation from the overall mean, and then summing all the 
        squares of the differences'''
        SS_total = (((df.iloc[:-1] - overall_mean)**2).sum()).sum()
        print("Sum of squares: " + str(SS_total))

        '''The sum of squares within is the sum of squared deviations 
        of scores around their group’s mean'''
        SS_within = (((df.iloc[:-1] - df.iloc[-1])**2).sum()).sum() 
        print("Sum of squares within the sum of squared deviation: " + str(SS_within))

        '''Next we calculate the sum of squares of the group means
        from the overall mean'''
        SS_between = (n * (df.iloc[-1] - overall_mean)**2).sum()
        print("Sum of squares of the group means from the overall mean: " + str(SS_between))

        '''SS_total = SS_between + SS_within'''
        mean_sq_between = SS_between / (k - 1)      
        mean_sq_within = \
            SS_within / (n*k - k)

        F = mean_sq_between / mean_sq_within
        print("F: " + str(F))

        fvalue, pvalue = stats.f_oneway(
            *df.iloc[:-1,0:17826].T.values
        )
        print("F-vlaue: " +str(fvalue))
        print("P-value: " +str(pvalue))
        
def higly_variable_genes(path):
    adata = sc.read(path)
    print(adata)
    sc.pp.scale(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    # sc.pl.umap(adata, color=["samples", "n_genes", "mt_frac", "doublet_score"], ncols=2)

    sc.pp.highly_variable_genes(adata, min_mean=0.01, max_mean=8, min_disp=1, n_top_genes=4000, flavor="cell_ranger", n_bins=20)
    sc.pl.highly_variable_genes(adata)
    print("Highly variable genes: ", np.sum(adata.var["highly_variable"]))
    print("Total genes:", adata.shape[1])
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
###################################### MAIN ########################################

path = "2023 Work/Data/Preprocessed_data/Gene Matrix/biogrid_included_genes_NIVOLUMAB.csv"
data, y = import_data(path)
# n_features = input("How many features do you want to select? (number/all) Maximum nmber of featues is: "+str(len(data.columns))+" ")
# anova_study(data, y, n_features)
higly_variable_genes(path)

#####################################################################################
