import pandas as pd
import numpy as np
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_classif
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from bioinfokit.analys import stat
import scipy.stats as stats
from sklearn.feature_selection import VarianceThreshold
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import ExtraTreesClassifier
def import_data(path):
    # load data file
    df = pd.read_csv(path)
    if 'Y' in df.columns:
        y = df['Y']
        df = df.iloc[:,1:17826]
    elif 'PFS' in df.columns:
        y = []
        for patient in range (len(df)):
            if df['PFS'][patient] < 3:
                y.append(0)
            else:
                y.append(1)
        df = df.iloc[:,122:44015]
    else:
        print('Target column has not been found.')

    print(df)
    return df, y


def anova_study(df, y, n_features):
    X = df
    # print(X)
    # print(y)
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
        included = compare_dataframe_with_biogris_genes(X_kbest)

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
        return included
    else:
        n = df.shape[0]   # 181 patients
        k = df.shape[1]   # 17826 genes

        # 1. Calculate the mean for allgenes:
        df.loc['Group Means'] = df.mean()
        # print(df)

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

def variance_threshold_study(data, threshold):
    sel = VarianceThreshold(threshold=(threshold * (1 - threshold)))
    result = sel.fit_transform(data)
    cols = sel.get_feature_names_out(data.columns.values)
    result = pd.DataFrame(result, columns=cols)
    # print(result)
    included = compare_dataframe_with_biogris_genes(result)
    return included

def tree_based_study(data, y, min_importance):
    clf = ExtraTreesClassifier(n_estimators=100)
    clf = clf.fit(data, y)
    importances = clf.feature_importances_
    elems_over_0 = np.fromiter((element for element in importances if element > min_importance), dtype = importances.dtype)
    # print('Features with relevance over 0: ', len(elems_over_0))
    model = SelectFromModel(clf, prefit=True)
    X_new = model.transform(data)
    # print(X_new.shape)
    cols = []
    for name, importance in zip(data, clf.feature_importances_):
        if importance > min_importance:
            cols.append(name)
    result = pd.DataFrame(X_new, columns=cols)
    included = compare_dataframe_with_biogris_genes(result)
    return included

def compare_dataframe_with_biogris_genes(df):
    path = "2023 Work/Graph networks/Biogrid/Data/genes_found.txt"
    infile = open(path, "r")
    rna_genes = infile.read().split("\n")
    infile.close()

    included = 0
    not_included = []
    for d in df:
        if d in rna_genes:
            included = included +1
            with open('2023 Work\Graph networks\Biogrid\Data\ANOVA\Anova_selected_'+str(n_features)+'_genes.txt',"a") as file2:
                file2.write(str(d) + "\n")
        else:
            not_included.append(d)
    return included
###################################### MAIN ########################################

# path = "2023 Work/Data/Preprocessed_data/Gene Matrix/biogrid_included_genes_NIVOLUMAB.csv"
path = "2023 Work/Data/Preprocessed_data/clinic_and_RNA_data_raw_NIVOLUMAB.csv"
data, y = import_data(path)
print(data.columns)
n_features = input("How many features do you want to select? (number/all) Maximum number of featues fo select: "+str(len(data.columns))+" ")
included = anova_study(data, y, n_features)
print('ANOVA: '+str(included))
# threshold = input("Remove features bellow: (0.0-1.0) ")
# included = variance_threshold_study(data, float(threshold))
# print("Variance threshold: "+str(included))
# included = tree_based_study(data, y, 0)
# print("Tree-based: "+str(included))

#####################################################################################

# Compare with RNA_total data
# path = "2023 Work\Data\Preprocessed_data\Feature selection\RNA_total.txt"
# path = "2023 Work/Graph networks/Biogrid/Data/genes_found.txt"
# infile = open(path, "r")
# rna_genes = infile.read().split("\n")
# infile.close()
# path = "2023 Work/Graph networks/Spectral_clustering+CNN/Feature selection/hvgs.txt"
# infile = open(path, "r")
# hvgs_genes = infile.read().split("\n")
# infile.close()
# included = 0
# not_included = []
# for k in hvgs_genes:
#     if k in rna_genes:
#         included = included +1
#     else:
#         not_included.append(k)
# print(not_included)
# print(included)
# print(len(not_included))
