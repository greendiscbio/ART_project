import pandas as pd
import numpy as np

from sklearn.metrics import  matthews_corrcoef, balanced_accuracy_score, f1_score, roc_auc_score, recall_score, precision_score, matthews_corrcoef, confusion_matrix, roc_curve

import sklearn.metrics._scorer as sc

from sklearn.model_selection import cross_val_score, StratifiedKFold

from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier

from functools import reduce
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import os
os.chdir('2023 Work/Graph networks/Graph2Vec')


def upload_data(graph2vec_file, rna_file):

    embeddings = pd.read_csv(graph2vec_file, index_col=0)
    metadata = pd.read_csv(rna_file, index_col=0)
    metadata['PFS_label'] = pd.cut(metadata['PFS'], bins=[0, 6, np.inf], labels=['NR', 'R'])
    results = pd.concat([embeddings, metadata[['PFS_label']]], axis=1)
    return results


def transform_data(results):

    x = results.drop(columns='PFS_label')
    y = []
    for i in results['PFS_label']:
        if i == 'NR':
            y.append(0)
        else:
            y.append(1)

    return x,y


def _get_model_name(model):
    """
            Returns a string with the name of a sklearn model
                model: Sklearn stimator class
    """
    name = str(model)[:str(model).find("(")]
    return name


def get_best_model(name, X_train , y_train, cv):

    if name == 'LogisticRegression':
        param_grid = {
            'C': [0.25,0.5,1,2,3,4,5,6,7,10],
            'max_iter':[1000,5000,10000],
            'random_state':[42],
            'solver': ['lbfgs', 'newton-cg', 'saga']
            }

        clf = GridSearchCV(LogisticRegression(), param_grid, cv = cv)

    elif name == 'SVC':
        param_grid = {
        'C': [0.25,0.5,1,2,3,4,5,6,7,10],
        'kernel':['linear', 'poly', 'rbf', 'sigmoid'],
        'degree':[1,2,3,4,5,6,7],
        'random_state' :[42],
        'gamma':['scale', 'auto']        
        }

        clf = GridSearchCV(SVC(), param_grid, cv = cv)

    elif name == 'RandomForestClassifier':
        param_grid = {'min_samples_leaf': [1, 2, 3],
            'min_samples_split': [2, 3, 4, 5],
            'random_state':[42],
            'n_estimators': [15, 20],
            'bootstrap': [True, False],
            'criterion': ['gini', 'entropy'],
            'max_depth':[None, 2, 5, 10,50]
            }

        clf = GridSearchCV(RandomForestClassifier(), param_grid, cv = cv)

    elif name == 'MLPClassifier':
        param_grid = {'max_iter': [200000, 50000,10000000],
            'activation': ['identity', 'logistic', 'tanh', 'relu'],
            'random_state': [42],
            'max_fun': [300, 500,1000, 5000, 10000, 15000, 20000],
            'hidden_layer_sizes': [3,5]}

        clf = GridSearchCV(MLPClassifier(), param_grid, cv = cv)

    clf.fit(X_train, y_train)
    model = clf.best_estimator_

    return model


def cv_score(X_train, X_test, y_train, y_test,  e, wl, dim, models_list, scoring_list = None, return_scores = True):

    train_mean_score = list()
    ldf = list(); ldf_std = list(); test_ldf = list()

    # cv = StratifiedKFold(n_splits=5)
    kf=StratifiedKFold(n_splits=10, random_state=None, shuffle=False)
    for i, model in enumerate(models_list):
        val_score = {'balanced_accuracy_score': [], 'f1_score':[], 'roc_auc_score': [], 'precision_score':[], 'recall_score':[], 'matthews_corrcoef':[]}
        val_mean_score = []; val_std_score = []

        # Get model
        name = _get_model_name(model)

        # Get best hyperparametres for current model
        model = get_best_model(name, X_train , y_train, kf)
                
        # Perform CV
        for train_index, val_index in kf.split(X_train, y_train):
             
            # Create train and validation datasets
            train_dataset = []; val_dataset = []; train_y_dataset = []; val_y_dataset = []
            # print("TRAIN: ", train_index, "VAL:", val_index)

            for i in train_index:
                train_dataset.append(X_train.iloc[i])
                train_y_dataset.append(y_train[i])
            for i in val_index:
                val_dataset.append(X_train.iloc[i])
                val_y_dataset.append(y_train[i])
            # print('Patients in train dataset: '+ str(len(train_dataset)))
            # print('Patients in val dataset: '+ str(len(val_dataset)))

            # Fit the model
            model.fit(train_dataset, train_y_dataset)

            # Get validation metrics
            val_pred = model.predict(val_dataset)

            val_score['balanced_accuracy_score'].append(balanced_accuracy_score(val_y_dataset, val_pred))
            val_score['f1_score'].append(f1_score(val_y_dataset, val_pred, average='weighted'))
            val_score['roc_auc_score'].append(roc_auc_score(val_y_dataset,val_pred, average = 'weighted'))
            val_score['precision_score'].append(precision_score(val_y_dataset, val_pred, average='weighted'))
            val_score['recall_score'].append(recall_score(val_y_dataset, val_pred, average='weighted'))
            val_score['matthews_corrcoef'].append(matthews_corrcoef(val_y_dataset, val_pred))

        # Save validation metrics
        for metric in val_score:
            val_mean_score.append(np.mean(val_score[metric]))
            val_std_score.append(np.std(val_score[metric]))


        # Save current model val scores
        tmp = pd.DataFrame({name: val_mean_score}, index = score_list)
        tmp_std = pd.DataFrame({name: val_std_score}, index = score_list)

        # Add current model val scores to final dataframe (all models and metrics included)
        ldf.append(tmp)
        ldf_std.append(tmp_std)

        # Refit model (train+val) for Test evaluation
        model.fit(X_train, y_train)

        # Test evaluation
        test_pred = model.predict(X_test)
        if name != 'SVC': # SVC fit method as it internally uses 5-fold cross-validation, and predict_proba may be inconsistent with predict
            test_proba = model.predict_proba(X_test)

        # Get test scores
        balance_accuracy = balanced_accuracy_score(y_test, test_pred)
        f1 = f1_score(y_test, test_pred, average='weighted')
        roc_auc = roc_auc_score(y_test,test_pred, average = 'weighted')
        precision = precision_score(y_test, test_pred, average='weighted')
        recall = recall_score(y_test, test_pred, average='weighted')
        matthews = matthews_corrcoef(y_test, test_pred)

        test_score = [balance_accuracy, f1, roc_auc, precision, recall, matthews]

        # Save test confussion matrix
        confussion_matrix(y_test, test_pred, e, wl, dim, name)

        # Plot ROC AUC
        if name != 'SVC': # SVC fit method as it internally uses 5-fold cross-validation, and predict_proba may be inconsistent with predict
            result_visualization(roc_auc, y_test, test_proba, name)

        # Save current model test scores
        test_tmp = pd.DataFrame({name: test_score}, index = score_list)

        # Add current model test scores to final dataframe (all models and metrics included)
        test_ldf.append(test_tmp)

    # Build final report with models and metrics
    train_scores = reduce(lambda x,y: pd.merge(x,y, left_index = True, right_index = True), ldf).T
    train_scores_std = reduce(lambda x,y: pd.merge(x,y, left_index = True, right_index = True), ldf_std).T
    test_scores = reduce(lambda x,y: pd.merge(x,y, left_index = True, right_index = True), test_ldf).T

    return train_scores, train_scores_std, test_scores


def plot_train_test_metrics(train_scores, test_scores):    

    fig, ax  = plt.subplots(1,1, figsize = (17,14))

    train_scores.plot.bar(ax = ax, cmap = 'RdYlBu', edgecolor = "black")
    ax.legend(loc = 'best')
    ax.set_xlabel("Score")
    ax.set_title(f'Cross validation for RNA_graph2vec_{e}e_{wl}WL_{dim}dim all_components: Train set')
    plt.savefig(f'Models/Results_level_2/images/Cross_validation_models_RNA_graph2vec_{e}e_{wl}WL_{dim}dim_2_all_components_train.png')
    plt.show()


    fig, ax  = plt.subplots(1,1, figsize = (17,14))

    test_scores.plot.bar(ax = ax, cmap = 'RdYlBu', edgecolor = "black")
    ax.legend(loc = 'best')
    ax.set_xlabel("Score")
    ax.set_title(f'Cross validation for RNA_graph2vec_{e}e_{wl}WL_{dim}dim all_components: Test set')
    plt.savefig(f'Models/Results_level_2/images/Cross_validation_models_RNA_graph2vec_{e}e_{wl}WL_{dim}dim_2__all_components_test.png')
    plt.show()



def confussion_matrix(y_test, test_pred, e, wl, dim, name):

    # with open('Models/Results_level_2/Cross_validation_models_RNA_graph2vec_2_all_components_confussion_matrix_train.txt', 'a') as file:
    #     file.write(f'RNA_graph2vec_{e}e_{wl}WL_{dim}dim_full train confussion matrix for {name}'+'\n')
    #     conf_matrix = confusion_matrix(y_train,val_pred)
    #     file.write(str(conf_matrix)+'\n')
    with open('Models/Results_level_2/Cross_validation_models_RNA_graph2vec_2_all_components_confussion_matrix_test.txt', 'a') as file:
        file.write(f'RNA_graph2vec_{e}e_{wl}WL_{dim}dim_full test confussion matrix for {name}'+'\n')
        conf_matrix = confusion_matrix(y_test, test_pred)
        file.write(str(conf_matrix)+'\n')

def result_visualization(roc_auc, y_test, test_proba, name):

    fpr, tpr, thresholds = roc_curve(y_test, test_proba[:,1])
    
    plt.figure()
    plt.plot(fpr, tpr, label='AUC Curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'RNA graph2vec {e}e {wl}WL {dim}dim 2 all components {name}')
    plt.legend(loc="lower right")
    plt.savefig(f'Models/Results_level_2/images/Auc_curve_RNA_graph2vec_{e}e_{wl}WL_{dim}dim_2_all_components_test_{name}.png')

    plt.show()

if __name__ == "__main__":
    for e in [100, 200, 300, 400, 500, 600, 700, 900, 1000, 1200, 1400]:
        for dim in [128,256,512]:
                for wl in [2,3]:
                    # Caso concreto para probar
                    # e = 500; wl = 3; dim = 512
                    results = upload_data(f'Results_2/RNA_graph2vec_{e}e_{wl}WL_{dim}dim_full.csv', 'Data/RNA_Disgenet_matrix_2_full.csv')
                    x, y = transform_data(results)
                    X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.2, stratify=y)
                    
                    pipe = make_pipeline(StandardScaler())
                    pipe.fit(X_train, y_train)  # apply scaling on training data

                    X_train.reset_index(drop = True, inplace = True)
                    # print(X_train)

                    models_list =[LogisticRegression(),
                                SVC(),
                                RandomForestClassifier(),
                                MLPClassifier()
                                ]

                    score_list = ["balanced_accuracy", "f1_weighted", "roc_auc", "precision_weighted", "recall_weighted", "matthews_corrcoef"]

                    train_scores, train_scores_std, test_scores = cv_score(X_train, X_test, y_train, y_test, e, wl, dim, models_list, score_list)
                    print(train_scores)
                    print(train_scores_std)
                    print(test_scores)

                    plot_train_test_metrics(train_scores, test_scores)

                    with open('Models/Results_level_2/Cross_validation_models_RNA_graph2vec_2_all_components_train_scores.txt', 'a') as file:
                        file.write(f'RNA_graph2vec_{e}e_{wl}WL_{dim}dim_all_components'+'\n')
                    train_scores.to_csv(r'Models/Results_level_2/Cross_validation_models_RNA_graph2vec_2_all_components_train_scores.txt', sep='\t', mode='a', encoding='utf-8')
                    with open('Models/Results_level_2/Cross_validation_models_RNA_graph2vec_2_all_components_train_scores.txt', 'a') as file:
                        file.write('\n')

                    with open('Models/Results_level_2/Cross_validation_models_RNA_graph2vec_2_all_components_train_scores_std.txt', 'a') as file:
                        file.write(f'RNA_graph2vec_{e}e_{wl}WL_{dim}dim_all_components'+'\n')
                    train_scores_std.to_csv(r'Models/Results_level_2/Cross_validation_models_RNA_graph2vec_2_all_components_train_scores_std.txt', sep='\t', mode='a', encoding='utf-8')
                    with open('Models/Results_level_2/Cross_validation_models_RNA_graph2vec_2_all_components_train_scores_std.txt', 'a') as file:
                        file.write('\n')

                    with open('Models/Results_level_2/Cross_validation_models_RNA_graph2vec_2_all_components_test_scores.txt', 'a') as file:
                        file.write(f'RNA_graph2vec_{e}e_{wl}WL_{dim}dim_all_components'+'\n')
                    test_scores.to_csv(r'Models/Results_level_2/Cross_validation_models_RNA_graph2vec_2_all_components_test_scores.txt', sep='\t', mode='a', encoding='utf-8')
                    with open('Models/Results_level_2/Cross_validation_models_RNA_graph2vec_2_all_components_test_scores.txt', 'a') as file:
                        file.write('\n')