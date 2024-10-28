import sys
import os
from importlib import reload

# Add the path of your package
package_path = os.path.abspath("/Genomics/pritykinlab/tamjeed/github_packages/Gloss/")
sys.path.insert(0, package_path)
package_path = os.path.abspath("/Genomics/pritykinlab/tamjeed/github_packages/GlossPath/")
sys.path.insert(0, package_path)

import scanpy as sc
import pandas as pd 
import numpy as np
import seaborn as sns

from sklearn.linear_model import ElasticNet, ElasticNetCV, Ridge, RidgeCV, Lasso, LassoCV, ElasticNet
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_validate, KFold
from sklearn.metrics import make_scorer, r2_score, mean_squared_error
import matplotlib.pyplot as plt

import pickle

import statistics
import skglm

from Gloss.prepdata import PrepData
from Gloss.regressor import Regressor

## function definitions

datapath = '../joint_notebooks/datasets/lcmv_sys_data_aug_11_2024.h5ad'

resolutions = {
    'annot v2' : ['macrophage']
}

prepped_data_hallmark = PrepData(datapath, 'hallmark', 'sample count', 'total_counts', 'biotin_interaction', donors_profiled=True)
prepped_data_kegg = PrepData(datapath, 'kegg', 'sample count', 'total_counts', 'biotin_interaction', donors_profiled=True)

with open('../model_ulipstic_lcmv/experiment_evaluations/lcmv_sys_new_gloss_tuning_kegg_oct24_2024.pickle', 'rb') as handle:
    kegg_g_cv = pickle.load(handle)
    
with open('../model_ulipstic_lcmv/experiment_evaluations/lcmv_sys_new_gloss_tuning_hallmark_oct24_2024.pickle', 'rb') as handle:
    hallmark_g_cv = pickle.load(handle)
    
with open('../model_ulipstic_lcmv/experiment_evaluations/gloss_lasso_tuning_lcmv_sys_oct23_2024.pickle', 'rb') as handle:
    lasso_dict = pickle.load(handle)
    
with open('../model_ulipstic_lcmv/experiment_evaluations/gloss_elasticnet_tuning_lcmv_sys_oct23_2024.pickle', 'rb') as handle:
    enet_dict = pickle.load(handle)
    
with open('../model_ulipstic_lcmv/experiment_evaluations/gloss_ridge_tuning_lcmv_sys_oct23_2024.pickle', 'rb') as handle:
    ridge_dict = pickle.load(handle)

optimal_hallmark_params = hallmark_g_cv.best_params

optimal_kegg_params = kegg_g_cv.best_params

optimal_ridge_params = {}
for res in ridge_dict:
    optimal_ridge_params[res] = {}
    for ctype in ridge_dict[res]:
        ctype_response = ridge_dict[res][ctype]
        best_alphas = [model.alpha_ for model in ctype_response['estimator']]
        optimal_ridge_params[res][ctype] = statistics.median(best_alphas)

optimal_lasso_params = {}
for res in lasso_dict:
    optimal_lasso_params[res] = {}
    for ctype in lasso_dict[res]:
        ctype_response = lasso_dict[res][ctype]
        best_alphas = [model.alpha_ for model in ctype_response['estimator']]
        optimal_lasso_params[res][ctype] = statistics.median(best_alphas)

optimal_enet_params = {}
for res in enet_dict:
    optimal_enet_params[res] = {}
    for ctype in enet_dict[res]:
        ctype_response = enet_dict[res][ctype]
        best_alphas = [model.alpha_ for model in ctype_response['estimator']]
        best_l1_ratios = [model.l1_ratio_ for model in ctype_response['estimator']]

        optimal_enet_params[res][ctype] = (statistics.median(best_l1_ratios), statistics.median(best_alphas))

r2 = make_scorer(r2_score)

gene_ad = prepped_data_hallmark.adata
rna_logcounts = gene_ad.obs['log_RNA_libsize']
hash_logcounts = gene_ad.obs['log_sample_hashtag']

matrix = gene_ad.layers['log_scaled'].copy()

features_matrix = np.concatenate((matrix, 
                                  rna_logcounts.to_numpy().reshape(-1,1),
                                  hash_logcounts.to_numpy().reshape(-1,1),
                                  ), axis=1)

comp_res = {}
for res in resolutions:
    comp_res[res] = {}
    for ctype in resolutions[res]:
        
        hallmark_test_res = []
        hallmark_train_res =[]
        
        kegg_test_res = []
        kegg_train_res =[]

        ridge_ftrain_res = []
        ridge_ftest_res = []

        elasticnet_ftrain_res = []
        elasticnet_ftest_res = []

        lasso_ftrain_res = []
        lasso_ftest_res = []
        
        print(ctype)
        index = gene_ad.obs[res].isin([ctype])
        adata = gene_ad[index].copy()
        y = adata.obs['new_biotin']

        my_hallmark_matrix = prepped_data_hallmark.subset(res, ctype)[0]
        my_kegg_matrix = prepped_data_kegg.subset(res, ctype)[0]
        my_features_matrix = features_matrix[index, :]
        
        for rand in range(50):
            print(rand)
            x_htrain, x_htest, y_htrain, y_htest = train_test_split(my_hallmark_matrix, y, test_size=0.2, random_state=rand)
            x_ktrain, x_ktest, y_ktrain, y_ktest = train_test_split(my_kegg_matrix, y, test_size=0.2, random_state=rand)
            x_ftrain, x_ftest, y_ftrain, y_ftest = train_test_split(my_features_matrix, y, test_size=0.2, random_state=rand)
            
            hallmark_gl = Regressor(prepped_data_hallmark.genes, 'hallmark', 
                                    group_reg=optimal_hallmark_params[res][ctype]['group_reg'],
                                    single_gene_reg=optimal_hallmark_params[res][ctype]['single_gene_reg'])
            hallmark_gl = hallmark_gl.fit(x_htrain, y_htrain)

            hallmark_train_res.append(r2(hallmark_gl, x_htrain, y_htrain))
            hallmark_test_res.append(r2(hallmark_gl, x_htest, y_htest))
            
            kegg_gl = Regressor(prepped_data_kegg.genes, 'kegg', 
                                group_reg=optimal_kegg_params[res][ctype]['group_reg'],
                                single_gene_reg=optimal_kegg_params[res][ctype]['single_gene_reg'])
            kegg_gl = kegg_gl.fit(x_ktrain, y_ktrain)

            kegg_train_res.append(r2(kegg_gl, x_ktrain, y_ktrain))
            kegg_test_res.append(r2(kegg_gl, x_ktest, y_ktest))
            
        methods_test = {
            'hallmark' : hallmark_test_res,
            'kegg' : kegg_test_res,
        }
        
        methods_train = {
            'hallmark' : hallmark_train_res,
            'kegg' : kegg_train_res,
        }
        
        comp_res[res][ctype] = [ methods_test, methods_train ]

with open('./experiment_evaluations/model_arena_new_oct24_2024.pickle', 'wb') as handle:
    pickle.dump(comp_res, handle, protocol=4)