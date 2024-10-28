import sys
import os
from importlib import reload

import scanpy as sc
import pandas as pd 
import numpy as np
import seaborn as sns

from scipy.stats import pearsonr

from sklearn.linear_model import ElasticNetCV, RidgeCV, LassoCV
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_validate, KFold
from sklearn.metrics import make_scorer, r2_score

package_path = os.path.abspath("/Genomics/pritykinlab/tamjeed/github_packages/Gloss/")
sys.path.insert(0, package_path)

from Gloss.prepdata import PrepData

datapath = '../joint_notebooks/datasets/lcmv_ln_data_aug_11_2024.h5ad'

prepped_data_hallmark = PrepData(datapath, 'hallmark', 'sample count', 'total_counts', 'biotin_interaction', donors_profiled=True)
gene_ad = prepped_data_hallmark.adata
rna_logcounts = gene_ad.obs['log_RNA_libsize']
hash_logcounts = gene_ad.obs['log_sample_hashtag']

matrix = gene_ad.layers['log_scaled'].copy()

features_matrix = np.concatenate((matrix, 
                                  rna_logcounts.to_numpy().reshape(-1,1),
                                  hash_logcounts.to_numpy().reshape(-1,1),
                                  ), axis=1)

y = gene_ad.obs['new_biotin']

resolutions = {
    'annot' : ['Cd4']
}

def pearson_corr(y, y_pred):
    return pearsonr(y, y_pred)[0]
def pearson_sig(y, y_pred):
    return pearsonr(y, y_pred)[1]

pearson = make_scorer(pearson_corr)
pearson_2 = make_scorer(pearson_sig)
r2 = make_scorer(r2_score)

score_dict = {
    'pearson' : pearson,
    'pearson_sig' : pearson_2,
    'r2_score' : r2
}

cell_type_res_dict = {}

r = 42
inner_cv = KFold(n_splits=5, shuffle=True, random_state=r)
outer_cv = KFold(n_splits=5, shuffle=True, random_state=r)

for res in resolutions:
    ctypes = resolutions[res]
    cell_type_res_dict[res] = {}
    
    for ctype in ctypes:
        print(ctype)
        index = gene_ad.obs[res].isin([ctype])
        adata = gene_ad[index].copy()
        y = adata.obs['new_biotin']
        my_matrix = features_matrix[index, :]

        print('*** {} {} ***'.format(res, ctype))
        print(' ')
        
        l1_ratios = [0.1, 0.5, 0.7, 0.9]
        alphas = [0.01, 0.03, 0.05, 0.07, 0.09, 0.2, 0.3, 0.4]
        model = ElasticNetCV(alphas=alphas, l1_ratio=l1_ratios,
                             cv=inner_cv, verbose=True, n_jobs=-1)
        
        regression_res = cross_validate(model, X=my_matrix, y=y, cv=outer_cv,
                    return_train_score=True,
                    return_estimator=True,
                    scoring=score_dict,
                    )
        cell_type_res_dict[res][ctype] = regression_res

cell_type_res_dict

import pickle
    
with open('./experiment_evaluations/gloss_elasticnet_tuning_lcmv_ln_oct23_2024.pickle', 'wb') as handle:
    pickle.dump(cell_type_res_dict, handle, protocol=4)