import sys
import os
from importlib import reload

import scanpy as sc
import pandas as pd 
import numpy as np
import seaborn as sns

from scipy.stats import pearsonr

from sklearn.linear_model import ElasticNetCV, RidgeCV, Ridge
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_validate, KFold
from sklearn.metrics import make_scorer, r2_score

from collections import Counter 
from sklearn.utils import resample 

import statistics

import pickle

package_path = os.path.abspath("/Genomics/pritykinlab/tamjeed/github_packages/Gloss/")
sys.path.insert(0, package_path)

from Gloss.prepdata import PrepData

datapath = '../joint_notebooks/datasets/lipstic_tumor_data_aug_11_2024.h5ad'

prepped_data_hallmark = PrepData(datapath, 'hallmark', 'hash_max', 'nCount_RNA', 'biotin_raw')
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
    'annotation_fine' : ['Mo/MF', 'cDC2', 'mRegDC2']
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

with open('./experiment_evaluations/gloss_ridge_tuning_l_tumor_oct23_2024.pickle', 'rb') as handle:
    cell_type_res_dict = pickle.load(handle)

optimal_ridge_params = {}
for res in cell_type_res_dict:
    optimal_ridge_params[res] = {}
    for ctype in cell_type_res_dict[res]:
        ctype_response = cell_type_res_dict[res][ctype]
        best_alphas = [model.alpha_ for model in ctype_response['estimator']]
        optimal_ridge_params[res][ctype] = statistics.median(best_alphas)

r = 42
inner_cv = KFold(n_splits=5, shuffle=True, random_state=r)
outer_cv = KFold(n_splits=5, shuffle=True, random_state=r)

bootstrap_res_dict = {}
gene_list = list(gene_ad.var.index)
gene_list.append('log_RNA_libsize')
gene_list.append('log_sample_hashtag')

for res in resolutions:
    ctypes = resolutions[res]
    bootstrap_res_dict[res] = {}
    
    for ctype in ctypes:
        print(ctype)
        index = gene_ad.obs[res].isin([ctype])
        adata = gene_ad[index].copy()
        y = adata.obs['new_biotin']
        my_matrix = features_matrix[index, :]

        print('*** {} {} ***'.format(res, ctype))
        print(' ')
        
        models = cell_type_res_dict[res][ctype]['estimator']

        alphas = []
        for i in range(len(models)):
            alphas.append(models[i].alpha_)
        alphas_counter = Counter(alphas)
        most_common_alpha = alphas_counter.most_common(1)[0][0]

        # Sample data
        X = my_matrix
        y = y

        # Elastic Net parameters
        alpha = optimal_ridge_params[res][ctype]
        n_bootstraps = 100

        # List to store bootstrap coefficients
        bootstrap_coefficients = []

        # Bootstrapping
        for i in range(n_bootstraps):
            print(i)
            X_resampled, y_resampled = resample(X, y, random_state=i)
            ridge = Ridge(alpha=alpha)
            ridge.fit(X_resampled, y_resampled)
            bootstrap_coefficients.append(ridge.coef_)

        # Convert to DataFrame for analysis
        bootstrap_coefficients = np.array(bootstrap_coefficients)
        coefficients_df = pd.DataFrame(bootstrap_coefficients)
        coefficients_df.columns = gene_list

        means = coefficients_df.mean()

        # Sort the columns by mean value
        sorted_columns = means.sort_values().index
        # Reorder the DataFrame columns
        sorted_df = coefficients_df[sorted_columns]

        bootstrap_res_dict[res][ctype] = sorted_df
    
with open('./experiment_evaluations/gloss_ridge_bootstrap_l_tumor_oct23_2024.pickle', 'wb') as handle:
    pickle.dump(bootstrap_res_dict, handle, protocol=4)