import sys
import os
from importlib import reload

# Add the path of your package
package_path = os.path.abspath("/Genomics/pritykinlab/tamjeed/github_packages/Gloss/")
sys.path.insert(0, package_path)

import scanpy as sc
import numpy as np
import pandas as pd
import pickle

gutpath = '../joint_notebooks/datasets/gut_data_aug_11_2024.h5ad'

resolutions = {
    'annotation' : ['CD4']
}

from Gloss.regresscv import RegressCV
from Gloss.regressbootstrap import RegressBootstrap

regcvpath = './experiment_evaluations/gut_new_gloss_tuning_hallmark_oct24_2024.pickle'
with open(regcvpath, 'rb') as handle:
    myregcv = pickle.load(handle)

regb = RegressBootstrap(gutpath, resolutions, 'hallmark', myregcv.best_params, 100,
                  'avg_sample_hto', 'total_counts', 'raw_biotin')

savepath = './experiment_evaluations/gut_new_gloss_bootstrap_hallmark_oct24_2024.pickle'

with open(savepath, 'wb') as handle:
    pickle.dump(regb, handle, protocol=4)