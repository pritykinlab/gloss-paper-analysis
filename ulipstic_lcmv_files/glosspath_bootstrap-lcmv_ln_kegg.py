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

datapath = '../joint_notebooks/datasets/lcmv_ln_data_aug_11_2024.h5ad'

resolutions = {
    'annot' : ['Cd4']
}

from Gloss.regresscv import RegressCV
from Gloss.regressbootstrap import RegressBootstrap

regcvpath = './experiment_evaluations/lcmv_ln_new_gloss_tuning_hallmark_oct24_2024.pickle'
with open(regcvpath, 'rb') as handle:
    myregcv = pickle.load(handle)

regb = RegressBootstrap(datapath, resolutions, 'kegg', myregcv.best_params, 100,
                  'sample count', 'total_counts', 'biotin_interaction', donors_profiled=True)

savepath = './experiment_evaluations/lcmv_ln_new_gloss_bootstrap_kegg_oct26_2024.pickle'

with open(savepath, 'wb') as handle:
    pickle.dump(regb, handle, protocol=4)