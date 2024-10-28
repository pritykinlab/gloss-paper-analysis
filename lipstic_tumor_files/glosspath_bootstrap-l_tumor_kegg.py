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

datapath = '../joint_notebooks/datasets/lipstic_tumor_data_aug_11_2024.h5ad'

resolutions = {
    'annotation_fine' : ['Mo/MF', 'cDC2', 'mRegDC2']
}

from Gloss.regresscv import RegressCV
from Gloss.regressbootstrap import RegressBootstrap

regcvpath = './experiment_evaluations/tumor_new_gloss_tuning_kegg_oct24_2024.pickle'
with open(regcvpath, 'rb') as handle:
    myregcv = pickle.load(handle)

regb = RegressBootstrap(datapath, resolutions, 'kegg', myregcv.best_params, 100,
                  'hash_max', 'nCount_RNA', 'biotin_raw')

savepath = './experiment_evaluations/tumor_new_gloss_bootstrap_kegg_oct24_2024.pickle'

with open(savepath, 'wb') as handle:
    pickle.dump(regb, handle, protocol=4)