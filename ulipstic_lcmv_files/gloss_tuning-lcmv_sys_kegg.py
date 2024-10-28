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

datapath = '../joint_notebooks/datasets/lcmv_sys_data_aug_11_2024.h5ad'

resolutions = {
    'annot v2' : ['macrophage']
}

from Gloss.regresscv import RegressCV

regcv = RegressCV(datapath, resolutions, 'kegg', 
                  'sample count', 'total_counts', 'biotin_interaction', donors_profiled=True)

savepath = './experiment_evaluations/lcmv_sys_new_gloss_tuning_kegg_oct24_2024.pickle'

with open(savepath, 'wb') as handle:
    pickle.dump(regcv, handle, protocol=4)