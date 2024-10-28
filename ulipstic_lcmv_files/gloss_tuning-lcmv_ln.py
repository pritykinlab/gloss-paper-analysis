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

regcv = RegressCV(datapath, resolutions, 'hallmark', 
                  'sample count', 'total_counts', 'biotin_interaction', 
                  group_regs = [0.004, 0.005, 0.006, 0.007], single_gene_regs=[1,2,3,4,5], donors_profiled=True)

savepath = './experiment_evaluations/lcmv_ln_new_gloss_tuning_hallmark_oct24_2024.pickle'

with open(savepath, 'wb') as handle:
    pickle.dump(regcv, handle, protocol=4)