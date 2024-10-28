import sys
import os
from importlib import reload

# Add the path of your package
package_path = os.path.abspath("/Genomics/pritykinlab/tamjeed/github_packages/GlossPath/")
sys.path.insert(0, package_path)

import scanpy as sc
import numpy as np
import pandas as pd
import pickle

datapath = '../joint_notebooks/datasets/0_lipstic_tumor_data_perturbed_v3_oct_19_2024.h5ad'

resolutions = {
    'annotation_fine' : ['Mo/MF', 'cDC1', 'cDC2', 'mRegDC1', 'mRegDC2']
}

from GlossPath.regresscv import RegressCV

regcv = RegressCV(datapath, resolutions, 'hallmark', 
                  'hash_max', 'nCount_RNA', 'biotin_raw_perturbed')

savepath = './experiment_evaluations/tumor_perturbed_glosspath_tuning_hallmark_v3_oct19_2024.pickle'

with open(savepath, 'wb') as handle:
    pickle.dump(regcv, handle, protocol=4)
