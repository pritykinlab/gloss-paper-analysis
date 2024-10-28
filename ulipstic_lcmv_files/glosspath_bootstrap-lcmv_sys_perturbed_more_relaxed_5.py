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

import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description="A script to demonstrate argument passing.")
    
# Add argument for parameter
parser.add_argument('--param', type=int, required=True, help='The parameter value to pass to the function')

# Parse the arguments from the command line
args = parser.parse_args()

i = args.param

datapath = '../joint_notebooks/datasets/{}_lcmv_sys_data_perturbed_combined_c_n_more_relaxed_5_v2_oct_26_2024.h5ad'.format(i)

resolutions = {
    'annot v2' : ['macrophage']
}

from Gloss.regresscv import RegressCV
from Gloss.regressbootstrap import RegressBootstrap

regcvpath = './experiment_evaluations/lcmv_sys_new_gloss_tuning_hallmark_oct24_2024.pickle'
with open(regcvpath, 'rb') as handle:
    myregcv = pickle.load(handle)

regb = RegressBootstrap(datapath, resolutions, 'hallmark', myregcv.best_params, 100,
                  'sample count', 'total_counts', 'biotin_raw_perturbed', donors_profiled=True)

savepath = './experiment_evaluations/{}_lcmv_sys_perturbed_more_relaxed_5_new_gloss_bootstrap_hallmark_oct26_2024.pickle'.format(i)

with open(savepath, 'wb') as handle:
    pickle.dump(regb, handle, protocol=4)