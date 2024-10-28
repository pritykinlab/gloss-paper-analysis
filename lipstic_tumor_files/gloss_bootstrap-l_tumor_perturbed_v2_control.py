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

datapath = '../joint_notebooks/datasets/{}_lipstic_tumor_data_perturbed_filtered_v2_oct_25_2024.h5ad'.format(i)

resolutions = {
    'annotation_fine' : ['Mo/MF', 'cDC2', 'mRegDC2']
}

from Gloss.regresscv import RegressCV
from Gloss.regressbootstrap import RegressBootstrap

regcvpath = './experiment_evaluations/tumor_new_gloss_tuning_hallmark_oct24_2024.pickle'
with open(regcvpath, 'rb') as handle:
    myregcv = pickle.load(handle)

regb = RegressBootstrap(datapath, resolutions, 'hallmark', myregcv.best_params, 100,
                  'hash_max', 'nCount_RNA', 'biotin_raw')

savepath = './experiment_evaluations/{}_tumor_perturbed_filtered_v2_control_new_glosspath_bootstrap_hallmark_oct26_2024.pickle'.format(i)

with open(savepath, 'wb') as handle:
    pickle.dump(regb, handle, protocol=4)