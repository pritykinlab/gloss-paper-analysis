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

gutpath = '../joint_notebooks/datasets/{}_gut_data_perturbed_itgae_v2_oct_26_2024.h5ad'.format(i)

resolutions = {
    'annotation' : ['CD4']
}

from Gloss.regresscv import RegressCV
from Gloss.regressbootstrap import RegressBootstrap

regcvpath = './experiment_evaluations/gut_new_gloss_tuning_hallmark_oct24_2024.pickle'
with open(regcvpath, 'rb') as handle:
    myregcv = pickle.load(handle)

regb = RegressBootstrap(gutpath, resolutions, 'hallmark', myregcv.best_params, 100,
                  'avg_sample_hto', 'total_counts', 'biotin_raw_perturbed')

savepath = './experiment_evaluations/{}_gut_perturbed_itgae_new_gloss_bootstrap_hallmark_oct26_2024.pickle'.format(i)

with open(savepath, 'wb') as handle:
    pickle.dump(regb, handle, protocol=4)