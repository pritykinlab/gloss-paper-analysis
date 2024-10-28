# gloss-paper-analysis
Companion to Gloss package, holding paper analysis.

The various folders hold code used to produce various results in the paper. The following primarily contain executable *.py files:

`ulipstic_gut_files` : For uLIPSTIC gut, code used for bootstrapping and tuning Gloss and baseline regression methods. Also includes all code for bootstrapping over all simulations made over this data subset.
`ulipstic_lcmv_files` : For uLIPSTIC LCMV Sys and LN, code used for bootstrapping and tuning Gloss and baseline regression methods. Also includes all code for bootstrapping over all simulations made over these subsets.
`lipstic_tumor_files` : For LIPSTIC v1 Tumor, code used for bootstrapping and tuning Gloss and baseline regression methods. Also includes all code for bootstrapping over all simulations made over these subsets.

The other folders primarily contain notebooks for various parts of the analysis.

`ulipstic_notebooks` : Notebooks used for perturbing the uLIPSTIC data as described in the paper.
`lipstic_tumor_notebooks` : Notebooks used for perturbing the LIPSTIC v1 data as described in the paper.
`processing_plotting_enrichment` : Where the experiments are primarily analyzed and plotted, holds all of the enrichment-based analysis for both non-perturbed and perturbed data, as well as performance and coefficient/pathway variation plots.

