import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.gridspec import GridSpec
import seaborn as sns
import sys
import os
import copy
from adjustText import adjust_text

sys.path.append("../../scripts/")
import utils
import analysis

from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr
from scipy.stats import entropy
import warnings
warnings.filterwarnings("ignore")

#path_HLCA = "../../data/HLCA_core_h5ads/HLCA_v2.h5ad"
path_HLCA = "../data/HLCA_v2_intermediates/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_log1p.h5ad"
dir_results = "../data/results/variance_explained_by_covariates/"
dir_figures = "../data/results/figures/"

adata = sc.read(path_HLCA)
covariates = [
    "dataset",
    "tissue_dissociation_protocol",
    "log10_total_counts",
    "mito_frac",
    "3'_or_5'",
    "BMI",
    "cell_ranger_version_short",
    "cell_viability_%",
    "fresh_or_frozen",
    "sample_type",
    "sequencing_platform_short",
    "sex",
    "single_cell_platform",
    "smoking_status_num",
    "subject_type",
    "subject_ID",
    "anatomical_region_ccf_score",
    "nose",
    "age",
]

for cov in covariates:
    if type(adata.obs[cov].values) == pd.Categorical:
        print(cov.upper())
        print(set(adata.obs[cov]))
        print("\n")
    else:
        print(cov.upper(), "not a categorical.\n")

int_types = ["integrated"]  # list needs to include "integrated" and/or "unintegrated" 
if "unintegrated" in int_types:
    n_pcs = 30

outdir_1 = os.path.join(dir_results, f"variance_explained_fractions/")
if not os.path.exists(outdir_1):
    print("Creating directory:", outdir_1)
    os.makedirs(outdir_1)
outdir_2 = os.path.join(dir_results, f"samples_included/")
if not os.path.exists(outdir_2):
    print("Creating directory:", outdir_2)
    os.mkdir(outdir_2)

# Initiate a dictionary, in which we will store which samples were included per single regression
samples_included = dict()
cts_to_skip = list()

cell_type_key = "leiden_3" # change to cell type annotation key
verbose = True
min_n_cells_total = 50  # in total
min_n_cells = 10  # per sample
min_n_samples = 2

for int_type in int_types:
    samples_included[int_type] = dict()

    for subset in sorted(adata.obs[cell_type_key].unique()) + ["whole_atlas"]:
        subset_no_space = subset.replace(" ", "_")
        samples_included[int_type][subset] = pd.DataFrame(index=adata.obs["sample"].unique(), columns=covariates)

        
        if os.path.isfile(os.path.join(dir_results, f"variance_explained_fractions/variance_explained_fractions_{subset_no_space}_{int_type}.csv",) ):
            continue
        print(f"Working on {int_type},  {subset}...")
        if subset == "whole_atlas":
            subadata = adata.copy()
            verbose = True
        else:
            subadata = adata[adata.obs[cell_type_key] == subset, :].copy()
            verbose = False

        if subadata.n_obs < min_n_cells_total:
            print(f"{subset} has fewer than {min_n_cells_total} cells! Skipping.")
            cts_to_skip.append(subset)
            continue

        if int_type == "unintegrated":
            emb_name = "X_pca"
            sc.tl.pca(subadata, n_comps=n_pcs, use_highly_variable=True)
        elif int_type == "integrated":
            emb_name = "X_scanvi_emb"
        else:
            raise ValueError("emb_name should be set either to 'integrated' or 'unintegrated'")
        

        # store the number of components in our embedding of choice
        n_comps = subadata.obsm[emb_name].shape[1]
        # initiate a dataframe in which we'll store the variance explained by each covariate, plus the total variance ("overall") observed
        var_explained = pd.DataFrame(index=range(n_comps), columns=covariates + ["overall"])


        # initiate a dataframe in which we will store the data for our linear regression (i.e. the PC/latent components, + covariates).
        # Rows are cells, but we will collapse this to samples below
        comp_sample_df = pd.DataFrame(index=subadata.obs.index)
        comp_sample_df["sample"] = subadata.obs["sample"]
        # prepare aggregation dictionary for collapsing into sample-wise observations
        agg_dict = {"sample": "count"}  # this will be number of cells
        for comp in range(n_comps):
            # store component scores per cell
            comp_sample_df[f"comp{comp}"] = subadata.obsm[emb_name][:, comp]
            # we will aggregate these later by taking the mean per sample
            agg_dict[f"comp{comp}"] = "mean"
        for cov in covariates:
            if cov in ["log10_total_counts", "mito_frac"]:
                # store values
                comp_sample_df[cov] = subadata.obs[cov]
                # we will aggregate by taking the mean
                agg_dict[cov] = "mean"
            else:
                # for all other covariates: these are sample-level covariates, so we will take the "first" overservation in the sample (which should be the only)
                comp_sample_df[cov] = subadata.obs[cov]
                agg_dict[cov] = "first"


        # now collapse into sample-level observations
        sample_df = comp_sample_df.groupby("sample").agg(agg_dict).rename(columns={"sample": "n_cells"})

        
        # filter out samples with fewer than min_n_cells cells of the cell type
        sample_df = sample_df.loc[sample_df.n_cells >= min_n_cells,].copy()
        # check number of samples left. If fewer than min_n_samples remain, we will skip the cell type
        if sample_df.shape[0] < min_n_samples:
            print(f"Only {sample_df.shape[0]} samples available for {subset}. Skipping.")
            cts_to_skip.append(subset)
            continue
        # Otherwise, move on to the linear regression: do a linear regression on each component, with the component scores as response variable...


        for comp in range(n_comps):
            # store the component values (for all samples i.e. unfiltered)
            y_true_unfiltered = sample_df.loc[:, f"comp{comp}"].values
            # and store variance of y_true as "overall" variance
            var_explained.loc[f"comp{comp}", "overall"] = np.var(y_true_unfiltered)
            # and the covariate as fixed variable
            for cov in covariates:
                # store covariate observations under x
                x = sample_df[cov].values.copy()
                # store samples to which they match
                x_samples = sample_df.index
                # check which of these samples have no observation (e.g.because BMI was unknown, or age, etc.) (the function used below checks for different kinds of nas, e.g. np.nan, "nan", None, "None" etc.)
                x_nans = np.vectorize(utils.check_if_nan)(x)
                # now keep only xs that have real observations
                x = x[~x_nans]
                # if only one or no observations are left, skip this covariate
                if len(x) < 2:
                    continue
                # filter samples according to x filtering
                x_samples = x_samples[~x_nans]
                # and store which samples were included in our samples_included dictionary, for later reference (this is our "n")
                samples_included[int_type][subset][cov] = samples_included[int_type][subset].index.isin(x_samples.tolist())
                # filter y_true according to x's filtering
                y_true = y_true_unfiltered[~x_nans].reshape(-1, 1)

                
                # prepare x for linear regression: if it is a float (e.g. BMI, age), all we need to do is reshape:
                if x.dtype in ["float32", "float", "float64"]:
                    x = x.reshape(-1, 1)
                    # print that we are treating as numerical (only for first comp, so that we don't print the same thing many times)
                    if comp == 0 and verbose:
                        print(f"treating {cov} as continuous variable")
                # otherwise we are dealing with a categorical...
                else:
                    # if it has only one category, there is 0 variance and we cannot perform linear regression. In that case, move on to the next covariate.
                    if len(set(x)) == 1:
                        var_explained.loc[comp, cov] = np.nan
                        continue
                    # Otherwise, convert x to dummied variable: print that we are converting to dummy (only do it for the first comp, otherwise we print the same thing many times)
                    if comp == 0 and verbose:
                        print(f"converting {cov} to dummy variable")
                    # drop_first means we ensure that are x is full rank, and we only encode all-1 categories
                    x = pd.get_dummies(x, drop_first=True)

                
                # now perform linear regression
                lrf = LinearRegression(fit_intercept=True).fit(x, y_true,)
                # predict y based on the fit linear model
                y_pred = lrf.predict(x)
                # and store the variance of the predicted y, this is the "variance explained" by the covariate, for this component
                var_explained.loc[comp, cov] = np.var(y_pred)

            
        # for each covariate, sum up how much variance it explains across the components (i.e. PCs or scANVI latent components)
        # Sort covariates from explaining most to explaining least
        total_variance_explained = np.sum(var_explained, axis=0).sort_values(ascending=False)

        # divide this by the total variance that was observed in the components, to get fraction of variance explained
        total_variance_explained_fractions = (total_variance_explained / total_variance_explained["overall"])

        # write to files:
        # 1) variance explained fractions, for this integration type and cell type
        total_variance_explained_fractions.to_csv(os.path.join(dir_results, f"variance_explained_fractions/variance_explained_fractions_{subset_no_space}_{int_type}.csv",))

        # 2) samples included, for this cell type
        samples_included[int_type][subset].fillna(False).to_csv( os.path.join(dir_results, f"samples_included/samples_included_{subset_no_space}.csv", ))
