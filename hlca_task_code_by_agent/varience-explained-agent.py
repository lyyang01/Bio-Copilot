# Step 1: Loading the Integrated Single-Cell Data File

# Environment Preparation
# Ensure Python 3.7 or higher is installed and accessible
# Install scanpy library if not already installed. This can be done outside this script using:
# pip install scanpy

# Import Required Libraries
import scanpy as sc
import os

# Set File Path
file_path = '../data/final_adata_with_subclusters.h5ad'

# Check if the file exists before proceeding
if os.path.exists(file_path):
    # Load the Data File
    adata = sc.read_h5ad(file_path)
    
    # Verify Data Loading
    print(adata)
    print(adata.var_names)
    
    # Initial Data Exploration (Optional)
    # View metadata associated with each cell
    print(adata.obs.head())
    # View the data matrix as a pandas DataFrame
    #print(adata.to_df().head())
else:
    print(f"File {file_path} does not exist. Please check the file path and permissions.")

import pandas as pd
import numpy as np

# Initialize an empty dictionary to hold DataFrames for each cluster, including a 'whole_atlas' cluster
cluster_dataframes = {}

# Define function for extracting the first value
def extract_first_value(data):
    return data.iloc[0]

# Add 'whole_atlas' to the list of unique clusters
unique_clusters = ['whole_atlas'] + list(adata.obs['cluster_level_3'].unique())
skip = []
# Loop through clusters including 'whole_atlas'
for cluster in unique_clusters:

    if cluster == "whole_atlas":
        subadata = adata.copy()
        verbose = True
    else:
        subadata = adata[adata.obs['cluster_level_3'] == cluster, :].copy()
        verbose = False
    
    if subadata.n_obs < 50:
        #print(f"{subset} has fewer than {min_n_cells_total} cells! Skipping.")
        skip.append(cluster)
        continue
    # Initialize an empty list to hold data for all samples in the current cluster
    cluster_samples_data = []

    # For 'whole_atlas', consider all samples; otherwise, filter by the current cluster
    if cluster == 'whole_atlas':
        cluster_samples = adata.obs['sample']
        subdata = adata.copy()
        
    else:
        cluster_samples = adata.obs[adata.obs['cluster_level_3'] == cluster]['sample']
        subdata = adata[adata.obs['cluster_level_3'] == cluster, :].copy()
        
    for sample in cluster_samples.unique():
        # Filter data for the current sample
        
        sample_data = subdata[subdata.obs['sample']==sample]
        # Use n_cells to filter out samples containing less than 10 cells
        
        if len(sample_data.obs['sample']) < 10:
            continue
        
        # Aggregate embedding representations (X_scANVI) using np.mean with axis=0
        X_scANVI_mean = np.mean(sample_data.obsm['X_scANVI'], axis=0)
        
        # Extract and aggregate metadata covariates
        log10_total_counts_mean = np.mean(sample_data.obs['log10_total_counts'])
        mito_frac_mean = np.mean(sample_data.obs['mito_frac'])
        
        # Extract first available values for other covariates
        covariates = {
            'sample': sample,  # Ensure 'sample' is included as a key
            #'log10_total_counts_mean': log10_total_counts_mean,
            'log10_total_counts': log10_total_counts_mean,
            #'mito_frac_mean': mito_frac_mean,
            'mito_frac': mito_frac_mean,
            'n_cells': len(sample_data),
            'dataset': extract_first_value(sample_data.obs['dataset']),
            'tissue_dissociation_protocol': extract_first_value(sample_data.obs['tissue_dissociation_protocol']),
            "3'_or_5'": extract_first_value(sample_data.obs["3'_or_5'"]),
            'BMI': extract_first_value(sample_data.obs['BMI']),
            'cell_ranger_version_short': extract_first_value(sample_data.obs['cell_ranger_version_short']),
            'cell_viability_%': extract_first_value(sample_data.obs['cell_viability_%']),
            'fresh_or_frozen': extract_first_value(sample_data.obs['fresh_or_frozen']),
            'sample_type': extract_first_value(sample_data.obs['sample_type']),
            'sequencing_platform_short': extract_first_value(sample_data.obs['sequencing_platform_short']),
            'sex': extract_first_value(sample_data.obs['sex']),
            'single_cell_platform': extract_first_value(sample_data.obs['single_cell_platform']),
            'smoking_status_num': extract_first_value(sample_data.obs['smoking_status_num']),
            'subject_type': extract_first_value(sample_data.obs['subject_type']),
            'subject_ID': extract_first_value(sample_data.obs['subject_ID']),
            'anatomical_region_ccf_score': extract_first_value(sample_data.obs['anatomical_region_ccf_score']),
            'nose': extract_first_value(sample_data.obs['nose']),
            'age': extract_first_value(sample_data.obs['age'])
        }
        
        # Initialize a row for aggregated data
        aggregated_row = covariates
        
        # Add X_scANVI components to the aggregated row
        for i, component in enumerate(X_scANVI_mean):
            aggregated_row[f'X_scANVI_component_{i}'] = component
        
        # Append the aggregated row to the cluster_samples_data list
        cluster_samples_data.append(aggregated_row)
    
    # Convert the list of dictionaries to a DataFrame
    if len(cluster_samples_data) >= 2:
    #if cluster_samples_data:  # Check if the list is not empty
        cluster_df = pd.DataFrame(cluster_samples_data)
        # Correctly set 'sample' as the index using the 'sample' key from each dictionary
        cluster_df.set_index('sample', inplace=True)
        # Store the DataFrame in the cluster_dataframes dictionary with the cluster as the key
        cluster_dataframes[cluster] = cluster_df
    else:
        skip.append(cluster)
        continue

# Now, you can access the aggregated data for each cluster from the cluster_dataframes dictionary
# For example, to access the data for the 'whole_atlas':
print(f"Data for 'whole_atlas':")
print(cluster_dataframes['whole_atlas'].head())


# Step 3: Final Modifications and Enhancements

# A. Setting Up for Variance Analysis with Linear Regression

# 1. Initialization:
# Import necessary libraries
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import numpy as np
import os

# Ensure the output directory exists
output_dir = './variance_explain/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


# 2. Prepare the Results DataFrame
num_components = 30  # Adjusted to 30 components in "X_scANVI_component_*"


covariates = ['log10_total_counts', 'mito_frac', 'n_cells', 'dataset', 
              'tissue_dissociation_protocol', "3'_or_5'", 'BMI', 'cell_ranger_version_short', 
              'cell_viability_%', 'fresh_or_frozen', 'sample_type', 'sequencing_platform_short', 
              'sex', 'single_cell_platform', 'smoking_status_num', 'subject_type', 
              'subject_ID', 'anatomical_region_ccf_score', 'nose', 'age']  # Actual covariate names


# B. Performing Variance Analysis at Cluster Level
# 1. Iterate Over Clusters
for cluster, df in cluster_dataframes.items():
    X_columns = covariates  # Directly use covariates
    Y_columns = ['X_scANVI_component_' + str(i) for i in range(num_components)]
    
    cluster_results = pd.DataFrame(index=Y_columns, columns=X_columns + ['overall_variance'])
    used_rows_tracker = pd.DataFrame(False, index=df.index, columns=covariates)
    
    for y_col in Y_columns:
        Y = df[y_col].dropna()  # Response: the current component of X_scANVI_component_*, dropping NAN values
        overall_variance = np.var(Y, ddof=1)  # Calculate overall variance of Y
        for x_col in X_columns:
            X = df[[x_col]].dropna()  # Predictor: each metadata covariate, dropping rows with NAN in X
            
            # Continue processing X's string values
            if X[x_col].dtype == 'object':
                X = pd.get_dummies(X[x_col])
            
            if X.shape[0] > 1 and Y.shape[0] > 0:  # Ensure both X and Y have data
                model = LinearRegression()                
                Y_ = Y.loc[X.index]
                model.fit(X, Y_)
                Y_pred = model.predict(X)
                
                variance_explained = np.var(Y_pred, ddof=1)
                cluster_results.loc[y_col, x_col] = variance_explained
                
                used_rows_tracker.loc[X.index, x_col] = True
            else:
                cluster_results.loc[y_col, x_col] = np.nan
                
        cluster_results.loc[y_col, 'overall_variance'] = overall_variance
    
    # 2. Sum up variance explained across components for each covariate and sort
    try:
        variance_sum = np.sum(cluster_results, axis=0).sort_values(ascending=False)
        total_variance_observed = cluster_results['overall_variance'].sum()
        fraction_variance_explained = variance_sum / total_variance_observed
    except:
        pass
    
    # 3. Save the fraction of variance explained and used_rows_tracker for each cluster
    fraction_variance_explained.to_csv(output_dir + f'{cluster}_fraction_variance_explained.csv')
    used_rows_tracker.to_csv(output_dir + f'{cluster}_used_rows_tracker.csv')

    print(f"Fraction of variance explained for cluster {cluster}:")
    #print(fraction_variance_explained)
    print(f"Rows used in regression for cluster {cluster}:")
    #print(used_rows_tracker)
