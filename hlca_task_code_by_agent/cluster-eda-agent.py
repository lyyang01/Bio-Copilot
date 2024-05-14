# Step 1: Loading the Integrated Single-Cell Dataset Using Scanpy

# Preparation
# Ensure scanpy is installed. If not, you can uncomment the next line to install it.
# !pip install scanpy

import scanpy as sc
import os

# Define the file path for the .h5ad file
file_path = '/data/yangliu/llm-agent/data/integrated_dataset_copy.h5ad'
# Check if the file exists before proceeding
if os.path.exists(file_path):
    # Loading the Dataset
    adata = sc.read_h5ad(file_path)
    # Verification of Successful Loading
    print(adata)
    # Check for the presence of the X_scANVI attribute within the adata object
    if 'X_scANVI' in adata.obsm.keys():
        print('X_scANVI is present')
    else:
        print('X_scANVI is missing')
else:
    print(f"The file {file_path} does not exist. Please check the file path.")

# Step 2: Initial Clustering and Visualization (Level 1) with Adjustments

# Check if cluster output result file exists at the start of the code
output_csv_path = 'cluster_level_1_assignments.csv'
output_h5ad_path = 'adata_with_clusters.h5ad'

if os.path.exists(output_csv_path) and os.path.exists(output_h5ad_path):
    print("Cluster output result file and AnnData object with clusters already exist.")
else:
    # A. Preparing Data and Computing Neighbors
    # Verify Data Integrity
    # Assuming 'adata' is loaded from previous steps and contains 'X_scANVI'
    if 'X_scANVI' not in adata.obsm.keys():
        raise ValueError("X_scANVI is missing from the AnnData object. Please ensure it's correctly loaded.")
    # Computing Neighbors
    # Use 'X_scANVI' for computing the neighborhood graph with specified n_neighbors
    sc.pp.neighbors(adata, use_rep='X_scANVI', n_neighbors=30, metric='euclidean')
    # B. Clustering with scanpy.tl.leiden
    # Using Leiden algorithm for clustering with a specified resolution parameter
    sc.tl.leiden(adata, resolution=0.01)
    # Finalizing and Saving Clustering Results
    # Saving the clustering results into 'adata.obs'
    adata.obs['cluster_level_1'] = adata.obs['leiden']
    # C. Saving the Clustering Results
    # Exporting Cluster Assignments
    if not os.path.exists(output_csv_path):
        adata.obs['cluster_level_1'].to_csv(output_csv_path)
    # Saving the AnnData Object
    if not os.path.exists(output_h5ad_path):
        adata.write(output_h5ad_path)
    # D. Visualization with UMAP
    # Visualizing the Clusters
    sc.tl.umap(adata)
    sc.pl.umap(adata, color='cluster_level_1', save='level_1_clustering.png')
    # Note: The 'save' parameter in 'sc.pl.umap' automatically saves the plot to a file.
    # The file will be saved in the current working directory with the specified filename.

#Step 3
import pandas as pd
import numpy as np

# Define the file path for the .h5ad file with initial clusters
file_path = 'adata_with_clusters.h5ad'
final_h5ad_path = 'final_adata_with_subclusters.h5ad'

# Check if the final h5ad file exists before proceeding with the loop
if os.path.exists(final_h5ad_path):
    print(f"The final h5ad file {final_h5ad_path} already exists. Skipping the clustering process.")
else:
    # Loading the Dataset with initial clusters
    adata = sc.read_h5ad(file_path)
    # Preparation for iterative sub-clustering
    current_level = 1  # Initialize current_level to track the clustering level
    for current_level in range(2, 5):  # Iterates through levels 2, 3, and 4
        # Determine unique clusters from the previous level
        previous_clusters = adata.obs[f'cluster_level_{current_level-1}'].unique()  
        # Placeholder for storing new sub-cluster labels
        new_sub_clusters = pd.Series(index=adata.obs.index, dtype='object') 
        for cluster in previous_clusters:        
            # Subset the data for cells in the current cluster
            subset_adata = adata[adata.obs[f'cluster_level_{current_level-1}'] == cluster]        
            # Recalculate the nearest-neighbor graph for the subset
            sc.pp.neighbors(subset_adata, use_rep='X_scANVI', n_neighbors=30, metric='euclidean')     
            # Apply the sub-clustering algorithm (Leiden) with adjusted resolution
            #resolution = 0.01 * current_level  # Example of adjusting resolution
            resolution = 0.1
            sc.tl.leiden(subset_adata, resolution=resolution)          
            # Store the sub-clustering results in new_sub_clusters
            new_labels = [f'{cluster}_{label}' for label in subset_adata.obs['leiden']]
            new_sub_clusters[subset_adata.obs.index] = new_labels   
        # Update adata.obs with new sub-cluster labels
        adata.obs[f'cluster_level_{current_level}'] = new_sub_clusters      
        # Save the cluster results to a csv file
        output_csv_path = f'cluster_level_{current_level}_assignments.csv'
        adata.obs[f'cluster_level_{current_level}'].to_csv(output_csv_path)     
        
    # After finishing the loop, save the final adata to a h5ad file
    adata.write(final_h5ad_path)
    print(f"Final adata with subclusters saved to {final_h5ad_path}.")


import os
import scanpy as sc
import pandas as pd

# Load the anndata file containing all clustering information
adata = sc.read_h5ad('final_adata_with_subclusters.h5ad')

# Identify all unique clustering levels and their corresponding cluster labels
clustering_levels = [col for col in adata.obs.columns if col.startswith('cluster_level_')]
unique_clusters = {level: adata.obs[level].unique() for level in clustering_levels}

# Initialize a dictionary to store consolidated marker genes results for each level
consolidated_marker_genes = {}
# Initialize a set to track processed fine-grained clusters
processed_clusters = set()

#consider all
marker_tpye = 'all'
# Calculate Marker Genes Considering Cluster Depth
for level, clusters in unique_clusters.items():
            consolidated_marker_genes[level] = pd.DataFrame()
            #if mark_ref == 'all' and level=='cluster_level_2':
            #    continue
            for cluster in clusters:
                # Determine the cluster depth
                cluster_depth = cluster.count('_') + 1
                
                # Identify fine-grained clusters if cluster depth is more than 1
                if cluster_depth > 1:
                    parent = "_".join(cluster.split("_")[: cluster_depth - 1])
                    fine_grained_clusters = [c for c in clusters if c.startswith(parent)]
                else:
                    fine_grained_clusters = [cluster]
                    #adata_ = adata
                    
                
                # Check if fine-grained clusters have already been processed
                if any(fg_cluster in processed_clusters for fg_cluster in fine_grained_clusters):
                    continue  # Skip the current loop if any of the fine-grained clusters have been processed
                
                # Mark fine-grained clusters as processed
                processed_clusters.update(fine_grained_clusters)
                
                # Calculate marker genes
                sc.tl.rank_genes_groups(adata, groupby=level, groups=fine_grained_clusters, method='t-test', n_genes=100)
                sc.tl.filter_rank_genes_groups(adata, min_in_group_fraction=0.25, max_out_group_fraction=0.5, min_fold_change=1)
               
                # Extract and merge marker genes results
                for fg_cluster in fine_grained_clusters:
                    result_df = pd.DataFrame({
                        'genes': adata.uns['rank_genes_groups']['names'][fg_cluster],
                        'scores': adata.uns['rank_genes_groups']['scores'][fg_cluster],
                        'pvals': adata.uns['rank_genes_groups']['pvals'][fg_cluster],
                        'pvals_adj': adata.uns['rank_genes_groups']['pvals_adj'][fg_cluster],
                        'logfoldchanges': adata.uns['rank_genes_groups']['logfoldchanges'][fg_cluster],
                        'cluster': fg_cluster
                    })
                    consolidated_marker_genes[level] = pd.concat([consolidated_marker_genes[level], result_df], ignore_index=True)
            
            # Create directory for each level of clustering
            output_dir = f'cluster_level_{level[-1]}'  # Extract level number from level string
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            
            # Save the consolidated marker gene list for the level
            output_file = f'{output_dir}/consolidated_marker_genes_{level}_{marker_tpye}.csv'
            consolidated_marker_genes[level].to_csv(output_file, index=False)
            print(f'Saved consolidated marker genes for {level} to {output_file}.')

# consider sister
marker_tpye = 'sister'
# Identify all unique clustering levels and their corresponding cluster labels
clustering_levels = [col for col in adata.obs.columns if col.startswith('cluster_level_')]
unique_clusters = {level: adata.obs[level].unique() for level in clustering_levels}
# Initialize a dictionary to store consolidated marker genes results for each level
consolidated_marker_genes = {}
# Initialize a set to track processed fine-grained clusters
processed_clusters = set()
# Calculate Marker Genes Considering Cluster Depth
for level, clusters in unique_clusters.items():
    consolidated_marker_genes[level] = pd.DataFrame()
    for cluster in clusters:
        # Determine the cluster depth
        cluster_depth = cluster.count('_') + 1
        # Identify fine-grained clusters if cluster depth is more than 1
        if cluster_depth > 1:
            parent = "_".join(cluster.split("_")[: cluster_depth - 1])
            fine_grained_clusters = [c for c in clusters if c.startswith(parent)]
            adata_ = adata[[cl.startswith(parent) for cl in adata.obs[level]],:].copy()
        else:
            fine_grained_clusters = [cluster]
            adata_ = adata        
        # Check if fine-grained clusters have already been processed
        if any(fg_cluster in processed_clusters for fg_cluster in fine_grained_clusters):
            continue  # Skip the current loop if any of the fine-grained clusters have been processed
        # Mark fine-grained clusters as processed
        processed_clusters.update(fine_grained_clusters)
        # Calculate marker genes
        try:
            sc.tl.rank_genes_groups(adata_, groupby=level, groups=fine_grained_clusters, method='t-test', n_genes=100)
        except:
            continue        
        sc.tl.filter_rank_genes_groups(adata_, min_in_group_fraction=0.25, max_out_group_fraction=0.5, min_fold_change=1)
        # Extract and merge marker genes results
        for fg_cluster in fine_grained_clusters:
            result_df = pd.DataFrame({
                'genes': adata_.uns['rank_genes_groups']['names'][fg_cluster],
                'scores': adata_.uns['rank_genes_groups']['scores'][fg_cluster],
                'pvals': adata_.uns['rank_genes_groups']['pvals'][fg_cluster],
                'pvals_adj': adata_.uns['rank_genes_groups']['pvals_adj'][fg_cluster],
                'logfoldchanges': adata_.uns['rank_genes_groups']['logfoldchanges'][fg_cluster],
                'cluster': fg_cluster
            })
            consolidated_marker_genes[level] = pd.concat([consolidated_marker_genes[level], result_df], ignore_index=True)   
    # Create directory for each level of clustering
    output_dir = f'cluster_level_{level[-1]}'  # Extract level number from level string
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)   
    # Save the consolidated marker gene list for the level
    output_file = f'{output_dir}/consolidated_marker_genes_{level}_{marker_tpye}.csv'
    consolidated_marker_genes[level].to_csv(output_file, index=False)
    print(f'Saved consolidated marker genes for {level} to {output_file}.')
# Perform necessary cleanup operations
del adata
# Conduct a quality check of the saved files
# This step is more manual and involves checking the files directly, so no code is provided here.
print("Step 4 completed. Please manually check the saved files for correctness.")
