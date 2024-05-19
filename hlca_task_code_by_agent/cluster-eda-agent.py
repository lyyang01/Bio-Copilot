import scanpy as sc

# Load the dataset
dataset = sc.read_h5ad('../data/integrated_dataset_copy.h5ad')

# Verify the X_scANVI attribute
if 'X_scANVI' in dataset.obsm.keys():
    print("X_scANVI attribute is present and correctly loaded.")
else:
    raise ValueError("X_scANVI attribute is missing from the dataset. Please check the dataset integrity.")

import scanpy as sc
from leidenalg import find_partition

# Assuming the dataset is already loaded and X_scANVI attribute verified
# dataset variable is 'dataset'

# Perform Leiden Clustering
# Access the X_scANVI representation
X_scANVI = dataset.obsm['X_scANVI']

# Convert the X_scANVI representation into a nearest neighbor graph
sc.pp.neighbors(dataset, use_rep='X_scANVI', n_neighbors=30)

# Apply the Leiden algorithm to the nearest neighbor graph to identify clusters
sc.tl.leiden(dataset, key_added='cluster_level_1', resolution=0.01)

# Save Cluster Labels
# The labels are automatically saved in dataset.obs under the specified key ('cluster_level_1' in this case)
print(dataset.obs['cluster_level_1'].head())

# Documentation and Cleanup
# Documenting parameters used for Leiden clustering: resolution=0.01, n_neighbors=30
# Cleanup not explicitly required in this code snippet as it's managed by Python's garbage collection

import pandas as pd
import scanpy as sc
import os

# Check if the final h5ad file exists before proceeding with the clustering
output_file = '../data/clustered_final.h5ad'
if os.path.exists(output_file):
    print("Final clustered dataset already exists. Skipping clustering process.")
else:
    # Initialization for Subsequent Levels of Clustering
    # Assuming the initial clustering has been completed and is accessible in 'dataset'

    # Iterative Clustering for Levels 2 to 4
    for level in range(2, 5):  # Levels 2 to 4
        previous_level = f'cluster_level_{level-1}'
        current_level = f'cluster_level_{level}'
        
        # Get unique clusters from the previous level
        unique_clusters = dataset.obs[previous_level].unique()
        
        for cluster in unique_clusters:
            # Data Subsetting
            subset = dataset[dataset.obs[previous_level] == cluster]
            
            # Performing Leiden Clustering with adjusted parameters
            sc.pp.neighbors(subset, use_rep='X_scANVI', n_neighbors=30)  # Set n_neighbors to 30
            sc.tl.leiden(subset, resolution=0.1, key_added=current_level)  # Set resolution to 0.1 at all levels
            
            # Generate new sub-clustering labels for each cluster
            new_sub_cluster_labels = f"{previous_level}_{cluster}_" + subset.obs[current_level].astype(str)
            
            # Update dataset with new sub-cluster labels
            dataset.obs.loc[subset.obs.index, current_level] = new_sub_cluster_labels

        # Save cluster results to a CSV file for each cluster level
        csv_output_file = f'../data/cluster_results_{current_level}.csv'
        if not os.path.exists(csv_output_file):
            dataset.obs[[current_level]].to_csv(csv_output_file)
            print(f"Cluster results for {current_level} saved to CSV.")
        else:
            print(f"CSV file for {current_level} already exists. Skipping save.")

    # Finalization and Saving the Clustered Dataset
    dataset.write_h5ad(output_file)
    print("Clustering process completed and dataset saved.")


import scanpy as sc
import pandas as pd
import os

# Load the final dataset
dataset = sc.read_h5ad('../data/clustered_final.h5ad')

# Ensure the output directory exists
output_dir = '../data/marker_genes/'
os.makedirs(output_dir, exist_ok=True)

# Iterate through each clustering level (1 to 4)
for level in range(1, 5):
    current_level = f'cluster_level_{level}'
    unique_clusters = dataset.obs[current_level].unique()
    
    # Initialize a DataFrame to store results for the current level
    level_results = pd.DataFrame()
    
    # Iterate through each cluster within the current level
    for cluster in unique_clusters:
        #TODO
        print(cluster)
        # Determine the cluster's depth and its parent cluster, if applicable
        cluster_depth = cluster.count('-') + 1
        parent_cluster = '-'.join(cluster.split('-')[:-1]) if cluster_depth > 1 else None
        
        # Differential Expression Analysis with Depth Consideration
        if cluster_depth > 1:
            # Identify sister clusters
            sister_clusters = [c for c in unique_clusters if c.startswith(parent_cluster) and c != cluster]
            groups = [cluster] + sister_clusters
            
            # Subset dataset for the current and sister clusters
            subset = dataset[dataset.obs[current_level].isin(groups)].copy()
            
            # Perform differential expression analysis considering 'all' view
            try:
              sc.tl.rank_genes_groups(subset, groupby=current_level, groups=[cluster], reference='rest', n_genes=100)
            # Filter marker genes
            except:
              continue
            sc.tl.filter_rank_genes_groups(subset, min_in_group_fraction=0.25, min_fold_change=2)
            # Save marker genes for 'all' view
            result_all = sc.get.rank_genes_groups_df(subset, group=cluster)
            result_all['view'] = 'all'
            result_all['cluster'] = cluster
            level_results = pd.concat([level_results, result_all])
            
            # Perform differential expression analysis considering 'sister' view
            try:
              sc.tl.rank_genes_groups(subset, groupby=current_level, groups=groups, reference='rest', n_genes=100)
            except:
              continue
            # Filter marker genes
            sc.tl.filter_rank_genes_groups(subset, min_in_group_fraction=0.25, min_fold_change=2)
            # Save marker genes for 'sister' view
            result_sister = sc.get.rank_genes_groups_df(subset, group=cluster)
            result_sister['view'] = 'sister'
            result_sister['cluster'] = cluster
            level_results = pd.concat([level_results, result_sister])
        else:
            # Perform differential expression analysis without depth consideration
            sc.tl.rank_genes_groups(dataset, groupby=current_level, groups=[cluster], reference='rest', n_genes=100)
            # Filter marker genes
            sc.tl.filter_rank_genes_groups(dataset, min_in_group_fraction=0.25, min_fold_change=2)
            # Save marker genes
            result = sc.get.rank_genes_groups_df(dataset, group=cluster)
            result['view'] = 'all'
            result['cluster'] = cluster
            level_results = pd.concat([level_results, result])
    
    # Save the results for the current level to a CSV file
    level_results.to_csv(f'{output_dir}{current_level}_marker_genes.csv', index=False)

print("Marker gene analysis completed for all clusters and levels, with results saved for each level.")



