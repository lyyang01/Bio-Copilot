
# The core code for cell types annotation using decoupler.

### Step 2: Load and Preprocess the Data
import os
import scanpy as sc
import pandas as pd

# 2.1 Load the `.h5ad` File
# 2.1.1 Import the necessary library by adding `import scanpy as sc` at the beginning of your script.
# This step is already completed above.

# 2.1.2 Use the `scanpy.read_h5ad('/home/zhoulu/hlca_test_data/final_adata_with_subclusters.h5ad')` function to load the `.h5ad` file.
h5ad_path = '../data/final_adata_with_subclusters.h5ad'
if os.path.exists(h5ad_path):
    adata = sc.read_h5ad(h5ad_path)
else:
    print(f"File {h5ad_path} does not exist.")

# 2.1.3 Verify the loaded data by printing the `adata` object's summary.
print(adata)

# 2.2 Load the Cell Marker Database
# 2.2.1 Import the pandas library by adding `import pandas as pd` at the beginning of your script.
# This step is already completed above.

# 2.2.2 Load the cell marker manual annotation database using `pandas.read_excel('/home/zhoulu/hlca_test_data/Cell_marker_Human.xlsx')`.
cell_markers_path = '../data/Cell_marker_Human.xlsx'
if os.path.exists(cell_markers_path):
    cell_markers = pd.read_excel(cell_markers_path)
else:
    print(f"File {cell_markers_path} does not exist.")

# 2.2.3 Inspect the first few rows of the DataFrame to ensure it has loaded correctly and to familiarize yourself with the structure of the data.
print(cell_markers.head())


### Step 3: Prepare the marker gene list for Decoupler
import pandas as pd
import scanpy as sc

# Step 3: Refined Plan for Data Preparation and Filtering

# Step 1: Update `adata.var_names` to Use 'gene_symbols'
adata.var_names = adata.var['gene_symbols']

# Step 2: Retain Highly Variable Genes for Cell Annotation
adata_filtered = adata[:, adata.var.highly_variable].copy()

# Step 3: Filter the Cell Marker DataFrame Based on Specific Criteria
# Define the lists of acceptable values for 'tissue_type' and 'tissue_class'.
tissue_types = ['Airway', 'Airway epithelium', 'Bronchi', 'Bronchus', 'Epithelium', 'Lung', 'Respiratory tract']
tissue_classes = ['Acinus', 'Airway', 'Airway epithelium', 'Alveolus', 'Bronchial vessel', 'Bronchiole', 'Bronchoalveolar lavage', 'Bronchus', 'Bronchoalveolar system', 'Distal airway', 'Epithelium', 'Lung', 'Nose', 'Respiratory tract', 'Trachea', 'Tracheal airway epithelium']

# Filter the DataFrame to include only rows where 'tissue_type' is in the `tissue_types` list OR 'tissue_class' is in the `tissue_classes` list, AND 'cancer_type' is 'Normal'.
filtered_df = cell_markers[
    (cell_markers['tissue_type'].isin(tissue_types) | cell_markers['tissue_class'].isin(tissue_classes)) &
    (cell_markers['cancer_type'] == 'Normal')
]

# Remove duplicates based on 'cell_name' and 'Symbol' to ensure each marker gene is uniquely associated with a cell type.
filtered_df = filtered_df.drop_duplicates(subset=['cell_name', 'Symbol'])

# Validation and Final Adjustments

# Inspect the first few rows of `adata_filtered` to ensure it contains only highly variable genes.
print(adata_filtered.var.head())

# Perform a similar inspection on `filtered_df` to verify that the filtering criteria have been correctly applied and duplicates have been removed.
print(filtered_df.head())

# If any discrepancies or unexpected results are observed, revisit the filtering criteria and ensure they are correctly implemented.



### Step 4: Run Decoupler for cell type annotation
import decoupler as dc
import numpy as np
import json
import os

# 4.1 Preparing the Data
# Assuming `adata_filtered` and `filtered_df` are already prepared from previous steps.

# 4.2 Performing Over-Representation Analysis (ORA)
dc.run_ora(mat=adata_filtered,
           net=filtered_df,
           source='cell_name',
           target='Symbol',
           min_n=3,
           use_raw=False)

# 4.3 Identifying Top Predicted Cell Types
acts = dc.get_acts(adata_filtered, obsm_key='ora_estimate')

# 4.4 Handling Infinite or NaN Values
finite_max = np.nanmax(acts.X[np.isfinite(acts.X)])
acts.X = np.nan_to_num(acts.X, nan=finite_max, posinf=finite_max, neginf=finite_max)

# 4.5 Annotating Clusters at Multiple Levels
for level in ['cluster_level_2', 'cluster_level_3', 'cluster_level_4']:
    result_df = dc.rank_sources_groups(acts,
                                       reference='rest',
                                       method='t-test_overestim_var',
                                       groupby=level)
    # 4.6 Saving the Results
    top_cell_types = result_df.groupby('group').first()['names'].to_dict()
    output_path = f'top_cell_types_{level}.json'
    if not os.path.exists(output_path):
        with open(output_path, 'w') as json_file:
            json.dump(top_cell_types, json_file)


