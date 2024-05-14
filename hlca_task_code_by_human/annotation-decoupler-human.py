
import scanpy as sc
import decoupler as dc
import numpy as np
import json

path_HLCA = "../data/final_adata_with_subclusters.h5ad"
adata = sc.read(path_HLCA)

adata.var['esmid'] = adata.var_names
adata.var_names = adata.var['gene_symbols']

# get human DB of cellmarker2.0
import pandas as pd
cellmarker_file = '../data/Cell_marker_Human.xlsx'
markers = pd.read_excel(cellmarker_file)

# Filter by canonical_marker and human
lung_related_class = ['Airway', 'Airway epithelium', 'Bronchi', 'Bronchus', 'Epithelium', 'Lung', 'Respiratory tract', ]
lung_related_tissue = ['Acinus', 'Airway', 'Airway epithelium', 'Alveolus', 'Bronchial vessel', 'Bronchiole', 'Bronchoalveolar lavage', 'Bronchus',
                     'Bronchoalveolar system', 'Distal airway', 'Epithelium', 'Lung', 'Nose', 'Respiratory tract', 'Trachea', 'Tracheal airway epithelium',
                    ]
markers = markers[(markers['tissue_type'].isin(lung_related_tissue) | markers['tissue_class'].isin(lung_related_class)) & (markers['cancer_type'] == 'Normal')]
# Remove duplicated entries
markers = markers[~markers.duplicated(['cell_name', 'Symbol'])]

adata_hvg = adata[:,adata.var.highly_variable].copy()

dc.run_ora(
    mat=adata_hvg,
    net=markers,
    source='cell_name', #cell_type
    target='Symbol', #genesymbol
    min_n=3,
    verbose=True,
    use_raw=False
)

acts = dc.get_acts(adata_hvg, obsm_key='ora_estimate')

# We need to remove inf and set them to the maximum value observed for pvals=0
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e

df = dc.rank_sources_groups(acts, groupby='cluster_level_2', reference='rest', method='t-test_overestim_var')
ctypes_dict_lev2 = df.groupby('group').head(1).groupby('group')['names'].apply(lambda x: list(x)).to_dict()

df = dc.rank_sources_groups(acts, groupby='cluster_level_2', reference='rest', method='t-test_overestim_var')
ctypes_dict_lev3 = df.groupby('group').head(1).groupby('group')['names'].apply(lambda x: list(x)).to_dict()

df = dc.rank_sources_groups(acts, groupby='cluster_level_2', reference='rest', method='t-test_overestim_var')
ctypes_dict_lev4 = df.groupby('group').head(1).groupby('group')['names'].apply(lambda x: list(x)).to_dict()

json.dump({**ctypes_dict_lev2,**ctypes_dict_lev3, **ctypes_dict_lev4},open("../result/annotation_data_decoupler.json","w"))