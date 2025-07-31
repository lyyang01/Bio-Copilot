# Cell type annotation from marker genes
In single-cell, we have no prior information of which cell type each cell belongs. To assign cell type labels, we first project all cells in a shared embedded space, then we find communities of cells that show a similar transcription profile and finally we check what cell type specific markers are expressed. If more than one marker gene is available, statistical methods can be used to test if a set of markers is enriched in a given cell population.


### Loading packages
```python
#First, we need to load the relevant packages, scanpy to handle scRNA-seq data and decoupler to use statistical methods.
import scanpy as sc
import decoupler as dc
import numpy as np
```

### Data preprocessing and Quality control
```python
### Basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

### Annotate the group of mitochondrial genes as 'mt'
adata.var['mt'] = adata.var_names.str.startswith('MT-')  
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

### Filter cells following standard QC criteria.
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

### Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['log_norm'] = adata.X.copy()
Then we group cells based on the similarity of their transcription profiles. To visualize the communities we perform UMAP reduction.

### Identify the highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

### Regress and scale the data
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
```

### Clustering and Visualization
```python
### Generate PCA features
sc.tl.pca(adata, svd_solver='arpack')

### Restore X to be norm counts
dc.swap_layer(adata, 'log_norm', X_layer_key=None, inplace=True)

### Compute distances in the PCA space, and find cell neighbors
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

### Generate UMAP features
sc.tl.umap(adata)

### Run leiden clustering algorithm
sc.tl.leiden(adata)

# Visualize
sc.pl.umap(adata, color='leiden', title='RNA UMAP', 
           frameon=False, legend_fontweight='normal', legend_fontsize=15)
```

### Calculate Marker genes
To annotate single cell clusters, we can use cell type specific marker genes. These are genes that are mainly expressed exclusively by a specific cell type, making them useful to distinguish heterogeneous groups of cells. Marker genes were discovered and annotated in previous studies and there are some resources that collect and curate them.

Omnipath is one of the largest available databases of curated prior knowledge. Among its resources, there is PanglaoDB, a database of cell type markers, which can be easily accessed using a wrapper to Omnipath from decoupler.

```python
# Query Omnipath and get PanglaoDB
markers = dc.get_resource('PanglaoDB')

# Filter by canonical_marker and human
markers = markers[markers['human'] & markers['canonical_marker'] & (markers['human_sensitivity'] > 0.5)]

# Remove duplicated entries
markers = markers[~markers.duplicated(['cell_type', 'genesymbol'])]

#Enrichment with Over Representation Analysis (ORA)
#To infer functional enrichment scores we will run the Over Representation Analysis (ora) method. As input data it accepts an expression matrix (decoupler.run_ora) or the results of differential expression analysis (decoupler.get_ora_df). For the former, by default the top 5% of expressed genes by sample are selected as the set of interest (S), and for the latter a user-defined significance filtering can be used. Once we have S, it builds a contingency table using set operations for each set stored in the gene set resource being used (net). Using the contingency table, ora performs a one-sided Fisher exact test to test for significance of overlap between sets. The final score is obtained by log-transforming the obtained p-values, meaning that higher values are more significant.

#We can run ora with a simple one-liner:

dc.run_ora(
    mat=adata,
    net=markers,
    source='cell_type',
    target='genesymbol',
    min_n=3,
    verbose=True,
    use_raw=False
)

#The obtained scores (-log10(p-value))(ora_estimate) and p-values (ora_pvals) are stored in the .obsm key:

adata.obsm['ora_estimate']
acts = dc.get_acts(adata, obsm_key='ora_estimate')

# We need to remove inf and set them to the maximum value observed for pvals=0
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e

dc.get_acts returns a new AnnData object which holds the obtained activities in its .X attribute, allowing us to re-use many scanpy functions, for example:

# Visualization
sc.pl.umap(acts, color=['NK cells', 'leiden'], cmap='RdBu_r')
sc.pl.violin(acts, keys=['NK cells'], groupby='leiden')
```

### Cell Annotation
```python
#With decoupler we can also identify which are the top predicted cell types per cluster using the function dc.rank_sources_groups. Here, it identifies "marker" cell types per cluster using same statistical tests available in scanpy's scanpy.tl.rank_genes_groups.

df = dc.rank_sources_groups(acts, groupby='leiden', reference='rest', method='t-test_overestim_var')

clusters = sorted(set(df.group))
df_new = pd.DataFrame()
for cluster in clusters:
    sub_df = df[df.group==cluster]
    stat_norm = (sub_df['statistic'] - sub_df['statistic'].min()) / (sub_df['statistic'].max() - sub_df['statistic'].min())
    meanchange_norm = (sub_df['meanchange'] - sub_df['meanchange'].min()) / (sub_df['meanchange'].max() - sub_df['meanchange'].min())
    sub_df['stat_norm'] = stat_norm
    sub_df['meanchange_norm'] = meanchange_norm
    sub_df['score'] = stat_norm + meanchange_norm
    df_new = pd.concat([df_new, sub_df.sort_values(by='score', ascending=False)], axis=0)

#We can then extract the top 3 predicted cell types per cluster:
n_ctypes = 3
ctypes_dict = df.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict()

#We can visualize the obtained top predicted cell types:
sc.pl.matrixplot(acts, ctypes_dict, 'leiden', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='RdBu_r')

#We can check individual cell types by plotting their distributions:
sc.pl.violin(acts, keys=['T cells', 'B cells', 'Platelets', 'Monocytes', 'NK cells'], groupby='leiden')

#The final annotation should be done manually based on the assessment of the enrichment results. However, an automatic prediction can be made by assigning the top predicted cell type per cluster. This approach does not require expertise in the tissue being studied but can be prone to errors. Nonetheless it can be useful to generate a first draft, let's try it:
annotation_dict = df.groupby('group').head(1).set_index('group')['names'].to_dict()

# Once we have selected the top cell type we can finally annotate:
# Add cell type column based on annotation
# the groups with statistic value less than 30 are assigned Unknown
annotation_list = []
for cluster in adata.obs['leiden']:
    if annotation_condidate.loc[cluster, 'statistic'] > 30:
        annotation_list.append(annotation_condidate.loc[cluster, 'names'])
    else:
        annotation_list.append('Unknown (%s)'%cluster)

adata.obs['cell_type'] = annotation_list

#adata.obs['cell_type'] = [annotation_dict[clust] for clust in adata.obs['leiden']]
```

