# cell2location
## The import form of the cell2location package
import cell2location
## Plot the gene filter given a set of cutoffs and return resulting list of genes.
cell2location.utils.filtering.filter_genes(adata, cell_count_cutoff=15, cell_percentage_cutoff2=0.05, nonz_mean_cutoff=1.12)
### the detail description of parameters:
1. adata – anndata object with single cell / nucleus data.
2. cell_count_cutoff – All genes detected in less than cell_count_cutoff cells will be excluded.
3. cell_percentage_cutoff2 – All genes detected in at least this percentage of cells will be included.
4. nonz_mean_cutoff – genes detected in the number of cells between the above mentioned cutoffs are selected only when their average expression in non-zero cells is above this cutoff.
### returns
a list of selected var_names
## Plot spatial abundance of cell types (regulatory programmes) with colour gradient and interpolation (from Visium anndata).
cell2location.plt.plot_spatial(adata, color, img_key='hires', show_img=True, **kwargs)
### the detail description of parameters:
1. adata – adata object with spatial coordinates in adata.obsm[‘spatial’]
2. color – list of adata.obs column names to be plotted
3. kwargs – arguments to plot_spatial_general
## Selects the data for one slide from the spatial anndata object.
cell2location.utils.select_slide(adata, s, batch_key='sample')
### the detail description of parameters:
1. adata – Anndata object with multiple spatial experiments
2. s – name of selected experiment
3. batch_key – column in adata.obs listing experiment name for each location
## Cell2location model for spatial mapping
cell2location.models.Cell2location(adata: anndata._core.anndata.AnnData, cell_state_df: pandas.core.frame.DataFrame, model_class: Optional[pyro.nn.module.PyroModule, None] = None, detection_mean_per_sample: bool = False, detection_mean_correction: float = 1.0, **model_kwargs)
### the detail description of parameters:
1. adata – spatial AnnData object that has been registered via setup_anndata().
2. cell_state_df – pd.DataFrame with reference expression signatures for each gene (rows) in each cell type/population (columns).
3. use_gpu – Use the GPU?
4. **model_kwargs – Keyword args for LocationModelLinearDependentWMultiExperimentLocationBackgroundNormLevelGeneAlphaPyroModel
## Sets up the AnnData object for Cell2location model.
cell2location.models.Cell2location.setup_anndata(adata: anndata._core.anndata.AnnData, layer: Optional[str, None] = None, batch_key: Optional[str, None] = None, labels_key: Optional[str, None] = None, categorical_covariate_keys: Optional[List[str], None] = None, continuous_covariate_keys: Optional[List[str], None] = None, **kwargs)
### the detail description of parameters:
1. layer – if not None, uses this as the key in adata.layers for raw count data.
2. batch_key – key in adata.obs for batch information. Categories will automatically be converted into integer categories and saved to adata.obs[‘_scvi_batch’]. If None, assigns the same batch to all the data.
3. labels_key – key in adata.obs for label information. Categories will automatically be converted into integer categories and saved to adata.obs[‘_scvi_labels’]. If None, assigns the same label to all the data.
4. categorical_covariate_keys – keys in adata.obs that correspond to categorical data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space). Thus, these should not be used for biologically-relevant factors that you do _not_ want to correct for.
5. continuous_covariate_keys – keys in adata.obs that correspond to continuous data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space). Thus, these should not be used for biologically-relevant factors that you do _not_ want to correct for.
## Train Cell2location model
cell2location.models.Cell2location.train(max_epochs: int = 30000, batch_size: int = None, train_size: float = 1, lr: float = 0.002, num_particles: int = 1, scale_elbo: float = 1.0, **kwargs)
### the detail description of parameters:
1. max_epochs – Number of passes through the dataset. If None, defaults to np.min([round((20000 / n_cells) * 400), 400])
2. train_size – Size of training set in the range [0.0, 1.0]. Use all data points in training because we need to estimate cell abundance at all locations.
3. batch_size – Minibatch size to use during training. If None, no minibatching occurs and all data is copied to device (e.g., GPU).
4. lr – Optimiser learning rate (default optimiser is ClippedAdam). Specifying optimiser via plan_kwargs overrides this choice of lr.
5. kwargs – Other arguments to scvi.model.base.PyroSviTrainMixin().train() method
## Summarise posterior distribution and export results (cell abundance) to anndata object of cell2location model
cell2location.models.Cell2location.export_posterior(adata, sample_kwargs: Optional[dict, None] = None, export_slot: str = 'mod', add_to_obsm: list = ['means', 'stds', 'q05', 'q95'], use_quantiles: bool = False)
### the detail description of parameters:
1. adata – anndata object where results should be saved
2. sample_kwargs – arguments for self.sample_posterior (generating and summarising posterior samples), namely:
    num_samples - number of samples to use (Default = 1000). batch_size - data batch size (keep low enough to fit on GPU, default 2048). use_gpu - use gpu for generating samples?
3. export_slot – adata.uns slot where to export results
4. add_to_obsm – posterior distribution summary to export in adata.obsm ([‘means’, ‘stds’, ‘q05’, ‘q95’]).
5. use_quantiles – compute quantiles directly (True, more memory efficient) or use samples (False, default). If True, means and stds cannot be computed so are not exported and returned.
## RegressionModel which estimates per cluster average mRNA count account for batch effects
cell2location.models.RegressionModel(adata: anndata._core.anndata.AnnData, model_class=None, use_average_as_initial: bool = True, **model_kwargs)
### the detail description of parameters:
1. adata – single-cell AnnData object that has been registered via setup_anndata().
2. use_gpu – Use the GPU?
3. **model_kwargs – Keyword args for LocationModelLinearDependentWMultiExperimentModel
## Sets up the AnnData object for RegressionModel.
cell2location.models.RegressionModel.setup_anndata(adata: anndata._core.anndata.AnnData, layer: Optional[str, None] = None, batch_key: Optional[str, None] = None, labels_key: Optional[str, None] = None, categorical_covariate_keys: Optional[List[str], None] = None, continuous_covariate_keys: Optional[List[str], None] = None, **kwargs)
### the detail description of parameters:
1. layer – if not None, uses this as the key in adata.layers for raw count data.
2. batch_key – key in adata.obs for batch information. Categories will automatically be converted into integer categories and saved to adata.obs[‘_scvi_batch’]. If None, assigns the same batch to all the data.
3. labels_key – key in adata.obs for label information. Categories will automatically be converted into integer categories and saved to adata.obs[‘_scvi_labels’]. If None, assigns the same label to all the data.
4. categorical_covariate_keys – keys in adata.obs that correspond to categorical data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space). Thus, these should not be used for biologically-relevant factors that you do _not_ want to correct for.
5. continuous_covariate_keys – keys in adata.obs that correspond to continuous data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space). Thus, these should not be used for biologically-relevant factors that you do _not_ want to correct for.
## Train the RegressionModel with useful defaults
cell2location.models.RegressionModel.train(max_epochs: Optional[int, None] = None, batch_size: int = 2500, train_size: float = 1, lr: float = 0.002, **kwargs)
### the detail description of parameters:
1. max_epochs – Number of passes through the dataset. If None, defaults to np.min([round((20000 / n_cells) * 400), 400])
2. train_size – Size of training set in the range [0.0, 1.0].
3. batch_size – Minibatch size to use during training. If None, no minibatching occurs and all data is copied to device (e.g., GPU).
4. lr – Optimiser learning rate (default optimiser is ClippedAdam). Specifying optimiser via plan_kwargs overrides this choice of lr.
5. kwargs – Other arguments to scvi.model.base.PyroSviTrainMixin().train() method
## Summarise posterior distribution and export results (cell abundance) to anndata object of RegressionModel
cell2location.models.RegressionModel.export_posterior(adata, sample_kwargs: Optional[dict, None] = None, export_slot: str = 'mod', add_to_varm: list = ['means', 'stds', 'q05', 'q95'], scale_average_detection: bool = True, use_quantiles: bool = False)
### the detail description of parameters:
1. adata – anndata object where results should be saved
2. sample_kwargs – arguments for self.sample_posterior (generating and summarising posterior samples), namely:
    num_samples - number of samples to use (Default = 1000). batch_size - data batch size (keep low enough to fit on GPU, default 2048). use_gpu - use gpu for generating samples?
3. export_slot – adata.uns slot where to export results
4. add_to_varm – posterior distribution summary to export in adata.varm ([‘means’, ‘stds’, ‘q05’, ‘q95’]).
5. use_quantiles – compute quantiles directly (True, more memory efficient) or use samples (False, default). If True, means and stds cannot be computed so are not exported and returned.
## Show quality control plots
plot_QC(summary_name: str = 'means', use_n_obs: int = 1000, scale_average_detection: bool = True)
