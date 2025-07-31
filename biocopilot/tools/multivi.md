# multiVI 
## MultiVI is used to integrate paired/unpaired multiome data, impute missing modality, normalizate other cell- and sample-level confounding factors.
## The import form of the multiVI package
import scvi
## setup data function of multiVI
scvi.model.MULTIVI.setup_anndata(adata, layer=None, batch_key=None, size_factor_key=None, categorical_covariate_keys=None, continuous_covariate_keys=None, protein_expression_obsm_key=None, protein_names_uns_key=None, **kwargs)
### the detail description of parameters:
1. adata (AnnData) – AnnData object. Rows represent cells, columns represent features.
2. layer (str | None (default: None)) – if not None, uses this as the key in adata.layers for raw count data.
3. batch_key (str | None (default: None)) – key in adata.obs for batch information. Categories will automatically be converted into integer categories and saved to adata.obs[‘_scvi_batch’]. If None, assigns the same batch to all the data.
4. size_factor_key (str | None (default: None)) – key in adata.obs for size factor information. Instead of using library size as a size factor, the provided size factor column will be used as offset in the mean of the likelihood. Assumed to be on linear scale.
5. categorical_covariate_keys (list[str] | None (default: None)) – keys in adata.obs that correspond to categorical data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space). Thus, these should not be used for biologically-relevant factors that you do _not_ want to correct for.
6. continuous_covariate_keys (list[str] | None (default: None)) – keys in adata.obs that correspond to continuous data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space). Thus, these should not be used for biologically-relevant factors that you do _not_ want to correct for.
7. protein_expression_obsm_key (str | None (default: None)) – key in adata.obsm for protein expression data.
8. protein_names_uns_key (str | None (default: None)) – key in adata.uns for protein names. If None, will use the column names of adata.obsm[protein_expression_obsm_key] if it is a DataFrame, else will assign sequential names to proteins.
## setup model function of multiVI
scvi.model.MULTIVI(adata, n_genes, n_regions, modality_weights='equal', modality_penalty='Jeffreys', n_hidden=None, n_latent=None, n_layers_encoder=2, n_layers_decoder=2, dropout_rate=0.1, region_factors=True, gene_likelihood='zinb', dispersion='gene', use_batch_norm='none', use_layer_norm='both', latent_distribution='normal', deeply_inject_covariates=False, encode_covariates=False, fully_paired=False, protein_dispersion='protein', **model_kwargs)
### the detail description of parameters:
1. adata (AnnData) – AnnData object that has been registered via setup_anndata().
2. n_genes (int) – The number of gene expression features (genes).
3. n_regions (int) – The number of accessibility features (genomic regions).
4. modality_weights (Literal['equal', 'cell', 'universal'] (default: 'equal')) – Weighting scheme across modalities. 
    One of the following:
        "equal": Equal weight in each modality
        "universal": Learn weights across modalities w_m.
        "cell": Learn weights across modalities and cells. w_{m,c}
5. modality_penalty (Literal['Jeffreys', 'MMD', 'None'] (default: 'Jeffreys')) – Training Penalty across modalities. 
    One of the following:
        "Jeffreys": Jeffreys penalty to align modalities
        "MMD": MMD penalty to align modalities
        "None": No penalty
6. n_hidden (int | None (default: None)) – Number of nodes per hidden layer. If None, defaults to square root of number of regions.
7. n_latent (int | None (default: None)) – Dimensionality of the latent space. If None, defaults to square root of n_hidden.
8. n_layers_encoder (int (default: 2)) – Number of hidden layers used for encoder NNs.
9. n_layers_decoder (int (default: 2)) – Number of hidden layers used for decoder NNs.
10. dropout_rate (float (default: 0.1)) – Dropout rate for neural networks.
11. model_depth – Model sequencing depth / library size.
12. region_factors (bool (default: True)) – Include region-specific factors in the model.
13. gene_dispersion – 
    One of the following
        'gene' - genes_dispersion parameter of NB is constant per gene across cells
        'gene-batch' - genes_dispersion can differ between different batches 
        'gene-label' - genes_dispersion can differ between different labels
14. protein_dispersion (Literal['protein', 'protein-batch', 'protein-label'] (default: 'protein')) – One of the following * 'protein' - protein_dispersion parameter is constant per protein across cells * 'protein-batch' - protein_dispersion can differ between different batches NOT TESTED * 'protein-label' - protein_dispersion can differ between different labels NOT TESTED
15. latent_distribution (Literal['normal', 'ln'] (default: 'normal')) – One of * 'normal' - Normal distribution * 'ln' - Logistic normal distribution (Normal(0, I) transformed by softmax)
16. deeply_inject_covariates (bool (default: False)) – Whether to deeply inject covariates into all layers of the decoder. If False, covariates will only be included in the input layer.
17. fully_paired (bool (default: False)) – allows the simplification of the model if the data is fully paired. Currently ignored.
18. **model_kwargs – Keyword args for MULTIVAE
## how to train the multiVI model
scvi.model.MULTIVI.train(max_epochs=500, lr=0.0001, accelerator='auto', devices='auto', train_size=0.9, validation_size=None, shuffle_set_split=True, batch_size=128, weight_decay=0.001, eps=1e-08, early_stopping=True, save_best=True, check_val_every_n_epoch=None, n_steps_kl_warmup=None, n_epochs_kl_warmup=50, adversarial_mixing=True, datasplitter_kwargs=None, plan_kwargs=None, **kwargs)
### the detail description of parameters:
1. max_epochs (int (default: 500)) – Number of passes through the dataset.
lr (float (default: 0.0001)) – Learning rate for optimization.
2. accelerator (str (default: 'auto')) – Supports passing different accelerator types (“cpu”, “gpu”, “tpu”, “ipu”, “hpu”, “mps, “auto”) as well as custom accelerator instances.
3. devices (int | list[int] | str (default: 'auto')) – The devices to use. Can be set to a non-negative index (int or str), a sequence of device indices (list or comma-separated str), the value -1 to indicate all available devices, or “auto” for automatic selection based on the chosen accelerator. If set to “auto” and accelerator is not determined to be “cpu”, then devices will be set to the first available device.
4. train_size (float (default: 0.9)) – Size of training set in the range [0.0, 1.0].
5. validation_size (float | None (default: None)) – Size of the test set. If None, defaults to 1 - train_size. If train_size + validation_size < 1, the remaining cells belong to a test set.
6. shuffle_set_split (bool (default: True)) – Whether to shuffle indices before splitting. If False, the val, train, and test set are split in the sequential order of the data according to validation_size and train_size percentages.
7. batch_size (int (default: 128)) – Minibatch size to use during training.
weight_decay (float (default: 0.001)) – weight decay regularization term for optimization
8. eps (float (default: 1e-08)) – Optimizer eps
9. early_stopping (bool (default: True)) – Whether to perform early stopping with respect to the validation set.
10. save_best (bool (default: True)) – DEPRECATED Save the best model state with respect to the validation loss, or use the final state in the training procedure.
11. check_val_every_n_epoch (int | None (default: None)) – Check val every n train epochs. By default, val is not checked, unless early_stopping is True. If so, val is checked every epoch.
12. n_steps_kl_warmup (int | None (default: None)) – Number of training steps (minibatches) to scale weight on KL divergences from 0 to 1. Only activated when n_epochs_kl_warmup is set to None. If None, defaults to floor(0.75 * adata.n_obs).
13. n_epochs_kl_warmup (int | None (default: 50)) – Number of epochs to scale weight on KL divergences from 0 to 1. Overrides n_steps_kl_warmup when both are not None.
14. adversarial_mixing (bool (default: True)) – Whether to use adversarial training to penalize the model for umbalanced mixing of modalities.
datasplitter_kwargs (dict | None (default: None)) – Additional keyword arguments passed into DataSplitter.
15. plan_kwargs (dict | None (default: None)) – Keyword args for TrainingPlan. Keyword arguments passed to train() will overwrite values present in plan_kwargs, when appropriate.
16. **kwargs – Other keyword args for Trainer.
## how to get results from trained model
model.get_latent_representation()
model.get_normalized_expression()