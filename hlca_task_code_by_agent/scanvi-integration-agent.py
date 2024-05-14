# Adjusted Step 2: Data Integration using scANVI with scArches on GPU with specified and additional essential parameters

# Import Required Libraries
import scanpy as sc
import scarches as sca
import os

# Set CUDA Device
os.environ['CUDA_VISIBLE_DEVICES'] = '4'

# Define early stopping parameters for unsupervised and semisupervised training
early_stopping_kwargs = {
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "patience": 10,
    "threshold": 0,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}

early_stopping_kwargs_scanvi = {
    "early_stopping_metric": "accuracy",
    "save_best_state_metric": "accuracy",
    "on": "full_dataset",
    "patience": 10,
    "threshold": 0.001,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}

# Check if the integrated dataset file already exists
output_file_path = '../data/integrated_dataset_scANVI.h5ad'
if not os.path.exists(output_file_path):
    # Load the dataset
    adata = sc.read_h5ad('../data/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_scanvi_label.h5ad')
    
    # Subset data to highly variable genes
    adata = adata[:, adata.var.highly_variable]
    
    # Set `adata.X` to counts and update `raw.X` with normalized counts for compatibility
    adata.raw = adata
    adata.X = adata.layers['counts']
    
    # scANVI Model Setup with scArches
    sca.dataset.setup_anndata(adata, labels_key='scanvi_label', batch_key='dataset')
    
    scanvi_model = sca.models.SCANVI(
        adata, 
        unlabeled_category='unlabeled', 
        use_cuda=True, 
        n_layers=2,
        n_latent=30,
        gene_likelihood="zinb"  # Use zero-inflated negative binomial distribution for gene expression likelihood
    )
    
    # Model Training with adjusted parameters
    scanvi_model.train(
        n_epochs_unsupervised=500,
        n_epochs_semisupervised=200,
        unsupervised_trainer_kwargs=dict(early_stopping_kwargs=early_stopping_kwargs),
        semisupervised_trainer_kwargs=dict(metrics_to_monitor=["elbo", "accuracy"], early_stopping_kwargs=early_stopping_kwargs_scanvi),
        frequency=1
    )
    
    # Data Integration and Correction
    adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()
    
    # Save Integrated Data before visualization
    adata.write_h5ad(output_file_path)
    
    # Visualization can be performed after saving the integrated data
else:
    print(f"Integrated dataset file already exists at {output_file_path}. Skipping data integration step.")

##
adata_ori = sc.read_h5ad('../data/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_scanvi_label.h5ad')
adata = sc.read_h5ad(output_file_path)
adata_ori.obsm['X_scANVI'] = adata.obsm['X_scANVI']
adata_ori.write('../data/integrated_dataset_copy.h5ad')
