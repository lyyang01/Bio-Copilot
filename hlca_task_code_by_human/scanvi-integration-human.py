import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

import scanpy as sc
import anndata
import scarches as sca

import torch
torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)
os.environ['CUDA_VISIBLE_DEVICES'] = '1'

path_HLCA_unintegrated = "../../data/output/HLCA_v1_intermediates/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_log1p.h5ad"
path_HLCA_unintegrated_prepped = "../../data/output/HLCA_v1_intermediates/HLCA_v1_scANVI_input.h5ad"
dir_out = "../../data/scANVI_integration"

adata = sc.read(path_HLCA_unintegrated)
adata_ori = adata.copy()
adata = adata[:,adata.var.highly_variable].copy()

adata.X = adata.layers['counts']
adata.raw = adata
raw = adata.raw.to_adata()
raw.X = adata.layers['counts']
adata.raw = raw

adata.obs.dataset.unique().tolist()

condition_key = 'dataset'
cell_type_key = 'scanvi_label'
unlabeled_category = "unlabeled"

vae_epochs = 500
scanvi_epochs = 200

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

sca.dataset.setup_anndata(adata, batch_key=condition_key, labels_key=cell_type_key)

vae = sca.models.SCANVI(
    adata,
    unlabeled_category,
    n_layers=2,
    n_latent = 30, # to allow for capturing more heterogeneity
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
    gene_likelihood="nb", # because we have UMI data
    use_cuda=True #to use GPU
)

print("Labelled Indices: ", len(vae._labeled_indices))
print("Unlabelled Indices: ", len(vae._unlabeled_indices))


vae.train(
    n_epochs_unsupervised=vae_epochs,
    n_epochs_semisupervised=scanvi_epochs,
    unsupervised_trainer_kwargs=dict(early_stopping_kwargs=early_stopping_kwargs),
    semisupervised_trainer_kwargs=dict(metrics_to_monitor=["elbo", "accuracy"],
                                       early_stopping_kwargs=early_stopping_kwargs_scanvi),
    frequency=1
)

model_dir = os.path.join(dir_out, "scanvi_model") # this is the directory name/path of the directory *to be created*
vae.save(model_dir, overwrite=False)
reference_latent = sc.AnnData(vae.get_latent_representation())
reference_latent.obs.index = adata.obs.index

reference_latent.write(os.path.join(dir_out, "scANVI_embedding.h5ad"))
adata_ori.obsm['X_scanvi_emb'] = vae.get_latent_representation()
adata_ori.write("../../data/output/HLCA_v1_intermediates/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_log1p_scanvi_embedding.h5ad")