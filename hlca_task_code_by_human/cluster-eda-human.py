import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import sys
import os

sys.path.append("../../scripts/")
import analysis
import warnings
warnings.filterwarnings("ignore")


def add_nested_clustering_blind(
    adata,
    cluster_label_previous,
    cluster_label_new,
    use_rep,
    cluster_alg="leiden",
    cluster_res=0.2,
    cluster_k=30,
    min_cluster_size=50,
    redo_pca=True,
    verbose=True,
):
    """Function that goes through one round of clustering of already existing
    clusters, based on the input cluster df. All clusters will be reclustered 
    individually. ("blind" because we don't take into account annotation
    purity of clusters.) 
    Args:
        adata - anndata object to be clustered
        cluster_label_previous - parent cluster label
        cluster_label_new - label for new clustering
        use_rep - name of .obsm object to be used for neighbor graph
        cluster_alg - <"leiden","phenograph">
        cluster_res - only applicable when using "leiden" as cluster_alg
        cluster_k - only applicable when using "phenograph" as cluster_alg.
        min_cluster_size - only applicable when using "phenograph" as cluster_alg
            Make sure that cluster_k < min_cluster_size
        redo_pca - boolean. whether to re-calculate PCA for subclusters
        verbose - boolean
    Returns adata with new clustering (under adata.obs[cluster_label_new].
    """
    # copy original clustering
    clusters_previous = adata.obs[cluster_label_previous].tolist()
    adata.obs[cluster_label_new] = clusters_previous
    if not redo_pca:
        print("Not re-doing pca before nested clustering iterations!")
    for cluster in sorted(set(clusters_previous)):
        if verbose:
            print("Cluster:", cluster)
        subadata = adata[adata.obs[cluster_label_previous] == cluster, :].copy()
        if subadata.shape[0] < min_cluster_size:
            if verbose:
                print("cluster size smaller than", min_cluster_size, "\n")
            continue
        if verbose:
            print("reclustering...\n")
        if redo_pca:
            if verbose:
                print("running pca...")
            sc.tl.pca(subadata)
        if cluster_alg == "leiden":
            if verbose:
                print("calculating 30 nearest neighbors")
                print("using rep:", use_rep)
            sc.pp.neighbors(subadata, n_neighbors=30, use_rep=use_rep)
            if verbose:
                print("clustering")
            sc.tl.leiden(subadata, resolution=cluster_res, key_added=cluster_label_new)
        elif cluster_alg == "phenograph":
            subadata.obs[cluster_label_new] = pd.Categorical(
                sce.tl.phenograph(subadata.obsm[use_rep], k=cluster_k)[0]
            )
        else:
            raise ValueError("Your cluster_alg argument is incorrect.")
        subadata.obs[cluster_label_new] = [
            "{}.{}".format(cluster, new_cluster)
            for new_cluster in subadata.obs[cluster_label_new]
        ]
        adata.obs.loc[subadata.obs.index, cluster_label_new] = subadata.obs[
            cluster_label_new
        ]
    # order categories "numerically" (so not 1, 10, 11 but 1, 2, 3... 10, 11):
    # convert all cluster names to strings, instead of a mix of strings and ints:
    adata.obs[cluster_label_new] = [
        str(clust) for clust in adata.obs[cluster_label_new]
    ]
    cluster_numbers = list(sorted(set(adata.obs[cluster_label_new])))
    prefix_cluster = [float(x.split(".")[0]) for x in cluster_numbers]
    cluster_numbers_ordered = [
        cluster_numbers[idx] for idx in np.argsort(prefix_cluster)
    ]
    adata.obs[cluster_label_new] = pd.Categorical(
        adata.obs[cluster_label_new], categories=cluster_numbers_ordered
    )

    return adata

def get_cluster_markers(adata, cluster_label, marker_ref, ngenes=100, verbose=True, filtered=False):
    """
    Calculates markers for every cluster, using either all other cells or 
    the parent cluster as a reference (i.e. for cluster 00.00.01, it 
    uses all clusters starting with 00.00 as reference. For cluster 
    00, it uses all cells as reference).
    sc.tl.rank_genes is used for marker gene calculation.

    Arguments:
    adata - AnnData object
    cluster_label - string
        label in adata.obs that contains nested-cluster names
    marker_ref - either "all" or "sisters". Which clusters to compare with.
    ngenes - number of marker genes to get per cluster
    
    Returns:
    cluster_markers - pd.DataFrame
        dataframe with, for each cluster, 100 highest scoring genes,
        plus matching logfc and adj pvalue
    """
    # input check:
    if marker_ref == "all":
        print("Doing one versus all differential expression analysis.")
    elif marker_ref == "sisters":
        print("Doing one versus sisters differential expression analysis.")
    else:
        raise ValueError("marker_ref argument should be set to either 'all' or 'sisters'.")
    # convert clusters to strings:
    adata.obs[cluster_label] = [str(cl) for cl in adata.obs[cluster_label]]
    # store cluster set
    clusters = sorted(set(adata.obs[cluster_label]))
    colnames_nested = [
        [clust + "_gene", clust + "_logfc", clust + "_pval_adj"] for clust in clusters
    ]
    colnames = [item for sublist in colnames_nested for item in sublist]
    cluster_markers = pd.DataFrame(index=range(100), columns=colnames)
    parents_tested = list()
    for clust in clusters:
        clust_depth = len(clust.split("."))
        if clust_depth == 1:
            parent = "all"
            if parent not in parents_tested:
                if verbose:
                    print("ranking genes for parent group", parent)
                parents_tested.append(parent)
                sc.tl.rank_genes_groups(adata, groupby=cluster_label, n_genes=ngenes)
                rank_genes = "rank_genes_groups"
                if filtered:
                    sc.tl.filter_rank_genes_groups(adata, min_in_group_fraction=0.25, max_out_group_fraction=0.5, min_fold_change=1)
                    rank_genes = "rank_genes_groups_filtered"
                # store results for all clusters from this parent
                # i.e. all clusters of depth 1
                for d1_cluster in [
                    clust for clust in clusters if len(clust.split(".")) == 1
                ]:
                    # create a subdf that will allow us to sort genes per cluster
                    submarker_df = pd.DataFrame(
                        index=range(ngenes),
                        columns=[
                            d1_cluster + "_gene",
                            d1_cluster + "_logfc",
                            d1_cluster + "_pval_adj",
                        ],
                    )
                    submarker_df[d1_cluster + "_gene"] = adata.uns[rank_genes]["names"][d1_cluster]
                    submarker_df[d1_cluster + "_logfc"] = adata.uns[rank_genes]["logfoldchanges"][d1_cluster]
                    submarker_df[d1_cluster + "_pval_adj"] = adata.uns[rank_genes]["pvals_adj"][d1_cluster]
                    # sort values:
                    submarker_df.sort_values(
                        by=[d1_cluster + "_pval_adj", d1_cluster + "_logfc"],
                        ascending=[True, False],
                        inplace=True,
                    )
                    submarker_df = submarker_df.reset_index().drop(columns="index")
                    # and add to big dataframe
                    cluster_markers.loc[
                        submarker_df.index, submarker_df.columns
                    ] = submarker_df.values
        else:
            parent = ".".join(clust.split(".")[: clust_depth - 1])
            if parent not in parents_tested:
                # depending on reference choice, use whole adata as reference
                # or only the parent cluster.
                if marker_ref == "all":
                    subadata = adata
                elif marker_ref == "sisters":
                    subadata = adata[[cl.startswith(parent) for cl in adata.obs[cluster_label]],:].copy()
                if verbose:
                    print("ranking genes for parent group", parent)
                parents_tested.append(parent)
                siblings = [c for c in clusters if c.startswith(parent)]
                if len(siblings) < 2 and marker_ref == "sisters":
                    print("Cluster {} has only one subcluster. Skipping DEA for this parent.".format(parent))
                else:
                    sc.tl.rank_genes_groups(subadata, groupby=cluster_label, groups=siblings, n_genes=ngenes)
                    rank_genes = "rank_genes_groups"
                    if filtered:
                        sc.tl.filter_rank_genes_groups(subadata, min_in_group_fraction=0.25, max_out_group_fraction=0.5, min_fold_change=1)
                        rank_genes = "rank_genes_groups_filtered"
                    for same_depth_sibling in [
                        sib for sib in siblings if len(clust.split(".")) == clust_depth
                    ]:
                        # create a subdf that will allow us to sort genes per cluster
                        submarker_df = pd.DataFrame(
                            index=range(ngenes),
                            columns=[
                                same_depth_sibling + "_gene",
                                same_depth_sibling + "_logfc",
                                same_depth_sibling + "_pval_adj",
                            ],
                        )
                        submarker_df[same_depth_sibling + "_gene"] = subadata.uns[rank_genes]["names"][same_depth_sibling]
                        submarker_df[same_depth_sibling + "_logfc"] = subadata.uns[rank_genes]["logfoldchanges"][same_depth_sibling]
                        submarker_df[same_depth_sibling + "_pval_adj"] = subadata.uns[rank_genes]["pvals_adj"][same_depth_sibling]
                        # sort values:
                        submarker_df.sort_values(
                            by=[
                                same_depth_sibling + "_pval_adj",
                                same_depth_sibling + "_logfc",
                            ],
                            ascending=[True, False],
                            inplace=True,
                        )
                        submarker_df = submarker_df.reset_index().drop(columns="index")
                        # add to big dataframe
                        cluster_markers.loc[
                            submarker_df.index, submarker_df.columns
                        ] = submarker_df.values
    return cluster_markers

dir_benchmarking_res = (
    "../data/results/integration_benchmarking/benchmarking_results/integration/"
)
dir_benchmarking_cluster_output = "../data/results/integration_benchmarking/clustering/"
path_HLCA = "../data/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_log1p_scanvi_embedding.h5ad"
dir_HLCA_cluster_output = "../data/results/DEAs/leiden_v3/"

dataset_name = "HLCA"
number_of_clust_levels = 4
print("Dataset name:", dataset_name)
if dataset_name == "scanvi":
    # load dataset
    adata = sc.read(os.path.join(dir_benchmarking_res, "unscaled/hvg/scanvi.h5ad"))
    # specify which obsm to use for calculating neighbor graph
    use_rep = "X_emb"
    # specify whether or not to re-calculate the PCA for subsets of the object
    redo_pca = False
elif dataset_name == "seuratrpca":
    adata = sc.read(
        os.path.join(dir_benchmarking_res, "unscaled/hvg/R/seuratrpca.h5ad")
    )
    sc.tl.pca(adata)
    use_rep = "X_pca"
    redo_pca = True
elif dataset_name == "harmony":
    adata = sc.read(os.path.join(dir_benchmarking_res, "scaled/hvg/R/harmony.h5ad"))
    adata.obsm["X_pca"] = adata.obsm["X_emb"]
    use_rep = "X_emb"
    redo_pca = False
elif dataset_name == "HLCA":
    adata = sc.read(path_HLCA)
    use_rep = "X_scanvi_emb"
    redo_pca = False

if "X_umap" in adata.obsm.keys():
    if "scanvi_labels" in adata.obs.columns:
        sc.pl.umap(adata, color="scanvi_labels")
    elif "scgen_labels" in adata.obs.columns:
        sc.pl.umap(adata, color="scgen_labels")

for clustering_level in range(1, number_of_clust_levels + 1):
    print("clustering level:", clustering_level, "...")
    if clustering_level == 1:
        # skip for re-run
        cluster_name = "leiden_1"
        # first clustering is not nested, so use normal function:
        sc.pp.neighbors(adata, n_neighbors=30, use_rep=use_rep)
        sc.tl.leiden(adata, resolution=0.01, key_added=cluster_name)
    else:
        previous_clustering = "leiden_" + str(clustering_level - 1)
        cluster_name = "leiden_" + str(clustering_level)
        #         perform nested clustering
        #         set parameters:
        res = 0.1
        if clustering_level == 2:
            k = 30
            min_cluster_size = 50
        elif clustering_level == 3:
            k = 15
            min_cluster_size = 30
        elif clustering_level == 4:
            k = 10
            min_cluster_size = 10

        adata = add_nested_clustering_blind(
            adata,
            previous_clustering,
            cluster_name,
            use_rep=use_rep,
            cluster_alg="leiden",
            cluster_k=k,
            cluster_res=res,
            min_cluster_size=min_cluster_size,
            redo_pca=redo_pca,  # SET THIS TO FALSE FOR SCANVI!!! OR OTHER EMBEDDING-OUTPUT METHODS!!!!!
        )
    # plot
    if "X_umap" in adata.obsm.keys():
        sc.pl.umap(adata, color=cluster_name)
    # store clustering:
    cluster_df = pd.DataFrame(adata.obs[cluster_name], index=adata.obs.index)
    # write to csv for benchmarking data:
    if dataset_name in ["harmony","scanvi","seuratrpca"]:
        cluster_df.to_csv(
            os.path.join(dir_benchmarking_cluster_output, f"{dataset_name}/{dataset_name}_{cluster_name}_cluster_assignment.csv")
        )
    # If wanted/needed, for final HLCA:
    if dataset_name == "HLCA":
        # store cluster assignments:
        cluster_df.to_csv(
            os.path.join(dir_HLCA_cluster_output, f"LCA_{cluster_name}_cluster_assignment.csv")
        )
        # calculate marker genes with respect to all other clusters, and with respect to sister clusters (i.e. other cluster from the same parent cluster):
        for marker_ref in ["sisters", "all"]:
            marker_gene_df = get_cluster_markers(
                adata=adata,
                cluster_label=cluster_name,
                marker_ref=marker_ref,
                ngenes=100,
                filtered=True,
            )
            # and store:
            marker_gene_df.to_csv(
                os.path.join(dir_HLCA_cluster_output, f"LCA_{cluster_name}_marker_genes_versus_{marker_ref}.csv")
            )

adata.write("../data/HLCA_v2_intermediates/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_log1p.h5ad")