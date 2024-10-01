### python scripts for transferring label to teratoma time series data

source anaconda3/ENTER/etc/profile.d/conda.sh
conda create -n scvi-env python=3.9
conda activate scvi-env
conda install scvi-tools -c conda-forge

import scanpy
import numpy as np
import pandas as pd
import scvi
import matplotlib.pyplot as plt
import seaborn as sns

# code has to be run in command line python because python version in jupyter notebook is outdated

import os
os.getcwd()
os.chdir("/media/Scratch_SSD_Voyager/sammi/Teratoma_TS")

# the genes in the anndata ref and anndata query need to match - currently no way to save gene names for some reason
# so just merged the query and kept the same order of genes in both ref & query

# read ref & prep
adata_ref = scanpy.read_h5ad("/media/Scratch_SSD_Voyager/sammi/Teratoma_TS/h5ad/"+"4H1_20230119"+".h5ad")
adata_ref.layers["counts"] = adata_ref.X.copy()

scanpy.pp.normalize_total(adata_ref, target_sum=1e4)
scanpy.pp.log1p(adata_ref)

#adata_ref.raw = adata_ref

scanpy.pp.highly_variable_genes(
    adata_ref,
    n_top_genes=2000,
    batch_key="teratoma",
    subset=True
)

# read query
adata_query = scanpy.read_h5ad("/media/Scratch_SSD_Voyager/sammi/Teratoma_TS/h5ad/"+"TS_all_RNA_int_20230119"+".h5ad")

adata_query.layers["counts"] = adata_query.X.copy()

scanpy.pp.normalize_total(adata_query, target_sum=1e4)
scanpy.pp.log1p(adata_query)

#adata_query.raw = adata_query

adata_query = adata_query[:, adata_ref.var_names].copy()


scvi.model.SCVI.setup_anndata(adata_ref, batch_key="teratoma", layer = "counts")

# define parm & run SCVI
arches_params = dict(
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=0.2,
    n_layers=2,
)

vae_ref = scvi.model.SCVI(
    adata_ref,
    **arches_params
)
vae_ref.train()

# save model
dir_path_scvi = "/media/Scratch_SSD_Voyager/sammi/Teratoma_TS/tera_model_scvi_20230119/"
vae_ref.save(dir_path_scvi, overwrite=True)

# save SCVI plotting
adata_ref.obsm["X_scVI"] = vae_ref.get_latent_representation()
scanpy.pp.neighbors(adata_ref, use_rep="X_scVI")
scanpy.tl.leiden(adata_ref)
scanpy.tl.umap(adata_ref)

scanpy.pl.umap(
    adata_ref,
    color=["teratoma", "cluster"],
    frameon=False,
    ncols=1,
    save = "terA_scvi_dimplot_20230119.pdf"
)

# run SCANVI
vae_ref_scan = scvi.model.SCANVI.from_scvi_model(
    vae_ref,
    unlabeled_category="Unknown",
    labels_key="cluster"
)
vae_ref_scan.train(max_epochs=20, n_samples_per_label=100)

# save model
dir_path_scan = "/media/Scratch_SSD_Voyager/sammi/Teratoma_TS/tera_model_scanvi_20230119/"
vae_ref_scan.save(dir_path_scan, overwrite=True)

vae_ref_scan = scvi.model.SCANVI.load(dir_path_scan, adata_ref)


# save SCANVI plotting
adata_ref.obsm["X_scANVI"] = vae_ref_scan.get_latent_representation()
scanpy.pp.neighbors(adata_ref, use_rep="X_scANVI")
scanpy.tl.leiden(adata_ref)
scanpy.tl.umap(adata_ref)

scanpy.pl.umap(
    adata_ref,
    color=["teratoma", "cluster"],
    frameon=False,
    ncols=1,
    save = "terA_scanvi_dimplot_20230119.pdf"
)

adata_query.obs['teratoma'] = adata_query.obs.mouse

# online update query
vae_q = scvi.model.SCANVI.load_query_data(
    adata_query,
    dir_path_scan,
)

# decreasing epoch rn for taking too long, max_epochs=100,check_val_every_n_epoch=10
vae_q.train(
    max_epochs=50,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=5,
)

adata_query.obsm["X_scANVI"] = vae_q.get_latent_representation()
adata_query.obs["scANVI_predictions"] = vae_q.predict()
adata_query.obs["scANVI_prediction_scores"] = vae_q.predict(soft = True).max(axis=1)

df = adata_query.obs.groupby(["predicted.celltype", "scANVI_predictions"]).size().unstack(fill_value=0)
norm_df = df / df.sum(axis=0)

plt.figure(figsize=(8, 8))
_ = plt.pcolor(norm_df)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("scANVI Predicted")
plt.ylabel("Seurat Label Transfer")
plt.savefig('figures/confusion_matrix_20230115.pdf', bbox_inches='tight')  

vae_q.predict(soft = True)
maxValues = vae_q.predict(soft = True).max(axis=1)
plt.figure(figsize=(8, 8))
plt.hist(maxValues)
plt.xlabel("scANVI prediction score")
plt.ylabel("# of cells")
plt.savefig('figures/scANVI_pred_hist_20230115.pdf', bbox_inches='tight')  

# count how many labels align
len(adata_query.obs['predicted.celltype'] == adata_query.obs['scANVI_predictions'])
sum(adata_query.obs['predicted.celltype'] == adata_query.obs['scANVI_predictions'])

# save prediction score
adata_query.obs['scANVI_prediction_score'] = maxValues

adata_query[adata_query.obs['predicted.celltype'] == adata_query.obs['scANVI_predictions']].obs.age.value_counts()

adata_ref.__dict__['_raw'].__dict__['_var'] = adata_ref.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
adata_ref.__dict__['_var'] = adata_ref.__dict__['_var'].rename(columns={'_index': 'features'})
adata_ref.write(filename = '/media/Scratch_SSD_Voyager/sammi/Teratoma_TS/h5ad/7H1_scanvi_20230121.h5ad')

adata_query.__dict__['_raw'].__dict__['_var'] = adata_query.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
adata_query.__dict__['_var'] = adata_query.__dict__['_var'].rename(columns={'_index': 'features'})
adata_query.write(filename = '/media/Scratch_SSD_Voyager/sammi/Teratoma_TS/h5ad/TS_RNA_scanvi_20230121.h5ad')

adata_query.obs.to_csv('anndata_anno_20230121.csv', sep = ",")

