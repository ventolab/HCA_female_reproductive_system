#!/usr/bin/env python3

import argparse
import os
import scanpy as sc
import numpy as np
import pandas as pd
import snapatac2 as snap
import anndata as ad
import matplotlib.pyplot as plt
import pyranges as pr #need to get pyranges for this
from matplotlib.backends.backend_pdf import PdfPages
import harmonypy as hm
ad.settings.allow_write_nullable_strings = True

def obsm_to_csv(adata,X_name,csv):
    X_out = adata.obsm[X_name]
    df_out = pd.DataFrame(
    	X_out,
    	index=adata.obs_names,
    	columns=[f"X_{i}" for i in range(X_out.shape[1])])
	# Save to CSV
    df_out.to_csv(csv)

#-------------------Args setup---------------------

parser = argparse.ArgumentParser(description="Convert h5ad to sparse")
parser.add_argument('--h5ad_dir', type=str, required=True, help="input anndata")
parser.add_argument('--meta', type=str, required=True, help="input metadata")
parser.add_argument('--latent_outdir', type=str, required=True, help="input anndata")
parser.add_argument('--plotdir', type=str, required=True, help="input anndata")
parser.add_argument('--n_features', type=int, required=False, help="number of features", default=500000)
parser.add_argument('--batch_correction', type=str, required=False, help="type of batch correction: None, Harmony") #None, Harmony
parser.add_argument('--write_h5ad', action='store_true', help="write h5ad at the end")

args = parser.parse_args()

h5ad_dir=args.h5ad_dir
meta_fname=args.meta
latent_outdir=args.latent_outdir
plotdir=args.plotdir

n_features=args.n_features
batch_correction=args.batch_correction

#results_plots="/nfs/team292/projects/PanTissue/data/temp/ATAC/ovary_sanger/processed/snapATAC2/plots/"
#results_dir="/nfs/team292/projects/PanTissue/data/temp/ATAC/ovary_sanger/processed/snapATAC2/save_latents_obs/"
#outh5ad_dir="/nfs/team292/projects/PanTissue/data/temp/ATAC/ovary_sanger/processed/snapATAC2/concatenated_h5ads/"

blacklist_bed="/lustre/scratch125/cellgen/vento/cc53/utils/hg38-blacklist.v2.bed"
#gtf_fname="/software/cellgen/cellgeni/refdata_10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"

comment=""

#for plotting and writing files
features_suffix=str(n_features)+"features"

if batch_correction is None:
    batchcor_suffix=""
else:
    batchcor_suffix=batch_correction+"_"


#-------------------Args setup---------------------

print("Loading data...")
adata=ad.read(h5ad_dir+'concatenated.h5ad')
print("Data loaded:")
print(adata)

#clean up any keys tha will cause errors
# clean up obsm keys
obsm_to_drop = ['X_spectral', 'X_spectral_harmony']
for key in obsm_to_drop:
    if key in adata.obsm:
        print(f"Dropping existing obsm key: {key}")
        del adata.obsm[key]

#add metadata
meta=pd.read_csv(meta_fname)
meta=meta[['Organ','Dataset','Sample_id','Library_id','Donor_id','Sex','Stage','Menstrual_stage']]
for col in meta.select_dtypes(include='object').columns:
    meta[col] = meta[col].astype(pd.StringDtype())
meta['sample'] = meta['Sample_id']
meta = meta.set_index('sample')
print(meta.head())

cols_to_drop = adata.obs.columns.intersection(meta.columns) #also cleanup
if len(cols_to_drop) > 0:
    print(f"Dropping existing meta obs columns: {list(cols_to_drop)}")
    adata.obs = adata.obs.drop(columns=cols_to_drop)
adata.obs = adata.obs.join(meta, on='sample')

print("Selecting features...")
snap.pp.select_features(adata, n_features=int(n_features), inplace=True, blacklist=blacklist_bed)

print("Features selected:")
print(adata)

print("Getting latent...")
X_spectral=pd.read_csv(latent_outdir+'X_spectral_'+features_suffix+comment+'.csv', index_col=0)
adata.obsm["X_spectral"]= X_spectral.loc[adata.obs_names].values

if batch_correction=="Harmony":
    X_spectral_harmony=pd.read_csv(latent_outdir+'X_spectral_Harmony_'+features_suffix+comment+'.csv',index_col=0)#reread in because was having issues with formatting
    adata.obsm["X_spectral_harmony"]= X_spectral_harmony.loc[adata.obs_names].values

    X_umap=pd.read_csv(latent_outdir+'X_spectral_Harmony_umap_'+features_suffix+comment+'.csv',index_col=0)
    adata.obsm["X_umap"]= X_umap.loc[adata.obs_names].values

    obs=pd.read_csv(latent_outdir+'leiden_obs_Harmony_'+features_suffix+comment+'.csv',index_col=0)
    cols_to_drop = adata.obs.columns.intersection(obs.columns) #also cleanup
    obs = obs.drop(columns=cols_to_drop)
    adata.obs = adata.obs.join(obs)

    #add scanpy embedding
    X_umap_scanpy=pd.read_csv(latent_outdir+'X_spectral_Harmony_umap_scanpy_'+features_suffix+comment+'.csv',index_col=0)
    adata.obsm["X_umap_scanpy"]= X_umap_scanpy.loc[adata.obs_names].values

    obs=pd.read_csv(latent_outdir+'leiden_obs_Harmony_wscanpy_'+features_suffix+comment+'.csv',index_col=0)
    cols_to_drop = adata.obs.columns.intersection(obs.columns) #also cleanup
    obs = obs.drop(columns=cols_to_drop)
    adata.obs = adata.obs.join(obs)

for col in meta.select_dtypes(include='object').columns:
    meta[col] = meta[col].astype(pd.StringDtype())

#save adata with umap and obs so can plot backed
if args.write_h5ad:
    adata.write(h5ad_dir+'concatenated_processed.h5ad',compression='gzip')

print("Plotting...")

with PdfPages(plotdir+'snapATAC2_clustering_umaps_'+batchcor_suffix+features_suffix+comment+'.pdf') as pdf:
    sc.set_figure_params(dpi=300, frameon=False, figsize=(3.5, 3.5))
    for var in ["leiden","leiden_res5","leiden_res10","n_fragment",'tsse','doublet_score','sample',"leiden_scanpy","leiden_scanpy_res5","leiden_scanpy_res10"]:
        fig = sc.pl.umap(adata, color=var, return_fig=True, use_raw=False,show=False)
        pdf.savefig(fig, bbox_inches='tight', dpi=300)
        plt.close(fig)

with PdfPages(plotdir+'scanpy_clustering_umaps_'+batchcor_suffix+features_suffix+comment+'.pdf') as pdf:
    sc.set_figure_params(dpi=300, frameon=False, figsize=(3.5, 3.5))
    for var in ["leiden","leiden_res5","leiden_res10","n_fragment",'tsse','doublet_score','sample',"leiden_scanpy","leiden_scanpy_res5","leiden_scanpy_res10"]:
        fig = sc.pl.embedding(adata, basis='X_umap_scanpy',color=var, return_fig=True, use_raw=False,show=False)
        pdf.savefig(fig, bbox_inches='tight', dpi=300)
        plt.close(fig)

print("All done!")



