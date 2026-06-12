#!/usr/bin/env python3

import argparse
import os
import scanpy as sc
import numpy as np
import snapatac2 as snap
import anndata as ad
from anndata.experimental import concat_on_disk
import scglue
import matplotlib as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import sparse

parser = argparse.ArgumentParser(description="Convert fragments to h5ad")
parser.add_argument('--h5ads_dir', type=str, required=True, help="h5ads directory")
parser.add_argument('--outdir', type=str, required=True, help="Output directory")

args = parser.parse_args()
h5ads_dir=args.h5ads_dir
outdir=args.outdir

#h5ads_dir="/nfs/team292/projects/PanTissue/data/temp/ATAC/ovary_sanger/processed/snapATAC2/persample_h5ads/"
#outh5ad_dir="/nfs/team292/projects/PanTissue/data/temp/ATAC/ovary_sanger/processed/snapATAC2/concatenated_h5ads/"

# Get only .tsv.gz files (not .tbi)
h5ad_files = sorted([
    os.path.join(h5ads_dir, f) 
    for f in os.listdir(h5ads_dir) 
    if f.endswith(".h5ad")])

sample_ids = [os.path.basename(f).split(".")[0] for f in h5ad_files]
in_files = dict(zip(sample_ids, h5ad_files))

#exclude = {'Hrv58multi', 'Hrv58multiBis','Hrv65multi','Hrv65multiBis','AUrv22','AUrv28'}
#in_files = {k: v for k, v in in_files.items() if k not in exclude}

#data = snap.AnnDataSet(
    #adatas=[(donor_id, adata) for donor_id, adata in zip(donor_ids, adatas)],)

ref = ad.read_h5ad(h5ad_files[0], backed='r')
saved_var = ref.var.copy()
saved_uns  = ref.uns.copy()
ref.file.close()

concat_on_disk(
        in_files = in_files,
        out_file = outdir+"concatenated.h5ad",
        join     = "inner",
        label    = "sample",
        uns_merge="same")

adata=ad.read(outdir+"concatenated.h5ad")
adata.var = saved_var.loc[adata.var_names]
adata.uns = saved_uns
adata.X = sparse.csr_matrix(adata.X)
adata.raw = None
adata.write(outdir + "concatenated.h5ad",compression='gzip')

#data.write(outh5ad_dir+"ovaries.h5ad", compression='gzip')
