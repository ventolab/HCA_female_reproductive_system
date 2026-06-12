#!/usr/bin/env python3

import argparse
import os
import scanpy as sc
import numpy as np
import snapatac2 as snap
import anndata as ad
import scglue
import matplotlib as plt
from matplotlib.backends.backend_pdf import PdfPages

#QCplots_dir="/nfs/team292/projects/PanTissue/data/temp/ATAC/ovary_sanger/processed/snapATAC2/QCplots/"
#h5ads_dir="/nfs/team292/projects/PanTissue/data/temp/ATAC/ovary_sanger/processed/snapATAC2/persample_h5ads/"

parser = argparse.ArgumentParser(description="Convert fragments to h5ad")
parser.add_argument('--fragfile', type=str, required=True, help="Fragments file")
parser.add_argument('--outdir', type=str, required=True, help="Output directory")
parser.add_argument('--plotdir', type=str, required=True, help="Plotting directory")
parser.add_argument('--sampleid', type=str, required=True, help="Sample ID")

args = parser.parse_args()
fragfile=args.fragfile
outdir=args.outdir
plotdir=args.plotdir
sampleid=args.sampleid

#donor_id = os.path.basename(fragfile).split("-")[0]
print('Processing '+sampleid+'...')
adata=snap.pp.import_fragments(fragfile,chrom_sizes=snap.genome.hg38,sorted_by_barcode=False)

snap.metrics.tsse(adata, snap.genome.hg38)
fig=snap.pl.frag_size_distr(adata, interactive=False,show=False)
fig.write_image(plotdir+sampleid+'-fragsize.pdf')  
fig=snap.pl.tsse(adata, interactive=False,show=False)
fig.write_image(plotdir+sampleid+'-TSSEbyNFrags.pdf')  
snap.pp.filter_cells(adata, min_counts=10000) #no min_tsse=4 but stringent min counts to include germ cells
snap.pp.add_tile_matrix(adata, bin_size=500)
snap.pp.select_features(adata, n_features=1000000)
snap.pp.scrublet(adata)

adata.obs_names = sampleid + '#' + adata.obs_names #to match ArchR

adata.write(outdir+sampleid+".h5ad")