#!/usr/bin/env python3

#do tfidf per cell type, get union of top 50 genes per cell type
#output: subset scRNA-seq object with only relevant genes

import os
import sys
import warnings
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import anndata as ad
import argparse

import pandas as pd
import matplotlib as mpl
import seaborn as sb
import scanpy.external as sce
import scipy
from matplotlib import rcParams
from matplotlib import colors
from matplotlib.backends.backend_pdf import PdfPages
import copy

#ad.settings.allow_write_nullable_strings = True

#import my helper functions
utils_path = '/nfs/team292/cc53/scripts/menstrual_fluid/scRNAseq_processing' 
sys.path.insert(0, utils_path)
import utils_noDEG as utils

parser = argparse.ArgumentParser(description="Get TFIDF top genes")
parser.add_argument('--h5ad', type=str, required=True, help="anndata path")
parser.add_argument('--outdir', type=str, required=True, help="output directory")
parser.add_argument('--cellannot', type=str, required=True, help="Cell type annotation")
args = parser.parse_args()

cellannotation=args.cellannot
outpath=args.outdir
h5ad=args.h5ad

adata=ad.read_h5ad(h5ad)
markers = utils.quick_markers(adata, cluster_key=cellannotation, n_markers=100, r_output=True)
keepgenes=markers['gene']

keepgenes.to_csv(outpath+'top100markers_'+cellannotation+'.csv', index=False)

