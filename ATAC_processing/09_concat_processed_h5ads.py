#!/usr/bin/env python3

import argparse
import os
import scanpy as sc
import pandas as pd
import numpy as np
import snapatac2 as snap
import anndata as ad
from anndata.experimental import concat_on_disk
import scglue
import matplotlib as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import sparse

#parser = argparse.ArgumentParser(description="Convert fragments to h5ad")
#parser.add_argument('--h5ads_dir', type=str, required=True, help="h5ads directory")
#parser.add_argument('--outdir', type=str, required=True, help="Output directory")

#args = parser.parse_args()
#h5ads_dir=args.h5ads_dir
#outdir=args.outdir

outdir='/nfs/team292/projects/PanTissue/data/temp/ATAC/processed/nonfetal_all/'
annot_dir='/nfs/team292/projects/PanTissue/results/freeze/ATAC/annotations/'
in_files={'nonfetal_fallopiantube':'/nfs/team292/projects/PanTissue/data/temp/ATAC/processed/nonfetal_fallopiantube/concatenated.h5ad',
            'nonfetal_ovary':'/nfs/team292/projects/PanTissue/data/temp/ATAC/processed/nonfetal_ovary/concatenated.h5ad',
            'nonfetal_uterus':'/nfs/team292/projects/PanTissue/data/temp/ATAC/processed/nonfetal_uterus/concatenated.h5ad'}

#data = snap.AnnDataSet(
    #adatas=[(donor_id, adata) for donor_id, adata in zip(donor_ids, adatas)],)

ref = ad.read_h5ad(in_files['nonfetal_fallopiantube'], backed='r')
saved_var = ref.var.copy()
saved_uns  = ref.uns.copy()
ref.file.close()

concat_on_disk(
        in_files = in_files,
        out_file = outdir+"concatenated.h5ad",
        join     = "inner",
        label    = "tissue",
        uns_merge="same")

#concat annotations
files = ['ovary.csv', 'fallopiantube.csv', 'uterus.csv']
annot = pd.concat([pd.read_csv(annot_dir+f, index_col=0) for f in files])
annot=annot[['Sample', 'TSSEnrichment', 'ReadsInTSS', 'ReadsInPromoter',
       'ReadsInBlacklist', 'PromoterRatio', 'PassQC', 'NucleosomeRatio',
       'nMultiFrags', 'nMonoFrags', 'nFrags', 'nDiFrags', 'DoubletScore',
       'DoubletEnrichment', 'BlacklistRatio', 'is_doublet',
       'predicted_cell_unconstrained', 'predicted_annotation_unconstrained',
       'predicted_score_unconstrained', 'n_fragment', 'frac_dup', 'frac_mito',
       'tsse', 'doublet_probability', 'doublet_score', 'sample', 'in_HCAv1',
       'exclusion-inclusion_criteria', 'Organ', 'Dataset', 'Sample_id',
       'Library_id', 'Donor_id', 'Sex', 'Ancestry', 'Stage',
       'Postnatal_age_years', 'Gestational_age_pcw', 'Developmental_stage',
       'Tanner Stage', 'Menstrual_stage', 'Disease', 'Clinical_diagnosis',
       'Organ_part', 'Tissue_ROI', 'Specimen_type', 'Sampled_site_condition',
       'Observed_pathology', 'Tissue_status', 'Dissociation_method',
       'Preservation_method', 'Cell_enrichment', 'Target_cell_population',
       'Sorting_method', 'Assay_type', 'Library_chemistry',
       'Sequencing_platform', 'Multiplexed', 'Dataset_id', 'Batch',
       'Collection_site', 'Sample_id_other', 'Cellranger_version', 'RNA_matched_sample_id', 'leiden',
       'leiden_res5', 'leiden_res10', 'leiden_scanpy', 'leiden_scanpy_res5',
       'leiden_scanpy_res10', 'coarse_predicted_annotation', 'final_tissue_annotation']]
annot=annot[~annot.final_tissue_annotation.isin(['doublet','lowQC'])]

adata=ad.read_h5ad(outdir+"concatenated.h5ad")
adata=adata[adata.obs_names.isin(annot.index)].copy()
adata.var = saved_var.loc[adata.var_names]
adata.uns = saved_uns
adata.X = sparse.csr_matrix(adata.X)
adata.raw = None

overlap_cols = adata.obs.columns.intersection(annot.columns).tolist()
adata.obs = adata.obs.drop(columns=overlap_cols)
adata.obs = adata.obs.join(annot)
adata = adata[adata.obs['final_tissue_annotation'].notna()].copy()

adata.write(outdir + "concatenated.h5ad",compression='gzip')

#data.write(outh5ad_dir+"ovaries.h5ad", compression='gzip')
