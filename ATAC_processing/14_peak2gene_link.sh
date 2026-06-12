#!/bin/bash

#BSUB -o /lustre/scratch125/cellgen/vento/cc53/trait2cell/logs/14_peak2gene_links.out
#BSUB -G team292
#BSUB -J 14_peak2gene_links
#BSUB -W180
#BSUB -n 8
#BSUB -q normal
#BSUB -R "select[mem>25000] rusage[mem=25000]"
#BSUB -M25000

SCRIPT=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/ATAC_processing/14_peak2gene_links.R

module load HGI/softpack/users/cc53/R_base

Rscript $SCRIPT