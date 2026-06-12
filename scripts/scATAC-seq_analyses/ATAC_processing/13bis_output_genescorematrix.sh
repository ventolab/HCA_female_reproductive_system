#!/bin/bash

#BSUB -o /lustre/scratch125/cellgen/vento/cc53/trait2cell/logs/13bis_output_genescorematrix.out
#BSUB -G team292
#BSUB -J 13bis_output_genescorematrix
#BSUB -W600
#BSUB -n 8
#BSUB -q normal
#BSUB -R "select[mem>400000] rusage[mem=400000]"
#BSUB -M400000

SCRIPT=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/ATAC_processing/13bis_output_genescorematrix.R

cd /lustre/scratch125/cellgen/vento/cc53/trait2cell/ovary/scATAC/data

module load HGI/softpack/users/cc53/R_base

Rscript $SCRIPT