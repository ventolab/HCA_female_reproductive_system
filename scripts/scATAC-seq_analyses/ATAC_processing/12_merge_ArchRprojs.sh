#!/bin/bash

#BSUB -o /lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/12_merge_ArchRprojs.out
#BSUB -G team292
#BSUB -J 12_merge_ArchRprojs_lowmem
#BSUB -W300
#BSUB -n 1
#BSUB -q normal
#BSUB -R "select[mem>300000] rusage[mem=300000]"
#BSUB -M300000

SCRIPT=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/ATAC_processing/12_merge_ArchRprojs.R

module load HGI/softpack/users/cc53/R_base

Rscript $SCRIPT