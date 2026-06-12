#!/bin/bash

#BSUB -o /lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/10_integrate_across_tissues.out
#BSUB -G team292
#BSUB -J 10_integrate_across_tissues
#BSUB -W600
#BSUB -n 1
#BSUB -q normal
#BSUB -R "select[mem>300000] rusage[mem=300000]"
#BSUB -M300000

SCRIPT=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/ATAC_processing/10_integrate_across_tissues.py

module load cellgen/conda
conda activate /software/cellgen/team292/cc53/my-conda-envs/snapATAC2_env

python $SCRIPT