#!/bin/bash

#BSUB -o /lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/09_concat_processed_h5ads.out
#BSUB -G team292
#BSUB -J 09_concat_processed_h5ads
#BSUB -W300
#BSUB -n 1
#BSUB -q normal
#BSUB -R "select[mem>500000] rusage[mem=500000]"
#BSUB -M500000

SCRIPT=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/ATAC_processing/09_concat_processed_h5ads.py

module load ISG/conda
conda activate /software/cellgen/team292/cc53/my-conda-envs/snapATAC2_env

python $SCRIPT