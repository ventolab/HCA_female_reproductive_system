#!/bin/bash

#BSUB -o /lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/03_concat_h5ads.out
#BSUB -G team292
#BSUB -J 03_concat_h5ads
#BSUB -W600
#BSUB -n 1
#BSUB -q normal
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -M2000

BASE_DIR=/lustre/scratch125/cellgen/vento/projects/PanTissue/data/temp/ATAC/processed
SCRIPT=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/ATAC_processing/03_concat_h5ads.py

LOGDIR=/lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/03_concat_h5ads
mkdir -p $LOGDIR

chmod a+x $SCRIPT 
module load cellgen/conda
conda activate /software/cellgen/team292/cc53/my-conda-envs/snapATAC2_env

for dir in fetal_ovary fetal_reptract nonfetal_fallopiantube nonfetal_ovary nonfetal_uterus; do
    DATA_DIR="$BASE_DIR/$dir/persample_h5ads/"
    OUTDIR="$BASE_DIR/$dir/"
    bsub -o $LOGDIR/03_concat_h5ads.$dir.out \
            -G team292 -J 03_concat_h5ads.$dir -W600 -n 2 -R "select[mem>500000] rusage[mem=500000]" -M500000 \
            $SCRIPT --h5ads_dir $DATA_DIR --outdir $OUTDIR
done