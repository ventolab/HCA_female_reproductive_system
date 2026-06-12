#!/bin/bash

#BSUB -o /lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/02_fragments_to_h5ad.out
#BSUB -G team292
#BSUB -J 02_fragments_to_h5ad
#BSUB -W600
#BSUB -n 1
#BSUB -q normal
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -M2000

BASE_DIR=/lustre/scratch125/cellgen/vento/projects/PanTissue/data/temp/ATAC/processed
SCRIPT=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/ATAC_processing/02_fragments_to_h5ad.py

LOGDIR=/lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/02_fragments_to_h5ad

chmod a+x $SCRIPT 
module load cellgen/conda
conda activate /software/cellgen/team292/cc53/my-conda-envs/snapATAC2_env

for dir in fetal_ovary fetal_reptract nonfetal_fallopiantube nonfetal_ovary nonfetal_uterus; do
    DATA_DIR="$BASE_DIR/$dir/filtered_fragments/"
    OUTDIR="$BASE_DIR/$dir/persample_h5ads/"
    PLOTDIR="$BASE_DIR/$dir/persample_h5ads/QCplots/"
    mkdir -p $OUTDIR
    mkdir -p $PLOTDIR
    for FILE in $DATA_DIR/*-fragments_filt.tsv.gz; do
        SAMPLE=$(basename $FILE -fragments_filt.tsv.gz)
        bsub -o $LOGDIR/02_fragments_to_h5ad.$dir.$SAMPLE.out \
            -G team292 -J 02_fragments_to_h5ad.$dir.$SAMPLE -W600 -n 2 -R "select[mem>50000] rusage[mem=50000]" -M50000 \
            $SCRIPT --sampleid $SAMPLE --fragfile $FILE --outdir $OUTDIR --plotdir $PLOTDIR
    done
done
