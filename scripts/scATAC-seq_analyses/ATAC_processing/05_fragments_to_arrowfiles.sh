#!/bin/bash

#BSUB -o /lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/05_fragments_to_arrowfiles.out
#BSUB -G team292
#BSUB -J 05_fragments_to_arrowfiles
#BSUB -W600
#BSUB -n 1
#BSUB -q normal
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -M2000

BASE_DIR=/lustre/scratch125/cellgen/vento/projects/PanTissue/data/temp/ATAC/processed
SCRIPT=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/ATAC_processing/05_fragments_to_arrowfiles.R

LOGDIR=/lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/05_fragments_to_arrowfiles
mkdir -p $LOGDIR

chmod a+x $SCRIPT 
module load HGI/softpack/users/cc53/R_base

for dir in fetal_ovary fetal_reptract nonfetal_fallopiantube nonfetal_ovary nonfetal_uterus; do
    DATA_DIR="$BASE_DIR/$dir/filtered_fragments/"
    for FILE in $DATA_DIR/*-fragments_filt.tsv.gz; do
        SAMPLE=$(basename $FILE -fragments_filt.tsv.gz)
        OUTDIR="$BASE_DIR/$dir/arrowfiles/$SAMPLE-arrow/"
        mkdir -p $OUTDIR
        bsub -o $LOGDIR/05_fragments_to_arrowfiles.$dir.$SAMPLE.out \
            -G team292 -J 05_fragments_to_arrowfiles.$dir.$SAMPLE -W600 -n 2 -R "select[mem>50000] rusage[mem=50000]" -M50000 \
            $SCRIPT $SAMPLE $FILE $OUTDIR
    done
done
