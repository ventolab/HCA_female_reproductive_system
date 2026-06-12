#!/bin/bash

#BSUB -o /lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/01_remove_empty_drops_loopsub.out
#BSUB -G team292
#BSUB -J 01_remove_empty_drops_loopsub
#BSUB -W600
#BSUB -n 1
#BSUB -q normal
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -M2000

BASE_INDIR=/nfs/team292/projects/PanTissue/data/freeze/ATAC/cellranger-atac200 #should contain peaks and fragments
BASE_OUTDIR=/lustre/scratch125/cellgen/vento/projects/PanTissue/data/temp/ATAC/processed
SCRIPT=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/ATAC_processing/01_remove_empty_drops.sh
LOGDIR=/lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs

chmod a+x $SCRIPT 

for dir in fetal_ovary fetal_reptract nonfetal_fallopiantube nonfetal_ovary nonfetal_uterus; do
    DATA_DIR="$BASE_INDIR/$dir/"
    OUTDIR="$BASE_OUTDIR/$dir/filtered_fragments/"
    mkdir -p $OUTDIR
    bsub -o $LOGDIR/01_remove_empty_drops.$dir.out \
            -G team292 -J 01_remove_empty_drops.$dir -W600 -n 2 -R "select[mem>50000] rusage[mem=50000]" -M50000 \
            $SCRIPT $DATA_DIR $OUTDIR
done

