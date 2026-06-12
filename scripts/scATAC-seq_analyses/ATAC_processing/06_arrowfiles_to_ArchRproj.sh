#!/bin/bash

#BSUB -o /lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/06_arrowfiles_to_ArchRproj.out
#BSUB -G team292
#BSUB -J 06_arrowfiles_to_ArchRproj
#BSUB -W600
#BSUB -n 1
#BSUB -q normal
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -M2000

BASE_DIR=/lustre/scratch125/cellgen/vento/projects/PanTissue/data/temp/ATAC/processed
SCRIPT=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/ATAC_processing/06_arrowfiles_to_ArchRproj.R

LOGDIR=/lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/06_arrowfiles_to_ArchRproj
mkdir -p $LOGDIR

chmod a+x $SCRIPT 
module load HGI/softpack/users/cc53/R_base

for dir in fetal_ovary nonfetal_fallopiantube nonfetal_ovary nonfetal_uterus; do #fetal_reptract
    ARROWS_DIR="$BASE_DIR/$dir/arrowfiles/"
    PROJ_DIR="$BASE_DIR/$dir/ArchRproj/"
    mkdir -p $PROJ_DIR
    bsub -o $LOGDIR/06_arrowfiles_to_ArchRproj.$dir.out \
            -G team292 -J 06_arrowfiles_to_ArchRproj.$dir -W600 -n 2 -R "select[mem>50000] rusage[mem=50000]" -M50000 \
            $SCRIPT $PROJ_DIR $ARROWS_DIR
    done
done
