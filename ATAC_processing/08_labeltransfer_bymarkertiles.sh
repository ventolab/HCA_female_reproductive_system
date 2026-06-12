#!/bin/bash

#BSUB -o /lustre/scratch125/cellgen/vento/cc53/trait2cell/logs/labeltransfer_bymarkertiles.out
#BSUB -G team292
#BSUB -J labeltransfer_bymarkertiles
#BSUB -W300
#BSUB -n 1
#BSUB -q normal
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -M2000

BASE_DIR=/lustre/scratch125/cellgen/vento/projects/PanTissue/data/temp/ATAC/processed
RNA_DIR=/nfs/team292/projects/PanTissue/data/temp/GWASxRNA/RNA_references
SCRIPT=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/ATAC_processing/08_labeltransfer_bymarkertiles.R

LOGDIR=/lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/08_labeltransfer_bymarkertiles
mkdir -p $LOGDIR

chmod a+x $SCRIPT 
module load HGI/softpack/users/cc53/R_base

for dir in nonfetal_ovary ; do #nonfetal_uterus nonfetal_fallopiantube
    ATAC_PROJ="$BASE_DIR/$dir/ArchRproj"
    ATAC_LATENT="$BASE_DIR/$dir/save_latent_obs/X_spectral_Harmony_500000features.csv"
    RNA_H5AD="$RNA_DIR/$dir.h5ad"
    RNA_MARKERS="$RNA_DIR/${dir}_top100markers_fine_celltype.csv"
    bsub -o $LOGDIR/08_labeltransfer_bymarkertiles.$dir.out -q normal \
            -G team292 -J 08_labeltransfer_bymarkertiles.$dir -W600 -n 4 -R "select[mem>300000] rusage[mem=300000]" -M300000 \
            $SCRIPT "fine_celltype" $ATAC_PROJ $ATAC_LATENT $RNA_H5AD $RNA_MARKERS
done