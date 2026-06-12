#!/bin/bash

#BSUB -o /lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/07_get_marker_genes.out
#BSUB -G team292
#BSUB -J 07_get_marker_genes
#BSUB -W30
#BSUB -n 1
#BSUB -q small
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -M2000

BASE_DIR=/nfs/team292/projects/PanTissue/data/temp/GWASxRNA/RNA_references/
SCRIPT=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/ATAC_processing/07_get_marker_genes.py
chmod a+x $SCRIPT 

LOGDIR=/lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/07_get_marker_genes
mkdir -p $LOGDIR

module load cellgen/conda
conda activate /software/cellgen/team292/cc53/my-conda-envs/singlecellatlas_cc53

for analysis in nonfetal_uterus nonfetal_ovary nonfetal_fallopiantube; do
    OUTDIR="$BASE_DIR/${analysis}_"
    H5AD="$BASE_DIR/$analysis.h5ad"
    bsub -o $LOGDIR/07_get_marker_genes.$analysis.out -q small -W30 \
            -G team292 -J 07_get_marker_genes.$analysis -W600 -n 2 -R "select[mem>150000] rusage[mem=150000]" -M150000 \
            $SCRIPT --h5ad $H5AD --outdir $OUTDIR --cellannot "fine_celltype"
done
