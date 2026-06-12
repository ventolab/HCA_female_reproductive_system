#!/bin/bash

#BSUB -o /lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/ENCODE_enrichment_permute.out
#BSUB -G team292
#BSUB -J ENCODE_enrichment_permute
#BSUB -W600
#BSUB -n 1
#BSUB -q normal
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -M2000

SCRIPT=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/benchmark_peaks/ENCODE_enrichment_permute.R

LOGDIR=/lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/ENCODE_enrichment_permute
mkdir -p $LOGDIR

module load HGI/softpack/users/cc53/R_base/19

for PEAKS in allpeaks linkedpeaks; do
    for ENCODE_ANNOT in pELS CA-CTCF CA-TF CA dELS TF CA-H3K4me3 PLS; do
        bsub -o $LOGDIR/$PEAKS.$ENCODE_ANNOT.out \
            -G team292 -J ENCODE_enrichment_permute.$PEAKS.$ENCODE_ANNOT -W600 -n 8 -R "select[mem>50000] rusage[mem=50000]" -M50000 \
            Rscript $SCRIPT $PEAKS $ENCODE_ANNOT intergenic
    done
done