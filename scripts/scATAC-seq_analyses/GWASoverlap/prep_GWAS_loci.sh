#!/bin/bash

#BSUB -o /lustre/scratch125/cellgen/vento/cc53/trait2cell/logs/get_GWAS_loci.out
#BSUB -G team292
#BSUB -J get_GWAS_loci
#BSUB -W180
#BSUB -n 2
#BSUB -q normal
#BSUB -R "select[mem>30000] rusage[mem=30000]"
#BSUB -M30000

SCRIPT=/nfs/team292/cc53/scripts/trait2cell/atac_processing/peaks2genes/get_GWAS_loci.R
INPUT_CSV=/nfs/team292/projects/PanTissue/data/freeze/ATAC/GWAS_loci/GWAS_indexSNPs.csv
LOCIDIR=/nfs/team292/projects/PanTissue/data/freeze/ATAC/GWAS_loci

module load HGI/softpack/users/cc53/R_base

LIFTOVER_PREP=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/GWASoverlap/prep4liftover.R
Rscript $LIFTOVER_PREP $INPUT_CSV \
    $LOCIDIR/GWASloci_snp_hg19.bed \
    $LOCIDIR/GWASloci_intervals_hg19.bed

CHAIN_HG19_TO_HG38=/lustre/scratch125/cellgen/vento/cc53/utils/hg19ToHg38.over.chain.gz
# Lift SNP positions
/software/team152/bh18/liftOver \
    $LOCIDIR/GWASloci_snp_hg19.bed \
    $CHAIN_HG19_TO_HG38 \
    $LOCIDIR/GWASloci_snp_hg38.bed \
    $LOCIDIR/GWASloci_snp_hg19.unmapped.bed

# Lift intervals
/software/team152/bh18/liftOver \
    $LOCIDIR/GWASloci_intervals_hg19.bed \
    $CHAIN_HG19_TO_HG38 \
    $LOCIDIR/GWASloci_intervals_hg38.bed \
    $LOCIDIR/GWASloci_intervals_hg19.unmapped.bed

LIFTOVER_PREP=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/GWASoverlap/postliftover.R
Rscript $LIFTOVER_PREP $INPUT_CSV \
    $LOCIDIR/GWASloci_snp_hg38.bed \
    $LOCIDIR/GWASloci_intervals_hg38.bed \
    $LOCIDIR/GWAS_indexSNPs_hg38.csv