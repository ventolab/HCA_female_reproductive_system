#!/bin/bash

#BSUB -o /lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/04_integrate_by_sample.out
#BSUB -G team292
#BSUB -J 04_integrate_by_sample
#BSUB -W600
#BSUB -n 1
#BSUB -q normal
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -M2000

BASE_DIR=/lustre/scratch125/cellgen/vento/projects/PanTissue/data/temp/ATAC/processed
#SCRIPT=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/ATAC_processing/04_integrate_by_sample.py
SCRIPT=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/ATAC_processing/04bis_integrate_by_sample_rerun.py
META=/nfs/team292/projects/PanTissue/data/freeze/ATAC/metadata/ATAC_samples_panrepro_metadata.csv

LOGDIR=/lustre/scratch125/cellgen/vento/projects/PanTissue/results/temp/ATAC/logs/04_integrate_by_sample
mkdir -p $LOGDIR

chmod a+x $SCRIPT 
module load cellgen/conda
conda activate /software/cellgen/team292/cc53/my-conda-envs/snapATAC2_env

#for dir in nonfetal_uterus ; do ##fetal_ovary fetal_reptract nonfetal_fallopiantube nonfetal_ovary
for dir in nonfetal_ovary; do #fetal_ovary fetal_reptract nonfetal_fallopiantube 
    #H5AD="$BASE_DIR/$dir/concatenated.h5ad"
    H5AD_DIR="$BASE_DIR/$dir/"
    PLOTDIR="$BASE_DIR/$dir/plots/"
    LATENT_OUTDIR="$BASE_DIR/$dir/save_latent_obs/"
    mkdir -p $PLOTDIR
    mkdir -p $LATENT_OUTDIR
    bsub -o $LOGDIR/04_integrate_by_sample.$dir.Harmony.out \
            -G team292 -J 04_integrate_by_sample.Harmony.$dir -W300 -n 1 -R "select[mem>200000] rusage[mem=200000]" -M200000 \
            $SCRIPT --h5ad_dir $H5AD_DIR --meta $META --latent_outdir $LATENT_OUTDIR --plotdir $PLOTDIR --batch_correction "Harmony" --write_h5ad
    #bsub -o $LOGDIR/04_integrate_by_sample.$dir.out \
            #-G team292 -J 04_integrate_by_sample.$dir -W600 -n 1 -R "select[mem>300000] rusage[mem=300000]" -M300000 \
            #$SCRIPT --h5ad $H5AD --meta $META --latent_outdir $LATENT_OUTDIR --plotdir $PLOTDIR
done