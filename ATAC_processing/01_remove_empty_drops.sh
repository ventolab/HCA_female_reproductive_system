#!/bin/bash

DATA_DIR=$1
OUTDIR=$2

SCRIPT=/nfs/team292/cc53/scripts/trait2cell/FINAL_SCRIPTS/ATAC/ATAC_processing/01_remove_empty_drops.R

#STEP 1 - identify real cells
echo "Identifying real cells..."
module load HGI/softpack/users/cc53/R_emptydrops/2

Rscript $SCRIPT $DATA_DIR $OUTDIR

#STEP 2 - filter fragments to real cells
echo "Filtering fragments..."
module load HGI/softpack/groups/team152/genetics_tools/1

for kc_file in ${OUTDIR}/*_keepcells.csv; do
    sample=$(basename "$kc_file" "_keepcells.csv")
    echo "Processing ${sample}..."
    
    # Extract barcodes (skip header, get first column)
    barcode_file=$(mktemp)
    tail -n +2 "$kc_file" | cut -d',' -f1 | tr -d '"' > "$barcode_file"
    echo "  Barcodes to keep: $(wc -l < $barcode_file)"
    
    # Find matching fragments file (anchor to sample- to avoid partial matches)
    frag_file=$(ls ${DATA_DIR}/${sample}-*fragments.tsv.gz 2>/dev/null | head -1)
    
    if [ -z "$frag_file" ]; then
        echo "  WARNING: No fragments file found for ${sample} - skipping"
        rm "$barcode_file"
        continue
    fi
    echo "  Using fragments file: $(basename $frag_file)"
    
    out_frag="${OUTDIR}/${sample}-fragments_filt.tsv"

    zcat "$frag_file" | awk -v barcode_file="$barcode_file" '
        BEGIN { while((getline line < barcode_file) > 0) barcodes[line]=1 }
        /^#/ { print; next }
        ($4 in barcodes) { print }
    ' > "$out_frag"
    
    echo "  Kept $(grep -v "^#" $out_frag | wc -l) fragments"
    
    # Cleanup temp file
    rm "$barcode_file"
    
    # Bgzip and index
    echo "  Compressing new fragments..."
    bgzip "$out_frag"
    echo "  Indexing new fragments..."
    tabix -p bed "${out_frag}.gz"
    
    echo "Done: ${sample}"
done
echo "All done!"