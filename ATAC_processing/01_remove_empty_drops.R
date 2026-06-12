library(Seurat)
library(EmptyDropsMultiome)

#peak_matrix_dir="/nfs/team292/projects/PanTissue/data/temp/ATAC/ovary_sanger/raw_peak_matrices_cellranger/"
#outdir="/nfs/team292/projects/PanTissue/data/temp/ATAC/ovary_sanger/EmptyDropsMultiome/"

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
outdir <- args[2]

peak_matrices <- list.files(data_dir, pattern = "matrix\\.h5", full.names = TRUE)
print(peak_matrices)

# Extract donor ID and create named vector
inputFiles <- setNames(peak_matrices, sub("-[^-]*$", "", basename(peak_matrices))) #extract text before the last dash
print(inputFiles)

for (i in seq_along(inputFiles)) {
    file_path <- inputFiles[i]        # full file path
    sample <- names(inputFiles)[i]  # extracted name
    message(paste0("Processing ",sample,"..."))
  
    sce <- Read10X_h5(file_path)
    rna = atac = sce
    
    eD.out <- emptydrops_multiome(count_matrix_rna=rna, count_matrix_atac=atac)
    realcells=eD.out[eD.out$FDR_chromatin<0.001 & ! is.na(eD.out$FDR_chromatin),]
    message(paste0("The number of ATAC cells detected is: ", nrow(realcells)))

    rownames(realcells) <- realcells$Row.names
    realcells$Row.names <- NULL
    write.csv(realcells,paste0(outdir,sample,"_keepcells.csv"),quote=F) #save as a checkpoint
}
message("Done!")