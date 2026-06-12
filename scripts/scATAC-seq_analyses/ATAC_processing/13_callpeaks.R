library(Seurat)
library(SingleCellExperiment)
library(schard)
library(dplyr)
library(ArchR)
library(Matrix)
library(GenomicRanges)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

proj_path="/nfs/team292/projects/PanTissue/data/temp/ATAC/processed/nonfetal_all/ArchRproj_subset"

latent_fname='/nfs/team292/projects/PanTissue/data/temp/ATAC/processed/nonfetal_all/save_latent_obs/X_spectral_Harmony_500000features.csv'
obs_fname='/nfs/team292/projects/PanTissue/results/freeze/ATAC/annotations/alltissues.csv'

setwd(proj_path)
gene_files_outdir=paste0(proj_path,"/genefiles/")
dir.create(gene_files_outdir, recursive = TRUE, showWarnings = FALSE)

addArchRGenome("hg38")
addArchRThreads(threads = 8)

proj <- loadArchRProject(path = proj_path)

# Read the annotation file
#cell_annot <- read.csv(paste0(out_dir,"peakcalling_annot.csv"), row.names = 1,header=T, fill = TRUE)

# Filter ArchRProject to only keep cells in this annotation
#cells_to_keep <- rownames(cell_annot)
#proj <- proj[cells_to_keep, ]

# Add annotation to ArchRProject
#cell_order <- match(getCellNames(proj), rownames(cell_annot))
#proj$peakcalling_annot <- cell_annot[cell_order, 'peakcalling_annot']

skip=TRUE
if(skip){
    'Skipping to peak calling'
} else {
# add my own reduced dims
spectral <- read.csv(latent_fname, header = TRUE, row.names = 1)
#get obs 
meta <- read.csv(obs_fname, header = TRUE,row.names = 1)
meta=meta %>% filter(!final_annotation %in% c("doublet", "lowQC"))
cells_to_keep=rownames(meta)

archr_cells <- getCellNames(proj)
meta <- meta[archr_cells, ]
for(col in colnames(meta)){
    proj <- addCellColData(
        ArchRProj = proj,
        data = meta[[col]],
        name = col,
        cells = rownames(meta),
        force = TRUE  # overwrites existing columns
    )
}

spectral <- spectral[archr_cells, ]
colnames(spectral) <- paste0("X_", seq_len(ncol(spectral)))
proj@reducedDims[["SnapATAC2_Spectral"]] <- NULL
proj@reducedDims[["SnapATAC2_Spectral"]] <- SimpleList(matDR = as.matrix(spectral),params = list(),date = Sys.time(),scaleDims = FALSE,corToDepth = list(none = rep(0, ncol(spectral))))
proj <- saveArchRProject(ArchRProj = proj, outputDirectory = proj_path)
}


# Call peaks by annotation
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "final_annotation", force=FALSE)
proj <- saveArchRProject(ArchRProj = proj, outputDirectory = proj_path)

proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "final_annotation",
    pathToMacs2 = findMacs2(), force=TRUE
)
proj <- saveArchRProject(ArchRProj = proj, outputDirectory = proj_path)

proj <- addPeakMatrix(ArchRProj = proj, force=TRUE)
proj <- saveArchRProject(ArchRProj = proj, outputDirectory = proj_path)

# Get peak matrix
peak_matrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
rowData(peak_matrix)$peak=paste0(seqnames(peak_matrix@rowRanges),":",start(peak_matrix@rowRanges@ranges),"-",end(peak_matrix@rowRanges@ranges))
rownames(peak_matrix)=rowData(peak_matrix)$peak

#RDS
saveRDS(peak_matrix, paste0(gene_files_outdir,"peak_matrix.rds"))

# Create matrix with peak IDs as rownames
writeMM(assay(peak_matrix), paste0(gene_files_outdir, "peak_matrix.mtx"))

# 2. Save peak information (genomic coordinates)
peak_info <- rowData(peak_matrix)
write.csv(as.data.frame(peak_info), 
          paste0(gene_files_outdir, "peaks.csv"), 
          row.names = FALSE, 
          quote = FALSE)

# 3. Save cell names
write.csv(colnames(peak_matrix), 
          paste0(gene_files_outdir, "cells.csv"), 
          row.names = FALSE, 
          quote = FALSE)

# 4. Save metadata (if needed)
metadata <- getCellColData(proj)
write.csv(metadata, 
          paste0(gene_files_outdir, "cell_metadata.csv"), 
          quote = FALSE)

#also save gene integration matrix so I can use it
gi_matrix <- getMatrixFromProject(proj, useMatrix = "GeneIntegrationMatrix")
saveRDS(gi_matrix, paste0(gene_files_outdir,"geneintegration_matrix.rds"))