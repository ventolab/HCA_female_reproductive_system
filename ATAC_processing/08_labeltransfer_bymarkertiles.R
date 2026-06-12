#!/usr/bin/env Rscript

library(Seurat)
library(SingleCellExperiment)
library(schard)
library(dplyr)
library(ArchR)
library(Matrix)
library(GenomicRanges)
library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
cellannotation=args[1]
ATAC_ArchRproj=args[2]
ATAC_latent=args[3]
RNA_h5ad=args[4]
RNA_markers=args[5]



skip_to_integration=TRUE
if(skip_to_integration){
    ATAC_ArchRproj_subset=paste0(ATAC_ArchRproj, "_subset")

    setwd(ATAC_ArchRproj_subset)
    addArchRGenome("hg38")
    addArchRThreads(threads = 8)
    proj=loadArchRProject(path = ATAC_ArchRproj_subset)

    sce=h5ad2sce(RNA_h5ad)
    markergenes=read.csv(RNA_markers,header=T)$gene

} else {

setwd(ATAC_ArchRproj)

addArchRGenome("hg38")
addArchRThreads(threads = 8)

#--------------- Do scRNA-seq integration -----------------

proj=loadArchRProject(path = ATAC_ArchRproj)

sce=h5ad2sce(RNA_h5ad)
markergenes=read.csv(RNA_markers,header=T)$gene

#-------------- Add reduced dimensions-----------------

# add my own reduced dims
spectral <- read.csv(ATAC_latent, header = TRUE, row.names = 1)

# After loading spectral and archr_cells
archr_cells <- getCellNames(proj)
message(paste("ATAC cells:", length(archr_cells)))
message(paste("Spectral rows:", nrow(spectral)))
message(paste("Cells missing from spectral:", sum(!archr_cells %in% rownames(spectral))))

# Compare formats
message(paste("ArchR cell example:", head(archr_cells, 3)))
message(paste("Spectral rowname example:", head(rownames(spectral), 3)))

rownames(spectral)=gsub(":","#",rownames(spectral)) #should not be needed if i made the snapatac2 files correctly

#subset project to cells in spectral
#archr_cells <- getCellNames(proj)
cells_to_keep <- archr_cells[archr_cells %in% rownames(spectral)]
message(paste("Keeping", length(cells_to_keep), "of", length(archr_cells), "ATAC cells"))

ATAC_ArchRproj_subset=paste0(ATAC_ArchRproj, "_subset")
proj <- subsetArchRProject(ArchRProj = proj,cells = cells_to_keep,outputDirectory = ATAC_ArchRproj_subset,force = TRUE)

archr_cells <- getCellNames(proj)
spectral <- spectral[archr_cells, ]
colnames(spectral) <- paste0("X_", seq_len(ncol(spectral)))
proj@reducedDims[["SnapATAC2_Spectral"]] <- NULL
proj@reducedDims[["SnapATAC2_Spectral"]] <- SimpleList(matDR = as.matrix(spectral),params = list(),date = Sys.time(),scaleDims = FALSE,corToDepth = list(none = rep(0, ncol(spectral))))
proj <- saveArchRProject(ArchRProj = proj, outputDirectory = ATAC_ArchRproj_subset)

}
#----------------- Integrate with Seurat -----------------

# Filter to marker genes that exist in the ArchR gene score matrix
genes_in_proj <- getFeatures(proj, useMatrix = "GeneScoreMatrix")
message(paste("Checking for bad gene names:", genes_in_proj[1:10]))
markergenes_filtered <- intersect(markergenes, genes_in_proj)
message(sprintf("Marker genes: %d total, %d found in ArchR project", length(markergenes), length(markergenes_filtered)))
if (length(markergenes_filtered) == 0) {
    stop("No marker genes found in ArchR project!")}

#remove any low qc cells
sce <- sce[, !is.na(colData(sce)[[cellannotation]])]

#sanity check
sce_genes_after_filter <- rownames(sce)
markergenes_still_present <- markergenes_filtered %in% sce_genes_after_filter
if(!all(markergenes_still_present)) {
    message("Some marker genes lost after cell filtering!")
    message(paste(markergenes_filtered[!markergenes_still_present], collapse = ", "))
    markergenes_filtered <- markergenes_filtered[markergenes_still_present]}

#integrate
proj <- addGeneIntegrationMatrix(ArchRProj = proj, useMatrix = "GeneScoreMatrix",matrixName = "GeneIntegrationMatrix",reducedDims = "SnapATAC2_Spectral",
                                  seRNA = sce,addToArrow = TRUE,force = TRUE,groupRNA = cellannotation, genesUse=markergenes_filtered,
                                  nameCell = "predicted_cell_unconstrained",nameGroup = "predicted_annotation_unconstrained",nameScore = "predicted_score_unconstrained") 
                                  #no longer doing with iterativeLSI_MarkerTiles

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = ATAC_ArchRproj_subset)


#---------------- Saving gene files -------------------

gene_files_outdir=paste0(ATAC_ArchRproj_subset,"/genefiles/")
dir.create(gene_files_outdir, recursive = TRUE, showWarnings = FALSE)

message("Writing gene files...")
# Save metadata
metadata <- getCellColData(proj)
write.csv(metadata, paste0(gene_files_outdir,"cell_metadata.csv"),quote=F)

# Save gene score matrix
gs_matrix <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
writeMM(assay(gs_matrix), paste0(gene_files_outdir,"genescore_matrix.mtx"))
write.csv(rowData(gs_matrix)$name, paste0(gene_files_outdir,"genes.csv"), row.names = FALSE,quote=F)
write.csv(colnames(gs_matrix), paste0(gene_files_outdir,"cells.csv"), row.names = FALSE,quote=F)

# Save gene integration matrix
gi_matrix <- getMatrixFromProject(proj, useMatrix = "GeneIntegrationMatrix")
writeMM(assay(gi_matrix), paste0(gene_files_outdir,"geneintegration_matrix.mtx"))

# Save the gene names for GeneIntegration too 
write.csv(rowData(gi_matrix)$name, paste0(gene_files_outdir,"genes_geneintegration.csv"),row.names = FALSE, quote = FALSE)