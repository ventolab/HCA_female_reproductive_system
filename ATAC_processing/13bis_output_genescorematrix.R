library(Seurat)
library(SingleCellExperiment)
library(schard)
library(dplyr)
library(ArchR)
library(Matrix)
library(GenomicRanges)

proj_path="/nfs/team292/projects/PanTissue/data/temp/ATAC/processed/nonfetal_all/ArchRproj_subset"

setwd(proj_path)
gene_files_outdir=paste0(proj_path,"/genefiles/")
dir.create(gene_files_outdir, recursive = TRUE, showWarnings = FALSE)

addArchRGenome("hg38")
addArchRThreads(threads = 8)

proj <- loadArchRProject(path = proj_path)

#---------------- Saving gene files -------------------

gene_files_outdir=paste0(proj_path,"/genefiles/")
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