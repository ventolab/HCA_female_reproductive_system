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
#cellannotation=args[1]
proj_path="/nfs/team292/projects/PanTissue/data/temp/ATAC/processed/nonfetal_all_copy/ArchRproj"
latent_fname='/nfs/team292/projects/PanTissue/data/temp/ATAC/processed/nonfetal_all_copy/save_latent_obs/X_spectral_Harmony_500000features.csv'
obs_fname='/nfs/team292/projects/PanTissue/results/freeze/ATAC/annotations/alltissues.csv'

tissue_proj_paths=c('/nfs/team292/projects/PanTissue/data/temp/ATAC/processed/nonfetal_fallopiantube/ArchRproj_subset',
                        '/nfs/team292/projects/PanTissue/data/temp/ATAC/processed/nonfetal_uterus/ArchRproj_subset',
                        '/nfs/team292/projects/PanTissue/data/temp/ATAC/processed/nonfetal_ovary/ArchRproj_subset')

setwd(proj_path)
geneAnnotation_fname='/nfs/team292/projects/PanTissue/data/temp/GWASxRNA/RNA_references/10x_ref_geneAnnotation_hg38_2020A.rds'

addArchRGenome("hg38")
addArchRThreads(threads = 8)

exclude_samples=c()
if(length(exclude_samples)>0){
    message(paste0("Excluding samples: ",exclude_samples))}


#--------------- Make ArchR project -----------------

# collect arrows from all project ArrowFiles subdirectories
ArrowFiles <- unlist(lapply(tissue_proj_paths, function(p) {
  list.files(file.path(p, 'ArrowFiles'), pattern = "\\.arrow$", 
             recursive = TRUE, full.names = TRUE)
}))
message("Checking arrow files:")
message(paste(ArrowFiles))

if(length(exclude_samples)>0){
    message("Excluding some samples...")
    pattern <- paste(exclude_samples, collapse = "|")
    ArrowFiles <- ArrowFiles[!grepl(pattern, basename(ArrowFiles))]}

#make archr project
proj <- ArchRProject(ArrowFiles = ArrowFiles,outputDirectory = proj_path,copyArrows = FALSE) #, geneAnnotation=geneAnnotation

#--------------- Get reduced dims -----------------

# add my own reduced dims
spectral <- read.csv(latent_fname, header = TRUE, row.names = 1)
#get obs 
meta <- read.csv(obs_fname, header = TRUE,row.names = 1)
meta=meta %>% filter(!final_annotation %in% c("doublet", "lowQC"))
cells_to_keep=rownames(meta)

archr_cells <- getCellNames(proj)
message(paste("ATAC cells:", length(archr_cells)))
message(paste("Obs cells:", length(cells_to_keep)))
message(paste("Spectral rows:", nrow(spectral)))
#message(paste("Cells missing from spectral:", sum(!archr_cells %in% rownames(spectral))))

proj_path_subset=paste0(proj_path, "_subset")
proj <- subsetArchRProject(ArchRProj = proj,cells = cells_to_keep,outputDirectory = proj_path_subset,force = TRUE)

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
proj <- saveArchRProject(ArchRProj = proj, outputDirectory = proj_path_subset)

