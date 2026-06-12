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
proj_path=args[1]
arrows_dir=args[2]

#setwd("/lustre/scratch125/cellgen/vento/cc53/trait2cell/ovary/scATAC/")
setwd(proj_path)

#proj_path=paste0("/nfs/team292/projects/PanTissue/data/temp/ATAC/ovary_sanger/processed/archRproj_crosslifespan_emptydrops")
geneAnnotation_fname='/nfs/team292/projects/PanTissue/data/temp/GWASxRNA/RNA_references/10x_ref_geneAnnotation_hg38_2020A.rds'
#adata_fname="/lustre/scratch125/cellgen/vento/cc53/trait2cell/ovary/scATAC/data/scRNAseq_forArchR_integration/ovaries_pediatric_rawcounts_allgenes.h5ad"
#markergenes_fname=paste0('/lustre/scratch125/cellgen/vento/cc53/trait2cell/ovary/scATAC/data/scRNAseq_forArchR_integration/ovaries_top50markers_',cellannotation,'.csv')
#out_dir="/lustre/scratch125/cellgen/vento/cc53/trait2cell/ovary/scATAC/data/archRproj_allpedsamples/" 

addArchRGenome("hg38")
addArchRThreads(threads = 8)

exclude_samples=c() #"AUrv22","AUrv28","Hrv58multi","Hrv58multiBis","Hrv65multi","Hrv65multiBis"
if(length(exclude_samples)>0){
    message(paste0("Excluding samples: ",exclude_samples))}


#--------------- Make ArchR project -----------------

ArrowFiles <- list.files(arrows_dir, pattern = "\\.arrow$", recursive = TRUE, full.names = TRUE)

message("Checking arrow files:")
message(paste(ArrowFiles))

if(length(exclude_samples)>0){
    message("Excluding some samples...")
    pattern <- paste(exclude_samples, collapse = "|")
    ArrowFiles <- ArrowFiles[!grepl(pattern, basename(ArrowFiles))]}

geneAnnotation=readRDS(geneAnnotation_fname)
message(paste0("Number of genes in geneAnnotation: ",length(geneAnnotation$genes))) # 36559 -> 22073

#make archr project
proj <- ArchRProject(ArrowFiles = ArrowFiles,outputDirectory = proj_path,copyArrows = FALSE) #, geneAnnotation=geneAnnotation
  

#--------------- Add gene score matrix -----------------

proj <- addGeneScoreMatrix(proj,genes = geneAnnotation$genes,matrixName = "GeneScoreMatrix",force=TRUE)
proj <- saveArchRProject(ArchRProj = proj, outputDirectory = proj_path)

#--------------- Get reduced dims -----------------

#proj <- addIterativeLSI(ArchRProj = proj,useMatrix = "TileMatrix",name = "IterativeLSI",iterations = 2,clusterParams = list(resolution = c(0.2),sampleCells = 10000,n.start = 10),varFeatures = 25000,force = TRUE)
#proj <- addHarmony(ArchRProj = proj,reducedDims = "IterativeLSI",name = "Harmony",groupBy = "Sample",force = TRUE)
#proj <- addUMAP(ArchRProj = proj,reducedDims = "IterativeLSI",name = "UMAP_IterativeLSI",force = TRUE)
#proj <- addUMAP(ArchRProj = proj,reducedDims = "Harmony",name = "UMAP_Harmony",force = TRUE)

#proj <- saveArchRProject(ArchRProj = proj, outputDirectory = proj_path)

#--------------- Cluster and flag doublets -----------------

#cluster
#for(reddim in c("IterativeLSI","Harmony")){
    #proj <- addClusters(input = proj,reducedDims = reddim,method = "Seurat",name = paste0("Clusters_",reddim),resolution = 1,force = TRUE)}

#flag doublets
doublet_flag <- rep(FALSE, nCells(proj))
for (s in unique(proj$Sample)) {
  cells <- which(proj$Sample == s)
  cutoff <- quantile(proj$DoubletEnrichment[cells],probs = 0.95,na.rm = TRUE)
  doublet_flag[cells] <- proj$DoubletEnrichment[cells] > cutoff} #remove top 5% of cells per sample
proj$is_doublet <- doublet_flag

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = proj_path)

#plotting donors
#plot_list <- list(
        #plotEmbedding(proj, colorBy = "cellColData", name = "Sample",embedding = 'UMAP_Harmony'),
        #plotEmbedding(proj, colorBy = "cellColData", name = "TSSEnrichment",embedding = 'UMAP_Harmony'),
        #plotEmbedding(proj, colorBy = "cellColData", name = "nFrags",embedding = 'UMAP_Harmony'),
        #plotEmbedding(proj, colorBy = "cellColData", name = "is_doublet",embedding ='UMAP_Harmony'),

        #plotEmbedding(proj, colorBy = "cellColData", name = "Sample",embedding = 'UMAP_IterativeLSI'),
        #plotEmbedding(proj, colorBy = "cellColData", name = "TSSEnrichment",embedding = 'UMAP_IterativeLSI'),
        #plotEmbedding(proj, colorBy = "cellColData", name = "nFrags",embedding = 'UMAP_IterativeLSI'),
        #plotEmbedding(proj, colorBy = "cellColData", name = "is_doublet",embedding ='UMAP_IterativeLSI'))
#plotPDF(plotList = plot_list, name = paste0("dimred_umaps.pdf"), ArchRProj = proj, addDOC = FALSE,width = 8, height = 8)
