library(ArchR)
library(SingleCellExperiment)
library(dplyr)

addArchRGenome("hg38")
addArchRThreads(threads = 8)

proj_path="/nfs/team292/projects/PanTissue/data/temp/ATAC/processed/nonfetal_all/ArchRproj_subset"

setwd(proj_path)

proj <- loadArchRProject(proj_path)

# ── 1. Downsample to max 1000 cells per final_annotation ──────────────────────
cell_meta <- getCellColData(proj, select = "final_annotation") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cellName")

set.seed(42)
cells_keep <- cell_meta %>%
  group_by(final_annotation) %>%
  group_modify(~ slice_sample(.x, n = min(1000, nrow(.x)))) %>%
  ungroup() %>%
  pull(cellName)

n_cells_ds    <- length(cells_keep)
knn_iter      <- max(50, round(n_cells_ds / 200 * 0.5))  # ~50 % coverage, min 50
message("Downsampled N = ", n_cells_ds, " | knnIteration = ", knn_iter)

proj <- addPeak2GeneLinks(ArchRProj = proj,reducedDims = "SnapATAC2_Spectral",
                            maxDist = 500000, #max distance between peak and gene to consider for correlation (default = 250000),
                            cellsToUse = cells_keep,
                            k = 100, #number of cells in a neighbourhood (default = 100)
                            knnIteration = knn_iter, #number of neighbourhoods (default = 500) #so aiming for 110,000 cells
                            overlapCutoff = 0.5) #proportion that two neighbourhoods can overlap (default = 0.8)
proj <- saveArchRProject(ArchRProj = proj, outputDirectory = proj_path)

p2g <- getPeak2GeneLinks(ArchRProj = proj,corCutOff = 0.3,resolution = 1,returnLoops = FALSE)

p2g$geneName <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
p2g$peakName <- (metadata(p2g)$peakSet %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2g$idxATAC]

gene_files_outdir=paste0(proj_path,"/genefiles/")
#dir.create(gene_files_outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(p2g, paste0(gene_files_outdir,"peak2gene_links.csv"),row.names = FALSE, quote = FALSE)