#!/usr/bin/env Rscript

library(ArchR)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)
sampleid<- args[1]
fragments_file <- args[2]
outdir <- args[3]

setwd(outdir)
addArchRThreads(threads = 2)
addArchRGenome("hg38")

inputFiles <- setNames(fragments_file, sampleid)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 0,
  minFrags = 0,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE
)

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10,
    knnMethod = "UMAP",
    LSIMethod = 1,
    force = TRUE
)
