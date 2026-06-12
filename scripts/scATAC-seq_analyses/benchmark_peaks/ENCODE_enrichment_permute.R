library(BiocManager)
library(regioneR)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

peak_type <- args[1]
category <- args[2]
if (length(args)>=3 && args[3]=='intergenic'){
    intergenic=TRUE
    plotcomment="_intergenic"
} else {
    intergenic=FALSE
    plotcomment=""
}

blacklist=read.table("/lustre/scratch125/cellgen/vento/cc53/utils/hg38-blacklist.v2.bed",sep='\t')
geneAnnotation=readRDS("/nfs/team292/projects/PanTissue/data/temp/GWASxRNA/RNA_references/10x_ref_geneAnnotation_hg38_2020A.rds")
encode=read.table('/nfs/team292/projects/PanTissue/data/freeze/ATAC/resources/GRCh38-cCREs.bed')
encode=encode[encode[[6]] == category, ]

if (peak_type=='allpeaks'){
    peaks=read.csv('/nfs/team292/projects/PanTissue/data/temp/ATAC/processed/nonfetal_all/ArchRproj_subset/genefiles/peaks.csv')
    peaks=peaks %>%
        separate(peak, into = c("chr", "start", "end"), sep = "[:-]") %>%
        mutate(start = as.integer(start), end = as.integer(end)) %>%
        select(chr, start, end)
}

if (peak_type=='linkedpeaks'){
    peaks=read.csv('/nfs/team292/projects/PanTissue/data/temp/ATAC/processed/nonfetal_all/ArchRproj_subset/genefiles/peak2gene_links.csv')
    peaks=peaks %>%
        filter(FDR<1e-5) %>%
	    distinct(peakName) %>%
        separate(peakName, into = c("chr", "start", "end"), sep = "[:-]") %>%
        mutate(start = as.integer(start), end = as.integer(end)) %>%
        select(chr, start, end)
}

# convert dfs to GRanges
query_gr       <- toGRanges(peaks,       format="BED")
blacklist_gr   <- toGRanges(blacklist,   format="BED")
encode_gr  <- toGRanges(encode, format="BED")

if (intergenic){
    #remove intergenic peaks
    gene_gr <- geneAnnotation$gene
    hits_query <- findOverlaps(query_gr, gene_gr)
    query_gr <- query_gr[-unique(queryHits(hits_query)), ]
}

# build genome mask
genome_gr <- getGenomeAndMask(genome="hg38", mask=NULL)$genome
genome_gr <- subtractRegions(genome_gr, blacklist_gr)

# run permutation test per category per peak set
pt <- overlapPermTest(
            A                  = query_gr,
            B                  = encode_gr,
            genome             = genome_gr,
            ntimes             = 1000,
            count.once         = TRUE,
            mc.cores           = 8,
            verbose            = FALSE)
        
obs      <- pt$numOverlaps$observed
exp_mean <- mean(pt$numOverlaps$permuted)
exp_sd   <- sd(pt$numOverlaps$permuted)
zscore   <- pt$numOverlaps$zscore
pval     <- pt$numOverlaps$pval
        
# fold enrichment and SE
enrichment <- obs / exp_mean
se <- exp_sd / exp_mean  # propagated SE on fold enrichment
        
res <- data.frame(
            category   = category,
            peak_type  = peak_type,
            observed   = obs,
            expected   = exp_mean,
            exp_sd     = exp_sd,
            enrichment = enrichment,
            se         = se,
            zscore     = zscore,
            pval       = pval)

write.csv(res,paste0("/nfs/team292/projects/PanTissue/results/freeze/ATAC/benchmark/ENCODE_enrichment_",category,"_",peak_type,plotcomment,".csv"))
saveRDS(pt, paste0("/nfs/team292/projects/PanTissue/results/freeze/ATAC/benchmark/permtest_", category, "_", peak_type, plotcomment,".rds"))

# forest plot
#ggplot(results_df, aes(x=enrichment, y=category, color=peak_type)) +
    #geom_vline(xintercept=1, linetype="dashed", color="grey50") +
    #geom_errorbarh(aes(xmin=enrichment-1.96*se, xmax=enrichment+1.96*se),
    #               height=0.2, position=position_dodge(width=0.5)) +
    #geom_point(size=3, position=position_dodge(width=0.5)) +
    #scale_color_manual(values=c("all"="#3a86ff", "linked"="#ff006e"),
    #                   labels=c("all"="All peaks", "linked"="Linked peaks")) +
    #labs(x="Fold enrichment over permuted background",
    #     y="ENCODE annotation category",
    #     color="Peak set") +
    #theme_classic() +
    #theme(legend.position="bottom")