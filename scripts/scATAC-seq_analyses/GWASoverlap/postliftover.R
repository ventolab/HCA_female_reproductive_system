library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
input_csv        <- args[1]
hg38_snp_bed     <- args[2]
hg38_int_bed     <- args[3]
output_csv       <- args[4]

df <- read_csv(input_csv)

hg38_snp <- read_tsv(hg38_snp_bed, col_names = c("chr", "start0", "end", "variantID")) %>%
    mutate(Position_hg38_lifted = end) %>%  # end = 1-based position for SNP
    select(variantID, Position_hg38_lifted)

hg38_int <- read_tsv(hg38_int_bed, col_names = c("chr", "start0", "end", "variantID")) %>%
    mutate(
        Interval_start_bp_hg38 = start0 + 1L,
        Interval_end_bp_hg38   = end
    ) %>%
    select(variantID, Interval_start_bp_hg38, Interval_end_bp_hg38)

df_out <- df %>%
    left_join(hg38_snp, by = "variantID") %>%
    left_join(hg38_int, by = "variantID")

cat(sprintf("SNPs lifted: %d  unmapped: %d\n",
    sum(!is.na(df_out$Position_hg38_lifted)),
    sum( is.na(df_out$Position_hg38_lifted))))
cat(sprintf("Intervals lifted: %d  unmapped: %d\n",
    sum(!is.na(df_out$Interval_start_bp_hg38)),
    sum( is.na(df_out$Interval_start_bp_hg38))))

write_csv(df_out, output_csv)