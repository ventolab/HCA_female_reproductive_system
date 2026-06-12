library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
input_csv      <- args[1]
output_bed_snp <- args[2]  # for index SNP positions
output_bed_int <- args[3]  # for intervals

df <- read_csv(input_csv)

# ── Index SNP positions (hg19) ────────────────────────────────────────────────
bed_snp <- df %>%
    filter(!is.na(Position_hg19)) %>%
    mutate(
        chr   = paste0("chr", recode(as.character(chromosome), "23" = "X")),
        start = as.integer(Position_hg19) - 1L,
        end   = as.integer(Position_hg19),
        name  = variantID
    ) %>%
    select(chr, start, end, name)

write_tsv(bed_snp, output_bed_snp, col_names = FALSE)
cat(sprintf("Wrote %d index SNPs to %s\n", nrow(bed_snp), output_bed_snp))

# ── Intervals (hg19) ─────────────────────────────────────────────────────────
bed_int <- df %>%
    filter(!is.na(Interval_start_bp_hg19) & !is.na(Interval_end_bp_hg19)) %>%
    mutate(
        chr   = paste0("chr", recode(as.character(chromosome), "23" = "X")),
        start = as.integer(Interval_start_bp_hg19) - 1L,
        end   = as.integer(Interval_end_bp_hg19),
        name  = variantID
    ) %>%
    select(chr, start, end, name)

write_tsv(bed_int, output_bed_int, col_names = FALSE)
cat(sprintf("Wrote %d intervals to %s\n", nrow(bed_int), output_bed_int))