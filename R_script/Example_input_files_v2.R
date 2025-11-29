# Example_input_files_v2
# Make example input files for genotype calling = reads and snps
library (data.table)
library (dplyr)
library(stringr)
library(tibble)
library(tidyr)

rm(list = ls(all = TRUE))

snp_num <- 100
samp_num <- 100

out_dir <- "../Input" # change / to \\ in windows os
dir.create(out_dir, recursive = T)

# Process the metaData file with SNP information
meta_dir <- "../MetaData" # change / to \\ in windows os
meta_file <- "metaData_2023-11-7.rds"

meta_snp <- readRDS(paste0(meta_dir, meta_file))[[1]] #SNP metadata
meta_snp <- dplyr::select(meta_snp, colnames(meta_snp)[!grepl("41602|41620", colnames(meta_snp))])
colnames(meta_snp)
meta_snp <- dplyr::select(meta_snp, rg.snp.name, fin.snp.slash, axiom.fA, axiom.fC, axiom.fG, axiom.fT)
colnames(meta_snp) <- c("snp_name", "snp_code", "fA", "fC", "fG", "fT")
snp_example <- dplyr::arrange(meta_snp, snp_name)[1:snp_num]
saveRDS(snp_example, paste0(out_dir, "snp_input_", snp_num, ".rds"))
saveRDS(meta_snp, paste0(out_dir, "snp_input_all.rds"))

# Process reads file
read.dir <- "../PrepReads" # change / to \\ in windows os
read.file <- "OP_BWA_CC_OP_reads_v3.rds"
reads <- readRDS(paste0(read.dir, read.file))

reads_out <- dplyr::rename(reads, "snp_name" = "rg.snp.name",
                               "snp_code" = "fin.snp",
                               "A" = "BWA.A",
                               "C" = "BWA.C",
                               "G" = "BWA.G",
                               "T" = "BWA.T")
reads_out <- dplyr::select(reads_out, snp_name, sample, snp_code, A, C, G, T)
reads_example <- dplyr::filter(reads_out, snp_name %in% snp_example$snp_name)
reads_example <- dplyr::filter(reads_example, sample %in% sort(unique(reads_example$sample))[1:samp_num])
saveRDS(reads_example, paste0(out_dir, "reads_input_", samp_num, ".rds"))
saveRDS(reads_out, paste0(out_dir, "reads_input_all.rds"))

