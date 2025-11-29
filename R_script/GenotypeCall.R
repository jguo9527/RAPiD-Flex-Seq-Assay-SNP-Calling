# GenotypeCall.R: Routine for calling biallelic genotypes based on read depths
# Open the script from the cloned directory

# FLAT_priors = 10 uniform genotype frequencies
# FLAT_SNP_priors = 3 uniform genotype frequencies based on the snp_code (e.g., C/T)

# HWE_priors = 10 genotype frequencies based on 4 uniform allele frequencies and HWE
# HWE_SNP_priors = 3 genotype frequencies based the snp_code and HWE

# FREQ_HWE_priors = 3 genotype frequencies based on pre-determined frequencies for 2 alleles
# EST_HWE_priors = 3 genotype frequencies based on frequencies for 2 alleles estimated from the input dataset


###################################################################################################
library(data.table)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(updog)

rm(list = ls(all = TRUE)) # Remove existing objects

###################################################################################################
## Set inputs to callGenotypes function
###################################################################################################
out_dir <- "../Output" # change / to \\ in windows os
dir.create(out_dir, recursive = T)

priors_list <- c("FLAT_priors", "FLAT_SNP_priors", "HWE_priors", "HWE_SNP_priors", "FREQ_HWE_priors", "EST_HWE_priors")
prior <- priors_list[4]

distr_list <- c("beta.binom", "binom") # Probability distribution
distr <- distr_list[1]
rho <- 0.33
  
BRE <- 0.005 # Expected read error (0.005 = 0.5% and 1e-06 = 1e-4%)

pErr.test <- 0.01 # Error value for binomial test of called genotype

# setwd(".")
source("FunctionPriors.R") #Function for calculating SNP priors
source("Function_callGenotypes.R") #Function for calling genotypes from read counts and input parameters

# Read two input files
snp_input <- readRDS("../Input/snp_input_100.rds") # Input file of SNP information
reads_input <- readRDS("../Input/reads_input_100.rds") # Input file of read counts by SNP and sample

snp_priors <- FunctionPriors(snp_input, reads_input) # Calculate genotype priors (6 options)
reads_priors <- left_join(reads_input, snp_priors) # Join reads and SNP priors

call_list <- apply(reads_priors, 1, function(x) callGenotypes(x, prior, distr, rho, BRE, pErr.test))
calls_out <- data.table(do.call(rbind, call_list))
calls_out
saveRDS(calls_out, file.path(out_dir, paste0("call_", prior, "_", distr, "_", rho, ".rds")))

