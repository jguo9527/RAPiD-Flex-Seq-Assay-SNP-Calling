FunctionPriors <- function(snp, reads) {
  ###################################################################################################
  ## HWE: Function to calculate genotype frequencies based on HWE
  ###################################################################################################
  priorsHWE <- function(alleles) {
    mat <- alleles%*%t(alleles)
    mat <- mat/sum(mat)
    priors <- c(alleles, diag(mat), 2*mat[lower.tri(mat)])
    names(priors) <- c("A","C","G","T", "AA","CC","GG","TT","AC","AG","AT","CG","CT","GT")
    return(priors)
  }
  
  ###################################################################################################
  ## priorSNP: Function for genotype frequencies based on HWE
  ###################################################################################################
  priorsSNP <- function (snp_code, priors) {
    snp_letters <- c(substr(snp_code, 1, 1), substr(snp_code, 3, 3))
    grid <- expand.grid(letter1 = snp_letters, letter2 = snp_letters)
    snp_names <- c(snp_letters, paste0(grid$letter1, grid$letter2))
    priors[names(priors)[!(names(priors) %in% snp_names)]] <- 0
    priors[1:4] <- priors[1:4]/sum(priors[1:4])
    priors[5:14] <- priors[5:14]/sum(priors[5:14])
    return(priors)
  }
  
  ###################################################################################################
  ## estimateFreq: Function to estimate allele frequencies from the input reads
  ###################################################################################################
  estimateFreq <- function(reads) {
    alleleCount <- function(x) { #== Function to identify two probable alleles for each sample=======
      top_two <- sort(as.numeric(x[4:7]), decreasing = TRUE)[1:2]
      suppressWarnings(x$count_A <- ifelse(as.numeric(x[4]) %in% top_two & as.numeric(x[4]) != 0, 1, 0))
        x$count_C <- ifelse(as.numeric(x[5]) %in% top_two & as.numeric(x[5]) != 0, 1, 0)
        x$count_G <- ifelse(as.numeric(x[6]) %in% top_two & as.numeric(x[6]) != 0, 1, 0)
        x$count_T <- ifelse(as.numeric(x[7]) %in% top_two & as.numeric(x[7]) != 0, 1, 0)
      return(x)
    }
    alleleFreq <- function(x) { #== Function to estimate 2 allele frequencies for each SNP ==========
      # 1 = probable biallelic allele; 0 = unlikely biallelic allele
      top_two <- sort(as.numeric(x[3:6]), decreasing = TRUE)[1:2]
      suppressWarnings(x$eA <- ifelse(as.numeric(x[3]) %in% top_two & as.numeric(x[3]) != 0, x[3], 0))
      x$eC <- ifelse(as.numeric(x[4]) %in% top_two & as.numeric(x[4]) != 0, x[4], 0)
      x$eG <- ifelse(as.numeric(x[5]) %in% top_two & as.numeric(x[5]) != 0, x[5], 0)
      x$eT <- ifelse(as.numeric(x[6]) %in% top_two & as.numeric(x[6]) != 0, x[6], 0)
      x$f2_sum <- as.numeric(sum(as.numeric(x$eA), as.numeric(x$eC), as.numeric(x$eG), as.numeric(x$eT), na.rm = T))
      
      # Calculate allele  frequence by SNP
      x$eA <- as.numeric(x$eA)/as.numeric(x$f2_sum)
      x$eC <- as.numeric(x$eC)/as.numeric(x$f2_sum)
      x$eG <- as.numeric(x$eG)/as.numeric(x$f2_sum)
      x$eT <- as.numeric(x$eT)/as.numeric(x$f2_sum)
      x$eSum <- as.numeric(sum(as.numeric(x$eA), as.numeric(x$eC), as.numeric(x$eG), as.numeric(x$eT), na.rm = T))
      return(x)
    }
    
    # Use the alleleCount function to identify probable biallelic alleles
    alleles_out <- data.table(do.call(rbind, apply(reads, 1, function(x) alleleCount(x))))
    # Convert list variables to numeric
    alleles_out <- data.table(alleles_out[,c(1,3)] %>% dplyr::mutate_if(is.list, as.character),
                              alleles_out[,4:ncol(alleles_out)] %>% dplyr::mutate_if(is.list, as.numeric))
    # print(alleles_out) # Check results
    
    # Obtain allelle frequencies by SNP. f1 indicates that these values include counts for non-likelyy alleles
    # Non-likely alleles are set to zero in the alleleFreq function
    freq1 <- alleles_out %>%
      group_by(snp_name, snp_code) %>%
      dplyr::summarise(
        f1A = sum(count_A, na.rm=TRUE),
        f1C = sum(count_C, na.rm=TRUE),
        f1G = sum(count_G, na.rm=TRUE),
        f1T = sum(count_T, na.rm=TRUE)
      ) %>%
      mutate(
        f1_sum = f1A + f1C + f1G + f1T,
        f1A = f1A / f1_sum,
        f1C = f1C / f1_sum,
        f1G = f1G / f1_sum,
        f1T = f1T / f1_sum
      )
    # print(freq1) # Check results
    
    # Use the alleleFreq function to remove unlikely biallelic alleles and calculate biallelic frequencies
    freq_out <- data.table(do.call(rbind, apply(freq1, 1, function(x) alleleFreq(x))))
    freq_out <- data.table(freq_out[,c(1,2)] %>% dplyr::mutate_if(is.list, as.character),
                           freq_out[,3:ncol(freq_out)] %>% dplyr::mutate_if(is.list, as.numeric))
    # print(freq_out) # Check results
    
    freq_est <- dplyr::select(freq_out, snp_name, snp_code, eA, eC, eG, eT)
    return(freq_est)
  }
  
  FLAT_priors <- c(rep(0.25, 4), rep(0.1, 10))
  names(FLAT_priors) <- c("A","C","G","T", "AA","CC","GG","TT","AC","AG","AT","CG","CT","GT")
  snp$FLAT_priors <- list(FLAT_priors)
  
  HWE_priors <- priorsHWE(FLAT_priors[1:4])
  snp$HWE_priors <- list(HWE_priors)
  
  FLAT_SNP_priors <- apply(snp, 1, function(x) list(priorsSNP(x[2], FLAT_priors)))
  snp <- data.table(snp, FLAT_SNP_priors = unlist(FLAT_SNP_priors, recursive = F))
  
  HWE_SNP_priors <- apply(snp, 1, function(x) list(priorsSNP(x[2], HWE_priors)))
  snp <- data.table(snp, HWE_SNP_priors = unlist(HWE_SNP_priors, recursive = F))
  
  FREQ_HWE_priors <- apply(snp, 1, function(x) list(priorsHWE(as.numeric(x[3:6]))))
  snp <- data.table(snp, FREQ_HWE_priors = unlist(FREQ_HWE_priors, recursive = F))
  
  freq_est <- estimateFreq(reads)
  snp <- left_join(snp, freq_est, by = c("snp_name", "snp_code"))
  
  EST_HWE_priors <- apply(snp, 1, function(x) list(priorsHWE(as.numeric(x[12:15]))))
  snp <- data.table(snp, EST_HWE_priors = unlist(EST_HWE_priors, recursive = F))
  
  snp <- dplyr::select(snp, snp_name, snp_code, fA, fC, fG, fT, eA, eC, eG, eT, everything())
  
  return(snp)
}
