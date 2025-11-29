callGenotypes <- function (x, prior, distr, rho, BRE, pErr.test) {
  # Initialize variables
  sample <- x$sample; snp_name <- x$snp_name; snp_code <- x$snp_code
  A <- x$A; C <- x$C; G <- x$G; T <- x$T
  prior_name <- call_letter <- pG <- pCall <- PO <- logPO <- FPR <- miss_read <- pErr <- pHet <- FPR <- NA
  
  D <- c(A, C, G, T) # Read counts
  names(D) <- c("A", "C", "G", "T")
  if(distr == "binom") {rho <- NA}
  
  ###################################################################################################
  ## Calculate pR_G_matrix: Probabilities of obtaining the indicated read (R) given the genotype
  ###################################################################################################
  pR_ho <-matrix(BRE/3,4,4)  # ho = homozygotes
  diag(pR_ho) <- (1-BRE)
  pR_he <- matrix(BRE/3, nrow=4, ncol=6) # he = heterozygotes
  pR_he <- mapply(function(x,y) replace(x, y, 0.5*(1-BRE)+BRE/6),
                  unname(split(pR_he,rep(1:ncol(pR_he)))),combn(1:4, 2, simplify=F))
  pR_G_matrix <- cbind(pR_ho, pR_he)
  names(pR_G_matrix) <- c("AA","CC","GG","TT","AC","AG","AT","CG","CT","GT")
   
  if(sum(D) > 0 ) {
    prior_name <- ifelse(!is.na(snp_code), prior, "FLAT_priors") # Use FLAT_priors if snp_code = NA
    prior_col <- which(names(x) %in% paste0(prior_name))
    pG <- unlist(x[[prior_col]])[5:14]
    
    #################################################################################################
    ## Use probabilities based on binomial (binom) or beta-binomial distribution (beta,binom)
    #################################################################################################
    # Initialize variables
    pR_G <- pd_G <- pD_G <- pGxpD_G <- pGxpD_G_sum <- pG_D <- logPO_vec <- FPR_vec <- NA
    pR_G <- pR_G_matrix[, 1:10] # Prob of obtaining the read (R) given the diploid genotype
    
    if(distr == "binom") {
      pd_G <- dbinom(D, sum(D), pR_G, log = TRUE) # Probability of d(k)|G (ln)
    }
    if(distr == "beta.binom") {
      pd_G <- rbind(
        updog::dbetabinom(rep(D[1], ncol(pR_G)), sum(D), mu = pR_G[1,], rho = rho, log = TRUE),
        updog::dbetabinom(rep(D[2], ncol(pR_G)), sum(D), mu = pR_G[2,], rho = rho, log = TRUE),
        updog::dbetabinom(rep(D[3], ncol(pR_G)), sum(D), mu = pR_G[3,], rho = rho, log = TRUE),
        updog::dbetabinom(rep(D[4], ncol(pR_G)), sum(D), mu = pR_G[4,], rho = rho, log = TRUE)
      )
    }
    #################################################################################################
    ## Calculate other statistics
    #################################################################################################
    pD_G <- colSums(pd_G) # P(D|Gj) = MULT [P(dk|Gj)] (ln)
  
    pGxpD_G <- pD_G + log(pG) # P(Gj)P(D|Gj) (ln)
    pGxpD_G_sum <- ifelse(log(sum(exp(pGxpD_G))) == -Inf,
                          max(pGxpD_G),
                          log(sum(exp(pGxpD_G))))
    pG_D <- pGxpD_G - pGxpD_G_sum # P(Gi|D) (ln)
    
    # Likelihood ratio (ln transformed). Handles the case where the top two calls have equal max values
    logPO_vec <- sapply(seq_along(pG_D), function(x) {
      pG_D[x] - max(pG_D[-x])
      })
    FPR_vec <- sum(sapply(pG_D, exp)) - exp(pG_D) # False positive risk
    FPR <- unname(FPR_vec[which.max(logPO_vec)]) # Find the minimum FPR
    logPO <- max(logPO_vec) # Find the maximum lnPO and PO
    PO <- exp(logPO)
    pCall <- unname(exp(pG_D[which.max(logPO_vec)])) # Calculate the multinomial prob of the genotype
    
    #################################################################################################
    ## Find the genotype with the highest ln = called genotype
    #################################################################################################
      call_letter <- names(which.max(logPO_vec))
    
    # Calculate error statistics
    # Calculate number of reads for two called alleles and combined error alleles
    reads_geno <- D[unique(unlist(strsplit(call_letter, split = "")))]
    reads_error <- sum(D) - sum(reads_geno)
    miss_read <- reads_error/sum(D) * 100
    
    # Probability miss_read < pErr.test
    pErr <- ifelse(sum(D) > 0,
                   binom.test(reads_error, sum(D), pErr.test, alternative="greater")$p.value,
                   NA)
    
    # Probability that there's an equal distribution of reads (50/50) between two alleles
    D.sort <- sort(D)
    pHet <- as.numeric(binom.test(D.sort[3], D.sort[3]+D.sort[4], 0.50, alternative="less")$p.value)
  }
  df <- data.frame(snp_name, sample, prior_name, distr, rho, BRE, pErr.test,
                   A, C, G, T, call = call_letter, pCall,
                   PO, logPO, miss_read, pErr, pHet, FPR)
  return(df)
}

