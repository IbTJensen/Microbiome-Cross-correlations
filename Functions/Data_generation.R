Sequencing_sim <- function(community.samples, n.seq, const.reads = F){
  n <- nrow(community.samples)
  k <- ncol(community.samples)
  
  A <- matrix(NA, n, k)
  colnames(A) <- colnames(community.samples)
  rownames(A) <- rownames(community.samples)
  if(const.reads == T){
    Sequences <- rep(n.seq, n)
  }
  if(const.reads == F){
    Sequences <- round(rlnorm(n, meanlog = log(n.seq)-0.5^2, sdlog = 0.5))
  }
  for(i in 1:nrow(A)){
    A[i,] <- rmultinom(1, Sequences[i], prob = community.samples[i,])
  }
  return(A)
}

Abs_abundance_sim_B <- function(N, OTU_template, OTU_vars, phen_means,
                                phen_var, cor_mat, zero_props, n.seq, const.reads = F){
  p <- length(OTU_template)
  q <- length(phen_means)
  # Simulating initial variable
  g <- mvtnorm::rmvnorm(n = N, sigma = cor_mat)
  
  # Defining the inverse cdf of zero-inflated log-normal distribution
  F_inv <- function(x, pi, mu, sd){
    y <- rep(NA, length(x))
    y[x<=pi] <- 0
    y[x>pi] <- qlnorm((x[x>pi] - pi)/(1 - pi), meanlog = mu, sdlog = sd)
    return(y)
  }
  
  # Computing standard deviations
  sds <- sqrt(OTU_vars)
  phen_sds <- sqrt(phen_var)
  
  A <- matrix(NA, N, p)
  for(i in 1:p){
    # Logical vector for when the absolute abundance should be set to zero
    idx <- g[,i] < qnorm(zero_props[i])
    # Computing absolute abundances based on the initial variable
    A[idx,i] <- 0
    A[!idx,i] <- F_inv(pnorm(g[!idx,i]),
                       pi = zero_props[i],
                       mu = OTU_template[i],
                       sd = sds[i])
  }
  # Computing phenotypic variables based on the initial variable
  B <- matrix(NA, N, q)
  for(i in 1:q){
    B[,i] <- F_inv(pnorm(g[,i+p]),
                   pi = 0,
                   mu = phen_means[i],
                   sd = phen_sds[i])
  }

  # Simulating the sequecning step
  OTU_reads <- Sequencing_sim(A, n.seq, const.reads = const.reads)
  Phenotype <- log(B)
  out <- list(OTU_reads = OTU_reads, Phenotype = Phenotype, Abs_OTU = A)
  return(out)
}

Abs_abundance_sim_C <- function(N, OTU_template, Gene_template, OTU_vars, Gene_vars,
                                cor_mat, zero_props_OTU, zero_props_Gene,
                                n.seq.OTU, n.seq.gene){
  p <- length(OTU_template)
  q <- length(Gene_template)
  # Simulating initial variable
  g <- mvtnorm::rmvnorm(n = N, sigma = cor_mat)
  
  # Defining the inverse cdf of zero-inflated log-normal distribution
  F_inv <- function(x, pi, mu, sd){
    y <- rep(NA, length(x))
    y[x<=pi] <- 0
    y[x>pi] <- qlnorm((x[x>pi] - pi)/(1 - pi), meanlog = mu, sdlog = sd)
    return(y)
  }
  
  # Computing standard deviations
  OTU_sds <- sqrt(OTU_vars)
  Gene_sds <- sqrt(Gene_vars)
  
  # Computing absolute abundances based on the initial variable
  A <- matrix(NA, N, p)
  for(i in 1:p){
    idx <- g[,i] < qnorm(zero_props_OTU[i])
    A[idx,i] <- 0
    A[!idx,i] <- F_inv(pnorm(g[!idx,i]),
                       pi = zero_props_OTU[i],
                       mu = OTU_template[i],
                       sd = OTU_sds[i])
  }
  
  # Computing absolute expression levels based on the initial variable
  B <- matrix(NA, N, q)
  for(i in 1:q){
    idx <- g[,i+p] < qnorm(zero_props_Gene[i])
    B[idx,i] <- 0
    B[!idx,i] <- F_inv(pnorm(g[!idx,i+p]),
                       pi = zero_props_Gene[i],
                       mu = Gene_template[i],
                       sd = Gene_sds[i])
  }
  
  # Simulating the sequecning step
  OTU_reads <- Sequencing_sim(A, n.seq.OTU)
  Gene_reads <- Sequencing_sim(B, n.seq.gene)
  out <- list(OTU_reads = OTU_reads, Gene_reads = Gene_reads, Abs_OTU = A, Abs_Gene = B)
  return(out)
}
