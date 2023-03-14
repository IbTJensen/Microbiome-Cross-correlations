# Fast SparCC-derived approximation of absolute abundance variances
Absolute_abn_var <- function(OTU.abn.log, var_min = 1e-5){
  p <- ncol(OTU.abn.log)
  x_var <- apply(OTU.abn.log, 2, var)
  log_sum <- rowSums(OTU.abn.log)
  t <- p*x_var + sum(x_var) - 2*apply(OTU.abn.log, 2, function(x) cov(x, log_sum))
  M <- matrix(1,p,p) + diag(p-2, p, p)
  alpha <- solve(M, t)
  alpha[alpha < var_min] <- var_min
  return(alpha)
}

SparCEV <- function(OTU, phenotype, pseudo_cout = 1, var_min = 1e-5){
  p <- ncol(OTU)
  # Estimating relative abundances
  OTU_TSS <- (OTU+pseudo_cout)/rowSums(OTU+pseudo_cout)
  log_OTU <- log(OTU_TSS)
  log_sum <- rowSums(log_OTU)
  # Approximating absolute abundance variances
  alpha <- Absolute_abn_var(log_OTU, var_min = var_min)
  # Approximation of correlations between abundances and phenotype
  log_ratio_cov <- apply(log_OTU, 2, function(x) cov(x - (log_sum-x)/(p-1), phenotype))
  phenotype_sd <- sd(phenotype)
  cors <- log_ratio_cov/(sqrt(alpha)*phenotype_sd)
  cors[abs(cors)>1] <- sign(cors[abs(cors)>1])
  return(cors)
}

# Dirichlet Monte-Carlo sampling
dir.sim <- function(X){
  A <- t(apply(X, 1, function(x) gtools::rdirichlet(1, alpha = x+1)))
  return(A)
}

# SparCEV utilizing Dirichlet Monte-Carlo sampling
SparCEV_MC <- function(OTU, phenotype, var_min = 1e-5, dir_it = 20){
  p <- ncol(OTU)
  # Apply Dirchlet Monte-Carlo sampling and obtain multiple data sets
  L <- lapply(1:dir_it, function(x) dir.sim(OTU))
  # Apply SparCEV to each Monte-Carlo sampled data set
  temp_res <- lapply(L, function(x) SparCEV(x, phenotype, pseudo_cout = 0))
  # Arrange in array
  res_arr <- array(unlist(temp_res), c(p,dir_it))
  # Correlation estimates are the median over the results obtained from the sampled data sets
  cors <- apply(res_arr, 1, median)
  return(cors)
}
