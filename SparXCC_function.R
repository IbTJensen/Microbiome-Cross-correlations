library(limma)

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

SparXCC <- function(OTU.abn, gene.expr, pseudo_count = 1, var_min = 1e-5){
  D1 <- t(apply(OTU.abn+pseudo_count, 1, function(x) x/sum(x)))
  D2 <- t(apply(gene.expr+pseudo_count, 1, function(x) x/sum(x)))
  
  log.OTU <- log(D1); p <- ncol(log.OTU)
  log.genes <- log(D2); q <- ncol(log.genes)
  
  # Below t_ik's are computed in an effiicent manner. See (8) in supplementary material of paper
  # Variances of log-ratios between OTUs and genes
  var_ik <- as.numeric(apply(log.genes, 2, function(y) apply(log.OTU, 2, function(x) var(x - y, na.rm = T)) ))# taxa varies first
  # Variances of OTU relative abundances
  var.OTU <- apply(log.OTU, 2, var, na.rm = T)
  var.all.OTU <- sum( var.OTU )
  # Variances of gene relative expression levels
  var.genes <- apply(log.genes, 2, var, na.rm = T)
  var.all.genes <- sum( var.genes )# - var.all
  # Covariances between OTU- and gene-sums
  cov.all <- cov(rowSums(log.OTU, na.rm = T), rowSums(log.genes, na.rm = T))
  # Covariances between log-ratios and sum-contrasts
  sum.contr <- q*rowSums(log.OTU, na.rm = T) - p*rowSums(log.genes, na.rm = T)
  cov_ik <- as.numeric(apply(log.genes, 2, function(y) apply(log.OTU, 2, function(x) cov(x - y, sum.contr)) )) # taxa varies first
  # Estimating t_ik
  t_ik <- p*q*var_ik + q*var.all.OTU + p*var.all.genes - 2*cov.all - 2*cov_ik
  
  # Basis variances
  alpha <- Absolute_abn_var(log.OTU)
  beta <- Absolute_abn_var(log.genes)
  
  alpha_rep <- rep(alpha, q)
  beta_rep <- rep(beta, each = p)
  
  # Terms in correlation approximation. See (6) in paper
  pqa <- (p-1)*q*alpha_rep
  pqb <- p*(q-1)*beta_rep
  qa_sum <- q*( sum(alpha)-alpha_rep )
  pb_sum <- p*( sum(beta)-beta_rep )
  OTU.sd <- sqrt(alpha_rep)
  gene.sd <- sqrt(beta_rep)
  
  # Computing correlation approximations
  rho <- ( pqa + pqb + qa_sum + pb_sum - t_ik )/( 2*(p-1)*(q-1)*OTU.sd*gene.sd )
  rho[abs(rho)>1] <- sign(rho[abs(rho)>1])

  cor.mat.est <- matrix(rho, nrow = p, ncol = q)
  rownames(cor.mat.est) <- colnames(OTU.abn)
  colnames(cor.mat.est) <- colnames(gene.expr)
  
  return(cor.mat.est)
}
