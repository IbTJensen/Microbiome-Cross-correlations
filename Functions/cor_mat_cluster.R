cor_mat_sim_cluster <- function(p, q, cor.prob = 0.1, cor, case = "B"){
  if(case == "C" & q == 1){
    stop("In case C, q>1 is required")
  }
  
  cor_amount_a <- max(1, floor(p*cor.prob))
  cor_amount_b <- max(1, floor(q*cor.prob))
  
  # Sampling the OTU-indices to be in the cluster
  p_clust <- sample(1:p, cor_amount_a)
  if(q == 1){
    q_clust <- p+1
  }
  if(q>1){
    q_clust <- sample(p+1:q, cor_amount_b)
  }
  
  c <- runif(1,0,1)
  p_clust_plus <- sample(p_clust, floor(length(p_clust)*c))
  p_clust_minus <- p_clust[!(p_clust %in% p_clust_plus)]
  
  # As a consequence of these lines, when q=1 in case B, b is always in the negative part of the cluster
  q_clust_plus <- sample(q_clust, floor(length(q_clust)*c))
  q_clust_minus <- q_clust[!(q_clust %in% q_clust_plus)]
  
  # Creating a "positive" and "negative" cluster
  pos_clust <- sort(c(p_clust_plus, q_clust_plus))
  neg_clust <- sort(c(p_clust_minus, q_clust_minus))
  
  # Filling out the correlation matrix
  cor.mat <- matrix(0, p+q, p+q)
  cor.mat[pos_clust, pos_clust] <- cor
  cor.mat[neg_clust, neg_clust] <- cor
  cor.mat[pos_clust, neg_clust] <- -cor
  cor.mat[neg_clust, pos_clust] <- -cor
  diag(cor.mat) <- 1
  
  # Assigning names to columns and rows
  colnames(cor.mat) <- c(paste0("OTU", 1:p), rep("", q))
  if(case == "B"){
    colnames(cor.mat)[p+1:q] <- paste0("Phenotype", 1:q)
  }
  if(case == "C"){
    colnames(cor.mat)[p+1:q] <- paste0("Gene", 1:q)
  }
  rownames(cor.mat) <- colnames(cor.mat)
  
  return(cor.mat)
}

cor_mat_sim_cluster_positive <- function(p, q, cor.prob = 0.1, cor, case = "B"){
  if(case == "C" & q == 1){
    stop("In case C, q>1 is required")
  }
  
  cor_amount_a <- max(1, floor(p*cor.prob))
  cor_amount_b <- max(1, floor(q*cor.prob))
  
  # Sampling the OTU-indices to be in the cluster
  p_clust <- sample(1:p, cor_amount_a)
  if(q == 1){
    q_clust <- p+1
  }
  if(q>1){
    q_clust <- sample(p+1:q, cor_amount_b)
  }
  
  cor.mat <- matrix(0, p+q, p+q)
  cor.mat[c(p_clust, q_clust), c(p_clust, q_clust)] <- cor
  diag(cor.mat) <- 1
  
  # Assigning names to columns and rows
  colnames(cor.mat) <- c(paste0("OTU", 1:p), rep("", q))
  if(case == "B"){
    colnames(cor.mat)[p+1:q] <- paste0("Phenotype", 1:q)
  }
  if(case == "C"){
    colnames(cor.mat)[p+1:q] <- paste0("Gene", 1:q)
  }
  rownames(cor.mat) <- colnames(cor.mat)
  
  return(cor.mat)
}
