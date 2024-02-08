# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pkg <- c("data.table", "Matrix", "hydroGOF", "ggplot2", "RColorBrewer",
         "compositions", "psych", "foreach", "parallel")
for(i in pkg){
  library(i, character.only = T)
}
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source(file = "../Data_generation.R")
source(file = "../SparCEV_function.R")

# Loading template
template <- fread("Templates/OTU_template.csv")

source(file = "../cor_mat_cluster.R")

SparCEV_cov <- function(OTU, phenotype, pseudo_cout = 1, var_min = 1e-5){
  p <- ncol(OTU)
  # Estimating relative abundances
  OTU_TSS <- (OTU+pseudo_cout)/rowSums(OTU+pseudo_cout)
  log_OTU <- log(OTU_TSS)
  log_sum <- rowSums(log_OTU)
  # Approximation of correlations between abundances and phenotype
  log_ratio_cov <- apply(log_OTU, 2, function(x) cov(x - (log_sum-x)/(p-1), phenotype))
  return(log_ratio_cov)
}

SparCEV_cov_single_OTU <- function(OTU, phenotype, pseudo_cout = 1, var_min = 1e-5, i){
  p <- ncol(OTU)
  # Estimating relative abundances
  OTU_TSS <- (OTU+pseudo_cout)/rowSums(OTU+pseudo_cout)
  log_OTU <- log(OTU_TSS)
  log_sum <- rowSums(log_OTU)
  # Approximation of correlations between abundances and phenotype
  # log_ratio_cov <- apply(log_OTU, 2, function(x) cov(x - (log_sum-x)/(p-1), phenotype))
  x <- log_OTU[,i]
  log_ratio_cov <- cov(x - (log_sum-x)/(p-1), phenotype)
  
  # cov(log_OTU[i,] - (log_sum-log_OTU[i,])/(p-1), phenotype)
  return(log_ratio_cov)
}

rarefy_single_sample <- function(x, N){
  p <- length(x)
  names(x) <- paste0("OTU", 1:p)
  rare_vec <- rep(0, p); names(rare_vec) <- names(x)
  long_vec <- rep(names(x), x)
  rarefy_OTUs <- sample(long_vec, replace = F, size = N)
  rarefy <- table(rarefy_OTUs)
  rare_vec[names(rarefy)] <- rarefy
  return(as.numeric(rare_vec))
}

SparCEV_test <- function(OTU_table, Phenotype, B = 1000, cores = 1, lib_cut = ncol(OTU_table)*10){
  lib_size <- rowSums(OTU_table)
  OTU_table <- OTU_table[lib_size > lib_cut,]
  Phenotype <- Phenotype[lib_size > lib_cut]
  N <- min(rowSums(OTU_table))
  OTU_table <- t(apply(OTU_table, 1 , rarefy_single_sample, N = N))
  C <- SparCEV_cov(OTU_table, Phenotype)
  p <- ncol(OTU_table)
  n <- nrow(OTU_table)
  p_val <- rep(NA, p)
  for(i in 1:p){
    cat(i, "\r")
    Permute_list <- lapply(1:B, 
                           function(x){
                             A2 <- OTU_table
                             A2[, i] <- A2[sample(1:n), i]
                             return(A2)
                           })
    if(cores == 1){
      Permute_cor <- lapply(Permute_list,
                            function(x) SparCEV_cov_single_OTU(x, Phenotype, i=i))
    }
    
    if(cores>1){
      Permute_cor <- mclapply(Permute_list,
                              function(x) SparCEV_cov_single_OTU(x, Phenotype, i=i),
                              mc.cores = cores)
    }
    cors <- unlist(Permute_cor)
    # b <- sum(abs(cors) > abs(C[i]))
    b <- sum(abs(cors-mean(cors)) > abs(C[i]-mean(cors)))
    p_val[i] <- (b+1)/(B+1)
  }
  return(p_val)
}

SparCEV_test_full_MC <- function(OTU_table, Phenotype, B = 1000, cores = 1, lib_cut = ncol(OTU_table)*10){
  lib_size <- rowSums(OTU_table)
  OTU_table <- OTU_table[lib_size > lib_cut,]
  Phenotype <- Phenotype[lib_size > lib_cut]
  N <- min(rowSums(OTU_table))
  OTU_table <- t(apply(OTU_table, 1 , rarefy_single_sample, N = N))
  C <- SparCEV_cov(OTU_table, Phenotype)
  p <- ncol(OTU_table)
  n <- nrow(OTU_table)
  p_val <- rep(NA, p)
  
  Permute_settings <- rep(1:p, B)
  
  Permute_list <- lapply(Permute_settings, 
                         function(i){
                           A2 <- OTU_table
                           A2[, i] <- A2[sample(1:n), i]
                           return(A2)
                         })
  
  Permute_list <- lapply(1:length(Permute_settings),
                         function(i) list(OTU_table = Permute_list[[i]],
                                          OTU_num = Permute_settings[i])
  )
  
  if(cores == 1){
    Permute_cor <- lapply(Permute_list,
                          function(x) SparCEV_cov_single_OTU(x$OTU_table,
                                                             Phenotype,
                                                             i=x$OTU_num))
  }
  
  if(cores>1){
    Permute_cor <- mclapply(Permute_list,
                            function(x) SparCEV_cov_single_OTU(x$OTU_table,
                                                               Phenotype,
                                                               i=x$OTU_num),
                            mc.cores = cores)
  }
  cors <- unlist(Permute_cor)
  Permute_cov_mat <- matrix(cors, p, B)
  Distance_from_NULL_permute <- abs( Permute_cov_mat - rowMeans(Permute_cov_mat))
  Distance_from_NULL_est <- abs( C - rowMeans(Permute_cov_mat) )
  
  b <- rowSums(Distance_from_NULL_permute > Distance_from_NULL_est)
  p_val <- (b+1)/(B+1)
  
  return(p_val)
}

SparCEV_test_full_MC2 <- function(OTU_table, Phenotype, B = 1000, cores = 1, lib_cut = ncol(OTU_table)*10){
  lib_size <- rowSums(OTU_table)
  OTU_table <- OTU_table[lib_size > lib_cut,]
  Phenotype <- Phenotype[lib_size > lib_cut]
  N <- min(rowSums(OTU_table))
  OTU_table <- t(apply(OTU_table, 1 , rarefy_single_sample, N = N))
  C <- SparCEV_cov(OTU_table, Phenotype)
  p <- ncol(OTU_table)
  n <- nrow(OTU_table)
  p_val <- rep(NA, p)
  
  Permute_settings <- rep(1:p, B)
  
  if(cores == 1){
    Permute_cov <- lapply(Permute_settings,
                          function(i){
                            A2 <- OTU_table
                            A2[, i] <- A2[sample(1:n), i]
                            cov <- SparCEV_cov_single_OTU(A2,
                                                          Phenotype,
                                                          i=i)
                            return(cov)
                          })
  }
  
  if(cores > 1){
    Permute_cov <- mclapply(Permute_settings,
                          function(i){
                            A2 <- OTU_table
                            A2[, i] <- A2[sample(1:n), i]
                            cov <- SparCEV_cov_single_OTU(A2,
                                                          Phenotype,
                                                          i=i)
                            return(cov)
                          },
                          mc.cores = cores)
  }
  
  Permute_cov_mat <- matrix(unlist(Permute_cov), p, B)
  Distance_from_NULL_permute <- abs( Permute_cov_mat - rowMeans(Permute_cov_mat))
  Distance_from_NULL_est <- abs( C - rowMeans(Permute_cov_mat) )
  
  b <- rowSums(Distance_from_NULL_permute > Distance_from_NULL_est)
  p_val <- (b+1)/(B+1)
  
  return(p_val)
}

# Constructing matrix with settings
p <- c(100, 250, 1000)
dens <- 0.1
cor <- 0.75
bio_zero <- "No"

# Defining other simulation settings
n.sim <- c(20, 50, 200, 1000)
phenotype_var <- 1

reps <- 25

Opt <- expand.grid(p = p, dens = dens, bio_zero = bio_zero, cor = cor, n.sim = n.sim)
Opt <- do.call("rbind", replicate(reps, Opt, simplify = FALSE))

Res <- data.table(Method = NA, Power = NA, FDR = NA, Opt[1,])[-1]

RNGkind("L'Ecuyer-CMRG")
set.seed(1707394915)
for(k in 1:nrow(Opt)){
  print(k)
  
  OTU.idx <- sample(1:1000, Opt[k, "p"])
  template_sim <- template[OTU.idx]
  if(Opt[k, "bio_zero"] == "No"){
    template_sim[,pi:=0]
  }
  n.seq.abn <- p*40
  
  # Constructing correlation matrix
  cor_mat <- cor_mat_sim_cluster(p = Opt[k, "p"],
                                 q = 1,
                                 cor.prob = Opt[k, "dens"],
                                 cor = Opt[k, "cor"],
                                 case = "B" )
  true.cor <- cor_mat[1:Opt[k, "p"], Opt[k, "p"]+1]
  
  # Simulating sequencing data 
  Abs_abundance_sim_B(N = Opt[k, "n.sim"], 
                      OTU_template = template_sim$mu,
                      OTU_vars = template_sim$sigma2,
                      phen_means = log(30) - log(1+1/30^2)/2,
                      phen_var = log(1+1/30^2),
                      cor_mat = cor_mat,
                      zero_props = template_sim$pi,
                      n.seq = n.seq.abn) -> A
  
  # true.cov <- true.cor*sqrt(log(30) - log(1+1/30^2)/2)*sqrt(template_sim$sigma2)
  
  # Data transformations
  OTU_TSS <- (A$OTU_reads+1)/rowSums(A$OTU_reads+1)
  OTU_CLR <- t(apply(OTU_TSS, 1, clr))
  
  # Test log-TSS
  p_val_log <- sapply(1:Opt[k,"p"], function(i) cor.test(log(A$OTU_reads[,i]+1), A$Phenotype)$p.value)
  p_val_log_adj <- p.adjust(p_val_log, method = "fdr")
  Power_log <- mean(p_val_log_adj[true.cor != 0] < 0.05)
  FDR_log <- mean(true.cor[p_val_log_adj < 0.05] == 0)
  print( paste0("log: Power: ", Power_log, " FDR: ", FDR_log) )
  
  Res <- rbind(Res, data.table(Method = "log",
                               Power = Power_log,
                               FDR = FDR_log,
                               Opt[k,]))
  
  # Test log-TSS
  p_val_TSS <- sapply(1:Opt[k,"p"], function(i) cor.test(log(OTU_TSS[,i]), A$Phenotype)$p.value)
  p_val_TSS_adj <- p.adjust(p_val_TSS, method = "fdr")
  Power_TSS <- mean(p_val_TSS_adj[true.cor != 0] < 0.05)
  FDR_TSS <- mean(true.cor[p_val_TSS_adj < 0.05] == 0)
  print( paste0("TSS: Power: ", Power_TSS, " FDR: ", FDR_TSS) )
  
  Res <- rbind(Res, data.table(Method = "log-TSS",
                               Power = Power_TSS,
                               FDR = FDR_TSS,
                               Opt[k,]))
  
  # Test CLR
  p_val_clr <- sapply(1:Opt[k,"p"], function(i) cor.test(OTU_CLR[,i], A$Phenotype)$p.value)
  p_val_clr_adj <- p.adjust(p_val_clr, method = "fdr")
  Power_CLR <- mean(p_val_clr_adj[true.cor != 0] < 0.05)
  FDR_CLR <- mean(true.cor[p_val_clr_adj < 0.05] == 0)
  print( paste0("CLR: Power: ", Power_CLR, " FDR: ", FDR_CLR) )
  
  Res <- rbind(Res, data.table(Method = "CLR",
                               Power = Power_CLR,
                               FDR = FDR_CLR,
                               Opt[k,]))
  
  # Test permutation testing scheme
  # a <- Sys.time()
  # p_val <- SparCEV_test(OTU_table = A$OTU_reads,
  #                       Phenotype = A$Phenotype,
  #                       B = 1000,
  #                       cores = 2)
  # Sys.time() - a
  # 
  # a <- Sys.time()
  # p_val2 <- SparCEV_test_full_MC(OTU_table = A$OTU_reads,
  #                       Phenotype = A$Phenotype,
  #                       B = 1000,
  #                       cores = 2)
  # Sys.time() - a
  
  # a <- Sys.time()
  p_val <- SparCEV_test_full_MC2(OTU_table = A$OTU_reads,
                                  Phenotype = A$Phenotype,
                                  B = 1000,
                                  cores = 40)
  # Sys.time() - a
  
  # p_val_adj <- p.adjust(p_val, method = "fdr")
  # Power_CEV[k] <- mean(p_val_adj[true.cor != 0] < 0.05)
  # FDR_CEV[k] <- mean(true.cor[p_val_adj < 0.05] == 0)
  # print( paste0("SparCEV: Power: ", Power_CEV[k], " FDR: ", FDR_CEV[k]) )
  
  p_val_adj <- p.adjust(p_val, method = "fdr")
  Power_CEV <- mean(p_val_adj[true.cor != 0] < 0.05)
  FDR_CEV <- mean(true.cor[p_val_adj < 0.05] == 0)
  print( paste0("SparCEV: Power: ", Power_CEV, " FDR: ", FDR_CEV) )
  
  Res <- rbind(Res, data.table(Method = "SparCEV",
                               Power = Power_CEV,
                               FDR = FDR_CEV,
                               Opt[k,]))
}

fwrite(Res, "Res.csv")

Res[is.na(FDR),FDR:=0]
Res[,.(FDR = mean(FDR), Power = mean(Power),
       FDR_lower = mean(FDR) - sd(FDR)*sqrt(reps),
       FDR_upper = mean(FDR) + sd(FDR)*sqrt(reps),
       Power_lower = mean(Power) - sd(Power)/sqrt(reps),
       Power_upper = mean(Power) + sd(Power)/sqrt(reps)),
    list(Method, p, dens, bio_zero, n.sim)] -> Res_summary

Res_temp1 <- Res_summary[,-c(7,10:11)]; Res_temp1[,Type:="FDR"]
Res_temp2 <- Res_summary[,-c(6,8:9)]; Res_temp2[,Type:="Power"]
colnames(Res_temp1)[6:8] <- c("Measure", "Lower", "Upper")
colnames(Res_temp2)[6:8] <- c("Measure", "Lower", "Upper")
Res_summary_long <- rbind(Res_temp1, Res_temp2)

fwrite(Res_summary_long, "Res_summary_long.csv")

Res_summary_long[,p:=paste0("p = ", dens)]

ggplot(data = Res_summary_long, aes(x = as.factor(n.sim),
                                    y = Measure,
                                    col = Method,
                                    fill = Method,
                                    group = Method))+
  geom_line(linewidth = 1)+
  geom_ribbon(alpha=0.35, aes(col = NULL, ymin = Lower, ymax = Upper))+
  # geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_zero_min, ymax = MAE_zero_max))+
  scale_color_brewer(" Method", palette = "Dark2")+
  scale_fill_brewer(" Method", palette = "Dark2")+
  # scale_linetype_manual("Correlations",values=c("Non-zero"=1,"Zero"=2))+
  labs(x = "Number of replicates", y = NULL)+
  facet_grid(p ~ Type, switch = "y")+
  guides(linetype = guide_legend(order = 1))+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = 10)) -> g

ggsave(filename = "Figures/Case_B_Testing.pdf", g, width = 190, height = 160, unit = "mm")
ggsave(filename = "Figures/Case_B_Testing.eps", g, width = 190, height = 160, dpi = 300, unit = "mm")

