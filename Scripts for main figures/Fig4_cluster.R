# Code by Ib Thorsgaard Jensen
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Calling other files
source(file = "Shared_settings_Fig4.R")
source(file = "../Functions/cor_mat_cluster.R")

# Constructing matrix with settings
p <- c(10,20,50,100,250,1000)
q <- c(10, 100, 1000)
dens <- c(0.1, 0.4, 0.7)
cor <- c(0.75)

Opt <- expand.grid(p = p, q = q, dens = dens, bio_zero = bio_zero, cor = cor)
Opt <- do.call("rbind", replicate(reps, Opt, simplify = FALSE))

# Selecting t via bootstrap for each combination of p and q.
# This is done once for each q and q rather than for every
# dataset to save time.
Opt_pq <- expand.grid(p = p, q = q)
OTU_thresh_vec <- rep(NA, nrow(Opt_pq))
Gene_thresh_vec <- rep(NA, nrow(Opt_pq))
set.seed(1676299396)
for(i in 1:nrow(Opt_pq)){
  cat(i, "\r")
  p <- Opt[i, "p"]
  q <- Opt[i, "q"]
  
  # Randomly selecting OTUs from the template
  OTU.idx <- sample(1:nrow(template_OTU), p)
  template_OTU_sim <- template_OTU[OTU.idx]
  if(Opt[i, "bio_zero"] == "No"){
    template_OTU_sim[,pi:=0]
  }
  
  # Randomly selecting genes from the template
  Gene.idx <- sample(1:nrow(template_gene), q)
  template_gene_sim <- template_gene[Gene.idx]
  if(Opt[i, "bio_zero"] == "No"){
    template_gene_sim[,pi:=0]
  }
  
  # Defining other simulation settings
  n.sim <- 50
  n.seq.abn <- p*40
  n.seq.expr <- q*200
  
  # Constructing correlation matrix
  cor_mat <- cor_mat_sim_cluster(p,
                                 q,
                                 cor.prob = 0.1,
                                 cor = 0.75,
                                 case = "C")
  true.cor.mat <- cor_mat[1:p,p+1:q]
  true.cor <- as.numeric(true.cor.mat)
  
  # Simulating sequencing data 
  Abs_abundance_sim_C(N = n.sim, 
                      OTU_template = template_OTU_sim$mu,
                      Gene_template = template_gene_sim$mu,
                      OTU_vars = template_OTU_sim$sigma2,
                      Gene_vars = template_gene_sim$sigma2,
                      cor_mat = cor_mat,
                      zero_props_OTU = template_OTU_sim$pi,
                      zero_props_Gene = template_gene_sim$pi,
                      n.seq.OTU = n.seq.abn,
                      n.seq.gene = n.seq.expr) -> A
  
  S <- mclapply(1:100,
                function(x){
                  S <- SparXCC_base(OTU.abn = A$OTU_reads[sample(1:n.sim),],
                                    gene.expr = A$Gene_reads[sample(1:n.sim),],
                                    Find_m = F)$cor
                  list(OTU = rowMeans(abs(S)), Gene = colMeans(abs(S)))
                  }, 
                mc.cores = 1)
  
  OTU_thresh_vec[i] <- lapply(S, function(x) quantile(x$OTU, 0.8) ) %>% unlist() %>% mean()
  Gene_thresh_vec[i] <- lapply(S, function(x) quantile(x$Gene, 0.8) ) %>% unlist() %>% mean()
}
data.table(t = OTU_thresh_vec, p = Opt_pq[,"p"], q = Opt_pq[,"q"])

OTU_thresh_mat <- matrix(OTU_thresh_vec, 6, 3)
rownames(OTU_thresh_mat) <- as.character(c(10,20,50,100,250,1000))
colnames(OTU_thresh_mat) <- as.character(c(10,100,1000))

Gene_thresh_mat <- matrix(Gene_thresh_vec, 6, 3)
rownames(Gene_thresh_mat) <- as.character(c(10,20,50,100,250,1000))
colnames(Gene_thresh_mat) <- as.character(c(10,100,1000))

# Setting up parallelization
my.cluster <- parallel::makeCluster(
  cores, 
  type = "FORK"
)

registerDoParallel(cl = my.cluster)
registerDoRNG(seed = 1676299396)

a <- Sys.time()
R <- foreach(
  i = 1:nrow(Opt),
  .combine = "rbind"
) %dopar% {
  p <- Opt[i, "p"]
  q <- Opt[i, "q"]
  
  # Randomly selecting OTUs from the template
  OTU.idx <- sample(1:nrow(template_OTU), p)
  template_OTU_sim <- template_OTU[OTU.idx]
  if(Opt[i, "bio_zero"] == "No"){
    template_OTU_sim[,pi:=0]
  }
  
  # Randomly selecting genes from the template
  Gene.idx <- sample(1:nrow(template_gene), q)
  template_gene_sim <- template_gene[Gene.idx]
  if(Opt[i, "bio_zero"] == "No"){
    template_gene_sim[,pi:=0]
  }
  
  # Defining other simulation settings
  n.sim <- 50
  n.seq.abn <- p*40
  n.seq.expr <- q*200
  
  # Constructing correlation matrix
  cor_mat <- cor_mat_sim_cluster(p,
                                 q,
                                 cor.prob = Opt[i,"dens"],
                                 cor = Opt[i,"cor"],
                                 case = "C")
  true.cor.mat <- cor_mat[1:p,p+1:q]
  true.cor <- as.numeric(true.cor.mat)
  
  # Simulating sequencing data 
  Abs_abundance_sim_C(N = n.sim, 
                      OTU_template = template_OTU_sim$mu,
                      Gene_template = template_gene_sim$mu,
                      OTU_vars = template_OTU_sim$sigma2,
                      Gene_vars = template_gene_sim$sigma2,
                      cor_mat = cor_mat,
                      zero_props_OTU = template_OTU_sim$pi,
                      zero_props_Gene = template_gene_sim$pi,
                      n.seq.OTU = n.seq.abn,
                      n.seq.gene = n.seq.expr) -> A
  
  colnames(A$OTU_reads) <- paste0("OTU", 1:Opt[i,"p"])
  colnames(A$Gene_reads) <- paste0("Gene", 1:Opt[i,"q"])
  
  # Data transformations
  OTU_TSS <- (A$OTU_reads+1)/rowSums(A$OTU_reads+1)
  Gene_TSS <- (A$Gene_reads+1)/rowSums(A$Gene_reads+1)
  OTU_CLR <- t(apply(OTU_TSS, 1, clr))
  Gene_CLR <- t(apply(Gene_TSS, 1, clr))
  
  # VST may sometimes fail, if there are not genes without zero-count replicates. In such cases, the estimation is skipped
  Gene_VST <- tryCatch(t(varianceStabilizingTransformation(t(A$Gene_reads), fitType = "local")), error = function(e) "Error")
  if(length(class(Gene_VST)) > 1){
    Mix <- as.numeric(cor(OTU_CLR, Gene_VST))
  } else{
    Mix <- NA
  }
  
  # Estimating correlations
  C <- as.numeric(cor(log(A$OTU_reads+1), log(A$Gene_reads+1)))
  CC <- as.numeric(cor(log(OTU_TSS), log(Gene_TSS)))
  CLR <- as.numeric(cor(OTU_CLR, Gene_CLR))
  XCC <- SparXCC_base(OTU.abn = A$OTU_reads,
                      gene.expr = A$Gene_reads,
                      Find_m = F)
  
  # Extracting bootstrap thresholds
  pc <- as.character(Opt[i,"p"])
  qc <- as.character(Opt[i,"q"])
  thr_OTU <- OTU_thresh_mat[pc, qc]
  thr_Gene <- Gene_thresh_mat[pc, qc]
  
  XCC2 <- SparXCC(OTU.abn = A$OTU_reads,
                  gene.expr = A$Gene_reads,
                  pseudo_count = 1,
                  iter = 10,
                  t1 = thr_OTU,
                  t2 = thr_Gene,
                  Find_m = F)
  
  XCC <- as.numeric(XCC$cor)
  XCC2 <- as.numeric(XCC2$cor)
  
  # Setting up tables with results
  mae_C <- data.table(Opt[i,],
                      MAE_nonzero=mae(true.cor[true.cor!=0], C[true.cor!=0]),
                      MAE_zero=mae(true.cor[true.cor==0], C[true.cor==0]),
                      Method = "log")
  mae_CC <- data.table(Opt[i,],
                       MAE_nonzero=mae(true.cor[true.cor!=0], CC[true.cor!=0]),
                       MAE_zero=mae(true.cor[true.cor==0], CC[true.cor==0]),
                       Method = "log-TSS")
  mae_CLR <- data.table(Opt[i,],
                        MAE_nonzero=mae(true.cor[true.cor!=0], CLR[true.cor!=0]),
                        MAE_zero=mae(true.cor[true.cor==0], CLR[true.cor==0]),
                        Method = "CLR")
  mae_mix <- data.table(Opt[i,],
                        MAE_nonzero=mae(true.cor[true.cor!=0], Mix[true.cor!=0]),
                        MAE_zero=mae(true.cor[true.cor==0], Mix[true.cor==0]),
                        Method = "CLR+VST")
  mae_XCC <- data.table(Opt[i,],
                        MAE_nonzero=mae(true.cor[true.cor!=0], XCC[true.cor!=0]),
                        MAE_zero=mae(true.cor[true.cor==0], XCC[true.cor==0]),
                        Method = "SparXCC base")
  mae_XCC2 <- data.table(Opt[i,],
                        MAE_nonzero=mae(true.cor[true.cor!=0], XCC2[true.cor!=0]),
                        MAE_zero=mae(true.cor[true.cor==0], XCC2[true.cor==0]),
                        Method = "SparXCC iterative")
  B <- rbind(mae_C, mae_CC, mae_CLR, mae_mix, mae_XCC, mae_XCC2); B
}
Sys.time() - a
# Saving results in a csv-file
fwrite(R, file = paste0("../Results/Results_caseC_cluster_", as.character(Sys.Date()), ".csv"))

# Summarizing accuracy separated by the selected simulation settings
R1 <- R[,.(MAE_nonzero = mean(MAE_nonzero, na.rm = T),
           MAE_nonzero_min = mean(MAE_nonzero, na.rm = T) - 1.96*sd(MAE_nonzero, na.rm = T)/sqrt(reps),
           MAE_nonzero_max = mean(MAE_nonzero, na.rm = T) + 1.96*sd(MAE_nonzero, na.rm = T)/sqrt(reps),
           MAE_zero = mean(MAE_zero, na.rm = T),
           MAE_zero_min = mean(MAE_zero, na.rm = T) - 1.96*sd(MAE_zero, na.rm = T)/sqrt(reps),
           MAE_zero_max = mean(MAE_zero, na.rm = T) + 1.96*sd(MAE_zero, na.rm = T)/sqrt(reps)),
        list(p, q, dens, Method, bio_zero)]

# R1[,dens:=paste0("c=", dens)]
# R1[,Cor_mat := "Cluster method"]
R1[,":="(Cor_mat = "Cluster method",
         q = factor(paste0("q=", q), levels = c("q=10", "q=100", "q=1000")))]

R1_yes_main <- R1[bio_zero == "Yes" & q == "q=1000"]
R1_yes_supp <- R1[bio_zero == "Yes" & q != "q=1000"]
R1_no_main <- R1[bio_zero == "No" & q == "q=1000"]
# R1_no_supp <- R1[bio_zero == "No" & dens != 0.1]
R1_yes_supp[,dens:=paste("c", dens, sep = "=")]
R1_yes_main[,dens:=paste("c", dens, sep = "=")]
R1_no_main[,dens:=paste("c", dens, sep = "=")]

# Constructing and saving figures
ggplot(data = R1_yes_main, aes(x = as.factor(p), col = Method, fill = Method, group = Method))+
  geom_line(linewidth = 1, aes(y = MAE_nonzero, linetype = "Non-zero"))+
  geom_line(linewidth = 1, aes(y = MAE_zero, linetype = "Zero"))+
  geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_nonzero_min, ymax = MAE_nonzero_max))+
  geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_zero_min, ymax = MAE_zero_max))+
  scale_color_manual(values = c("#1B9E77", "#E6AB02", "#D95F02", "#7570B3", "#E7298A", "#66A61E"),
                     breaks = c("CLR", "CLR+VST", "log", "log-TSS", "SparXCC iterative", "SparXCC base"))+
  scale_fill_manual(values = c("#1B9E77", "#E6AB02", "#D95F02", "#7570B3", "#E7298A", "#66A61E"),
                    breaks = c("CLR", "CLR+VST", "log", "log-TSS", "SparXCC iterative", "SparXCC base"))+
  scale_linetype_manual("Correlations",values=c("Non-zero"=1,"Zero"=2))+
  labs(x = "Number of OTUs", y = "Mean Absolute Error")+
  facet_grid(dens ~ Cor_mat)+
  guides(linetype = guide_legend(order = 1))+
  theme(text = element_text(size = 10)) -> g_yes_main

# Supplementary figure without zero-inflation
ggplot(data = R1_no_main, aes(x = as.factor(p), col = Method, fill = Method, group = Method))+
  geom_line(linewidth = 1, aes(y = MAE_nonzero, linetype = "Non-zero"))+
  geom_line(linewidth = 1, aes(y = MAE_zero, linetype = "Zero"))+
  geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_nonzero_min, ymax = MAE_nonzero_max))+
  geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_zero_min, ymax = MAE_zero_max))+
  scale_color_manual(values = c("#1B9E77", "#E6AB02", "#D95F02", "#7570B3", "#E7298A", "#66A61E"),
                     breaks = c("CLR", "CLR+VST", "log", "log-TSS", "SparXCC iterative", "SparXCC base"))+
  scale_fill_manual(values = c("#1B9E77", "#E6AB02", "#D95F02", "#7570B3", "#E7298A", "#66A61E"),
                    breaks = c("CLR", "CLR+VST", "log", "log-TSS", "SparXCC iterative", "SparXCC base"))+
  scale_linetype_manual("Correlations",values=c("Non-zero"=1,"Zero"=2))+
  labs(x = "Number of OTUs", y = "Mean Absolute Error")+
  facet_grid(dens ~ Cor_mat)+
  guides(linetype = guide_legend(order = 1))+
  theme(text = element_text(size = 10)) -> g_no_main

# Supplementary figure with all densities
ggplot(data = R1_yes_supp, aes(x = as.factor(p), col = Method, fill = Method, group = Method))+
  geom_line(linewidth = 1, aes(y = MAE_nonzero, linetype = "Non-zero"))+
  geom_line(linewidth = 1, aes(y = MAE_zero, linetype = "Zero"))+
  # geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_nonzero_min, ymax = MAE_nonzero_max))+
  # geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_zero_min, ymax = MAE_zero_max))+
  scale_color_manual(values = c("#1B9E77", "#E6AB02", "#D95F02", "#7570B3", "#E7298A", "#66A61E"),
                     breaks = c("CLR", "CLR+VST", "log", "log-TSS", "SparXCC iterative", "SparXCC base"))+
  scale_fill_manual(values = c("#1B9E77", "#E6AB02", "#D95F02", "#7570B3", "#E7298A", "#66A61E"),
                    breaks = c("CLR", "CLR+VST", "log", "log-TSS", "SparXCC iterative", "SparXCC base"))+
  scale_linetype_manual("Correlations",values=c("Non-zero"=1,"Zero"=2))+
  labs(x = "Number of OTUs", y = "Mean Absolute Error")+
  facet_grid(dens  ~ q)+
  guides(linetype = guide_legend(order = 1))+
  theme(text = element_text(size = 10), legend.position = "bottom") -> g_supp

save(file = "../Temp/g_yes.txt", g_yes_main)
save(file = "../Temp/g_no.txt", g_no_main)

ggsave(filename = "../Supplementary figures/SFig4.pdf", g_supp, width = 190, height = 160, unit = "mm")
