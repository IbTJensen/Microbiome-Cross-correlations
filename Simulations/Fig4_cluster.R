# Code by Ib Thorsgaard Jensen

# Calling other files
source(file = "Shared_settings_Fig4.R")
source(file = "../cor_mat_cluster.R")

# Constructing matrix with settings
p <- c(10,20,50,100,250,500,1000)
q <- c(1000)
dens <- c(0.1, 0.4, 0.7)
cor <- c(0.75)

Opt <- expand.grid(p = p, q = q, dens = dens, bio_zero = bio_zero, cor = cor)
Opt <- do.call("rbind", replicate(reps, Opt, simplify = FALSE))

# Setting up parallelization
my.cluster <- parallel::makeCluster(
  cores, 
  type = "FORK"
)

registerDoParallel(cl = my.cluster)
registerDoRNG(seed = 1676299396)

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
  true.cor <- as.numeric(cor_mat[1:p,p+1:q])
  
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
  XCC <- as.numeric(SparXCC(A$OTU_reads, A$Gene_reads))
  
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
                        Method = "SparXCC")
  B <- rbind(mae_C, mae_CC, mae_CLR, mae_mix, mae_XCC)
}
# Saving results in a csv-file
fwrite(R, file = paste0("Results/Results_caseC_cluster_dense_", as.character(Sys.Date()), ".csv"))

# Summarizing accuracy separated by the selected simulation settings
R1 <- R[,.(MAE_nonzero = mean(MAE_nonzero, na.rm = T),
           MAE_nonzero_min = mean(MAE_nonzero, na.rm = T) - 1.96*sd(MAE_nonzero, na.rm = T)/sqrt(reps),
           MAE_nonzero_max = mean(MAE_nonzero, na.rm = T) + 1.96*sd(MAE_nonzero, na.rm = T)/sqrt(reps),
           MAE_zero = mean(MAE_zero, na.rm = T),
           MAE_zero_min = mean(MAE_zero, na.rm = T) - 1.96*sd(MAE_zero, na.rm = T)/sqrt(reps),
           MAE_zero_max = mean(MAE_zero, na.rm = T) + 1.96*sd(MAE_zero, na.rm = T)/sqrt(reps)),
        list(p, q, dens, Method, bio_zero)]

R1[,dens:=paste0("c=", dens)]
R1[,Cor_mat := "Cluster method"]

# Constructing and saving figures
ggplot(data = R1[bio_zero == "Yes"], aes(x = as.factor(p), col = Method, fill = Method, group = Method))+
  geom_line(size = 1, aes(y = MAE_nonzero, linetype = "Non-zero"))+
  geom_line(size = 1, aes(y = MAE_zero, linetype = "Zero"))+
  geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_nonzero_min, ymax = MAE_nonzero_max))+
  geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_zero_min, ymax = MAE_zero_max))+
  scale_color_manual("Method" ,values = cols)+
  scale_fill_manual("Method", values = cols)+
  scale_linetype_manual("Correlations",values=c("Non-zero"=1,"Zero"=2))+
  labs(x = "Number of OTUs", y = "Mean Absolute Error")+
  facet_grid(dens ~ Cor_mat)+
  guides(linetype = guide_legend(order = 1))+
  theme(text = element_text(size = 10)) -> g_yes

ggplot(data = R1[bio_zero == "No"], aes(x = as.factor(p), col = Method, fill = Method, group = Method))+
  geom_line(size = 1, aes(y = MAE_nonzero, linetype = "Non-zero"))+
  geom_line(size = 1, aes(y = MAE_zero, linetype = "Zero"))+
  geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_nonzero_min, ymax = MAE_nonzero_max))+
  geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_zero_min, ymax = MAE_zero_max))+
  scale_color_manual("Method", values = cols)+
  scale_fill_manual("Method", values = cols)+
  scale_linetype_manual("Correlations",values=c("Non-zero"=1,"Zero"=2))+
  labs(x = "Number of OTUs", y = "Mean Absolute Error")+
  facet_grid(dens ~ Cor_mat)+
  guides(linetype = guide_legend(order = 1))+
  theme(text = element_text(size = 10)) -> g_no

save(file = "g_yes.txt", g_yes)
save(file = "g_no.txt", g_no)

# R2 <- R[,.(MAE_nonzero = mean(MAE_nonzero, na.rm = T),
#            MAE_nonzero_min = mean(MAE_nonzero, na.rm = T) - 1.96*sd(MAE_nonzero, na.rm = T)/sqrt(reps),
#            MAE_nonzero_max = mean(MAE_nonzero, na.rm = T) + 1.96*sd(MAE_nonzero, na.rm = T)/sqrt(reps),
#            MAE_zero = mean(MAE_zero, na.rm = T),
#            MAE_zero_min = mean(MAE_zero, na.rm = T) - 1.96*sd(MAE_zero, na.rm = T)/sqrt(reps),
#            MAE_zero_max = mean(MAE_zero, na.rm = T) + 1.96*sd(MAE_zero, na.rm = T)/sqrt(reps)),
#         list(p, q, dens, Method, bio_zero)]
# 
# R3 <- rbind(R2[,c("p", "q", "dens", "Method")], R2[,c("p", "q", "dens", "Method")])
# R3[,":="(MAE = c(R2$MAE_nonzero, R2$MAE_zero),
#          MAE_min = c(R2$MAE_nonzero_min, R2$MAE_zero_min),
#          MAE_max = c(R2$MAE_nonzero_max, R2$MAE_zero_max),
#          Cor = rep(c("Non-zero correlation", "Zero correlations"), each = nrow(R2))
# )]
# 
# R3[,dens:=paste0("c=",dens)]
# R3[,Method:=factor(Method, levels = c("CLR","None","SparXCC","TSS", "CLR+VST"))]
# cols <- brewer.pal(n = 8, "Dark2")[-5]
# 
# ggplot(data = R3, aes(x = as.factor(p), y = MAE, ymin=MAE_min, ymax=MAE_max, col = Method, fill = Method, group = Method))+
#   geom_line(size = 1)+
#   geom_ribbon(alpha=0.35, aes(col = NULL))+
#   scale_color_manual(values = cols)+
#   scale_fill_manual(values = cols)+
#   labs(x = "Number of OTUs", y = "Mean Absolute Error")+
#   facet_grid(Cor ~ dens)+
#   theme(text = element_text(size = 20))-> g
# 
# ggsave(filename = paste0("CaseC_cluster_sparse_", as.character(Sys.Date()),".pdf"), g, width = 300, height = 180, units = "mm")
