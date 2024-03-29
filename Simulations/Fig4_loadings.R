# Code by Ib Thorsgaard Jensen

# Calling other files
source(file = "Shared_settings_Fig4.R")
source(file = "../cor_mat_loadings.R")

# Constructing matrix with settings
p <- c(10, 100, 1000)
q <- c(10, 100, 1000)

Opt <- expand.grid(p = p, q = q, bio_zero = bio_zero)
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
  n.sim <- 50
  n.seq.abn <- p*40
  
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
  cor_mat <- random_cor_mat(p+q, k = 5)
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
                      MAE = sapply(seq(0,0.8,0.1), function(x) mae(true.cor[abs(true.cor)>x], C[abs(true.cor)>x])),
                      Threshold = seq(0,0.8,0.1),
                      Method = "log")
  mae_CC <- data.table(Opt[i,],
                       MAE = sapply(seq(0,0.8,0.1), function(x) mae(true.cor[abs(true.cor)>x], CC[abs(true.cor)>x])),
                       Threshold = seq(0,0.8,0.1),
                       Method = "log-TSS")
  mae_CLR <- data.table(Opt[i,],
                        MAE = sapply(seq(0,0.8,0.1), function(x) mae(true.cor[abs(true.cor)>x], CLR[abs(true.cor)>x])),
                        Threshold = seq(0,0.8,0.1),
                        Method = "CLR")
  mae_mix <- data.table(Opt[i,],
                        MAE = sapply(seq(0,0.8,0.1), function(x) mae(true.cor[abs(true.cor)>x], Mix[abs(true.cor)>x])),
                        Threshold = seq(0,0.8,0.1),
                        Method = "CLR+VST")
  mae_XCC <- data.table(Opt[i,],
                        MAE = sapply(seq(0,0.8,0.1), function(x) mae(true.cor[abs(true.cor)>x], XCC[abs(true.cor)>x])),
                        Threshold = seq(0,0.8,0.1),
                        Method = "SparXCC")
  B <- rbind(mae_C, mae_CC, mae_CLR, mae_mix, mae_XCC)
}
# Saving results in a csv-file
fwrite(R, file = paste0("Results/Results_caseC_loadings_", as.character(Sys.Date()), ".csv"))

R3 <- R[,.(MAE = mean(MAE, na.rm = T),
           MAE_min = mean(MAE, na.rm = T) - 1.96*sd(MAE, na.rm = T)/sqrt(reps),
           MAE_max = mean(MAE, na.rm = T) + 1.96*sd(MAE, na.rm = T)/sqrt(reps)),
        list(p, q, Threshold, Method, bio_zero)]
# R3[,p:=factor(paste0("p=",p), levels = c("p=10", "p=100", "p=1000"))]
# R3[,Cor_mat := "Loadings method"]

R3_yes_main <- R3[bio_zero == "Yes" & p == q]
R3_supp <- R3[bio_zero == "Yes"]
R3_no_main <- R3[bio_zero == "No" & p == q]
# R1_no_supp <- R1[bio_zero == "No" & dens != 0.1]

R3_yes_main[,":="(Cor_mat = "Loadings method",
                  p=factor(paste0("p=q=",p), levels = c("p=q=10", "p=q=100", "p=q=1000")),
                  q=factor(paste0("p=q=",q), levels = c("p=q=10", "p=q=100", "p=q=1000")))]

R3_no_main[,":="(Cor_mat = "Loadings method",
                 p=factor(paste0("p=q=",p), levels = c("p=q=10", "p=q=100", "p=q=1000")),
                 q=factor(paste0("p=q=",q), levels = c("p=q=10", "p=q=100", "p=q=1000")))]

R3_supp[,":="(p=factor(paste0("p=",p), levels = c("p=10", "p=100", "p=1000")),
              q=factor(paste0("q=",q), levels = c("q=10", "q=100", "q=1000")))]

# Summarizing accuracy separated by the selected simulation settings
ggplot(data = R3_yes_main, aes(x = Threshold, y = MAE, ymin=MAE_min, ymax=MAE_max, col = Method, fill = Method, group = Method))+
  geom_line(size = 1, aes(linetype = "Non-zero"))+
  geom_line(size = 1, aes(linetype = "Zero"))+
  geom_ribbon(alpha=0.35, aes(col = NULL))+
  scale_color_manual("Method", values = cols)+
  scale_fill_manual("Method", values = cols)+
  labs(x = "Correlation threshold", y = "Mean Absolute Error")+
  facet_grid(p ~ Cor_mat)+
  scale_linetype_manual("Correlations",values=c("Non-zero"=1,"Zero"=2))+
  guides(linetype = guide_legend(order = 1))+
  theme(text = element_text(size = 10)) -> g2_yes_main

# Supplemntary figure without zero-inflation
ggplot(data = R3_no_main, aes(x = Threshold, y = MAE, ymin=MAE_min, ymax=MAE_max, col = Method, fill = Method, group = Method))+
  geom_line(size = 1, aes(linetype = "Non-zero"))+
  geom_line(size = 1, aes(linetype = "Zero"))+
  geom_ribbon(alpha=0.35, aes(col = NULL))+
  scale_color_manual("Method", values = cols)+
  scale_fill_manual("Method", values = cols)+
  labs(x = "Correlation threshold", y = "Mean Absolute Error")+
  facet_grid(p ~ Cor_mat)+
  scale_linetype_manual("Correlations",values=c("Non-zero"=1,"Zero"=2))+
  guides(linetype = guide_legend(order = 1))+
  theme(text = element_text(size = 10)) -> g2_no_main

# Supplementary figure with all tested combinations of p and q
ggplot(data = R3_supp, aes(x = Threshold, y = MAE, ymin=MAE_min, ymax=MAE_max, col = Method, fill = Method, group = Method))+
  geom_line(size = 1)+
  geom_ribbon(alpha=0.35, aes(col = NULL))+
  scale_color_manual("Method", values = cols)+
  scale_fill_manual("Method", values = cols)+
  labs(x = "Correlation threshold", y = "Mean Absolute Error")+
  facet_grid(p ~ q)+
  theme(text = element_text(size = 10), legend.position="bottom") -> g2_supp

save(file = "g2_yes.txt", g2_yes)
save(file = "g2_no.txt", g2_no)

ggsave(filename = "Figures/Supp_all_pq.pdf", g2_supp, width = 190, height = 160, unit = "mm")
