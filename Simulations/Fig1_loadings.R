# Code by Ib Thorsgaard Jensen

# Calling other files
source(file = "Shared_settings_Fig1.R")
source(file = "../cor_mat_loadings.R")

# Constructing matrix with settings
p <- c(10,100,1000)

Opt <- expand.grid(p = p, bio_zero = bio_zero)
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
  # Randomly selecting OTUs from the template
  OTU.idx <- sample(1:1000, Opt[i, "p"])
  template_sim <- template[OTU.idx]
  if(Opt[i, "bio_zero"] == "No bioloigcal zeros"){
    template_sim[,pi:=0]
  }
  # Defining other simulation settings
  n.sim <- 50
  n.seq.abn <- Opt[i, "p"]*40
  phenotype_var <- 1
  
  # Constructing correlation matrix
  cor_mat <- random_cor_mat(Opt[i, "p"]+1, 5)
  true.cor <- cor_mat[1:Opt[i, "p"],Opt[i, "p"]+1]
  
  # Simulating sequencing data 
  Abs_abundance_sim_B(N = n.sim, 
                      OTU_template = template_sim$mu,
                      OTU_vars = template_sim$sigma2,
                      phen_means = 30,
                      phen_var = phenotype_var,
                      cor_mat = cor_mat,
                      zero_props = template_sim$pi,
                      n.seq = n.seq.abn) -> A
  
  # Data transformations
  OTU_TSS <- (A$OTU_reads+1)/rowSums(A$OTU_reads+1)
  OTU_CLR <- t(apply(OTU_TSS, 1, clr))
  
  # Estimating correlations
  C <- cor(log(A$OTU_reads+1), A$Phenotype)[,1]
  CC <- cor(log(OTU_TSS), A$Phenotype)[,1]
  CLR <- cor(OTU_CLR, A$Phenotype)[,1]
  CEV <- SparCEV(A$OTU_reads, A$Phenotype); mae(CEV, true.cor)
  
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
  mae_CEV <- data.table(Opt[i,],
                        MAE = sapply(seq(0,0.8,0.1), function(x) mae(true.cor[abs(true.cor)>x], CEV[abs(true.cor)>x])),
                        Threshold = seq(0,0.8,0.1),
                        Method = "SparCEV")
  
  B <- rbind(mae_C, mae_CC, mae_CLR, mae_CEV)
}
# Saving results in a csv-file
fwrite(R, file = paste0("Results/Results_sim_loadings_", as.character(Sys.Date()), ".csv"))

# Summarizing accuracy separated by the selected simulation settings
R3 <- R[,.(MAE = mean(MAE, na.rm = T),
           MAE_min = mean(MAE, na.rm = T) - 1.96*sd(MAE, na.rm = T)/sqrt(reps),
           MAE_max = mean(MAE, na.rm = T) + 1.96*sd(MAE, na.rm = T)/sqrt(reps)),
        list(p, Threshold, Method, bio_zero)]
R3[,p:=factor(paste0("p=",p), levels = c("p=10", "p=100", "p=1000"))]
R3[,Cor_mat := "Loadings method"]

# Constructing and saving figures
ggplot(data = R3[bio_zero == "Yes"], aes(x = Threshold, y = MAE, ymin=MAE_min, ymax=MAE_max, col = Method, fill = Method, group = Method))+
  geom_line(size = 1, aes(linetype = "Non-zero"))+
  geom_line(size = 1, aes(linetype = "Zero"))+
  geom_ribbon(alpha=0.35, aes(col = NULL))+
  scale_color_brewer(" Method", palette = "Dark2")+
  scale_fill_brewer(" Method", palette = "Dark2")+
  labs(x = "Correlation threshold", y = "Mean Absolute Error")+
  facet_grid(p ~ Cor_mat)+
  scale_linetype_manual("Correlations",values=c("Non-zero"=1,"Zero"=2))+
  guides(linetype = guide_legend(order = 1))+
  theme(text = element_text(size = 10)) -> g2_yes

ggplot(data = R3[bio_zero == "No"], aes(x = Threshold, y = MAE, ymin=MAE_min, ymax=MAE_max, col = Method, fill = Method, group = Method))+
  geom_line(size = 1, aes(linetype = "Non-zero"))+
  geom_line(size = 1, aes(linetype = "Zero"))+
  geom_ribbon(alpha=0.35, aes(col = NULL))+
  scale_color_brewer(" Method", palette = "Dark2")+
  scale_fill_brewer(" Method", palette = "Dark2")+
  labs(x = "Correlation threshold", y = "Mean Absolute Error")+
  facet_grid(p ~ Cor_mat)+
  scale_linetype_manual("Correlations",values=c("Non-zero"=1,"Zero"=2))+
  guides(linetype = guide_legend(order = 1))+
  theme(text = element_text(size = 10)) -> g2_no

save(file = "g2_yes.txt", g2_yes)
save(file = "g2_no.txt", g2_no)
