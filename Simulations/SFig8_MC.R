# Code by Ib Thorsgaard Jensen

# Calling other files
source(file = "Shared_settings_Fig1.R") # This run uses the same settings as Fig1
source(file = "../cor_mat_loadings.R")

# Constructing matrix with settings
p <- c(10, 100, 1000)
bio_zero <- c("Yes")

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
  if(Opt[i, "bio_zero"] == "No"){
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
  
  # Estimating correlations
  CEV <- SparCEV(A$OTU_reads, A$Phenotype)
  CEV_MC <- SparCEV_MC(A$OTU_reads, A$Phenotype)
  
  # Setting up tables with results
  mae_CEV <- data.table(Opt[i,],
                        MAE = sapply(seq(0,0.8,0.1), function(x) mae(true.cor[abs(true.cor)>x], CEV[abs(true.cor)>x])),
                        Threshold = seq(0,0.8,0.1),
                        Method = "Pseudo-count")
  mae_CEV_MC <- data.table(Opt[i,],
                           MAE = sapply(seq(0,0.8,0.1), function(x) mae(true.cor[abs(true.cor)>x], CEV_MC[abs(true.cor)>x])),
                           Threshold = seq(0,0.8,0.1),
                           Method = "Monte-carlo")
  
  B <- rbind(mae_CEV, mae_CEV_MC)
}
# Saving results in a csv-file
fwrite(R, file = paste0("Results/Results_sim_MC_", as.character(Sys.Date()), ".csv"))

# Summarizing accuracy separated by the selected simulation settings
R2 <- R[,.(MAE = mean(MAE, na.rm = T),
           MAE_min = mean(MAE, na.rm = T) - 1.96*sd(MAE, na.rm = T)/sqrt(reps),
           MAE_max = mean(MAE, na.rm = T) + 1.96*sd(MAE, na.rm = T)/sqrt(reps)),
        list(p, Threshold, Method, bio_zero)]
R2[,p:=factor(paste0("p=",p), levels = c("p=10", "p=100", "p=1000"))]

# Constructing and saving figures
ggplot(data = R2, aes(x = Threshold, y = MAE, ymin=MAE_min, ymax=MAE_max, col = Method, fill = Method, group = Method))+
  geom_line(size = 1)+
  geom_ribbon(alpha=0.35, aes(col = NULL))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  labs(x = "Threshold", y = "MAE over threshold")+
  facet_grid(p~bio_zero)+
  theme(text = element_text(size = 10)) -> g

ggsave(filename = "Figures/SFig4.pdf", g, width = 150, height = 160, units = "mm")
