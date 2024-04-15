# Code by Ib Thorsgaard Jensen
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Calling other files
source(file = "Shared_settings_Fig1.R")
source(file = "../Functions/cor_mat_cluster.R")

reps <- 1000
cores <- 12

# Constructing matrix with settings
p <- c(10,20,50,100,250,1000)
dens <- c(0.1, 0.4, 0.7)
cor <- c(0.75)

Opt <- expand.grid(p = p, dens = dens, bio_zero = bio_zero, cor = cor)
Opt <- do.call("rbind", replicate(reps, Opt, simplify = FALSE))

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
  cor_mat <- cor_mat_sim_cluster(p = Opt[i, "p"],
                                 q = 1,
                                 cor.prob = Opt[i, "dens"],
                                 cor = Opt[i, "cor"],
                                 case = "B" )
  true.cor <- cor_mat[1:Opt[i, "p"],Opt[i, "p"]+1]
  
  # Simulating sequencing data 
  Abs_abundance_sim_B(N = n.sim, 
                      OTU_template = template_sim$mu,
                      OTU_vars = template_sim$sigma2,
                      phen_means = log(30) - log(1+1/30^2)/2,
                      phen_var = log(1+1/30^2),
                      cor_mat = cor_mat,
                      zero_props = template_sim$pi,
                      n.seq = n.seq.abn) -> A
  
  colnames(A$OTU_reads) <- paste0("OTU", 1:Opt[i,"p"])
  # true.cov <- true.cor*sqrt(log(30) - log(1+1/30^2)/2)*sqrt(template_sim$sigma2)
  
  # Data transformations
  OTU_TSS <- (A$OTU_reads+1)/rowSums(A$OTU_reads+1)
  OTU_CLR <- t(apply(OTU_TSS, 1, clr))
  
  # Estimating correlations
  C <- cor(log(A$OTU_reads+1), A$Phenotype)[,1]
  CC <- cor(log(OTU_TSS), A$Phenotype)[,1]
  CLR <- cor(OTU_CLR, A$Phenotype)[,1]
  CEV <- SparCEV_base(OTU.abn = A$OTU_reads,
                      phenotype = A$Phenotype,
                      Find_m = F)
  
  CEV2 <- SparCEV(OTU.abn = A$OTU_reads, 
                  phenotype = A$Phenotype,
                  pseudo_count = 1,
                  B_t = 20,
                  t = NULL,
                  cores = 1,
                  t_quant = 0.8,
                  Find_m = F)
  
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
  mae_CEV <- data.table(Opt[i,],
                        MAE_nonzero=mae(true.cor[true.cor!=0], CEV$cor[true.cor!=0]),
                        MAE_zero=mae(true.cor[true.cor==0], CEV$cor[true.cor==0]),
                        Method = "SparCEV base")
  if(all(is.na(CEV2$cor))){
    mae_CEV2 <- data.table(Opt[i,],
                           MAE_nonzero=NA,
                           MAE_zero=NA,
                           Method = "SparCEV iterative")
  }
  if(!all(is.na(CEV2$cor))){
    mae_CEV2 <- data.table(Opt[i,],
                           MAE_nonzero=mae(true.cor[true.cor!=0], CEV2$cor[true.cor!=0]),
                           MAE_zero=mae(true.cor[true.cor==0], CEV2$cor[true.cor==0]),
                           Method = "SparCEV iterative")
  }
  
  B <- rbind(mae_C, mae_CC, mae_CLR, mae_CEV, mae_CEV2)
}
Sys.time() - a
# Saving results in a csv-file
fwrite(R, file = paste0("../Results/Results_sim_cluster_", as.character(Sys.Date()), ".csv"))

# Summarizing accuracy separated by the selected simulation settings
R1 <- R[,.(MAE_nonzero = mean(MAE_nonzero, na.rm = T),
           MAE_nonzero_min = mean(MAE_nonzero, na.rm = T) - 1.96*sd(MAE_nonzero, na.rm = T)/sqrt(reps),
           MAE_nonzero_max = mean(MAE_nonzero, na.rm = T) + 1.96*sd(MAE_nonzero, na.rm = T)/sqrt(reps),
           MAE_zero = mean(MAE_zero, na.rm = T),
           MAE_zero_min = mean(MAE_zero, na.rm = T) - 1.96*sd(MAE_zero, na.rm = T)/sqrt(reps),
           MAE_zero_max = mean(MAE_zero, na.rm = T) + 1.96*sd(MAE_zero, na.rm = T)/sqrt(reps)),
        list(p, dens, Method, bio_zero)]

R1[,dens:=paste0("c=", dens)]
R1[,Cor_mat := "Cluster method"]

# Constructing and saving figures
ggplot(data = R1[bio_zero == "Yes"], aes(x = as.factor(p), col = Method, fill = Method, group = Method))+
  geom_line(linewidth = 1, aes(y = MAE_nonzero, linetype = "Non-zero"))+
  geom_line(linewidth = 1, aes(y = MAE_zero, linetype = "Zero"))+
  geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_nonzero_min, ymax = MAE_nonzero_max))+
  geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_zero_min, ymax = MAE_zero_max))+
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"),
                     breaks = c("CLR", "log", "log-TSS", "SparCEV iterative", "SparCEV base"))+
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"),
                    breaks = c("CLR", "log", "log-TSS", "SparCEV iterative", "SparCEV base"))+
  # scale_color_brewer(" Method", palette = "Dark2")+
  # scale_fill_brewer(" Method", palette = "Dark2")+
  scale_linetype_manual("Correlations",values=c("Non-zero"=1,"Zero"=2))+
  labs(x = "Number of OTUs", y = "Mean Absolute Error")+
  facet_grid(dens ~ Cor_mat)+
  guides(linetype = guide_legend(order = 1))+
  theme(text = element_text(size = 10)) -> g_yes

ggplot(data = R1[bio_zero == "No"], aes(x = as.factor(p), col = Method, fill = Method, group = Method))+
  geom_line(linewidth = 1, aes(y = MAE_nonzero, linetype = "Non-zero"))+
  geom_line(linewidth = 1, aes(y = MAE_zero, linetype = "Zero"))+
  geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_nonzero_min, ymax = MAE_nonzero_max))+
  geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_zero_min, ymax = MAE_zero_max))+
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"),
                     breaks = c("CLR", "log", "log-TSS", "SparCEV iterative", "SparCEV base"))+
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"),
                    breaks = c("CLR", "log", "log-TSS", "SparCEV iterative", "SparCEV base"))+
  # scale_color_brewer(" Method", palette = "Dark2")+
  # scale_fill_brewer(" Method", palette = "Dark2")+
  scale_linetype_manual("Correlations",values=c("Non-zero"=1,"Zero"=2))+
  labs(x = "Number of OTUs", y = "Mean Absolute Error")+
  facet_grid(dens ~ Cor_mat)+
  guides(linetype = guide_legend(order = 1))+
  theme(text = element_text(size = 10)) -> g_no

save(file = "../Temp/g_yes.txt", g_yes)
save(file = "../Temp/g_no.txt", g_no)
