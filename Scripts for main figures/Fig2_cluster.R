# Code by Ib Thorsgaard Jensen
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Calling other files
source(file = "Shared_settings_Fig2.R")
source(file = "../Functions/cor_mat_cluster.R")

reps <- 1000
cores <- 10

# Constructing matrix with settings
p_eff <- c(2,5,10,25,50,100)
dens <- c(0.1, 0.4, 0.7)
cor <- c(0.75)

Opt <- expand.grid(p_eff = p_eff, dens = dens, cor = cor)
Opt <- do.call("rbind", replicate(reps, Opt, simplify = FALSE))

# Constructs templates with effective number of species given by the vector p_eff
template_list <- lapply(as.list(p_eff), function(x) construct_template(sp_eff = x, p = p, bac.load = 1000))
names(template_list) <- p_eff

# Function that extracts a template from template_list given an input from p_eff
template_table <- function(p){
  idx <- which(names(template_list) == as.character(p))
  if(!(as.character(p) %in% names(template_list))){
    stop("Incorrect input")
  }
  return(template_list[[idx]])
}

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
  # Constructing mean-templates
  template_mean <- template_table(Opt[i, "p_eff"])
  
  # Defining other simulation settings
  n.sim <- 50
  n.seq.abn <- 4000
  phenotype_var <- 1

  # Constructing correlation matrix
  cor_mat <- cor_mat_sim_cluster(p = p,
                                 q = 1,
                                 cor.prob = Opt[i, "dens"],
                                 cor = Opt[i, "cor"],
                                 case = "B" )
  true.cor <- cor_mat[1:p,p+1]
  
  # Simulating sequencing data 
  Abs_abundance_sim_B(N = n.sim, 
                      OTU_template = log(template_mean) - 1/2,
                      OTU_vars = rep(1, p),
                      phen_means = log(30) - log(1+1/30^2)/2,
                      phen_var = log(1+1/30^2),
                      cor_mat = cor_mat,
                      zero_props = rep(0, p),
                      n.seq = n.seq.abn) -> A

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
# Saving results in a csv-file
fwrite(R, file = paste0("../Results/Results_eff_spc_cluster_", as.character(Sys.Date()), ".csv"))

# Summarizing accuracy separated by the selected simulation settings
R1 <- R[,.(MAE_nonzero = mean(MAE_nonzero, na.rm = T),
           MAE_nonzero_min = mean(MAE_nonzero, na.rm = T) - 1.96*sd(MAE_nonzero, na.rm = T)/sqrt(reps),
           MAE_nonzero_max = mean(MAE_nonzero, na.rm = T) + 1.96*sd(MAE_nonzero, na.rm = T)/sqrt(reps),
           MAE_zero = mean(MAE_zero, na.rm = T),
           MAE_zero_min = mean(MAE_zero, na.rm = T) - 1.96*sd(MAE_zero, na.rm = T)/sqrt(reps),
           MAE_zero_max = mean(MAE_zero, na.rm = T) + 1.96*sd(MAE_zero, na.rm = T)/sqrt(reps)),
        list(p_eff, dens, Method)]

R1[,dens:=paste0("c=", dens)]
R1[,Cor_mat := "Cluster method"]

# Constructing and saving figures
ggplot(data = R1, aes(x = as.factor(p_eff), y = MAE_nonzero, col = Method, fill = Method, group = Method))+
  geom_line(linewidth = 1)+
  geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_nonzero_min, ymax = MAE_nonzero_max))+
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"),
                     breaks = c("CLR", "log", "log-TSS", "SparCEV iterative", "SparCEV base"))+
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"),
                    breaks = c("CLR", "log", "log-TSS", "SparCEV iterative", "SparCEV base"))+
  # scale_color_brewer(" Method", palette = "Dark2")+
  # scale_fill_brewer(" Method", palette = "Dark2")+
  labs(x = "Effective number of OTUs", y = "Mean Absolute Error")+
  facet_grid(dens ~ Cor_mat, scales = "free")+
  theme(text = element_text(size = 10)) -> g

save(g, file = "../Temp/g.txt")

ggplot(data = R1, aes(x = as.factor(p_eff), y = MAE_zero, col = Method, fill = Method, group = Method))+
  geom_line(linewidth = 1)+
  geom_ribbon(alpha=0.35, aes(col = NULL, ymin = MAE_zero_min, ymax = MAE_zero_max))+
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"),
                     breaks = c("CLR", "log", "log-TSS", "SparCEV iterative", "SparCEV base"))+
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"),
                    breaks = c("CLR", "log", "log-TSS", "SparCEV iterative", "SparCEV base"))+
  # scale_color_brewer(" Method", palette = "Dark2")+
  # scale_fill_brewer(" Method", palette = "Dark2")+
  labs(x = "Effective number of OTUs", y = "Mean Absolute Error")+
  facet_grid(dens ~ Cor_mat, scales = "free")+
  theme(text = element_text(size = 10)) -> g_zero

ggsave(filename = "../Supplementary figures/SFig2_div.pdf", g_zero, width = 150, height = 160, unit = "mm")
