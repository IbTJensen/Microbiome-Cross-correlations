# Code by Ib Thorsgaard Jensen
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Loading packages
pkg <- c("data.table", "Matrix", "hydroGOF", "ggplot2", "ggpubr", "RColorBrewer",
         "compositions", "psych", "foreach", "parallel", "doParallel", "doRNG", "CompoCor")
for(i in pkg){
  library(i, character.only = T)
}

# Calling other files
source(file = "../Functions/Data_generation.R")
source(file = "../Functions/cor_mat_cluster.R")

# Loading template
template <- fread("../Templates/OTU_template.csv")

# Constructing matrix with settings
p <- c(100, 250, 1000)
n.sim = c(20,50,200,1000)
dens <- c(0.1)
bio_zero <- c("Yes")
cor <- c(0.75)

Opt <- expand.grid(p = p, n.sim = n.sim, dens = dens, bio_zero = bio_zero, cor = cor)
reps <- 200
Opt <- do.call("rbind", replicate(reps, Opt, simplify = FALSE))

# Setting up parallelization
my.cluster <- parallel::makeCluster(
  10, 
  type = "FORK"
)

registerDoParallel(cl = my.cluster)
registerDoRNG(seed = 1676299396)

a <- Sys.time()
set.seed(1658388669)
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
  n.sim <- Opt[i, "n.sim"]
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
  
  # Data transformations
  OTU_CLR <- t(apply(A$OTU_reads+1, 1, clr))
  
  # Estimating correlations and carrying out a t-test
  Pclr <- corr.test(x = OTU_CLR, y = A$Phenotype, adjust = "fdr", ci = FALSE)
  cor_est <- SparCEV(OTU.abn = A$OTU_reads, 
                     phenotype = A$Phenotype,
                     pseudo_count = 1,
                     iter = 20,
                     t = NULL,
                     B_t = 20,
                     t_quant = 0.8,
                     Find_m = T,
                     B_m = 20,
                     cores = 1)
  
  fdr_cor_thresh <- mean(true.cor[cor_est$Correlated] == 0)
  sens_cor_thresh <- mean(cor_est$Correlated[true.cor != 0])
    
  p_adj <- as.numeric(Pclr$p.adj)
  fdr_cor <- mean(true.cor[p_adj<0.05] == 0)
  sens_cor <- mean(p_adj[true.cor != 0] < 0.05)
  
  # Setting up tables with results
  mae_Cor <- data.table(Opt[i,],
                        Sensitivity=sens_cor,
                        FDR=fdr_cor,
                        Method = "CLR (t-test)")
  mae_cons <- data.table(Opt[i,],
                         Sensitivity=sens_cor_thresh,
                         FDR=fdr_cor_thresh,
                         Method = "SparCEV (permutation)")
  B <- rbind(mae_Cor, mae_cons)
}
Sys.time() - a
# Saving results in a csv-file
R[is.na(FDR), FDR:=0]
fwrite(R, file = paste0("../Results/Results_caseB_testing_", as.character(Sys.Date()), ".csv"))

# Summarizing accuracy separated by the selected simulation settings
R2 <- R[,.(Sensitivity = mean(Sensitivity, na.rm = T),
           Sensitivity_min = mean(Sensitivity, na.rm = T) - 1.96*sd(Sensitivity, na.rm = T)/sqrt(reps),
           Sensitivity_max = mean(Sensitivity, na.rm = T) + 1.96*sd(Sensitivity, na.rm = T)/sqrt(reps),
           FDR = mean(FDR, na.rm = T),
           FDR_min = mean(FDR, na.rm = T) - 1.96*sd(FDR, na.rm = T)/sqrt(reps),
           FDR_max = mean(FDR, na.rm = T) + 1.96*sd(FDR, na.rm = T)/sqrt(reps)),
        list(p, dens, Method, n.sim)]

R4 <- rbind(R2[,c("p", "n.sim", "Method")], R2[,c("p",  "n.sim", "Method")])

R4[,":="(MAE = c(R2$Sensitivity, R2$FDR),
         MAE_min = c(R2$Sensitivity_min, R2$FDR_min),
         MAE_max = c(R2$Sensitivity_max, R2$FDR_max),
         Metric = rep(c("Sensitivity", "False Discovery Rate"), each = nrow(R2))
)]
# R4[,p:=paste0("p=", p) %>% factor(levels = c("p=100", "p=250", "p=1000"))]

R4[,line:=ifelse(Metric == "False Discovery Rate", 0.05, 0)]
R4[,linetype:=ifelse(Metric == "False Discovery Rate", "solid", "dashed")]

R4[,p:=paste("p =", p)]
R4[,p:=factor(p, levels = unique(p))]

# Constructing and saving figures
ggplot(data = R4, aes(x = as.factor(n.sim), y = MAE, ymin=MAE_min, ymax=MAE_max, col = Method, fill = Method, group = Method))+
  geom_line(size = 1)+
  geom_ribbon(alpha=0.35, aes(col = NULL))+
  scale_colour_manual(breaks = c("CLR (t-test)", "SparCEV (permutation)"),
                      values = c("#1B9E77", "#E7298A"))+
  scale_fill_manual(breaks = c("CLR (t-test)", "SparCEV (permutation)"),
                    values = c("#1B9E77", "#E7298A"))+
  labs(x = "Number of replicates", y = "Metric")+
  geom_hline(data = R4, aes(yintercept = line, linetype = linetype))+
  guides(color=guide_legend("Method"), linetype = "none")+
  facet_grid(p~Metric)+
  # ggh4x::facet_grid2(cor~Metric, labeller = label_bquote(rho == .(cor)), scales = "free_y", independent = "y")+
  # ggh4x::facet_grid2(p ~ Metric, scales = "free_y", independent = "y")+
  theme(text = element_text(size = 10),
        legend.position="bottom") -> g

ggsave(filename = "../Supplementary figures/SFig12.pdf", g, width = 160, height = 160, unit = "mm")
