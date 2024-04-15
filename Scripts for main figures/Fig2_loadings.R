# Code by Ib Thorsgaard Jensen
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Calling other files
source(file = "../Functions/cor_mat_loadings.R")
source(file = "Shared_settings_Fig2.R")

reps <- 1000
cores <- 10

# Constructing matrix with settings
p_eff <- c(2,20,100)
Opt <- expand.grid(p_eff = p_eff)
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
  template_means <- template_table(Opt[i, "p_eff"])
  
  # Defining other simulation settings
  n.sim <- 50
  n.seq.abn <- p*20
  phenotype_var <- 1
  
  # Constructing correlation matrix
  cor_mat <- random_cor_mat(p+1, 2)
  true.cor <- cor_mat[1:p,p+1]
  
  # Simulating sequencing data 
  Abs_abundance_sim_B(N = n.sim, 
                      OTU_template = log(template_means) - 1/2,
                      OTU_vars = rep(1, p),
                      phen_means = 30,
                      phen_var = phenotype_var,
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
  mae_C <- data.table(p_eff = Opt[i,],
                      MAE = sapply(seq(0,0.8,0.1), function(x) mae(true.cor[abs(true.cor)>x], C[abs(true.cor)>x])),
                      Threshold = seq(0,0.8,0.1),
                      Method = "log")
  mae_CC <- data.table(p_eff = Opt[i,],
                       MAE = sapply(seq(0,0.8,0.1), function(x) mae(true.cor[abs(true.cor)>x], CC[abs(true.cor)>x])),
                       Threshold = seq(0,0.8,0.1),
                       Method = "log-TSS")
  mae_CLR <- data.table(p_eff = Opt[i,],
                        MAE = sapply(seq(0,0.8,0.1), function(x) mae(true.cor[abs(true.cor)>x], CLR[abs(true.cor)>x])),
                        Threshold = seq(0,0.8,0.1),
                        Method = "CLR")
  mae_CEV <- data.table(p_eff = Opt[i,],
                        MAE = sapply(seq(0,0.8,0.1), function(x) mae(true.cor[abs(true.cor)>x], CEV$cor[abs(true.cor)>x])),
                        Threshold = seq(0,0.8,0.1),
                        Method = "SparCEV base")
  if(all(is.na(CEV2$cor))){
    mae_CEV2 <- data.table(p_eff = Opt[i,],
                           MAE = rep(NA, 9),
                           Threshold = seq(0,0.8,0.1),
                           Method = "SparCEV iterative")
  }
  if(!all(is.na(CEV2$cor))){
    mae_CEV2 <- data.table(p_eff = Opt[i,],
                           MAE = sapply(seq(0,0.8,0.1),
                                        function(x) mae(true.cor[abs(true.cor)>x],
                                                        CEV2$cor[abs(true.cor)>x])),
                           Threshold = seq(0,0.8,0.1),
                           Method = "SparCEV iterative")
  }
  
  B <- rbind(mae_C, mae_CC, mae_CLR, mae_CEV, mae_CEV2)
}
# Saving results in a csv-file
fwrite(R, file = paste0("../Results/Results_eff_spc_loadings_", as.character(Sys.Date()), ".csv"))

# Summarizing accuracy separated by the selected simulation settings
R3 <- R[,.(MAE = mean(MAE, na.rm = T),
           MAE_min = mean(MAE, na.rm = T) - 1.96*sd(MAE, na.rm = T)/sqrt(reps),
           MAE_max = mean(MAE, na.rm = T) + 1.96*sd(MAE, na.rm = T)/sqrt(reps)),
        list(p_eff, Threshold, Method)]
R3[,p_eff:=factor(paste0("p_eff=",p_eff), levels = c("p_eff=2", "p_eff=20", "p_eff=100"))]
R3[,Cor_mat := "Loadings method"]

# Constructing and saving figures
ggplot(data = R3, aes(x = Threshold, y = MAE, ymin=MAE_min, ymax=MAE_max, col = Method, fill = Method, group = Method))+
  geom_line(linewidth = 1)+
  geom_ribbon(alpha=0.35, aes(col = NULL))+
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"),
                     breaks = c("CLR", "log", "log-TSS", "SparCEV iterative", "SparCEV base"))+
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"),
                    breaks = c("CLR", "log", "log-TSS", "SparCEV iterative", "SparCEV base"))+
  # scale_color_brewer(" Method", palette = "Dark2")+
  # scale_fill_brewer(" Method", palette = "Dark2")+
  labs(x = "Correlation threshold", y = "Mean Absolute Error")+
  facet_grid(p_eff ~ Cor_mat, scales = "free")+
  theme(text = element_text(size = 10)) -> g2

save(g2, file = "../Temp/g2.txt")
