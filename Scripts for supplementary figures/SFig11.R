# Code by Ib Thorsgaard Jensen
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Loading packages
pkg <- c("data.table", "Matrix", "hydroGOF", "ggplot2", "ggpubr", "RColorBrewer", "magrittr",
         "compositions", "psych", "foreach", "parallel", "doParallel", "doRNG", "SpiecEasi", "CompoCor")
for(i in pkg){
  library(i, character.only = T)
}

# Calling other files
source(file = "../Functions/Data_generation.R")
source(file = "../Functions/cor_mat_cluster.R")

# Loading template
template_OTU <- fread("../Templates/OTU_template.csv")
template_gene <- fread("../Templates/Gene_template.csv")

# Constructing matrix with settings
# p <- c(100,250,1000)
p <- 100
q <- 100
n.sim = c(20, 50, 200, 1000)
dens <- c(0.1)
bio_zero <- c("No", "Yes")
cor <- c(0.75)

Opt <- expand.grid(p = p, q = q, n.sim = n.sim, dens = dens, bio_zero = bio_zero, cor = cor)
reps <- 50
Opt <- do.call("rbind", replicate(reps, Opt, simplify = FALSE))

a <- Sys.time()
R <- data.table(p = NA,
                q = NA,
                n.sim = NA,
                dens = NA,
                bio_zero = NA,
                cor = NA,
                Sensitivity = NA,
                FDR = NA,
                Method = NA)[-1]
set.seed(1658388669)
for(i in 1:nrow(Opt)){
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
  n.sim <- Opt[i, "n.sim"]
  n.seq.abn <- p*40
  n.seq.expr <- q*200
  
  # Constructing correlation matrix
  cor_mat <- cor_mat_sim_cluster(p,
                                 q,
                                 cor.prob = Opt[i,"dens"],
                                 cor = Opt[i,"cor"],
                                 case = "C")
  true.cor.mat <- cor_mat[1:p,p+1:q]
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
  OTU_CLR <- t(apply(A$OTU_reads+1, 1, clr))
  Gene_CLR <- t(apply(A$Gene_reads+1, 1, clr))
  
  # Estimating correlations and carrying out a t-test
  Pclr <- corr.test(x = OTU_CLR, y = Gene_CLR, adjust = "fdr", ci = FALSE)
  fdr_cor <- mean(true.cor[Pclr$p.adj<0.05] == 0)
  sens_cor <- mean(Pclr$p.adj[true.cor != 0] < 0.05)
  
  # Estimating correlations with SparXCC
  XCC <- SparXCC(OTU.abn = A$OTU_reads,
                 gene.expr = A$Gene_reads,
                 pseudo_count = 1,
                 iter = 10,
                 t1 = NULL,
                 t2 = NULL,
                 B_t = 20,
                 t1_quant = 0.8,
                 t2_quant = 0.8,
                 Find_m = T,
                 B_m = 20,
                 cores = 1)

  fdr_cor_thresh <- mean(true.cor.mat[XCC$Correlated] == 0)
  sens_cor_thresh <- mean(XCC$Correlated[true.cor.mat != 0])
  
  # Applying SPIEC-EASI
  S <- multi.spiec.easi(list(A$OTU_reads, A$Gene_reads))
  Con <- as.matrix(getRefit(S)[1:p, q+1:p]) == 1
  fdr_spiec <- mean(true.cor[Con == 1] == 0)
  sens_spiec <- mean(Con[true.cor != 0] == 1)
  
  # Setting up tables with results
  mae_Cor <- data.table(Opt[i,],
                        Sensitivity=sens_cor,
                        FDR=fdr_cor,
                        Method = "CLR (test)")
  mae_perm <- data.table(Opt[i,],
                         Sensitivity=sens_cor_thresh,
                         FDR=fdr_cor_thresh,
                         Method = "SparXCC (permute)")
  mae_spiec <- data.table(Opt[i,],
                          Sensitivity=sens_spiec,
                          FDR=fdr_spiec,
                          Method = "SPIEC-EASI")
  B <- rbind(mae_Cor,
             mae_spiec,
             mae_perm)
  
  R <- rbind(R, B)
  fwrite(R, "SFig9_res_temp.csv")
  
  rm(S, perm_thresh, XCC); gc()
}
Sys.time() - a
# Saving results in a csv-file
fwrite(R, file = "../Results/SparXCC_SPIEC_EASI.csv")

# Summarizing accuracy separated by the selected simulation settings
R2 <- R[,.(Sensitivity = mean(Sensitivity, na.rm = T),
           Sensitivity_min = mean(Sensitivity, na.rm = T) - 1.96*sd(Sensitivity, na.rm = T)/sqrt(reps),
           Sensitivity_max = mean(Sensitivity, na.rm = T) + 1.96*sd(Sensitivity, na.rm = T)/sqrt(reps),
           FDR = mean(FDR, na.rm = T),
           FDR_min = mean(FDR, na.rm = T) - 1.96*sd(FDR, na.rm = T)/sqrt(reps),
           FDR_max = mean(FDR, na.rm = T) + 1.96*sd(FDR, na.rm = T)/sqrt(reps)),
        list(p, q, dens, Method, n.sim, bio_zero)]

R4 <- rbind(R2[,c("p", "q", "n.sim", "Method", "bio_zero")], R2[,c("p", "q", "n.sim", "Method", "bio_zero")])

R4[,":="(MAE = c(R2$Sensitivity, R2$FDR),
         MAE_min = c(R2$Sensitivity_min, R2$FDR_min),
         MAE_max = c(R2$Sensitivity_max, R2$FDR_max),
         Metric = rep(c("Sensitivity", "False Discovery Rate"), each = nrow(R2))
)]
R4[,p:=paste0("p=", p) %>% factor(levels = c("p=100", "p=250", "p=1000"))]

R4[,line:=ifelse(Metric == "False Discovery Rate", 0.05, 0)]
R4[,linetype:=ifelse(Metric == "False Discovery Rate", "solid", "dashed")]

# Constructing and saving figures
ggplot(data = R4, aes(x = as.factor(n.sim), y = MAE, ymin=MAE_min, ymax=MAE_max, col = Method, fill = Method, group = Method))+
  geom_line(linewidth = 1)+
  geom_ribbon(alpha=0.35, aes(col = NULL))+
  # scale_color_brewer(palette = "Dark2")+
  scale_color_manual(limits = c("CLR (test)", "SparXCC (permute)","SPIEC-EASI"),
                     values = c("#1B9E77", "#E7298A", "#B15928"))+
  # scale_fill_brewer(palette = "Dark2")+
  scale_fill_manual(limits = c("CLR (test)", "SparXCC (permute)","SPIEC-EASI"),
                     values = c("#1B9E77", "#E7298A", "#B15928"))+
  labs(x = "Number of replicates", y = "Metric")+
  geom_hline(data = R4, aes(yintercept = line, linetype = linetype))+
  guides(color=guide_legend("Method"), linetype = "none")+
  # facet_wrap(~Metric, scales = "free")+
  # facet_grid(bio_zero~Metric, scales = "free")+
  ggh4x::facet_grid2(bio_zero~Metric, scales = "free_y", independent = "y")+
  theme(text = element_text(size = 10),
        legend.position="bottom") -> g; g

ggsave(filename = "../Supplementary figures/SFig11.pdf", g, width = 200, height = 160, unit = "mm")
