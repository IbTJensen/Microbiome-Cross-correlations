# Code by Ib Thorsgaard Jensen

# Loading packages
pkg <- c("data.table", "ggplot2", "ggpubr", "RColorBrewer", "magrittr",
         "limma", "network", "GGally", "compositions", "psych", "CompoCor")
for(i in pkg){
  library(i, character.only = T)
}
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading data
Bac <- fread("../Data/otu_table_bacteria.txt")
Fung <- fread("../Data/otu_table_fungi.txt")
Fung[,OTUid:=gsub("OTU", "OTUF", OTUid)]

# Loading metadata and taxonomy
meta_data <- fread("../Data/meta_data.csv")
taxonomy_Bac <- fread("../Data/taxonomy_bacteria.txt", fill = T)
taxonomy_Fung <- fread("../Data/taxonomy_fungi.txt", fill = T)
taxonomy_Fung[,V1:=gsub("OTU", "OTUF", V1)]

# Renaming samples, so they match between datasets
num <- paste0("_", 1:6)
num2 <- as.character(1:6)
for(i in 1:6){
  colnames(Fung)[grepl(num[i], colnames(Fung))] <- gsub(num[i],num2[i],colnames(Fung)[grepl(num[i], colnames(Fung))])
}

merge(meta_data,
      transpose(Bac, make.names = "OTUid", keep.names = "SampleID"),
      by = "SampleID") -> Bac_full

merge(meta_data,
      transpose(Fung, make.names = "OTUid", keep.names = "SampleID"),
      by = "SampleID") -> Fung_full

# Constructing subtable with only root samples from Gifu
Bac_gifu <- Bac_full[Compartment == "Root" & Genotype == "Gifu"]
Fung_gifu <- Fung_full[Compartment == "Root" & Genotype == "Gifu"]

# Computing the number of non-zero samples for each OTU
Bac_gifu[,lapply(.SD, function(x) sum(x!= 0)),
         .SDcols = colnames(Bac_full)[-(1:4)]] %>% 
  as.matrix() -> Bac_gifu_num

Fung_gifu[,lapply(.SD, function(x) sum(x!= 0)),
          .SDcols = colnames(Fung_gifu)[-(1:4)]] %>% 
  as.matrix() -> Fung_gifu_num

# Calculating relative abundances (i.e applying the TSS transformation)
RA_bac <- Bac_gifu[,-(1:4)]/rowSums(Bac_gifu[,-(1:4)])
RA_fung <- Fung_gifu[,-(1:4)]/rowSums(Fung_gifu[,-(1:4)])

# Keeping only OTUs with more than 0.1% RA and that are present in at least 10 samples in Gifu
Bac_keep <- colnames(Bac_gifu_num)[Bac_gifu_num>=10 & colMeans(RA_bac)>1e-3]
Fung_keep <- colnames(Fung_gifu_num)[Fung_gifu_num>=10 & colMeans(RA_fung)>1e-3]

Bac <- Bac[OTUid %in% Bac_keep]
Fung <- Fung[OTUid %in% Fung_keep]

# Genotypes
gt <- c("Gifu", "ccamk", "ram1", "symrk", "nfr5")

# Subsetting data for each genotype
Bac_subset <- list()
Fung_subset <- list()
meta_sub <- list()
for(i in 1:5){
  col_keep <- c("OTUid", meta_data[Compartment == "Root" & Genotype == gt[i], SampleID])
  meta_sub[[i]] <- meta_data[Compartment == "Root" & Genotype == gt[i]]
  Bac_subset[[i]] <- Bac[,..col_keep]
  Fung_subset[[i]] <- Fung[,..col_keep]
}
N <- unlist(lapply(meta_sub, nrow))

# Examining level of sparsity for each subset of the data
lapply(Bac_subset, function(x) mean(x==0))
lapply(Fung_subset, function(x) mean(x==0))

p <- unlist(lapply(Bac_subset, nrow))
q <- unlist(lapply(Fung_subset, nrow))

colorBlindBlack8  <- c("Actinobacteria" = "#D55E00",
                       "Ascomycota" = "#009E73",
                       "Proteobacteria" = "#0072B2",
                       "Glomeromycota" = "#56B4E9",
                       "Bacteroidetes" = "#E69F00",
                       "Other fungi" = "#000000",
                       "Other bacteria" = "#CC79A7")

# Constructing graph over cross-correlations for each genotype
chosen_thresholds <- c(0.8, 0.9, 1, 1, 1)

fig_dt_list <- list()
g <- list()
fun_taxa <- list()
bac_taxa <- list()
cors <- list()
set.seed(1671522881)
for(i in 1:5){
  cat(i, "\r")
  Bac_log <- log(Bac_subset[[i]][,-1]+1)
  Fung_log <- log(Fung_subset[[i]][,-1]+1)
  
  Bac_log <- removeBatchEffect(Bac_log, batch = meta_sub[[i]]$Experiment)
  Fung_log <- removeBatchEffect(Fung_log, batch = meta_sub[[i]]$Experiment)
  rownames(Bac_log) <- Bac_subset[[i]]$OTUid
  rownames(Fung_log) <- Fung_subset[[i]]$OTUid
  Bac_log <- exp(Bac_log); Fung_log <- exp(Fung_log)
  
  C <- SparXCC(OTU.abn = t(Bac_log),
               gene.expr = t(Fung_log),
               pseudo_count = 0,
               iter = 10,
               t1 = NULL,
               t2 = NULL,
               B_t = 100,
               t1_quant = 0.8,
               t2_quant = 0.8,
               Find_m = F,
               cores = 1)
  
  C2 <- SparXCC(OTU.abn = t(Bac_log),
                gene.expr = t(Fung_log),
                pseudo_count = 0,
                iter = 10,
                t1 = NULL,
                t2 = NULL,
                B_t = 100,
                t1_quant = chosen_thresholds[i],
                t2_quant = chosen_thresholds[i],
                Find_m = F,
                cores = 1)

  C3 <- SparXCC_base(OTU.abn = t(Bac_log), 
                     gene.expr = t(Fung_log),
                     pseudo_count = 0,
                     Find_m = F)
  
  dt <- data.table(SparXCC_iter = as.numeric(C$cor),
                   SparXCC_base = as.numeric(C3$cor),
                   Genotype = gt[i], 
                   Threshold = "Bootstrap 80th percentile",
                   a = mean(C$cor - C3$cor))
  
  dt2 <- data.table(SparXCC_iter = as.numeric(C2$cor),
                   SparXCC_base = as.numeric(C3$cor),
                   Genotype = gt[i], 
                   Threshold = "Manually chosen percentile",
                   a = mean(C2$cor - C3$cor))
  
  fig_dt_list[[i]] <- dt_comb <- rbind(dt, dt2)

}
fig_dt <- rbindlist(fig_dt_list)
fig_dt[,Genotype:=factor(Genotype, levels = gt)]

ggplot(data = fig_dt, aes(x = SparXCC_base, y = SparXCC_iter))+
  geom_point(shape = 21, color = "black", fill = "blue")+
  facet_grid(Genotype~Threshold)+
  geom_abline(aes(intercept = a, slope = 1), color = "red", linewidth = 1)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1)+
  labs(x = "SparXCC base", y = "SparXCC iterative") -> gg

ggsave("../Supplementary figures/SFig7.pdf", gg, width = 190, height = 240, units = "mm")
