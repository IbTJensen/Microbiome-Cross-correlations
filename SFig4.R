# Code by Ib Thorsgaard Jensen

# Loading packages
pkg <- c("data.table", "ggplot2", "ggpubr", "RColorBrewer", "magrittr",
         "limma", "network", "GGally", "compositions", "psych")
for(i in pkg){
  library(i, character.only = T)
}
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading data
Bac <- fread("Data/otu_table_bacteria.txt")
Fung <- fread("Data/otu_table_fungi.txt")
Fung[,OTUid:=gsub("OTU", "OTUF", OTUid)]
meta_data <- fread("Data/meta_data.csv")

# Renaming data so they match between datasets
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

# Keeping only root samples from Gifu
Fung_gifu <- Fung_full[Compartment == "Root" & Genotype == "Gifu"]
Bac_gifu <- Bac_full[Compartment == "Root" & Genotype == "Gifu"]

# Computing the number of non-zero samples for each OTU
Bac_full[Compartment == "Root" & Genotype == "Gifu",
         lapply(.SD, function(x) sum(x!= 0)),
         .SDcols = colnames(Bac_full)[-(1:4)]] %>% 
  as.matrix() -> Bac_gifu_num

Fung_full[Compartment == "Root" & Genotype == "Gifu",
          lapply(.SD, function(x) sum(x!= 0)),
          .SDcols = colnames(Fung_full)[-(1:4)]] %>% 
  as.matrix() -> Fung_gifu_num

# Calculating relative abundances (i.e applying the TSS transformation)
Lib_size_bac <- rowSums(Bac_full[,-(1:4)])
Lib_size_fung <- rowSums(Fung_full[,-(1:4)])

RA_bac <- Bac_full[,-(1:4)]/Lib_size_bac
RA_fung <- Fung_full[,-(1:4)]/Lib_size_fung

# Keeping only root samples from Gifu
RA_bac_gifu <- RA_bac[Bac_full$Compartment == "Root" & Fung_full$Genotype == "Gifu"]
RA_fung_gifu <- RA_fung[Fung_full$Compartment == "Root" & Fung_full$Genotype == "Gifu"]

# Keeping only OTUs with more than 0.1% RA and that are present in at least 10 samples
common_bac <- colnames(RA_bac_gifu)[colMeans(RA_bac_gifu)>1e-3 & colSums(RA_bac_gifu != 0)>= 10]
common_fung <- colnames(RA_fung_gifu)[colMeans(RA_fung_gifu)>1e-3  & colSums(RA_fung_gifu != 0)>= 10]

Bac <- Bac[OTUid %in% common_bac]
Fung <- Fung[OTUid %in% common_fung]

RA_bac_gifu_filtered <- RA_bac_gifu[,..common_bac]
RA_fung_gifu_filtered <- RA_fung_gifu[,..common_fung]

# CLR transformation
Bac_CLR <- apply(Bac[,-1]+1, 1, clr); colnames(Bac_CLR) <- Bac$OTUid
Fung_CLR <- apply(Fung[,-1]+1, 1, clr); colnames(Fung_CLR) <- Fung$OTUid

# Keeping only the root-Gifu samples for CLR transformed data
gifu_smp <- rownames(Bac_CLR)[grepl("GR", rownames(Bac_CLR))]
gifu_smp <- gifu_smp[!grepl("Rh", gifu_smp)]
Bac_CLR <- Bac_CLR[gifu_smp,]
Fung_CLR <- Fung_CLR[gifu_smp,]

# Recreating correlations from Thiergart et al.
CC_spear <- corr.test(x = as.matrix(RA_bac_gifu_filtered),
                      y = as.matrix(RA_fung_gifu_filtered),
                      method = "spearman",
                      ci = F,
                      adjust = "none",
                      alpha = 0.001)
sum(CC_spear$p<0.001)

# CLR-Pearson correlations
CC_CLR <- corr.test(x = as.matrix(Bac_CLR),
                    y = as.matrix(Fung_CLR),
                    method = "pearson",
                    ci = F,
                    adjust = "none",
                    alpha = 0.001)
sum(CC_CLR$p<0.001)

# CLR-Spearman correlations
CC_CLR_spear <- corr.test(x = as.matrix(Bac_CLR),
                          y = as.matrix(Fung_CLR),
                          method = "spearman",
                          ci = F,
                          adjust = "none",
                          alpha = 0.001)
sum(CC_CLR_spear$p<0.001)

# Proportion of Thiergart et. al correlations also found found by CLR-Pearson
sum(CC_CLR_spear$p<0.001 & CC_CLR$p<0.001)/sum(CC_CLR_spear$p<0.001)
# Proportion of Thiergart et. al correlations also found found by CLR-Spearman
sum(CC_CLR_spear$p<0.001 & CC_spear$p<0.001)/sum(CC_CLR_spear$p<0.001)

# Permutation threshold
N <- nrow(Bac_CLR)
p <- ncol(Bac_CLR)
q <- ncol(Fung_CLR)
set.seed(1678191958)
cor_thresh <- sapply(1:50, function(x) cor(Bac_CLR[sample(1:N),], Fung_CLR[sample(1:N),]) %>% abs() %>% quantile(1 - 1/(p*q))) %>% mean()
sum(abs(CC_CLR$r) > cor_thresh)
# Proportion of Thiergart et. al correlations also found found by CLR-Pearson with permutation thresholding
sum(CC_CLR_spear$p<0.001 & abs(CC_CLR$r) > cor_thresh)/sum(CC_CLR_spear$p<0.001)

# Collecting results from Thiergart et al. and CLR-Pearson into one data table
# NOTE: This may only work with data.table version 1.14.6 or earlier
long_CLR <- melt(CC_CLR$r, value.name = "CLR_cor")
long_CLR_p <- melt(CC_CLR$p<0.001, value.name = "CLR_cor")
long_spear <- melt(CC_spear$r, value.name = "spear_cor")
long_spear_p <- melt(CC_spear$p<0.001, value.name = "spear_cor")
long_CLR <- as.data.table(long_CLR)
long_CLR[,":="(spear_cor = long_spear$spear_cor,
               CLR_sig = long_CLR_p$CLR_cor,
               Spearman_sig = long_spear_p$spear_cor)]
long_CLR[,Correlated:=fcase(
  CLR_sig == T & Spearman_sig == T, "Both",
  CLR_sig == T & Spearman_sig == F, "CLR",
  CLR_sig == F & Spearman_sig == T, "Spearman",
  CLR_sig == F & Spearman_sig == F, "None"
)]

# Constructing figure
ggplot(data = long_CLR, aes(x = spear_cor, y = CLR_cor, col = Correlated))+
  geom_point()+
  scale_color_manual(limits = c("Spearman", "CLR", "Both", "None"),
                     values = c("#6A3D9A", "#33A02C", "#1F78B4", "#A6CEE3"),
                     name = "Significant")+
  xlim(-1,1)+
  ylim(-1,1)+
  labs(x = "TSS-Spearman", y = "CLR-Pearson")+
  geom_abline(slope = 1, intercept = 0)+
  NULL -> g

ggsave("Results and Figures/Spearman_vs_CLR.pdf", g, width = 200, height = 120, units = "mm")

# Maintaining only relevant information and saving
long_CLR <- long_CLR[Correlated != "None"]
long_CLR[,Correlated:=NULL]
colnames(long_CLR)[1:2] <- c("Bacterial OTU", "Fungal OTU")

fwrite(long_CLR, "Results and Figures/S4 Table.csv")
