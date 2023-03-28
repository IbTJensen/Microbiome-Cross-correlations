# Code by Ib Thorsgaard Jensen

# Loading packages
pkg <- c("data.table", "ggplot2", "ggpubr", "RColorBrewer", "magrittr",
         "limma", "network", "GGally", "compositions", "psych")
for(i in pkg){
  library(i, character.only = T)
}
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("SparXCC_function.R")

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

# RA with a pseudo-count
RA_bac_pc <- t(apply(Bac_full[Compartment == "Root" & Genotype == "Gifu",-(1:4)]+1, 1, function(x) x/sum(x)))[,common_bac]
RA_fung_pc <- t(apply(Fung_full[Compartment == "Root" & Genotype == "Gifu",-(1:4)]+1, 1, function(x) x/sum(x)))[,common_fung]

# RA without pseudo-count
RA_bac_gifu_filtered <- RA_bac_gifu[,..common_bac]
RA_fung_gifu_filtered <- RA_fung_gifu[,..common_fung]

# Keeping only the root-Gifu samples for CLR transformed data
gifu_smp <- colnames(Bac)[grepl("GR", colnames(Bac))]
gifu_smp <- gifu_smp[!grepl("Rh", gifu_smp)]
Bac1 <- Bac[,..gifu_smp]
Fung1 <- Fung[,..gifu_smp]

# Recreating correlations from Thiergart et al. (Note that log-transformation will 
# not have any effect on Spearman correlations, since it is does not affect ranks)
CC_spear <- corr.test(x = as.matrix(RA_bac_gifu_filtered),
                      y = as.matrix(RA_fung_gifu_filtered),
                      method = "spearman",
                      ci = F,
                      adjust = "none",
                      alpha = 0.001)
sum(CC_spear$p<0.001)

# Pearson correlations on RA
CC_pears <- corr.test(x = as.matrix(log(RA_bac_pc)),
                      y = as.matrix(log(RA_fung_pc)),
                      method = "pearson",
                      ci = F,
                      adjust = "none",
                      alpha = 0.001)
sum(CC_pears$p<0.001)

# SparXCC correlations
CC_XCC <- SparXCC(t(Bac1), t(Fung1))
colnames(CC_XCC) <- Fung$OTUid
rownames(CC_XCC) <- Bac$OTUid

# Permutation threshold
N <- ncol(Bac1)
p <- nrow(Bac1)
q <- nrow(Fung1)
set.seed(1678191958)
cor_thresh <- sapply(1:20, function(x) SparXCC(t(Bac1)[sample(1:N),],
                                               t(Fung1)[sample(1:N),] ) %>% abs() %>% quantile(1 - 1/(p*q))) %>% mean()
sum(abs(CC_XCC) > cor_thresh)

# log-TSS Pearson vs TSS Spearman
sum(CC_spear$p<0.001 & CC_pears$p<0.001)/sum(CC_spear$p<0.001)
# Proportion of Thiergart et. al correlations also found found by SparXCC
sum(CC_spear$p<0.001 & abs(CC_XCC) > cor_thresh)/sum(CC_spear$p<0.001)
# Proportion of TSS-pearson correlations also found found by CLR-Spearman
sum(CC_pears$p<0.001 & abs(CC_XCC) > cor_thresh)/sum(abs(CC_XCC) > cor_thresh)

# Collecting results from Thiergart et al. and CLR-Pearson into one data table
# NOTE: This may only work with data.table version 1.14.6 or earlier
long_CLR <- melt(CC_XCC, value.name = "SparXCC")
long_CLR_p <- melt(abs(CC_XCC)>cor_thresh, value.name = "SparXCC")
long_spear <- melt(CC_spear$r, value.name = "spear_cor")
long_spear_p <- melt(CC_spear$p<0.001, value.name = "spear_cor")
long_CLR <- as.data.table(long_CLR)
long_CLR[,":="(spear_cor = long_spear$spear_cor,
               SparXCC_above_threshold = long_CLR_p$SparXCC,
               Spearman_sig = long_spear_p$spear_cor)]
long_CLR[,Correlated:=fcase(
  SparXCC_above_threshold == T & Spearman_sig == T, "Both",
  SparXCC_above_threshold == T & Spearman_sig == F, "SparXCC",
  SparXCC_above_threshold == F & Spearman_sig == T, "Spearman",
  SparXCC_above_threshold == F & Spearman_sig == F, "None"
)]

# Constructing figure
ggplot(data = long_CLR, aes(x = spear_cor, y = SparXCC, col = Correlated))+
  geom_point()+
  scale_color_manual(limits = c("Spearman", "SparXCC", "Both", "None"),
                     values = c("#6A3D9A", "#33A02C", "#1F78B4", "#A6CEE3"),
                     name = "Correlated")+
  xlim(-1,1)+
  ylim(-1,1)+
  labs(x = "TSS-Spearman", y = "SparXCC")+
  geom_abline(slope = 1, intercept = 0)+
  NULL -> g

ggsave("Results and Figures/SFig6.pdf", g, width = 200, height = 120, units = "mm")

# Maintaining only relevant information and saving
long_CLR <- long_CLR[Correlated != "None"]
long_CLR[,Correlated:=NULL]
colnames(long_CLR)[1:2] <- c("Bacterial OTU", "Fungal OTU")

fwrite(long_CLR, "Results and Figures/S4 Table.csv")

# Examining the correlations where Spearman and SparXCC disagree on the sign
w <- which(sign(CC_XCC) != sign(CC_spear$r) & abs(CC_XCC) > cor_thresh, arr.ind = TRUE)
data.table(Bac_OTU = rownames(CC_XCC)[w[,1]],
           Fung_OTU = colnames(CC_XCC)[w[,2]]) -> D
quantile(as.matrix(Fung[OTUid == "OTUF_11", ..gifu_smp]))
quantile(as.matrix(Fung[OTUid == "OTUF_5", ..gifu_smp]))

# Examining the correlations where Pearson and SparXCC disagree on the sign
w <- which(sign(CC_XCC) != sign(CC_pears$r) & abs(CC_XCC) > cor_thresh, arr.ind = TRUE)
data.table(Bac_OTU = rownames(CC_XCC)[w[,1]],
           Fung_OTU = colnames(CC_XCC)[w[,2]]) -> D2

# Examining pair whose sign differ for SparXCC vs Spearman but not SparXCC vs Pearson
apply(D, 1, paste, collapse = "_") %in% apply(D2, 1, paste, collapse = "_")
D[1]
CC_pears$r["OTU_10", "OTUF_11"]; CC_spear$r["OTU_10", "OTUF_11"]
CC_XCC["OTU_10", "OTUF_11"]

# Examining pair whose sign differ for SparXCC vs Spearman but not SparXCC vs Pearson
apply(D2, 1, paste, collapse = "_") %in% apply(D, 1, paste, collapse = "_")
D2[c(17, 22)]
CC_pears$r["OTU_533", "OTUF_5"]; CC_spear$r["OTU_533", "OTUF_5"]
CC_XCC["OTU_533", "OTUF_5"]
CC_pears$r["OTU_37", "OTUF_5"]; CC_spear$r["OTU_37", "OTUF_5"]
CC_XCC["OTU_37", "OTUF_5"]
