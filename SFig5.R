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

# Loading metadata and taxonomy
meta_data <- fread("Data/meta_data.csv")
taxonomy_Bac <- fread("Data/taxonomy_bacteria.txt", fill = T)
taxonomy_Fung <- fread("Data/taxonomy_fungi.txt", fill = T)
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

# Constructing subtable with only rhizosphere samples from Gifu
Bac_gifu <- Bac_full[Compartment == "Rhizosphere" & Genotype == "Gifu"]
Fung_gifu <- Fung_full[Compartment == "Rhizosphere" & Genotype == "Gifu"]

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
  col_keep <- c("OTUid", meta_data[Compartment == "Rhizosphere" & Genotype == gt[i], SampleID])
  meta_sub[[i]] <- meta_data[Compartment == "Rhizosphere" & Genotype == gt[i]]
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
g <- list()
fun_taxa <- list()
bac_taxa <- list()
cors <- list()
set.seed(1671522881)
for(i in 1:5){
  # CLR transformation
  Bac_clr <- t(apply(Bac_subset[[i]][,-1]+1, 2, clr))
  Fung_clr <- t(apply(Fung_subset[[i]][,-1]+1, 2, clr))
  colnames(Bac_clr) <- Bac_subset[[i]]$OTUid
  colnames(Fung_clr) <- Fung_subset[[i]]$OTUid
  
  # Removing batch-effect
  Bac_clr_adj <- t(removeBatchEffect(t(Bac_clr), meta_sub[[i]]$Experiment))
  Fung_clr_adj <- t(removeBatchEffect(t(Fung_clr), meta_sub[[i]]$Experiment))
  
  # Estimating correlations
  C <- cor(Bac_clr_adj, Fung_clr_adj)
  # Permutation thresholding
  cor_thresh <- sapply(1:50, function(x) cor(Bac_clr_adj[sample(1:N[i]),], Bac_clr_adj[sample(1:N[i]),]) %>% abs() %>% quantile(1 - 1/(p[i]*q[i]))) %>% mean()
  sum(abs(C) > cor_thresh)
  
  # Constructing matrix with 1 when correlation is above threshold and 0 otherwise
  # Rows and columns with all zeros are excluded 
  C2 <- ifelse(abs(C) > cor_thresh & C != 1, 1, 0)
  C3 <- C2[apply(C2, 1, function(x) any(x!=0)), apply(C2, 2, function(x) any(x!=0)), drop = F]
  
  ## saving Results
  A <- which(C2 == 1, arr.ind = TRUE)
  Res <- data.table(BOTU = rownames(C)[A[,1]],
                    FOTU = colnames(C)[A[,2]],
                    Correlation = as.numeric(C)[C2 == 1],
                    Genotype = gt[i])
  
  cors[[i]] <- Res
  
  # Taxonomy for OTUs in Res
  bac_taxa <- taxonomy_Bac[V1 %in% Res$BOTU]
  fun_taxa <- taxonomy_Fung[V1 %in% Res$FOTU]
  
  ### Phylum-level nodes
  # Taxonomy for the OTUs in the dataset after filtering
  tax_bac <- taxonomy_Bac[V1 %in% colnames(Bac_clr)]
  tax_fun <- taxonomy_Fung[V1 %in% colnames(Fung_clr)]
  
  # Making nicer names for the phyla
  tax_bac[,V3:=substr(V3, 4, 1000)]
  tax_bac[!(V3 %in% c("Bacteroidetes", "Actinobacteria", "Proteobacteria")),
          V3 := "Other bacteria"]
  tax_fun[,V3:=substr(V3, 4, 1000)]
  tax_fun[!(V3 %in% c("Ascomycota", "Glomeromycota")),
          V3 := "Other fungi"]
  tax <- rbind(tax_bac, tax_fun, fill = T)
  tax[,V2:=substr(V2, 4,1000)]
  tax[!grepl("F", V1), V2:="Bacteria"]
  
  # Saving only Kingdom and Phylum level taxonomy
  df <- data.frame(tax[,c("V2","V3")], row.names = tax$V1)
  
  # Construct fill adjacency matrix
  adj_mat <- matrix(0, sum(dim(C3)), sum(dim(C3)))
  p2 <- dim(C3)[1]; q2 <- dim(C3)[2]
  adj_mat[1:p2, p2+1:q2] <- C3
  adj_mat <- adj_mat + t(adj_mat)
  colnames(adj_mat) <- c( rownames(C3), colnames(C3) )
  rownames(adj_mat) <- c( rownames(C3), colnames(C3) )
  
  # Constructing network and adding Phylum and Kingdom as node-information
  nn <- network(adj_mat)
  nn %v% "Phylum" = df[rownames(adj_mat),2]
  nn %v% "Kingdom" = df[rownames(adj_mat),1]
  
  # Plotting network
  ggnet2(nn,
         color = "Phylum",
         palette = colorBlindBlack8,
         shape = "Kingdom",
         shape.palette = c("Bacteria" = 19, "Fungi" = 15),
         node.size = 3,
         legend.position = "bottom",
         mode = "fruchtermanreingold") + 
    scale_color_manual(name = "Phylum", limits = names(colorBlindBlack8), values = colorBlindBlack8)+
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
    guides(color = guide_legend(override.aes = list(size=3)))+
    NULL -> g[[i]]
}

# Constructing plot
g1 <- ggarrange(g[[2]] + theme(legend.position = "none"), g[[3]] + theme(legend.position = "none"), nrow = 1)
gg <- ggarrange(g[[1]], g[[4]], g[[5]], g1, common.legend = T, legend = "bottom")

ggsave("Results and Figures/Cor_networks_permute_rhizosphere.pdf", gg, width = 200, height = 160, units = "mm")

# Saving correlation results as table
Cors_filt <- rbindlist(cors)
fwrite(Cors_filt, "Results and Figures/S3 Table.csv")
