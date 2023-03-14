# Code by Ib Thorsgaard Jensen

# Loading packages
pkg <- c("data.table", "ggplot2", "ggpubr", "RColorBrewer", "magrittr", "phyloseq")
for(i in pkg){
  library(i, character.only = T)
}
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("SparCEV_function.R")

# Loading data files
OTU_data <- fread("Data/Byrd_OTU_table.csv")
taxonomy <- fread("Data/Byrd_taxonomy.csv")
meta <- fread("Data/byrd_metadata.txt")

taxonomy <- taxonomy[family != ""]
keep_OTU <- c("SampleID",taxonomy$OTUid)
OTU_data <- OTU_data[,..keep_OTU]

# Keeping data from only one body-site
keep_sample <- meta[Site == "Pc-R", Sample_ID]
OTU_data <- OTU_data[SampleID %in% keep_sample]
meta_sub <- meta[Sample_ID %in% keep_sample]

# Collapsing the OTU-table into the family level
OTU_table <- transpose(OTU_data, keep.names = "OTUid", make.names = "SampleID")
Taxon_list <- rowsum(OTU_table[,-1], taxonomy$family)

# Omitting families present in less than 10% of the samples
Taxon_list <- Taxon_list[apply(Taxon_list, 1, function(x) mean(x!=0)>0.1 & mean(x) > 10),]
Taxon_data <- data.table(SampleID = colnames(Taxon_list), t(Taxon_list))

# Running SparCEV
S <- SparCEV(Taxon_data[,-1], meta_sub$`Objective SCORAD`)

# Comparison between SparCEV and other methods
Taxon_data_TSS <- (Taxon_data[,-1]+1)/rowSums((Taxon_data[,-1]+1))
Taxon_data_clr <- t(apply(Taxon_data[,-1]+1, 1, compositions::clr))
A <- data.table(Family = colnames(Taxon_data)[-1],
                SparCEV = S,
                TSS = cor(log(Taxon_data_TSS), meta_sub$`Objective SCORAD`)[,1],
                CLR = cor(Taxon_data_clr, meta_sub$`Objective SCORAD`)[,1])

# Computes diversity for all samples
phs <- phyloseq(otu_table(Taxon_data[,-1], taxa_are_rows = F))
set.seed(1673967505)
RR <- rarefy_even_depth(phs)
M <- as.matrix(otu_table(RR))

e <- apply(M, 1, entropy::entropy)
C <- data.table(p_eff = exp(e),
                SCORAD = meta_sub$`Objective SCORAD`,
                Staphylococcaceae = Taxon_data_clr[,354])

# Constructing right column of Fig 3
ggplot(data = C, aes(x = SCORAD, y = p_eff))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x)+
  labs(x = "Objective SCORAD", y = "Effective number of families")+
  theme(text = element_text(size = 10))+
  NULL -> p1

ggplot(data = C, aes(x = SCORAD, y = Staphylococcaceae))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x)+
  labs(x = "Objective SCORAD", y = "CLR(Staphylococcaceae)")+
  theme(text = element_text(size = 10))+
  NULL -> p2

ggplot(data = A, aes(x = SparCEV, y = TSS))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)+
  theme(text = element_text(size = 10))+
  NULL -> p3

gg <- ggarrange(p2, p1, p3, ncol = 1)

# Performing permutation thresholding
set.seed(1670509317)
cor_thresh <- mean(sapply(1:16, function(x) S_perm <- SparCEV(Taxon_data[,-1], meta_sub$`Objective SCORAD`[sample(1:nrow(meta_sub))]) %>% abs() %>% quantile(1-1/407)))

# Collecting results in data table
D_temp <- data.table(Taxon = colnames(Taxon_data)[-1],
                     Correlation = S,
                     above_thresh = ifelse(abs(S)>cor_thresh, "Yes", "No"))
D <- D_temp[order(Correlation, decreasing = T)]

# Keeping only OTUs with correlation above the permutation threshold
D <- D[above_thresh == "Yes"]
D[,Sign:=ifelse(sign(Correlation) == 1, "Positive", "Negative")]
D <- D[order(abs(Correlation))]
D[,Correlation:=abs(Correlation)]
D[,Taxon:=factor(Taxon, levels = Taxon)]

# Confidence intervals
Conf_int <- matrix(NA, length(S), 3); colnames(Conf_int) <- c("Estimate", "Lower 2.5%", "Upper 2.5%")
Conf_int[,"Estimate"] <- S
set.seed(1670244819)
B <- 1000
S_est <- matrix(NA, length(S), B)
# Bootstrapping the correlation estimates
for(i in 1:B){
  cat(i, "\r")
  resample <- sample(1:nrow(meta_sub), replace = T)
  S_perm <- SparCEV(Taxon_data[resample,-1], meta_sub$`Objective SCORAD`[resample])
  S_est[,i] <- S_perm
}
# Calculating empirical confidence intervals based on bootstraps
z0 <- sapply(rowSums(S_est<S)/B, qnorm)
theta_i <- sapply(1:nrow(meta_sub), function(x) SparCEV(Taxon_data[-x,-1], meta_sub$`Objective SCORAD`[-x]))
U <- (nrow(meta_sub)-1)*(rowMeans(theta_i) - theta_i)
a <- 1/6*((rowSums(U^3))/(rowSums(U^2))^1.5)
qnt_low <- pnorm(z0 + (z0 + qnorm(0.025))/(1 + a*(z0 + qnorm(0.025))))
qnt_high <- pnorm(z0 + (z0 + qnorm(0.975))/(1 + a*(z0 + qnorm(0.975))))
for(i in 1:nrow(Conf_int)){
  Conf_int[i,2] <- quantile(S_est[i,], qnt_low[i])
  Conf_int[i,3] <- quantile(S_est[i,], qnt_high[i])
}
rownames(Conf_int) <- colnames(Taxon_data)[-1]
fwrite(file = "Conf_int.csv", data.table(Taxon = rownames(Conf_int), Conf_int))

# Saving Results for all families
All_res <- cbind(D_temp, Conf_int)
All_res[,Estimate:=NULL]
colnames(All_res) <- c("Family", "Correlation", "Above threshold", "Lower 2.5% CI", "Upper  2.5% CI")
setcolorder(All_res, c("Family", "Correlation",  "Lower 2.5% CI", "Upper  2.5% CI", "Above threshold"))
fwrite(All_res, "S1 Table.csv")

# Keeping only families with confidence-intervals that do not cross 0
CC <- Conf_int[as.character(D$Taxon),]
D <- D[which(sign(CC[,2]) == sign(CC[,3]))]
CC <- CC[sign(CC[,2]) == sign(CC[,3]),]
D[,":="(Lower = abs(CC)[,"Lower 2.5%"], Upper = abs(CC)[,"Upper 2.5%"])]

# Constructing figure
ggplot(data = D, aes(x = Correlation, y = Taxon, fill = Sign))+
  geom_bar(stat="identity")+
  theme(text = element_text(size = 10))+
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8),
                     minor_breaks = seq(0,0.9,0.1))+
  scale_fill_brewer(palette = "Dark2")+
  labs(x = "Absolute Correlation", y = "Family")+
  geom_errorbar(aes(xmin=Lower, xmax=Upper), width=.2,
                position=position_dodge(.9))+
  NULL -> g

ggg <- ggarrange(g, gg, ncol = 2, widths = c(1,0.5), legend = "right")

ggsave("Simulations/Figures/Fig3.pdf", ggg, width = 200, height = 190, units = "mm")
