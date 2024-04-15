# Code by Ib Thorsgaard Jensen

# Loading packages
pkg <- c("data.table", "ggplot2", "ggpubr", "RColorBrewer", "magrittr",
         "phyloseq", "entropy", "parallel", "CompoCor")
for(i in pkg){
  library(i, character.only = T)
}
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading data files
OTU_data <- fread("../Data/Byrd_OTU_table.csv")
taxonomy <- fread("../Data/Byrd_taxonomy.csv")
meta <- fread("../Data/byrd_metadata.txt")

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

RNGkind("L'Ecuyer-CMRG")
set.seed(1708331668)
S <- SparCEV(OTU.abn = as.data.frame(Taxon_data[,-1]),
             phenotype = meta_sub$`Objective SCORAD`,
             pseudo_count = 1,
             iter = 20,
             t = NULL,
             B_t = 100,
             t_quant = 0.8,
             Find_m = T,
             B_m = 100,
             cores = 10)
S_cor <- S$cor
cor_fam <- names(S$Correlated)[S$Correlated]

S_cor_base <- SparCEV_base(OTU.abn = as.data.frame(Taxon_data[,-1]),
                           phenotype = meta_sub$`Objective SCORAD`,
                           Find_m = F)

# Comparison between SparCEV and other methods
Taxon_data_TSS <- (Taxon_data[,-1]+1)/rowSums((Taxon_data[,-1]+1))
Taxon_data_clr <- t(apply(Taxon_data[,-1]+1, 1, compositions::clr))
Taxon_data_CEV <- t(apply(log(Taxon_data[,-1]+1), 1, function(x) x - mean(x[S$Final_idx]) ))

pclr <- psych::corr.test(Taxon_data_clr, meta_sub$`Objective SCORAD`, adjust = "fdr")
cor_fam_clr <- rownames(pclr$p.adj)[pclr$p.adj<0.05]
pTSS <- psych::corr.test(log(Taxon_data_TSS), meta_sub$`Objective SCORAD`, adjust = "fdr")
cor_fam_TSS <- rownames(pTSS$p.adj)[pTSS$p.adj<0.05]

setdiff(cor_fam, cor_fam_clr)
setdiff(cor_fam, cor_fam_TSS)

intersect(cor_fam, cor_fam_clr)
intersect(cor_fam, cor_fam_TSS)

set.seed(1708331668)
lapply(1:100,
       function(x){
         cor(Taxon_data_clr, meta_sub$`Objective SCORAD`[sample(1:27)])
       }) -> b_clr
clr_m <- b_clr %>% lapply(abs) %>% lapply(max) %>% unlist() %>% mean()
cor_fam_clr_perm <- rownames(pclr$r)[abs(pclr$r) > clr_m]

B <- data.table(Family = names(S_cor),
                cor_TSS = pTSS$r[,1],
                cor_clr = pclr$r[,1],
                cor_CEV = S_cor,
                cor_CEV_base = S_cor_base$cor,
                sig_TSS = names(S_cor) %in% cor_fam_TSS,
                sig_clr = names(S_cor) %in% cor_fam_clr,
                sig_CEV = S$Correlated)

ggplot(data = B, aes(x = cor_clr, cor_CEV))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)+
  labs(x = "CLR", y = "SparCEV") -> CEV_v_CLR

ggsave("../Supplementary figures/SFig3.pdf", CEV_v_CLR)
fwrite(B, "../Tables/STable1.csv")

A <- data.table(Family = colnames(Taxon_data)[-1],
                SparCEV = S_cor,
                SparCEV_base = S_cor_base$cor,
                TSS = cor(log(Taxon_data_TSS), meta_sub$`Objective SCORAD`)[,1],
                CLR = pclr$r[,1])

# Computes diversity for all samples
phs <- phyloseq(otu_table(Taxon_data[,-1], taxa_are_rows = F))
set.seed(1673967505)
RR <- rarefy_even_depth(phs)
M <- as.matrix(otu_table(RR))

e <- apply(M, 1, entropy::entropy)
C <- data.table(p_eff = exp(e),
                SCORAD = meta_sub$`Objective SCORAD`,
                Staphylococcaceae = Taxon_data_CEV[,"Staphylococcaceae"])

# Constructing right column of Fig 3
ggplot(data = C, aes(x = SCORAD, y = p_eff))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x)+
  labs(x = "Objective SCORAD", y = "Effective number of families")+
  theme(text = element_text(size = 10),
        axis.title = element_text(size = 9))+
  NULL -> p1

ggplot(data = A, aes(x = SparCEV, y = TSS))+
  geom_point()+
  xlim(-0.8, 0.8)+
  ylim(-0.8, 0.8)+
  geom_abline(slope = 1, intercept = 0)+
  theme(text = element_text(size = 10),
        axis.title = element_text(size = 9))+
  labs(x = "SparCEV iterative", y = "log-TSS")+
  NULL -> p3

ggplot(data = A, aes(x = SparCEV, y = SparCEV_base))+
  geom_point()+
  xlim(-0.8, 0.8)+
  ylim(-0.8, 0.8)+
  geom_abline(slope = 1, intercept = 0)+
  theme(text = element_text(size = 10),
        axis.title = element_text(size = 9))+
  labs(y = "SparCEV base", x = "SparCEV iterative")+
  NULL -> p4

gg <- ggarrange(p1, p3, p4, ncol = 1, labels = c("B", "C", "D"))

intercept_OLS <- mean(A$SparCEV_base - A$SparCEV)

# Finding Confidence intervals with bootstrap
B <- 1000
set.seed(1711015815)
bt_for_CI <- mclapply(1:B,
                      function(x){
                        idx <- sample(1:nrow(Taxon_data), replace = T)
                        SparCEV(OTU.abn = as.data.frame(Taxon_data[idx,-1]),
                                phenotype = meta_sub$`Objective SCORAD`[idx],
                                pseudo_count = 1,
                                iter = 20,
                                t = S$t,
                                Find_m = F)$cor
                      }, mc.cores = 10)
bt_cor_CI_mat <- matrix(unlist(bt_for_CI), ncol(Taxon_data)-1, B)

# S_cor_f <- Fisher_transformation(S_cor)
G_est <- rowMeans(bt_cor_CI_mat<S_cor)
z0 <- sapply(G_est, qnorm)
theta_i <- sapply(1:nrow(meta_sub), function(x) SparCEV(OTU.abn = as.data.frame(Taxon_data[-x,-1]),
                                                        phenotype = meta_sub$`Objective SCORAD`[-x],
                                                        pseudo_count = 1,
                                                        iter = 20,
                                                        t = S$t,
                                                        t_quant = 0.8,
                                                        Find_m = F)$cor)
# theta_i <- Fisher_transformation(theta_i)
U <- (nrow(meta_sub)-1)*(rowMeans(theta_i) - theta_i)
a <- 1/6*((rowSums(U^3))/(rowSums(U^2))^1.5)
z_norm_025 <- qnorm(0.025)
z_025 <- z0 + (z0+z_norm_025)/( 1 - a*(z0+z_norm_025) )
qnt_low <- pnorm(z_025)
z_norm_975 <- qnorm(0.975)
z_975 <- z0 + (z0+z_norm_975)/( 1 - a*(z0+z_norm_975) )
qnt_high <- pnorm(z_975)
CI_lower <- rep(NA, ncol(Taxon_data)-1)
CI_upper <- rep(NA, ncol(Taxon_data)-1)
for(i in 1:length(CI_lower)){
  CI_lower[i] <- quantile(bt_cor_CI_mat[i,], qnt_low[i])
  CI_upper[i] <- quantile(bt_cor_CI_mat[i,], qnt_high[i])
}

CI_lower[198]; CI_upper[198]
CI_lower[301]; CI_upper[301]

# Collecting results in data table
D_temp <- data.table(Taxon = colnames(Taxon_data)[-1],
                     Correlation = S_cor,
                     Significance = S$Correlated,
                     Lower = CI_lower,
                     Upper = CI_upper)
D <- D_temp[order(Correlation, decreasing = T)]

# Keeping only OTUs with correlation above the permutation threshold
D <- D[Significance == T]
D <- D[sign(Lower) == sign(Upper)]
D[,Sign:=ifelse(sign(Correlation) == 1, "Positive", "Negative")]
D <- D[order(abs(Correlation))]
D[,":="(Correlation = abs(Correlation), Lower = abs(Lower), Upper = abs(Upper))]
D[,Taxon:=factor(Taxon, levels = Taxon)]

# Constructing figure
ggplot(data = D, aes(x = Correlation, y = Taxon, fill = Sign))+
  geom_bar(stat="identity")+
  theme(text = element_text(size = 10),
        axis.title = element_text(size = 9))+
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),
                     minor_breaks = seq(0,0.9,0.1))+
  scale_fill_brewer(palette = "Dark2")+
  labs(x = "Absolute Correlation", y = "Family")+
  geom_errorbar(aes(xmin=Lower, xmax=Upper), width=.2,
                position=position_dodge(.9))+
  expand_limits(x = 1)+
  NULL -> g

ggg <- ggarrange(g, gg, ncol = 2, widths = c(1,0.5), legend = "bottom", common.legend = T, labels = c("A", NULL))

ggsave("../Main figures/Fig3.pdf", ggg, width = 190, height = 200, units = "mm")

