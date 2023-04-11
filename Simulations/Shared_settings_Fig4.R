# Code by Ib Thorsgaard Jensen
# Auxilliary file called by the files Fig1_cluster.R and Fig1_loadings.R
pkg <- c("data.table", "Matrix", "hydroGOF", "ggplot2", "ggpubr", "RColorBrewer",
         "compositions", "DESeq2", "foreach", "parallel", "doParallel", "doRNG")
for(i in pkg){
  library(i, character.only = T)
}
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source(file = "../Data_generation.R")
source(file = "../SparXCC_function.R")

# Loading template
template_OTU <- fread("Templates/OTU_template.csv")
template_gene <- fread("Templates/Gene_template.csv")
reps <- 200
bio_zero <- c("Yes", "No")
cores <- 20
cols <- c("#1B9E77", "#E6AB02", "#D95F02", "#7570B3", "#E7298A")
