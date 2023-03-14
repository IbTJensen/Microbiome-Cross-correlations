# Code by Ib Thorsgaard Jensen
# Auxilliary file called by the files Fig1_cluster.R and Fig1_loadings.R
pkg <- c("data.table", "Matrix", "hydroGOF", "ggplot2", "ggpubr", "RColorBrewer",
         "compositions", "psych", "foreach", "parallel", "doParallel", "doRNG")
for(i in pkg){
  library(i, character.only = T)
}
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source(file = "../Data_generation.R")
source(file = "../SparCEV_function.R")

# Loading template
template <- fread("Templates/OTU_template.csv")
reps <- 1000
bio_zero <- c("Yes", "No")
cores <- 20
