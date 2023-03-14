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

# Function that returns the right hand-side of (5) in the paper as function of nu_1
fun_gen <- function(p, sp_eff){
  f <- function(x){
    return(-x*log(x)-(1-x)*log((1-x)/(p-1)) - log(sp_eff))
  }
  return(f)
}

# A function that constructs an abundance-template for a given effective number of speicies 
construct_template <- function(sp_eff, p, bac.load){
  if(sp_eff > p){
    stop("sp_eff <= p is required")
  }
  if(sp_eff == p){
    return(bac.load*rep(1/p, p))
  }
  fun <- fun_gen(p, sp_eff)
  r1 <- uniroot(fun, c(1e-15, 1-1e-15))
  rel.abn <- rep(NA, p)
  rel.abn[1] <- r1$root
  rel.abn[-1] <- (1-r1$root)/(p-1)
  return(rel.abn*bac.load)
}

p <- 100
reps <- 1000
cores <- 20
