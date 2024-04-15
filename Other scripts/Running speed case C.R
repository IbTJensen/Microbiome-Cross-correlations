pkg <- c("data.table", "Matrix", "hydroGOF", "ggplot2", "ggpubr", "RColorBrewer",
         "compositions", "DESeq2", "microbenchmark")
for(i in pkg){
  library(i, character.only = T)
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source(file = "../Functions/Data_generation.R")
source(file = "../Functions/SparXCC_function.R")
source(file = "../Functions/cor_mat_cluster.R")

# Loading template
template_OTU <- fread("../Templates/OTU_template.csv")
template_gene <- fread("../Templates/Gene_template.csv")

p <- 1000
q <- c(100, 1000, 10000)

Opt <- expand.grid(p = p, q = q, dens = 0.1, bio_zero = "Yes", cor = 0.75)

Res <- data.table(Time = NA, Pairs = NA, Method = NA)[-1]

set.seed(1711034849)
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
  n.sim <- 50
  n.seq.abn <- p*40
  n.seq.expr <- q*200
  
  # Constructing correlation matrix
  cor_mat <- cor_mat_sim_cluster(p,
                                 q,
                                 cor.prob = 0.1,
                                 cor = 0.75,
                                 case = "C")
  true.cor.mat <- cor_mat[1:p,p+1:q]
  true.cor <- as.numeric(true.cor.mat)
  
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
  
  # CLR
  microbenchmark( cor(t(apply(A$OTU_reads+1, 1, clr)), 
                      t(apply(A$Gene_reads+1, 1, clr)))
                  ) -> m_CLR
  
  # CLR+VST
  CLR_VST <- function(x){
    Gene_VST <- tryCatch(t(varianceStabilizingTransformation(t(A$Gene_reads), fitType = "local")), error = function(e) "Error")
    if(length(class(Gene_VST)) > 1){
      OTU_CLR <- t(apply(A$OTU_reads+1, 1, clr))
      Mix <- cor(OTU_CLR, Gene_VST)
    } else{
      Mix <- NA
    }
    return(Mix)
  }
  microbenchmark(CLR_VST(0)) -> m_VST
  mean(m_VST$time/(1e+9))
  
  # SparXCC base
  microbenchmark(SparXCC_base(A$OTU_reads, A$Gene_reads)) -> m_XCC
  mean(m_XCC$time/(1e+9))
  
  # # SparXCC iterative
  # microbenchmark(XCC2 <- SparXCC(OTU.abn = A$OTU_reads,
  #                                gene.expr = A$Gene_reads,
  #                                pseudo_count = 1,
  #                                var_min = 1e-5,
  #                                iter = 10,
  #                                cutoff = 0.12,
  #                                cutoff2 = 0.12)) -> m_XCC2
  # mean(m_XCC2$time/(1e+9))
  
  data.table(Time = c(mean(m_CLR$time/(1e+9)),
                      mean(m_VST$time/(1e+9)),
                      mean(m_XCC$time/(1e+9))),
             Pairs = p*q,
             Method = c("CLR", "CLR+VST", "SparXCC")
             ) -> B
  
  Res <- rbind(Res, B)
}

Res_table <- data.frame(p = 1000, q = Opt[,"q"], matrix(Res$Time, 3, 3, byrow = T))
colnames(Res_table)[3:5] <- c("CLR", "CLR+VST", "SparXCC (base)")

xtable(x = Res_table, align = "cccccc", digits = 4,
       caption = "Average running time for cross-correlation estimation methods in seconds.
       Time benchmark was carried out with the R package microbenchmark [cite] on a 
       Lenovo X1 Carbon labptop equipped with a 13th Gen Intel(R) Core(TM) i7-1365U processor.")
# 
# ggplot(data = Res, aes(x = Pairs, y = Time, col = Method, fill = Method, group = Method))+
#   geom_line()+
#   labs(y = "Time [seconds]")+
#   scale_x_continuous(trans='log10')+
#   scale_y_continuous(trans='log10')+
#   # geom_ribbon(alpha=0.35, aes(col = NULL))+
#   scale_colour_manual(values = c("#1B9E77", "#E6AB02", "#E7298A"),
#                       breaks = c("CLR", "CLR+VST", "SparXCC"))+
#   scale_fill_manual(values = c("#1B9E77", "#E6AB02", "#E7298A"),
#                     breaks = c("CLR", "CLR+VST", "SparXCC"))+
#   # expand_limits(y=10)+
#   theme(legend.position = "bottom")+
#   NULL -> g
# 
# ggsave("Case_C_Running_time.pdf", g, width = 190, height = 160, units = "mm")
