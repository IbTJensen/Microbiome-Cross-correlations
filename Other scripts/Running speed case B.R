pkg <- c("data.table", "Matrix", "hydroGOF", "ggplot2", "ggpubr", "RColorBrewer",
         "compositions", "DESeq2", "microbenchmark")
for(i in pkg){
  library(i, character.only = T)
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source(file = "../Functions/Data_generation.R")
source(file = "../Functions/SparCEV_function.R")
source(file = "../Functions/cor_mat_cluster.R")

# Loading template
template <- fread("../Templates/Gene_template.csv")

p <- c(100, 1000, 10000)

Opt <- expand.grid(p = p, dens = 0.1, bio_zero = "Yes", cor = 0.75)

Res <- data.table(Time = NA, p = NA, Method = NA)[-1]

set.seed(1711034849)
for(i in 1:nrow(Opt)){
  cat(i, "\r")
  
  p <- Opt[i, "p"]
  
  # Randomly selecting genes from the template
  idx <- sample(1:nrow(template), p)
  template_sim <- template[idx]
  if(Opt[i, "bio_zero"] == "No"){
    template_sim[,pi:=0]
  }
  
  # Defining other simulation settings
  n.sim <- 50
  n.seq.abn <- p*40
  
  # Constructing correlation matrix
  cor_mat <- cor_mat_sim_cluster(p,
                                 q = 1,
                                 cor.prob = 0.1,
                                 cor = 0.75,
                                 case = "B")
  true.cor.mat <- cor_mat[1:p,p+1]
  
  # Simulating sequencing data 
  Abs_abundance_sim_B(N = n.sim, 
                      OTU_template = template_sim$mu,
                      OTU_vars = template_sim$sigma2,
                      phen_means = log(30) - log(1+1/30^2)/2,
                      phen_var = log(1+1/30^2),
                      cor_mat = cor_mat,
                      zero_props = template_sim$pi,
                      n.seq = n.seq.abn) -> A
  
  # CLR
  microbenchmark( cor(t(apply(A$OTU_reads+1, 1, clr)), 
                      A$Phenotype) ) -> m_CLR
  
  # SparXCC base
  microbenchmark(SparCEV_base(A$OTU_reads, A$Phenotype)) -> m_CEV
  mean(m_CEV$time/(1e+9))
  
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
                      mean(m_CEV$time/(1e+9))),
             p = p,
             Method = c("CLR", "SparCEV")
  ) -> B
  
  Res <- rbind(Res, B)
}

Res_table <- data.frame(p = Opt[,"p"], matrix(Res$Time, 3, 2, byrow = T))
colnames(Res_table)[2:3] <- c("CLR", "SparXCC (base)")

xtable(x = Res_table, align = "cccc", digits = 4,
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
