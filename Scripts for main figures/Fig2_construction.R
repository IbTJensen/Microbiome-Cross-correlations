# Code by Ib Thorsgaard Jensen
# This file should be run after the files Fig2_cluster.R and Fig2_loadings.R
pkg <- c("ggplot2", "ggpubr", "RColorBrewer")
for(i in pkg){
  library(i, character.only = T)
}

load("../Temp/g.txt")
load("../Temp/g2.txt")

gg <- ggarrange(g + rremove("legend"), g2 + rremove("ylab"), common.legend = T, legend = "bottom", ncol = 2, label.y = "Mean Absolute Error")
ggsave(filename = "../Main figures/Fig2_CaseB.pdf", gg, width = 190, height = 160, unit = "mm")
