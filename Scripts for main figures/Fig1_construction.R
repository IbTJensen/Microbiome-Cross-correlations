# Code by Ib Thorsgaard Jensen
# This file should be run after the files Fig1_cluster.R and Fig1_loadings.R
pkg <- c("ggplot2", "ggpubr", "RColorBrewer")
for(i in pkg){
  library(i, character.only = T)
}

load("../Temp/g_yes.txt")
load("../Temp/g_no.txt")
load("../Temp/g2_yes.txt")
load("../Temp/g2_no.txt")

gg_no <- ggarrange(g_no, g2_no + rremove("ylab") + rremove("legend"), common.legend = T, legend = "bottom", ncol = 2, label.y = "Mean Absolute Error")
ggsave(filename = "../Supplementary figures/SFig1_CaseB.pdf", gg_no, width = 190, height = 160, unit = "mm")

gg_yes <- ggarrange(g_yes, g2_yes + rremove("ylab") + rremove("legend"), common.legend = T, legend = "bottom", ncol = 2, label.y = "Mean Absolute Error")
ggsave(filename = "../Main figures/Fig1_CaseB.pdf", gg_yes, width = 190, height = 160, unit = "mm")

