# Code by Ib Thorsgaard Jensen
# This file should be run after the files Fig4_cluster.R and Fig4_loadings.R
pkg <- c("ggplot2", "ggpubr", "RColorBrewer")
for(i in pkg){
  library(i, character.only = T)
}

load("g_yes.txt")
load("g_no.txt")
load("g2_yes.txt")
load("g2_no.txt")

gg_no <- ggarrange(g_no_main, g2_no_main + rremove("legend") + rremove("ylab"), common.legend = T, legend = "bottom", ncol = 2, label.y = "Mean Absolute Error")
ggsave(filename = "Figures/SFig3_CaseC.pdf", gg_no, width = 190, height = 160, unit = "mm")

gg_yes <- ggarrange(g_yes_main, g2_yes_main + rremove("legend") + rremove("ylab"), common.legend = T, legend = "bottom", ncol = 2, label.y = "Mean Absolute Error")
ggsave(filename = "Figures/Fig4_CaseC.pdf", gg_yes, width = 190, height = 160, unit = "mm")
ggsave(filename = "Figures/Fig4.eps", gg_yes, width = 190, height = 160, dpi = 300, unit = "mm")
