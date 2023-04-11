#script to plot volin plots for Fig 3b (Lakshmikanth et al, BioRxiv, 2023)
# Author: Petter Brodin (petter.brodin@ki.se)

library(tidyverse)
library(ggthemes)
#setwd("enter source data folder")

#Focus on Hallmark TNFa genes in monocytes
tnfa <- read.delim("violin_df_MonosLPS_Hallmark TNF.csv", header=T, sep=",")
tnfa <- tnfa[tnfa$gene %in% c("TNF", "IL1B", "IL1A", "IL6","SOCS3", "GOS2","F3", "CCL20", "CD44", "CXCL2", "SERPINB2"), ]
tnfa.2 <- spread(tnfa, V, value)

# create plot with ggplot2
p <- ggplot(tnfa.2) +
  # Top violin
  geom_density(aes(x = log(V1), y = ..density..), color = NA, fill="#999999" ) +
  # Bottom violin
  geom_density(aes(x = log(V2), y = -..density..), color = NA, fill= "#E2703C") +
  theme_bw() + facet_wrap(~gene, scales="free_x") +
  ylim(-1, 1) +
  theme(axis.line = element_blank(),
        aspect.ratio = 1, 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
p


