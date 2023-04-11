#script to plot volin plot for Fig 2h - pDC R848 stim
# By Petter Brodin (petter.brodin@ki.se)

#script to plot volin plot for Fig 2h - pDC R848 stim
# By Petter Brodin (petter.brodin@ki.se)

library(tidyverse)
library(ggthemes)

#setwd("enter source data folder")

#Filter source data on key Hallmark IFNa genes
ifna <- read.delim("violin_df_pdcR848HallmarkIFNa.csv", header=T, sep=",")
ifna <- ifna[ifna$gene %in% c("EIF2AK2", "ISG20", "MX1", "PARP14", "SP110"), ]
infa.2 <- spread(ifna, V, value)

#plot each BTM genes
data <- infa.2 # change according to btm

#code for plotting:
p <- ggplot(data) +
  # Top
  geom_density(aes(x = V1, y = ..density..), color = NA, fill="#999999" ) +
  # Bottom
  geom_density(aes(x = V2, y = -..density..), color = NA, fill= "#E2703C") +
  theme_bw() + facet_wrap(~gene, ncol=1, scales="free_x") +
  xlim(-10, 80) +
  theme(axis.line = element_blank(),
        aspect.ratio = .5, 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
p
#ggsave()