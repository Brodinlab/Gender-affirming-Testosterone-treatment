# plot of pre-DC frequency during GAHT

library(tidyverse)
library(ggpubr)
library(Hmisc)

setwd("./data/Figure2/")
d <- read.csv("GAHT_Samples_CyTOF_DC subsets_v3.csv", header=T)

#plot %preDC by visit grouped by subject from manually gated CyTOF data
plot.preDC <- ggplot(d, aes(x=(Month), y=(Pre_DC), group=Subject_ID)) + 
  geom_hline(yintercept = 0, col = "gray60") + 
  geom_line() + geom_point(size=5) +
  theme_bw()  + xlim(0,12.5) + ggtitle("pre-DC%") 
plot.preDC

aov.preDC <- aov(formula = Pre_DC ~ Month , data = d)
summary(aov2.preDC)

