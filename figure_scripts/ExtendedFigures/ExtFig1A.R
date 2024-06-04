# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## Extended Fig1A ##
library(ggplot2)
library(dplyr)
library(tidyverse)
library(factoextra)
library(cowplot)


# meta
setwd('./data/ExtFig1/')
meta <- read.csv("Fig1_metadata2024.csv", row.names = 1)

# KW pvalues
pval_df <- sapply(meta[,c(2:9,11)], function(i) kruskal.test(i ~ meta$VisitMonths))[3,] %>% data.frame() %>% pivot_longer(1:9)
colnames(pval_df) <- c("name", "KWpval")
pval_df$KWpval <- round(pval_df$KWpval, 9)
pval_df$KWpval <- paste("KWpval=", pval_df$KWpval, sep="")
pval_df$name <- gsub("\\.", " ", pval_df$name)
pval_df$name <-  gsub("mol L", "mol/L", pval_df$name)
pval_df$name <- str_remove(pval_df$name, "X")
pval_df$name <-  gsub("alfa Ohp", "alfa-Ohp", pval_df$name)

# plot these with ranges for each hormone
# Reference intervals of 20-40 year age group  from Bae et al 2019: https://doi.org/10.1016/j.jsbmb.2019.105409
# estradiol <1109 pmol/L for women and <158 pmol/L for men; 
# testosterone 0.3-2 nmol/L for women and 8-29 nmol/L for men. 
# NOTE: DHEAS in the publication is conjugated DHEA, not free as we measured.

# 17A OHP
meta$X17alfa.Ohp.nmol.L
A <- meta %>% 
  ggplot(aes(x=VisitMonths, y=X17alfa.Ohp.nmol.L)) + 
  #male range
  geom_rect(aes(ymin = 0.85, ymax = 6.28, xmin = -Inf, xmax = Inf), fill = "midnightblue", alpha = 0.006) +
  #female range
  geom_rect(aes(ymin = 0, ymax = 8.02, xmin = -Inf, xmax = Inf), fill = "deeppink4", alpha = 0.006) +
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)") + 
  geom_text(aes(x = -Inf, y = -Inf, label=KWpval), data=pval_df[pval_df$name %in% "17alfa-Ohp nmol/L",], hjust   = -0.1, vjust   = -0.1, size=3)

#Androstendion
meta$Androstendion.nmol.L
B <- meta %>% 
  ggplot(aes(x=VisitMonths, y=Androstendion.nmol.L)) + 
  #male range
  geom_rect(aes(ymin = 1.41, ymax = 6.20, xmin = -Inf, xmax = Inf), fill = "midnightblue", alpha = 0.006) +
  #female range
  geom_rect(aes(ymin = 1.32, ymax = 7.37, xmin = -Inf, xmax = Inf), fill = "deeppink4", alpha = 0.006) +
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)") + 
  geom_text(aes(x = -Inf, y = -Inf, label=KWpval), data=pval_df[pval_df$name %in% "Androstendion nmol/L",], hjust   = -0.1, vjust   = -0.1, size=3)

#DHEA
meta$DHEA.nmol.L
C <- meta %>% 
  ggplot(aes(x=VisitMonths, y=DHEA.nmol.L)) + 
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)") + 
  geom_text(aes(x = -Inf, y = -Inf, label=KWpval), data=pval_df[pval_df$name %in% "DHEA nmol/L",], hjust   = -0.1, vjust   = -0.1, size=3)

#DHT
meta$DHT.nmol.L
D <- meta %>% 
  ggplot(aes(x=VisitMonths, y=DHT.nmol.L)) + 
  #male range
  #female range
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)") + 
  geom_text(aes(x = -Inf, y = -Inf, label=KWpval), data=pval_df[pval_df$name %in% "DHT nmol/L",], hjust   = -0.1, vjust   = -0.1, size=3)

#Testosterone
meta$Testosteron.nmol.L
E <- meta %>% 
  ggplot(aes(x=VisitMonths, y=Testosteron.nmol.L)) + 
  #male range
  geom_rect(aes(ymin = 7.98, ymax = 29.14, xmin = -Inf, xmax = Inf), fill = "midnightblue", alpha = 0.006) +
  #female range
  geom_rect(aes(ymin = 0.33, ymax = 2.05, xmin = -Inf, xmax = Inf), fill = "deeppink4", alpha = 0.006) +
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)") + 
  geom_text(aes(x = -Inf, y = -Inf, label=KWpval), data=pval_df[pval_df$name %in% "Testosteron nmol/L",], hjust   = -0.1, vjust   = -0.1, size=3)

#Estrone
meta$Estron.pmol.L
F1 <- meta %>% 
  ggplot(aes(x=VisitMonths, y=Estron.pmol.L)) + 
  #male range
  #female range
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)") + 
  geom_text(aes(x = -Inf, y = -Inf, label=KWpval), data=pval_df[pval_df$name %in% "Estron pmol/L",], hjust   = -0.1, vjust   = -0.1, size=3)

# ExtFig 1A
plot_grid(A,B,C,D,E,F1,ncol=3)
