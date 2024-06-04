# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## FIGURES 1B, 1C, 1D, 1E ##
library(ggplot2)
library(dplyr)
library(tidyverse)
library(factoextra)


# meta
setwd("path_to_data/Figure1")
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

# Fig 1B: Bioavailable Testosterone
meta$NewBioavailableTestosterone
meta %>% 
  ggplot(aes(x=VisitMonths, y=NewBioavailableTestosterone)) + 
  #male range
  geom_rect(aes(ymin = 4.54, ymax = 23.64, xmin = -Inf, xmax = Inf), fill = "midnightblue", alpha = 0.006) +
  #female range
  geom_rect(aes(ymin = 0, ymax = 1.8, xmin = -Inf, xmax = Inf), fill = "deeppink4", alpha = 0.006) + 
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)") + 
  geom_text(aes(x = -Inf, y = -Inf, label=KWpval), data=pval_df[pval_df$name %in% "NewBioavailableTestosterone",], hjust   = -0.1, vjust   = -0.1, size=3)

# Fig 1C: Estradiol
meta$`Estradiol pmol/L`
meta %>% 
  ggplot(aes(x=VisitMonths, y=`Estradiol pmol/L`)) + 
  #male range
  geom_rect(aes(ymin = 0, ymax = 158.4, xmin = -Inf, xmax = Inf), fill = "midnightblue", alpha = 0.006) +
  #female range
  geom_rect(aes(ymin = 236, ymax = 1109, xmin = -Inf, xmax = Inf), fill = "deeppink4", alpha = 0.006) +
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)") + 
  geom_text(aes(x = -Inf, y = -Inf, label=KWpval), data=pval_df[pval_df$name %in% "Estradiol pmol/L",], hjust   = -0.1, vjust   = -0.1, size=3)

# Fig 1D: Progesterone
meta$`Progesteron nmol/L`
meta %>% 
  ggplot(aes(x=VisitMonths, y=`Progesteron nmol/L`)) + 
  #male range
  geom_rect(aes(ymin = 0, ymax = 0.65, xmin = -Inf, xmax = Inf), fill = "midnightblue", alpha = 0.006) +
  #female range
  geom_rect(aes(ymin = 0, ymax = 61.61, xmin = -Inf, xmax = Inf), fill = "deeppink4", alpha = 0.006) +
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)") + 
  geom_text(aes(x = -Inf, y = -Inf, label=KWpval), data=pval_df[pval_df$name %in% "Progesteron nmol/L",], hjust   = -0.1, vjust   = -0.1, size=3)

# Calculation of PCA for Fig 1E
rownames(meta) <- meta$ID
colnames(meta)
pca <- prcomp(meta[,c(2:9,11)], scale. = T)
pca_x <- as.data.frame(pca$x) %>% rownames_to_column() %>% as.data.frame()
pca_x <- merge(meta, pca_x, by.x="ID", by.y="rowname", sort=FALSE)  %>% as.data.frame()

scree <- fviz_eig(pca)
scree
screedata <- scree$data
sum(screedata$eig[1:5]) 

#Fig 1E: PC1 vs PC2
ggplot(pca_x, aes(y = PC1, x = PC2)) + theme_classic() +
    geom_path(aes(group=SubjectID), color="gray") +
  geom_point(aes(color=factor(VisitMonths))) + 
    labs(x=paste(
      "PC2 ", round(screedata$eig[screedata$dim %in% 2], 2), "%", sep=""), 
      y=paste(
        "PC1 ", round(screedata$eig[screedata$dim %in% 1], 2), "%", sep="")) + coord_flip() +
  scale_color_manual(values = c("#F078AE", "#9290C6", "#273B91")) + labs(color="Visit")

