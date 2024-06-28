# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## FIGURES 1B, 1C, 1D, 1E ##
library(ggplot2)
library(dplyr)
library(tidyverse)
library(factoextra)


# meta
meta <- read.csv("Fig1_metadata2024.csv", row.names = 1)

# KW pvalues
pval_df <- sapply(meta[,c(2:9,11)], function(i) kruskal.test(i ~ meta$VisitMonths))[3,] %>% data.frame() %>% pivot_longer(1:9)
colnames(pval_df) <- c("name", "KWpval")
pval_df$KWfdr <- p.adjust(pval_df$KWpval, method = "fdr", n = 9)

# plot these with ranges for each hormone
# Reference intervals of 20-40 year age group  from Bae et al 2019: https://doi.org/10.1016/j.jsbmb.2019.105409
# estradiol <1109 pmol/L for women and <158 pmol/L for men; 
# testosterone 0.3-2 nmol/L for women and 8-29 nmol/L for men. 
# NOTE: DHEAS in the publication is conjugated DHEA, not free as we measured.

# Fig 1B: Bioavailable Testosterone
meta %>% 
  ggplot(aes(x=VisitMonths, y=NewBioavailableTestosterone)) + 
  #male range
  geom_rect(aes(ymin = 4.54, ymax = 23.64, xmin = -Inf, xmax = Inf), fill = "#1F8F89", alpha = 0.006) +
  #female range
  geom_rect(aes(ymin = 0, ymax = 1.8, xmin = -Inf, xmax = Inf), fill = "#EE5A45", alpha = 0.006) + 
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)", y="Bioavailable testosterone") 

# Fig 1C: Estradiol
meta %>% 
  ggplot(aes(x=VisitMonths, y=Estradiol.pmol.L)) + 
  #male range
  geom_rect(aes(ymin = 0, ymax = 158.4, xmin = -Inf, xmax = Inf), fill = "#1F8F89", alpha = 0.006) +
  #female range
  geom_rect(aes(ymin = 236, ymax = 1109, xmin = -Inf, xmax = Inf), fill = "#EE5A45", alpha = 0.006) +
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)", y="Estradiol") 


# Fig 1D: Progesterone
meta %>% 
  ggplot(aes(x=VisitMonths, y=Progesteron.nmol.L)) + 
  #male range
  geom_rect(aes(ymin = 0, ymax = 0.65, xmin = -Inf, xmax = Inf), fill = "#1F8F89", alpha = 0.006) +
  #female range
  geom_rect(aes(ymin = 0, ymax = 61.61, xmin = -Inf, xmax = Inf), fill = "#EE5A45", alpha = 0.006) +
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)", y="Progesterone") + 
  geom_text(aes(x = -Inf, y = -Inf, label=KWpval), data=pval_df[pval_df$name %in% "Progesteron nmol/L",], hjust   = -0.1, vjust   = -0.1, size=3)

# Calculation of PCA for Fig 1E
rownames(meta) <- meta$ID
pca <- prcomp(meta[,c(2:9,11)], scale. = T)
pca_x <- as.data.frame(pca$x) %>% rownames_to_column() %>% as.data.frame()
pca_x <- merge(meta, pca_x, by.x="ID", by.y="rowname", sort=FALSE)  %>% as.data.frame()

#scree <- fviz_eig(pca)
#screedata <- scree$data

#Fig 1E PC1 vs PC2
ggplot(pca_x, aes(y = PC1, x = PC2)) + theme_classic() +
  geom_path(aes(group=SubjectID), color="gray") +
  geom_point(aes(color=factor(VisitMonths))) + 
  labs(x=paste(
    "PC2 ", round(screedata$eig[screedata$dim %in% 2], 2), "%", sep=""), 
    y=paste(
      "PC1 ", round(screedata$eig[screedata$dim %in% 1], 2), "%", sep="")) + coord_flip() +
  scale_color_manual(values = c("#EE5A45", "#9EC2C0", "#1F8F89")) + labs(color="Visit")

sessionInfo()
#R version 4.2.1 (2022-06-23)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Monterey 12.5

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

#locale:
#[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] factoextra_1.0.7 lubridate_1.9.2  forcats_1.0.0    stringr_1.5.0    purrr_1.0.1     
#[6] readr_2.1.4      tidyr_1.3.0      tibble_3.2.1     tidyverse_2.0.0  dplyr_1.1.1     
#[11] ggplot2_3.4.2   

#loaded via a namespace (and not attached):
#[1] Rcpp_1.0.10      rstudioapi_0.14  magrittr_2.0.3   hms_1.1.3        tidyselect_1.2.0
#[6] munsell_0.5.0    timechange_0.2.0 colorspace_2.1-0 R6_2.5.1         rlang_1.1.0     
#[11] fansi_1.0.4      tools_4.2.1      grid_4.2.1       gtable_0.3.3     utf8_1.2.3      
#[16] cli_3.6.1        withr_2.5.0      lifecycle_1.0.3  crayon_1.5.2     farver_2.1.1    
#[21] tzdb_0.4.0       vctrs_0.6.1      ggrepel_0.9.3    glue_1.6.2       labeling_0.4.2  
#[26] stringi_1.7.12   compiler_4.2.1   pillar_1.9.0     generics_0.1.3   scales_1.2.1    
#[31] pkgconfig_2.0.3 
