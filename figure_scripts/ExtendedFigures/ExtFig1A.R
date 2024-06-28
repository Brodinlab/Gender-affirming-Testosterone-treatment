# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## Extended Fig1A ##
library(ggplot2)
library(dplyr)
library(tidyverse)
library(factoextra)
library(cowplot)


# meta
setwd("/Users/camilaconsiglio/Library/CloudStorage/OneDrive-LundUniversity/BrodinLab/SRT/Data/Metadata/")
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

# 17A OHP
A <- meta %>% 
  ggplot(aes(x=VisitMonths, y=X17alfa.Ohp.nmol.L)) + 
  #male range
  geom_rect(aes(ymin = 0.85, ymax = 6.28, xmin = -Inf, xmax = Inf), fill = "#1F8F89", alpha = 0.006) +
  #female range
  geom_rect(aes(ymin = 0, ymax = 8.02, xmin = -Inf, xmax = Inf), fill = "#EE5A45", alpha = 0.006) +
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)") 

#Androstendion
B <- meta %>% 
  ggplot(aes(x=VisitMonths, y=Androstendion.nmol.L)) + 
  #male range
  geom_rect(aes(ymin = 1.41, ymax = 6.20, xmin = -Inf, xmax = Inf), fill = "#1F8F89", alpha = 0.006) +
  #female range
  geom_rect(aes(ymin = 1.32, ymax = 7.37, xmin = -Inf, xmax = Inf), fill = "#EE5A45", alpha = 0.006) +
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)") 

#DHEA
C <- meta %>% 
  ggplot(aes(x=VisitMonths, y=DHEA.nmol.L)) + 
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)") 

#DHT
D <- meta %>% 
  ggplot(aes(x=VisitMonths, y=DHT.nmol.L)) + 
  #male range
  #female range
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)") 

#Testosterone
E <- meta %>% 
  ggplot(aes(x=VisitMonths, y=Testosteron.nmol.L)) + 
  #male range
  geom_rect(aes(ymin = 7.98, ymax = 29.14, xmin = -Inf, xmax = Inf), fill = "#1F8F89", alpha = 0.006) +
  #female range
  geom_rect(aes(ymin = 0.33, ymax = 2.05, xmin = -Inf, xmax = Inf), fill = "#EE5A45", alpha = 0.006) +
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)") 

#Estrone
F1 <- meta %>% 
  ggplot(aes(x=VisitMonths, y=Estron.pmol.L)) + 
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)") 

# ExtFig 1A
plot_grid(A,B,C,D,E,F1,ncol=3)

sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.5

# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] cowplot_1.1.1    factoextra_1.0.7 lubridate_1.9.2  forcats_1.0.0    stringr_1.5.0   
# [6] purrr_1.0.1      readr_2.1.4      tidyr_1.3.0      tibble_3.2.1     tidyverse_2.0.0 
# [11] dplyr_1.1.1      ggplot2_3.4.2   

# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.10      pillar_1.9.0     compiler_4.2.1   tools_4.2.1      lifecycle_1.0.3 
# [6] gtable_0.3.3     timechange_0.2.0 pkgconfig_2.0.3  rlang_1.1.0      cli_3.6.1       
# [11] rstudioapi_0.14  ggrepel_0.9.3    withr_2.5.0      generics_0.1.3   vctrs_0.6.1     
# [16] hms_1.1.3        grid_4.2.1       tidyselect_1.2.0 glue_1.6.2       R6_2.5.1        
# [21] fansi_1.0.4      tzdb_0.4.0       farver_2.1.1     magrittr_2.0.3   scales_1.2.1    
# [26] colorspace_2.1-0 labeling_0.4.2   utf8_1.2.3       stringi_1.7.12   munsell_0.5.0   
# [31] crayon_1.5.2
