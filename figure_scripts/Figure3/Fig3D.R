# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## FIGURE 3D ##
library(OlinkAnalyze)
library(stringr)
library(reshape2)
library(tidyverse)
library(factoextra)
library(ggplot2)
library(cowplot)
library(MuMIn)
library(lme4)
library(ggthemes)
library(dplyr)

setwd("Please insert working directory here to point at directory '/data/Figure3/'")

# Load df
olinkX <- read.csv("220718_olink_fourbatches_FtM_renamed.csv", row.names = 1)
meta <- read.csv("Fig1_metadata.csv", row.names = 1) 
m_olink <- merge(meta, olinkX, by.x="Subject", by.y="row.names") 
colnames(m_olink)

# prep df
sapply(m_olink[,12:ncol(m_olink)], var) %>% sort()
m_olink <- m_olink %>% select(!"IL4")
m_olink <- m_olink %>% filter(!is.na(Age)) 
m_olink$Visit <- as.factor(m_olink$Visit)

# Mixed effects model
varlist = colnames(m_olink)[12:ncol(m_olink)]

MEM_PP <- lapply(varlist, function(x) {
  mod2 = try(lmer(substitute(i ~ Visit + Age + OlinkBatch +
                               (1|SubjectID), list(i = as.name(x))), 
                  data = m_olink, na.action=na.exclude))
  if(isTRUE(class(mod2)=='try-error')) {return(NULL)} else{return(mod2)}
})

# Model check
MEM_PP_1 = as.list(MEM_PP)
MEM_PP_1[sapply(MEM_PP_1, is.null)] <- NULL 

# Model info extraction

# Mixed-effect model expression - extract PP ID
Avarlist1 <- lapply(MEM_PP_1, function(f) summary(f)$call[2])
Avarlist1 <- sapply(str_split(string = Avarlist1, pattern = "~"), `[`, 1)
Avarlist1 <- gsub(pattern = "\\(", replacement = "", Avarlist1)
Avarlist1 <- gsub(pattern = "\\ ", replacement = "", Avarlist1)

## Bcoefs for fixed effects
Aestimate_V2 <- lapply(MEM_PP_1, function(f) summary(f)$coefficients[2,1])
Aestimate_V3 <- lapply(MEM_PP_1, function(f) summary(f)$coefficients[3,1])
estimate_age <- lapply(MEM_PP_1, function(f) summary(f)$coefficients[4,1])

## p-values for covariates (fixed-effect)
# V2
Atest_V2 <- lapply(MEM_PP_1, function(f) parameters::p_value(f, method = "wald",)[2,2])
# V3
Atest_V3 <- lapply(MEM_PP_1, function(f) parameters::p_value(f, method = "wald",)[3,2])
# Age
Atest_pAge <-  lapply(MEM_PP_1, function(f) parameters::p_value(f, method = "wald",)[4,2])


# Median residuals
Amed <- lapply(MEM_PP_1, function(f) summary(f)$residuals) #get residuals
Amed_na <- lapply(Amed, function(f) na.exclude(f)) #exclude NAs from residuals
Amed_calc <- lapply(Amed_na, function(f) median(f)) #calculate median

# R^2 
r2 <- lapply(MEM_PP_1, function(f) r.squaredGLMM(f)[2]) 

# Prepare dataframe with extracted info. for downstream use
Atest_data = list(
  Avarlist1, Amed_calc, r2,
  Aestimate_V2, Aestimate_V3, estimate_age, 
  Atest_V2, Atest_V3, Atest_pAge)
names(Atest_data) <- c('proteins', 'residuals', 'Rsq',
                       'Estimate_V2', 'Estimate_V3', 'Estimate_Age', 
                       'pValueV2', 'pValueV3', 'pValueAge')

Atest_final <- as.data.frame(do.call(rbind, Atest_data))
Adf <- data.frame(matrix(unlist(Atest_final), nrow=length(Atest_final), byrow=T), stringsAsFactors = F)
colnames(Adf) <- c('proteins', 'residuals', 'Rsq',
                   'Estimate_V2', 'Estimate_V3', 'Estimate_Age', 
                   'pValueV2', 'pValueV3', 'pValueAge')
Adf[,-1] <- sapply(Adf[,-1], as.numeric)
Adf$proteins <- factor(Adf$proteins, levels = Adf$proteins[order(Adf$Estimate_V2)])

plot_v2 <- Adf %>%
  filter(proteins %in% c(Adf$proteins[Adf$pValueV2 < 0.05], Adf$proteins[Adf$pValueV3 < 0.05])) %>%
  ggplot( aes(x=Estimate_V2, y=proteins)) + geom_point(size=4, aes(color= pValueV2 < 0.05)) + # 
  geom_vline(xintercept = 0, alpha=0.3, linetype="dotted") + scale_color_manual(values = c("grey", "black")) + scale_x_continuous(limits = c(-1, 1)) +
  labs(x="B coefficient for 3 months", y=NULL, size=NULL)  + theme_classic() + theme(legend.position = "none")
plot_v2 #Fig 3A

sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.5
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] ggthemes_4.2.4     lme4_1.1-32        Matrix_1.5-4       MuMIn_1.47.5       cowplot_1.1.1     
# [6] factoextra_1.0.7   lubridate_1.9.2    forcats_1.0.0      dplyr_1.1.1        purrr_1.0.1       
# [11] readr_2.1.4        tidyr_1.3.0        tibble_3.2.1       ggplot2_3.4.2      tidyverse_2.0.0   
# [16] reshape2_1.4.4     stringr_1.5.0      OlinkAnalyze_3.6.2
# 
# loaded via a namespace (and not attached):
# [1] ggrepel_0.9.3       Rcpp_1.0.10         mvtnorm_1.1-3       lattice_0.21-8     
# [5] utf8_1.2.3          R6_2.5.1            cellranger_1.1.0    plyr_1.8.8         
# [9] backports_1.4.1     stats4_4.2.1        pillar_1.9.0        rlang_1.1.0        
# [13] readxl_1.4.2        rstudioapi_0.14     minqa_1.2.5         car_3.1-2          
# [17] nloptr_2.0.3        labeling_0.4.2      splines_4.2.1       munsell_0.5.0      
# [21] broom_1.0.4         compiler_4.2.1      numDeriv_2016.8-1.1 pkgconfig_2.0.3    
# [25] parameters_0.20.3   lmerTest_3.1-3      insight_0.19.1      tidyselect_1.2.0   
# [29] fansi_1.0.4         tzdb_0.4.0          withr_2.5.0         MASS_7.3-58.3      
# [33] grid_4.2.1          nlme_3.1-162        xtable_1.8-4        gtable_0.3.3       
# [37] lifecycle_1.0.3     magrittr_2.0.3      bayestestR_0.13.1   scales_1.2.1       
# [41] datawizard_0.7.1    zip_2.2.2           estimability_1.4.1  cli_3.6.1          
# [45] stringi_1.7.12      carData_3.0-5       farver_2.1.1        generics_0.1.3     
# [49] vctrs_0.6.1         boot_1.3-28.1       tools_4.2.1         glue_1.6.2         
# [53] hms_1.1.3           emmeans_1.8.5       abind_1.4-5         timechange_0.2.0   
# [57] colorspace_2.1-0    rstatix_0.7.2 
