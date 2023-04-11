# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## FIGURE 3A ##
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
olinkX <- read.csv("220718_olink_fourbatches_FtM.csv", row.names = 1)
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

