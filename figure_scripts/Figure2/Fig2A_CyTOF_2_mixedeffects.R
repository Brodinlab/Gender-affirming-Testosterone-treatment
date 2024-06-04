# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## CyTOF population mixed effects modelling, for FIGURES 2A, 2B, 2C ##
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

# 4 batches
# EXP20CL8796
# EXP20DG3622
# EXP21DG3628
# EXP22DG3665

setwd("Please insert working directory here to point at directory './data/Figure2/'")

# Load frequency data generated with flowSOM 
imm.freq <-  read.csv('221003_flowsom_frequency_level2_populations_FtM.csv')

# load metadata
meta <- read.csv("SRT_allmetadata.csv", row.names = 1) 

# Merge
m_imm.freq <- merge(meta, imm.freq, by="FCSfile") 

# Filter FtM
m_imm.freq <- m_imm.freq %>% filter(Transition %in% "FtM") 

# Filter NA for age
m_imm.freq <- m_imm.freq %>% filter(!is.na(Age)) 

# Visit as factors
m_imm.freq$Visit <- as.factor(m_imm.freq$Visit)

# Filter out dead and unknown
m_imm.freq <- m_imm.freq %>% select(!c("dead", "unknown_35_CD99hi"))

# MEM lmer
varlist = colnames(m_imm.freq)[c(14:48)]

MEM_A <- lapply(varlist, function(x) {
  mod2 = try(lmer(substitute(i ~ Visit + Age + 
                               (1|SubjectID), list(i = as.name(x))), 
                  data = m_imm.freq, na.action=na.exclude))
  if(isTRUE(class(mod2)=='try-error')) {return(NULL)} else{return(mod2)} })

MEM_Cells_A1 = as.list(MEM_A)
MEM_Cells_A1[sapply(MEM_Cells_A1, is.null)] <- NULL # Remove NULL models that failed to converge

# Extract info from model

# Cell
Avarlist1 <- lapply(MEM_Cells_A1, function(f) summary(f)$call[2])
Avarlist1 <- sapply(str_split(string = Avarlist1, pattern = "~"), `[`, 1)
Avarlist1 <- gsub(pattern = "\\(", replacement = "", Avarlist1)
Avarlist1 <- sub("\\s+$", "", Avarlist1) %>% as.character()

# Correlation coefficients for fixed effects
# Corr coef V2
estimate_visit2 <- lapply(MEM_Cells_A1, function(f) summary(f)$coefficients[2,1])
# Corr coef V3
estimate_visit3 <- lapply(MEM_Cells_A1, function(f) summary(f)$coefficients[3,1])
# Corr coef Age
estimate_age <- lapply(MEM_Cells_A1, function(f) summary(f)$coefficients[4,1]) 

# pValues
#  V2
Atest_V2 <- lapply(MEM_Cells_A1, function(f) parameters::p_value(f, method = "wald",)[2,2])
#  V3
Atest_V3 <- lapply(MEM_Cells_A1, function(f) parameters::p_value(f, method = "wald",)[3,2])
# Age
Atest_age <- lapply(MEM_Cells_A1, function(f) parameters::p_value(f, method = "wald",)[4,2])

# Residuals
Amed <- lapply(MEM_Cells_A1, function(f) summary(f)$residuals) 
Amed_na <- lapply(Amed, function(f) na.exclude(f)) 
Amed_calc <- lapply(Amed_na, function(f) median(f)) 
#hist(unlist(Amed_na))

# R^2 and adjusted
r2 <- lapply(MEM_Cells_A1, function(f) r.squaredGLMM(f)[2]) 

# prepare df
Atest_data = list(Avarlist1, r2, Amed_calc, estimate_visit2, estimate_visit3, estimate_age, 
                  Atest_V2, Atest_V3, Atest_age)
names(Atest_data) <- c('Cells', 'r2','residuals','Estimate_V2', 'Estimate_V3','Estimate_Age', 
                       'pValueV2',  'pValueV3','pValueAge')
Atest_final <- as.data.frame(do.call(rbind, Atest_data))
Adf <- data.frame(matrix(unlist(Atest_final), nrow=length(Atest_final), byrow=T), stringsAsFactors = F)
colnames(Adf) <- c('Cells', 'r2','residuals','Estimate_V2', 'Estimate_V3','Estimate_Age', 
                   'pValueV2',  'pValueV3','pValueAge')
Adf[,-1] <- sapply(Adf[,-1], as.numeric)

Adf$Cells <- factor(Adf$Cells, levels = Adf$Cells[order(Adf$Estimate_V3)])

# Save df for plotting
write.csv(Adf, "221003_CyTOF_flowSOMlevel2_MEM_FtM_visit_age.csv")

# Plots
plot_v2 <- Adf %>%
  ggplot( aes(x=Estimate_V2, y=Cells)) + geom_point(aes(size=-log(pValueV2),color= pValueV2 < 0.05)) + # 
  geom_vline(xintercept = 0, alpha=0.3, linetype="dotted") + scale_color_manual(values = c("grey", "#e11f28")) + 
  scale_x_continuous(limits = c(-4.2, 4.2)) +
  labs(x="Beta coefficient for V2", y=NULL, size=NULL)  + theme_bw() + theme(legend.position = "none")

plot_v3 <- Adf %>%
  ggplot( aes(x=Estimate_V3, y=Cells)) + geom_point(aes(size=-log(pValueV3),color= pValueV3 < 0.05)) + # 
  geom_vline(xintercept = 0, alpha=0.3, linetype="dotted") + scale_color_manual(values = c("grey", "#e11f28")) + 
  scale_x_continuous(limits = c(-4.2, 4.2)) +
  labs(x="Beta coefficient for V3", y=NULL, size=NULL)  + theme_bw() + theme(legend.position = "none")

plot_grid(plot_v2, plot_v3, ncol=2)
