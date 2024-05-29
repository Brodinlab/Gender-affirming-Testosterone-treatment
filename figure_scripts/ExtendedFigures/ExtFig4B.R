# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## EXT FIG 3B ##
library(OlinkAnalyze)
library(stringr)
library(reshape2)
library(tidyverse)
library(factoextra)
library(ggplot2)
library(cowplot)
library(lme4)
library(MuMIn)
library(cowplot)
library(ggthemes)

# Load data
setwd("./ExtFig4")
Jul22 <- read_NPX(filename = "SexChange_NPX_LOD.xlsx")
Jul22.meta <- read.csv("20220714_OlinkPlateLayout.csv")
Jul22.meta <- Jul22.meta %>% filter(Project %in% "FemaleWBstim")
Jul22 <- Jul22 %>% filter(SampleID %in% Jul22.meta$BarcodeID) 
Jul22.meta$Subject_Group <- paste(Jul22.meta$Subject, Jul22.meta$Group, sep = "__")
Jul22$SampleID<- Jul22.meta$Subject_Group[match(Jul22$SampleID, Jul22.meta$BarcodeID)]
olink <- Jul22[,c(1,5,12)]
olink$Assay <- gsub(" ", "", olink$Assay)
olink$Assay <- gsub("-", "", olink$Assay)
olink <- dcast(olink, SampleID ~ factor(Assay, levels = unique(olink$Assay)), fun.aggregate=function(i) mean(i, na.rm=TRUE))
rownames(olink) <- olink$SampleID

olink <- olink %>% 
  mutate(Subject= sapply(str_split(SampleID, "__"), "[", 1),
         Group= sapply(str_split(SampleID, "__"), "[", 2),
         Stim=sapply(str_split(Group, " "), "[", 1),
         Hormone=paste(sapply(str_split(Group, " "), "[", 2), sapply(str_split(Group, " "), "[", 3))  )

olink$Hormone <- gsub(" NA", "", olink$Hormone)

# Filtering 
olink <- olink %>% filter(Stim %in% "unstim")
olink_sups <- olink %>% filter(!Group %in% "x")
PPvar <- sapply(olink_sups[,2:181], var) %>% sort()
olink_var <- olink %>% 
  select(Subject, Group, Stim, Hormone, names(PPvar[PPvar > 0.02])) # 150 pp

# mixed effects model
olink_var$Hormone <- factor(olink_var$Hormone, levels = c("x", "test", "test enz"))
varlist = colnames(olink_var)[5:(ncol(olink_var))]

MEM_PP <- lapply(varlist, function(x) {
  mod2 = try(lmer(substitute(i ~ Hormone + 
                               (1|Subject), list(i = as.name(x))), 
                  data = olink_var, na.action=na.exclude))
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

# Estimates for fixed effects
Aestimate_test <- lapply(MEM_PP_1, function(f) summary(f)$coefficients[2,1])
Aestimate_testEnz <- lapply(MEM_PP_1, function(f) summary(f)$coefficients[3,1])

## p-values for covariates (fixed-effect)
Atest_test <- lapply(MEM_PP_1, function(f) parameters::p_value(f, method = "wald",)[2,2])
Atest_testEnz <-  lapply(MEM_PP_1, function(f) parameters::p_value(f, method = "wald",)[3,2])

# Median
Amed <- lapply(MEM_PP_1, function(f) summary(f)$residuals) #get residuals
Amed_na <- lapply(Amed, function(f) na.exclude(f)) #exclude NAs from residuals
Amed_calc <- lapply(Amed_na, function(f) median(f)) #calculate median
#hist(as.numeric(Amed_calc))

# R^2 and adjusted
r2 <- lapply(MEM_PP_1, function(f) r.squaredGLMM(f)[2]) 

# Prepare dataframe with extracted info. for downstream use
Atest_data = list(
  Avarlist1, Amed_calc, r2,
  Aestimate_test, Aestimate_testEnz, 
  Atest_test, Atest_testEnz)

names(Atest_data) <- c('proteins', 'residuals', 'Rsq',
                       'Estimate_test', 'Estimate_testEnz',
                       'pValue_test', 'pValue_testEnz')

Atest_final <- as.data.frame(do.call(rbind, Atest_data))
Adf <- data.frame(matrix(unlist(Atest_final), nrow=length(Atest_final), byrow=T), stringsAsFactors = F)
colnames(Adf) <- c('proteins', 'residuals', 'Rsq',
                   'Estimate_test', 'Estimate_testEnz',
                   'pValue_test', 'pValue_testEnz')
Adf[,-1] <- sapply(Adf[,-1], as.numeric)

# RANKL
Adf[Adf$proteins %in% "TRANCE",] #pval test=.006304067
olink %>% 
  mutate(Hormone=factor(Hormone, levels = c("x", "test", "test enz"))) %>%
  pivot_longer(cols = 2:181) %>%
  filter(name %in% c("TRANCE")) %>%
  ggplot(aes(x=Hormone, y=value)) + geom_line(aes(group=Subject)) + 
  geom_point(aes(fill=Hormone), size=6, shape = 21, color="black") + 
  facet_wrap(~name, scales = "free", ncol=3) + theme_base() + labs(y="NPX", x=NULL) + 
  scale_fill_manual(values=c("white", "black", "gray"))

# TNF
Adf[Adf$proteins %in% "TNF",] #ns
olink %>% 
  mutate(Hormone=factor(Hormone, levels = c("x", "test", "test enz"))) %>%
  pivot_longer(cols = 2:181) %>%
  filter(name %in% c("TNF")) %>%
  ggplot(aes(x=Hormone, y=value)) + geom_line(aes(group=Subject)) + 
  geom_point(aes(fill=Hormone), size=6, shape = 21, color="black") + 
  facet_wrap(~name, scales = "free", ncol=3) + theme_base() + labs(y="NPX", x=NULL) + 
  scale_fill_manual(values=c("white", "black", "gray"))

# FGF23, IL17C, OSM, CXCL9
Adf[Adf$proteins %in% "FGF23",] #ns, not a variable feature in this dataset
Adf[Adf$proteins %in% "IL17C",] #ns
Adf[Adf$proteins %in% "OSM",] #ns
Adf[Adf$proteins %in% "CXCL9",] #ns
olink %>% 
  mutate(Hormone=factor(Hormone, levels = c("x", "test", "test enz"))) %>%
  pivot_longer(cols = 2:181) %>%
  filter(name %in% c("FGF23", "IL17C", "OSM", "CXCL9")) %>%
  mutate(name=factor(name, levels = c("FGF23", "IL17C", "OSM", "CXCL9"))) %>%
  ggplot(aes(x=Hormone, y=value)) + geom_line(aes(group=Subject)) + 
  geom_point(aes(fill=Hormone), size=6, shape = 21, color="black") + 
  facet_wrap(~name, scales = "free", ncol=1) + theme_base() + labs(y="NPX", x=NULL) + 
  scale_fill_manual(values=c("white", "black", "gray"))
