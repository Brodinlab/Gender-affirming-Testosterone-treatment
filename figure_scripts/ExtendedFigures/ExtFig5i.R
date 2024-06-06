# script for plotting FACS IFN-g levels following blood stimulation
# after pretreatment of female blood with DHT, Enz, Fulvestrant.
library(tidyverse)
library(ggpubr)
setwd("./data/ExtFig5")

# stimulated blood in vitro and TNFa, IFN-I analyzed by SIMOA (Pasteur)
d.stim <-read.delim("FACS_GAHT_PMA.csv", sep=";", header=T)

#remove Donor 5 diagnosed with PCOS after this study.
d.stim <- d.stim[d.stim$donor!="Healthy_Donor5", ]
d.stim$pretreatment <- factor(d.stim$pretreatment, levels = c("untreated", "DHT", "DHT_ENZ", "Fulvestrant")) # reorder X-label

# %IFNg+ by pre-treatment condition after PMA stim
p1 <- ggplot(d.stim, aes(x=pretreatment, y=X.IFNg, color=donor)) + 
  geom_point(size=6) + geom_line(aes(group = donor)) + 
  facet_wrap(~cell) + theme_bw() + ggtitle("%IFNg")
p1

#running 2-way ANOVA on all cell types
res.aov2.perc <- aov(X.IFNg ~ pretreatment + donor + cell + pretreatment:cell, data = d.stim)
summary(res.aov2.perc)
TukeyHSD(res.aov2.perc)

#subset on NK cells and run 2-way ANOVA:
d.stim.nk <- d.stim[d.stim$cell=="NK", ]
res.aov2.perc.nk <- aov(X.IFNg ~ pretreatment + donor, data=d.stim.nk)
summary(res.aov2.perc.nk)
TukeyHSD(res.aov2.perc.nk)

#subset on CD8 T cells and run 2-way ANOVA:
d.stim.cd8 <- d.stim[d.stim$cell=="CD8T", ]
res.aov2.perc.cd8 <- aov(X.IFNg ~ pretreatment + donor, data=d.stim.cd8)
summary(res.aov2.perc.cd8)
TukeyHSD(res.aov2.perc.cd8)

#subset on CD4 T cells and run 2-way ANOVA:
d.stim.cd4 <- d.stim[d.stim$cell=="CD4T", ]
res.aov2.perc.cd4 <- aov(X.IFNg ~ pretreatment + donor, data=d.stim.cd4)
summary(res.aov2.perc.cd4)
TukeyHSD(res.aov2.perc.cd4)


#--------------- Plot MFI IFNg apart from % positive cells------
# IFNg MFI by pre-treatment condition after PMA stim
p2 <- ggplot(d.stim, aes(x=pretreatment, y=log10(IFNg_MFI), color=donor)) + 
  geom_point(size=6) + geom_line(aes(group = donor)) + 
  facet_wrap(~cell) + theme_bw() + ggtitle("IFNg MFI")
p2
