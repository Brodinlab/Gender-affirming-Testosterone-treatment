# script for plotting SIMOA IFNa levels following stimulation
# Author: Petter Brodin (petter.brodin@ki.se)

library(tidyverse)
setwd("./data/Figure2/")

# stimulated PBMCs in vitro and IFN-I analyzed
d.stim <-read.delim("Simoa_GAHT.csv", sep=",", header=T)
stim <- d.stim[d.stim$Visit %in% c("V1", "V2") & d.stim$Sample_type=="PBMC stim supernatant", ] #Focus on Baseline vs 3mo due to few samples at 12mo

stim$Subject.ID <- as.factor(stim$Subject.ID)
stim$Visit <-  as.factor(stim$Visit)

# stimulated panIFNa levels
p1 <- ggplot(stim, aes(x=Visit, y=log10(IFNa.ratio))) + 
  geom_point(size=6) + geom_line(aes(group = Subject.ID)) +
  theme_bw() + ggtitle("IFNa (Stim/unstim (NTC = Non treated ctrl))") +
  facet_wrap(~Group)
p1

# stimulated panIFNb levels
p2 <- ggplot(stim, aes(x=Visit, y=log10(IFNb.ratio))) + 
  geom_point(size=6) + geom_line(aes(group = Subject.ID)) +
  theme_bw() + ggtitle("IFNb (Stim/unstim (NTC = Non treated ctrl))") +
  facet_wrap(~Group)
p2

#-----------------------------------------------------------
#For R848/NTC (unstim) ratios run 1-way t-tests
#IFNa ratio:
pairwise.t.test(stim$IFNa.ratio, stim$Visit, p.adjust.method = "none", paired = T, alternative = "less")

#IFNb ratio:
pairwise.t.test(stim$IFNb.ratio, stim$Visit, p.adjust.method = "none", paired = T, alternative = "less")

