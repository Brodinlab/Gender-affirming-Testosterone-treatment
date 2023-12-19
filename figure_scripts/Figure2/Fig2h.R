# script for plotting SIMOA IFNa levels following stimulation
# Author: Petter Brodin (petter.brodin@ki.se)

library(tidyverse)
setwd("Please insert working directory here to point at directory '../data/Figure2/'")
# stimulated PBMCs in vitro and IFN-I analyzed
d.stim <-read.delim("Simoa_stimulated_cultures_gender_affirming_project_V2.csv", sep=",", header=T)

# run paired t-test Visit 1/2:
r848stim <- d.stim[d.stim$Group=="R848", ]


# stimulated panIFNa levels
p1 <- ggplot(d.stim, aes(x=Visit, y=log10(IFNa.ratio), color=Subject.ID)) + 
  geom_point(size=6) + geom_line(aes(group = Subject.ID)) +
  theme_bw() + ggtitle("IFNa (Stim/unstim (NTC = Non treated ctrl))") +
  facet_wrap(~Group)
p1
