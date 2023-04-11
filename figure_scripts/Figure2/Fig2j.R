# script for plotting SIMOA IFNa results
# Author: Petter Brodin (petter.brodin@ki.se)
library(tidyverse)

setwd("Insert working directory here")

# ex vivo plasma sample SIMOA analysis
d.in.vivo <-read.delim("./data/SIMOA/Simoa_plasma_gender_affirming_project.csv", sep=",", header=T)[, -c(7:10)]

#panIFNa levels
p <- ggplot(d.in.vivo, aes(x=Visit, y=log(panIFNa))) + 
  geom_point(size=6) + geom_line(aes(group = Subject.ID)) +
  theme_bw() + ggtitle("pan IFNa") + xlab("")
p
