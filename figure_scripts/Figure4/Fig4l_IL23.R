# script for plotting IL23A mRNA etc in CD4 T cells
# Author: Petter Brodin (petter.brodin@ki.se)
library(tidyverse)
library(ggthemes)

# IL23A and STAT3 mRNA in V1 and V2 CD4 T cells (unstimulated and not Tregs). 
d_V1_IL23A_STAT3 <-read.delim("CD4T_NTCv1_IL23A_STAT3.csv", sep=",", header=T)[, -1]

#plot counts with jitter noise
p1 <- ggplot(d_V1_IL23A_STAT3, aes(x=STAT3, y=IL23A)) + 
  geom_jitter(size=3, width = 0.5, height = 0.5, alpha=0.4) + xlim(-1, 7) + ylim(-1, 7) +
  theme_tufte() + ggtitle("IL23A vs STAT3 mRNA expression in CD4 T at baseline")
p1

d_V2_IL23A_STAT3 <-read.delim("CD4T_NTCv2_IL23A_STAT3.csv", sep=",", header=T)[, -1]


#plot counts with jitter noise - 3 months
p2 <- ggplot(d_V2_IL23A_STAT3, aes(x=STAT3, y=IL23A)) + 
  geom_jitter(size=3, width = 0.5, height = 0.5, alpha=0.4) + xlim(-1, 7) + ylim(-1, 7) +
  theme_tufte() + ggtitle("IL23A vs STAT3 mRNA expression in CD4 T at 3mo.")
p2
