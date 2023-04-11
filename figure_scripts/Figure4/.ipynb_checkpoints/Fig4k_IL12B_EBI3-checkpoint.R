# script for plotting IL23, IL12B and EIB3 mRNA counts in CD4 T cells before/after testosterone treatment
# Author: Petter Brodin (petter.brodin@ki.se)
library(tidyverse)
library(ggthemes)

# IL23A and STAT3 mRNA in V1 and V2 CD4 T cells (unstimulated and not Tregs). 
d_V1 <-read.delim("CD4T_NTCv1_EBI3_IL12B.csv", sep=",", header=T)
d_V2 <-read.delim("CD4T_NTCv2_EBI3_IL12B.csv", sep=",", header=T)

d_V1$name = 'Baseline'
d_V2$name = '3 Months'

d_V1 <- subset(d_V1, select = -X)
d_V1 <- d_V1 %>% gather(-name, key = "var", value = "value")

d_V2 <- subset(d_V2, select = -X)
d_V2 <- d_V2 %>% gather(-name, key = "var", value = "value")

d_V1$name <- paste(d_V1$name,d_V1$var)
d_V2$name <- paste(d_V2$name,d_V2$var)

plotframe <- full_join(d_V1, d_V2)
plotframe <- subset(plotframe, select = -var)

#plot counts with jitter noise
p1 <- ggplot(plotframe, aes(x=name, y=value)) + 
  geom_jitter(size=.8, width = 0.4, height = 0.1) +
  theme_tufte() + ggtitle("mRNA expression in CD4 T") + coord_flip()
p1
