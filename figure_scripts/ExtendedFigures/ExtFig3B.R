# script for plotting SIMOA IFNa results
library(tidyverse)
setwd("./data/ExtFigure3")

# ex vivo plasma sample SIMOA analysis
d <-read.delim("Simoa_GAHT.csv", sep=",", header=T)
d$Visit[d$Visit=="V1"] <- 0
d$Visit[d$Visit=="V2"] <- 3
d$Visit[d$Visit=="V3"] <- 12

d$Subject.ID <- as.factor(d$Subject.ID)
d$Visit <- as.numeric(d$Visit)
d$Sample_type <- as.factor(d$Sample_type)

d.in.vivo <- d[d$Sample_type=="Plasma", ]

#panIFNa levels in vivo
p1 <- ggplot(d.in.vivo, aes(x=Visit, y=log10(panIFNa..fg.ml.), group = Subject.ID)) + 
  geom_point(size=6) + geom_line(aes(group = Subject.ID)) +
  theme_bw() + ggtitle("pan IFNa") + xlab("Month") + ylab("pan-IFNa fg/ml")
p1

aov.pan.IFNa <- aov(formula = panIFNa..fg.ml. ~ Visit , data = d.in.vivo)
summary(aov.pan.IFNa)

#IFNb levels in vivo
p2 <- ggplot(d.in.vivo, aes(x=Visit, y=log10(IFNb..pg.ml.), group = Subject.ID)) + 
  geom_point(size=6) + geom_line(aes(group = Subject.ID)) +
  theme_bw() + ggtitle("IFNb") + xlab("Month") + ylab("IFNb pg/ml")
p2

aov.IFNb <- aov(formula = IFNb..pg.ml. ~ Visit , data = d.in.vivo)
summary(aov.IFNb)
