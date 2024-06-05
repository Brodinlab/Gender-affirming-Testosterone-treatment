#plot wbc values for Extended Data Fig 2a.
library(tidyverse)
library(rstatix)

setwd("./data/Ext_Figure2")
d <- read.delim("GAHT_wbc.csv", header=T, sep=",")
str(d)

p1 <- ggplot(d, aes(x=month, y=wbc, group=studyID)) + ylim(0,10) +
  geom_point(size=4) + theme_bw() + geom_line()#+ ylim(-1,1) + theme_bw()
p1

#run within-subjects ANOVA
res.aov <- anova_test(data = d, dv = wbc, wid = studyID, within = month)
get_anova_table(res.aov)

