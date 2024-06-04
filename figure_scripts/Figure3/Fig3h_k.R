# plots of Nanostring data after DHT, Fulvestrant pretreatment and stimulation of blood. 

library(tidyverse)
library(ggpubr)
library(Hmisc)
library(stringr)
library(ggrepel)

#load data and preprocess
setwd("./data/Figure3")
d <- as.data.frame(read.csv("Nanostring_all.csv", header=T))

#Switch to long format for plotting
d <- d %>% separate(`X`, c('stim', 'hormone.pretreatment','batch'), sep="_")

pivoted.dat <- pivot_longer(d, cols=!c('stim', 'hormone.pretreatment','batch'), names_to = "gene", values_to = "count")

#Change gene formats to factors and numbers
pivoted.dat$stim <- as.factor(pivoted.dat$stim) 
pivoted.dat$hormone.pretreatment <- as.factor(pivoted.dat$hormone.pretreatment) 
pivoted.dat$batch <- as.factor(pivoted.dat$batch)
pivoted.dat$gene <- as.factor(pivoted.dat$gene)

#reorder levels for plotting
levels(d$hormone.pretreatment) <- c("untreated","dht", "dht-enze","fulve")

#Z-score transform gene counts per batch and gene.
dat_z = pivoted.dat %>%
  group_by(batch, gene) %>%
  mutate(z_score = scale(count))
pivoted.dat <- dat_z

#------------------------------------------------
#For paper fig individual gene plots are selected
plot.gene <- dat_z[dat_z$gene %in% c("SOCS3"), ]

plot.single.genes <- ggplot(plot.gene, aes(x=as.factor(hormone.pretreatment), y=(z_score), group=hormone.pretreatment)) + 
  geom_boxplot(fill="grey", alpha=0.5, na.rm = T, outliers=F) + geom_point(size=1, alpha=0.8) + 
  theme_bw() + xlab("") + ylab("Gene count (Z-score)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(gene ~ stim, scales="fixed") +
  scale_colour_brewer(type="qual", palette = 2)
plot.single.genes

#--------------------------------------------------------------------------------------
#calculate Anova for selected gene above and LPS-stim to compare hormone pre-treatments
lps <- plot.gene[plot.gene$stim=="lps", ]
pairwise.t.test(lps$z_score, lps$hormone.pretreatment, p.adjust.method = "none", alternative = "less")

#-----------------
#calculate Anova for selected gene above and R848-stim to compare hormone pre-treatments
r848 <- plot.gene[plot.gene$stim=="r848", ]
pairwise.t.test(r848$z_score, r848$hormone.pretreatment, p.adjust.method = "none", alternative = "less")

#-----------------
#calculate Anova for selected gene above and R848-stim to compare hormone pre-treatments
unstim <- plot.gene[plot.gene$stim=="ntc", ]
pairwise.t.test(unstim$z_score, unstim$hormone.pretreatment, p.adjust.method = "none", alternative = "less")
 
