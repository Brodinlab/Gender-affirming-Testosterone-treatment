#Script to generate Fig. 2D of AR/ESR expression across immune cells (Blood atlas, Uhlen et al Science, 2019)
#Author Petter Brodin (petter.brodin@ki.se)

library(tidyverse)
library(ggplot2)
library(paletteer)
library(scales)

setwd("enter source data folder")

d <- read.delim("data/rna_blood_cell.tsv", header=T, sep="\t")
d <- d[d$Blood.cell != "total PBMC", ]

#subset the data on hormone receptors only
horm.rec <- d[d$Gene.name %in% c("AR", "ESR1", "ESR2"), ]

# plot heatmap of hormone receptors in all cell subsets

#define plot theme
spot.theme <- list(
  theme_classic(),
  theme(axis.ticks.x=element_blank(), axis.text.x=element_text(colour="black",size = 20, angle = 45, hjust = 0)),
  theme(axis.ticks.y=element_blank(), axis.text.y=element_text(colour="black", size = 20)),
  theme(axis.title.y = element_blank()),
  theme(axis.title.x = element_blank()),
  theme(axis.line=element_blank()),
  theme(text = element_text(size = 22)),
  theme(legend.position = "none"),
  theme(plot.margin = unit(c(10,10,10,10), "mm")),
  scale_size_continuous(range = c(-0.3, 15)),
  scale_x_discrete(position = "top"))

#create plot with ggplot2
p1 <- ggplot(horm.rec, aes(x=Blood.cell, y=Gene.name, color=NX)) + geom_point(size=10)+ geom_point(aes(fill=NX), colour="black", pch=21, size=15) + 
  spot.theme + scale_fill_gradient(low = "#FFFFFF", high = "#A32F33")
p1