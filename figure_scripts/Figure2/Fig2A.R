# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## Force-directed graph of FlowSOM cell clusters, FIGURES 2A, 2B, 2C ##

# 0. Load packages ----------
library(vite) 
library(tidyverse)
library(ggraph) 
library(igraph) 
library(flowCore) 
library(tidyverse)
library(magrittr)
library(scales)
library(ggthemes)
library(RColorBrewer)
library(paletteer)

setwd("Please insert working directory here to point at directory '../data/Figure2/'")

# Load mixed effects results
res_df = read.csv('221003_CyTOF_flowSOMlevel2_MEM_FtM_visit_age.csv')
colnames(res_df)[2] <- "flowSOM_level2"
res_df$flowSOM_level2[res_df$flowSOM_level2 %in% "Monocytes_NCM"] <- "debris"

# load flowSOM cluster info
res_bind = read.table('221013_flowsom_clustered_FtM_updatedclusternumbers.txt', header = TRUE, sep = "\t", check.names = FALSE)
res_bind$flowSOM_level2[res_bind$flowSOM_level2 %in% "Monocytes_NCM"] <- "debris"

res_bind <- res_bind %>% dplyr::filter(!flowSOM_level2 %in% c("dead", "unknown_35_CD99hi"))
res_bind_scaled <- res_bind
res_bind_scaled[,1:48] <- sapply(res_bind_scaled[,1:48], scale)

# 1. Create unsupervised graph ----------
common_marker = colnames(res_bind)[1:48]

### after scaling by marker
set.seed(824)

# run vite::get_unsupervised_graph function step by step, from https://rdrr.io/github/ParkerICI/scgraphs/src/R/unsupervised.R#sym-get_unsupervised_graph
G <- NULL
G <- build_graph(res_bind_scaled, common_marker, filtering_T = 5)
for (i in names(res_bind_scaled)) {
  G <- igraph::set.vertex.attribute(G, name = i, value = res_bind_scaled[, i])
}
cc <- igraph::multilevel.community(G)
V(G)$community_id <- as.character(cc$membership)
V(G)$name <- seq_along(V(G))
V(G)$type <- "cluster"
V(G)$Label <- paste("c", V(G)$cellType, sep = "")

l_g = create_layout(G, layout="fr")
l_g$x <- V(G)$x
l_g$y <- V(G)$y
l_g %<>% left_join(res_df)
ann_col <- c("pos"= "#F8766D", "neg" = "#619CFF")
title = 'FtM cell populations - scaled marker'

# FIG 2A
ggraph(l_g) +
  geom_edge_link(alpha=.1) + 
  geom_node_point(aes(size= percentage, fill = flowSOM_level2), shape=21,  alpha = 1) +
  paletteer::scale_fill_paletteer_d("ggsci::default_igv") +
  scale_size_continuous(range = c(10, 20)) +
  #geom_node_text(aes(label = flowSOM_level2), size = 3, repel = T, max.overlaps = 1000) +
  geom_node_text(aes(label = clusters), size = 4, color = 'white', show.legend = FALSE, repel = F) +
  theme_graph(base_family = 'Helvetica') +
  labs(title = title)


