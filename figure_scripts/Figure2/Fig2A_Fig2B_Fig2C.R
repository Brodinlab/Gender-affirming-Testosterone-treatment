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


# EXT FIG 2A, Heatmap of all cell clusters
myplotd <- pheatmap::pheatmap(as.matrix(res_bind_scaled[,1:48]),
                              color = paletteer::paletteer_d("rcartocolor::PurpOr", n=100, type='continuous'),
                              scale = 'none',
                              labels_row = paste(res_bind_scaled$clusters, "_", res_bind_scaled$flowSOM_level2,' (', round(res_bind_scaled$percentage,1), '%', ')',sep = ''),
                              display_numbers = FALSE,
                              angle_col = 45,
                              breaks = seq(0, 2, 2/90),
                              main = "Median marker expression per cluster")
myplotd 


## Plot results from mixed effects, color by up or down

# FIG2B, 3 months vs baseline
l_g$log2fc_total_directionV2 <- ifelse(l_g$Estimate_V2 > 0, "pos", "neg")
ann_col <- c("pos" = "#F38E6E", "neg" = "#6E51A1")

plot_3m <- ggraph(l_g) +
  geom_edge_link(alpha=.1) + 
  geom_node_point(data = l_g, aes(size= percentage), fill = "#999999", shape=21, alpha = 0.9) +
  geom_node_point(data = subset(l_g, pValueV2 < 0.05), aes(size= percentage, fill = log2fc_total_directionV2), shape=21, alpha = 0.9) +
  scale_fill_manual(values=ann_col)+
  scale_size_continuous(range = c(10, 20)) +
  #geom_node_text(data = l_g, aes(label = flowSOM_level2), size = 3, repel = T, max.overlaps = 1000) +
  geom_node_text(aes(label = clusters), size = 4, color = 'white', show.legend = FALSE, repel = F) +
  theme_graph(base_family = 'Helvetica') +
  labs(title = "FtM: 3m vs baseline")
plot_3m

# FIG2C, 12 months vs baseline
l_g$log2fc_total_directionV3 <- ifelse(l_g$Estimate_V3 > 0, "pos", "neg")

plot_12m <- ggraph(l_g) +
  geom_edge_link(alpha=.1) + 
  geom_node_point(data = l_g, aes(size= percentage), fill = "#999999", shape=21, alpha = 0.9) +
  geom_node_point(data = subset(l_g, pValueV3 < 0.05), aes(size= percentage, fill = log2fc_total_directionV3), shape=21, alpha = 0.9) +
  scale_fill_manual(values=ann_col)+
  scale_size_continuous(range = c(10, 20)) +
  #geom_node_text(data = l_g, aes(label = flowSOM_level2), size = 3, repel = T, max.overlaps = 1000) +
  geom_node_text(aes(label = clusters), size = 4, color = 'white', show.legend = FALSE, repel = F) +
  theme_graph(base_family = 'Helvetica') +
  labs(title = 'FtM: 12m vs baseline')
plot_12m





