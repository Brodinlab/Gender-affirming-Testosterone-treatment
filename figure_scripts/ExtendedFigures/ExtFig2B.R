# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## Heatmap marker expression per cluster, EXT FIG 2A ##

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

# Load mixed effects results
setwd("./data/ExtFig2/")
res_df = read.csv('221003_CyTOF_flowSOMlevel2_MEM_FtM_visit_age.csv')
colnames(res_df)[2] <- "flowSOM_level2"
res_df$flowSOM_level2[res_df$flowSOM_level2 %in% "Monocytes_NCM"] <- "debris"

# load flowSOM cluster info
setwd("./data/ExtFig2/")
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

# EXT FIG 2B, Heatmap of all cell clusters
myplotd <- pheatmap::pheatmap(as.matrix(res_bind_scaled[,1:48]),
                              color = paletteer::paletteer_d("rcartocolor::PurpOr", n=100, type='continuous'),
                              scale = 'none',
                              labels_row = paste(res_bind_scaled$clusters, "_", res_bind_scaled$flowSOM_level2,' (', round(res_bind_scaled$percentage,1), '%', ')',sep = ''),
                              display_numbers = FALSE,
                              angle_col = 45,
                              breaks = seq(0, 2, 2/90),
                              main = "Median marker expression per cluster")
myplotd 

sessionInfo()
# Platform: x86_64-apple-darwin17.0 (64-bit)
# R version 4.2.1 (2022-06-23)
# Running under: macOS Monterey 12.5
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] paletteer_1.5.0    RColorBrewer_1.1-3 ggthemes_4.2.4     scales_1.2.1       magrittr_2.0.3    
# [6] flowCore_2.10.0    igraph_1.4.1       ggraph_2.1.0       lubridate_1.9.2    forcats_1.0.0     
# [11] stringr_1.5.0      dplyr_1.1.1        purrr_1.0.1        readr_2.1.4        tidyr_1.3.0       
# [16] tibble_3.2.1       ggplot2_3.4.2      tidyverse_2.0.0    vite_0.4.10       
# 
# loaded via a namespace (and not attached):
# [1] ggrepel_0.9.3       Rcpp_1.0.10         mvtnorm_1.1-3       RProtoBufLib_2.10.0
# [5] digest_0.6.31       utf8_1.2.3          ggforce_0.4.1       R6_2.5.1           
# [9] stats4_4.2.1        pillar_1.9.0        rlang_1.1.0         rstudioapi_0.14    
# [13] S4Vectors_0.36.2    pheatmap_1.0.12     polyclip_1.10-4     munsell_0.5.0      
# [17] compiler_4.2.1      pkgconfig_2.0.3     BiocGenerics_0.44.0 parameters_0.20.3  
# [21] insight_0.19.1      tidyselect_1.2.0    gridExtra_2.3       matrixStats_0.63.0 
# [25] graphlayouts_0.8.4  fansi_1.0.4         viridisLite_0.4.1   tzdb_0.4.0         
# [29] withr_2.5.0         prismatic_1.1.1     MASS_7.3-58.3       grid_4.2.1         
# [33] xtable_1.8-4        gtable_0.3.3        lifecycle_1.0.3     bayestestR_0.13.1  
# [37] datawizard_0.7.1    estimability_1.4.1  cli_3.6.1           stringi_1.7.12     
# [41] farver_2.1.1        viridis_0.6.2       generics_0.1.3      vctrs_0.6.1        
# [45] rematch2_2.1.2      tools_4.2.1         Biobase_2.58.0      glue_1.6.2         
# [49] tweenr_2.0.2        hms_1.1.3           emmeans_1.8.5       timechange_0.2.0   
# [53] colorspace_2.1-0    cytolib_2.10.1      tidygraph_1.2.3 

