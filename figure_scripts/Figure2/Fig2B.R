# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## FIGURE 2B ##
library(tidyverse)

# Load mixed effects results
setwd('./data/Figure2')
res_df = read.csv('221003_CyTOF_flowSOMlevel2_MEM_FtM_visit_age.csv')
colnames(res_df)[2] <- "flowSOM_level2"
res_df$flowSOM_level2[res_df$flowSOM_level2 %in% "Monocytes_NCM"] <- "debris"

# Load cell frequencies
#m_imm.freq <- read.csv("CYTOF_cellfrequencies.csv")
cells_to_plot <- c("pDC", "CD8Tcells_MAIT", "CD8Tcells_TCM_CD24")

# Fig 2B 
m_imm.freq %>% 
  pivot_longer(cols=10:46) %>%
  ggplot(aes(x=VisitMonths, y=value, fill=factor(VisitMonths))) +
  geom_path(aes(group=SubjectID), alpha=0.4) + theme_classic() +
  scale_x_continuous(breaks = c(0,3,12)) +
  facet_wrap(~name, scales = "free") +
  geom_boxplot(outlier.shape = NA) + 
  geom_point() +
  scale_fill_manual(values = c("#EE5A45", "#9EC2C0", "#1F8F89")) + labs(y="%") + 
  theme(legend.position = "none")

sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-apple-darwin17.0 (64-bit)
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
# [1] lubridate_1.9.2 forcats_1.0.0   stringr_1.5.0   dplyr_1.1.1     purrr_1.0.1     readr_2.1.4    
# [7] tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.2   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
# [1] rstudioapi_0.14  magrittr_2.0.3   hms_1.1.3        tidyselect_1.2.0 munsell_0.5.0   
# [6] timechange_0.2.0 colorspace_2.1-0 R6_2.5.1         rlang_1.1.0      fansi_1.0.4     
# [11] tools_4.2.1      grid_4.2.1       gtable_0.3.3     utf8_1.2.3       cli_3.6.1       
# [16] withr_2.5.0      lifecycle_1.0.3  farver_2.1.1     tzdb_0.4.0       vctrs_0.6.1     
# [21] glue_1.6.2       labeling_0.4.2   stringi_1.7.12   compiler_4.2.1   pillar_1.9.0    
# [26] generics_0.1.3   scales_1.2.1     pkgconfig_2.0.3 
