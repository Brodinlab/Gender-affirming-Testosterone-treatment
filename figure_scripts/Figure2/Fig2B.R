# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## FIGURE 2B ##
library(tidyverse)

# Load mixed effects results
res_df = read.csv('221003_CyTOF_flowSOMlevel2_MEM_FtM_visit_age.csv')
colnames(res_df)[2] <- "flowSOM_level2"
res_df$flowSOM_level2[res_df$flowSOM_level2 %in% "Monocytes_NCM"] <- "debris"

# Load cell frequencies
m_imm.freq <- read.csv("CYTOF_cellfrequencies.csv")
cells_to_plot <- c("pDC", "DC_CD141DC", "MDSC_M", "CD8Tcells_TEMRA_CD24", "Monocytes_CM")

# Fig 2B
m_imm.freq %>% 
  pivot_longer(cols=10:46) %>% 
  filter(name %in% cells_to_plot) %>%
  mutate(name=factor(name, levels = cells_to_plot)) %>%
  ggplot(aes(x=VisitMonths, y=value, fill=factor(VisitMonths))) +
  geom_path(aes(group=SubjectID), alpha=0.4) + theme_classic() +
  scale_x_continuous(breaks = c(0,3,12)) +
  facet_wrap(~name, scales = "free", ncol = 2) +
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(values = c("#F38E6E", "#D1BED6", "#6E51A1")) + labs(y="Fraction of cells (%)", x="Time (months)") + 
  theme(legend.position = "none")
