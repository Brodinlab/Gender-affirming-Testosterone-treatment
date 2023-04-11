# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## CD52 gMFI, EXT FIG 2A ##
library(tidyverse)
library(data.table)
library(stringr)
library(ggthemes)
library(ggridges)

# dfs
all_bind_FtM <- fread('221003_CyTOF_1234batches_ComBat_corrected_FlowSOM_FtM.csv') %>% column_to_rownames(var = "V1")
all_bind_FtM$Subject__level2 <- paste(all_bind_FtM$Subject, all_bind_FtM$flowSOM_level2, sep="__")

gm_mean = function(x, na.rm=TRUE) {
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

summary_markers <- all_bind_FtM %>% 
  select(1:48,60) %>%
  group_by(Subject__level2) %>%
  summarise_all(list(mean=mean, med=median, gMFI=gm_mean))
summary_markers$Subject <- sapply(str_split(summary_markers$Subject__level2, "__"), "[", 1)
summary_markers$level2 <- sapply(str_split(summary_markers$Subject__level2, "__"), "[", 2)
summary_markers$SubjectID <- sapply(str_split(summary_markers$Subject, "_"), "[", 1)
summary_markers$Visit <- sapply(str_split(summary_markers$Subject, "_"), "[", 2)

summary_markers_long <- summary_markers %>% 
  pivot_longer(cols = 2:145, names_to = "marker", values_to = "value") 
summary_markers_long$statistic <- sapply(str_split(summary_markers_long$marker, "_"), "[", 2)
summary_markers_long$marker <- sapply(str_split(summary_markers_long$marker, "_"), "[", 1)
summary_markers_long$Cell <- sapply(str_split(summary_markers_long$Subject__level2, "__"), "[", 2)

unique(summary_markers_long$Cell) %>% sort()
cells_interest <- c("CD4Tcells_Naive", "CD4Tcells_TCM", "CD4Tcells_EM",
                    "CD8Tcells_Naive", "CD8Tcells_TCM", "CD8Tcells_EM")

summary_markers_long %>%
  filter(marker %in% "CD52") %>%
  filter(statistic %in% "gMFI") %>%
  filter(Cell %in% cells_interest) %>%
  mutate(Cell=factor(Cell, levels = cells_interest)) %>%
  ggplot(aes(x=Visit, y=value)) + geom_line(aes(group=SubjectID), alpha=0.5) + 
  geom_point(aes(fill=Visit), shape = 21, color="black") + 
  facet_wrap(~Cell, scales = "free") + theme_base() + labs(y="CD52 gMFI", x=NULL) + 
  scale_fill_manual(values=c("#E58D9E", "#7E81B9", "#2B3C8B"))

