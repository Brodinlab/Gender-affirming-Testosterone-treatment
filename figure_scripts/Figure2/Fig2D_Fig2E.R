library(zellkonverter)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(tibble)
library(wesanderson)

setwd("Please insert working directory here to point at directory './data/Figure2/'")

sub_name <- 'pDC'
file_path = 'pDC_sample5000.h5ad'

datah5 = readH5AD(file = file_path)
data = t(assay(datah5)) %>% as_tibble()
meta = colData(datah5) %>% 
  as_tibble() %>% 
  select(Visit, SubjectID, Subject, leiden)
data_anno <- data %>% 
  add_column(leiden = meta$leiden,
             Visit = meta$Visit,
             SubjectID = meta$SubjectID)
# CD81 for pDC
lapply(unique(data_anno$SubjectID), function(x) {
  d <- data_anno %>% filter(SubjectID == x) %>% select(CD81, Visit)
  ggplot(d, 
         aes(x = CD81, y = Visit, fill = Visit)) +
    geom_density_ridges(scale = 4, rel_min_height=.01) +
    scale_fill_manual(values = wes_palette("Rushmore1", nlevels(meta$Visit), type = "continuous")) +
    theme_ridges() + 
    theme(legend.position = "none") +
    xlim(quantile(d$marker_expression, 0.01),
         quantile(d$marker_expression, 0.99)) +
    labs(title = x)
}) %>%
  ggarrange(plotlist = ., ncol = 3, nrow = 6) %>%
  ggexport(filename = file.path('./figure_scripts/Figure2/output figures', paste0(sub_name, '_CD81_expression_by_SubjectID.pdf')),
           width = 10,
           height = 20)

lapply(unique(data_anno$SubjectID), function(x) {
  d <- data_anno %>% filter(SubjectID == x) %>% select(CD81, Visit)
  ggplot(d, aes(x=Visit, y=CD81, fill=Visit)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size=0.3, alpha=0.1) + 
    scale_fill_manual(values = wes_palette("Rushmore1", nlevels(meta$Visit), type = "continuous")) +
    stat_compare_means(comparisons = list(c("Visit1", "Visit2"), c("Visit1", "Visit3"), c("Visit2", "Visit3")),) +
    labs(title = x)
}) %>% ggarrange(plotlist = ., ncol = 3, nrow = 6) %>% 
  ggexport(filename = file.path('./figure_scripts/Figure2/output figures', paste0(sub_name, '_CD81_expression_by_SubjectID_stats.pdf')),
           width = 10,
           height = 20) 


ggplot(data_anno, aes(x=Visit, y=CD81, fill=Visit)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.3, alpha=0.1) + 
  scale_fill_manual(values = wes_palette("Rushmore1", nlevels(meta$Visit), type = "continuous")) +
  stat_compare_means(comparisons = list(c("Visit1", "Visit2"), c("Visit1", "Visit3"), c("Visit2", "Visit3")),) +
  labs(title = 'pDC from all donors') +
  scale_x_discrete(labels=c('Baseline', '3 mo.', '12 mo.')) +
  theme_bw() +
  theme(legend.position="none") +
  scale_y_continuous(trans='log1p')
ggsave(filename = file.path('./figure_scripts/Figure2/output figures', 'pDC_CD81_expression_all_donors_log1p.pdf'),
       width=5, height = 15)





