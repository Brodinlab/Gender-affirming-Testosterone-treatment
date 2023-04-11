library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(rstatix)

setwd("Set working directory here")
info <- read_csv('./data/frequency table/info.csv')

# 4a

data_freq <- read_csv('./data/frequency table/CyTOF flowSOM level1.csv')

d <- data_freq %>% 
  select(-n) %>% 
  tidyr::pivot_wider(names_from = flowSOM_level1, values_from = freq, values_fill = 0) %>%
  mutate(CD4CD8_ratio = CD4Tcells/CD8Tcells) %>%
  tidyr::pivot_longer(cols = !Subject, names_to = 'flowSOM_level1', values_to = 'freq') %>%
  left_join(info, by='Subject')

tb <- d %>% filter(flowSOM_level1 == 'CD4CD8_ratio') %>% mutate(timepoint_group = case_when(
  Visit == 1 ~ 'baseline',
  Visit == 2 ~ '3 mo',
  Visit == 3 ~ '12 mo'
)) %>% ungroup() %>% pairwise_t_test(freq ~ timepoint_group, paired = T)

ggplot(d %>% filter(flowSOM_level1 == 'CD4CD8_ratio'), aes(x=timepoint, y=freq, group=SubjectID)) + 
  geom_point() + geom_path() +
  ylab('CD4T to CD8T ratio') + 
  theme_bw() +
  scale_x_continuous(breaks=c(0,3, 6, 9, 12))

# 4bcg
data_freq <- read_csv('./data/frequency table/CyTOF flowSOM level2.csv')

d <- data_freq %>% filter(flowSOM_level1 == 'CD4Tcells') %>%
  mutate(custom_subpop = case_when(
    flowSOM_level2 == 'CD4Tcells_Naive' ~ 'Naive CD4',
    flowSOM_level2 %in% c('CD4Tcells_EM', 'CD4Tcells_TCM', 'CD4Tcells_TCM_CD127', 'CD4Tcells_TEMRA') ~ 'Memory CD4',
    flowSOM_level2 == 'CD4Tcells_Treg' ~ 'Treg'
  )) %>% 
  group_by(Subject, custom_subpop) %>% 
  summarise(freq = sum(freq)) %>%
  ungroup() %>% 
  tidyr::complete(Subject, custom_subpop, fill = list(freq=0)) %>%
  group_by(Subject) %>%
  mutate(freq_of_prev = freq/sum(freq)) %>%
  left_join(info, by='Subject')

tb <- d %>% mutate(timepoint_group = case_when(
  Visit == 1 ~ 'baseline',
  Visit == 2 ~ '3 mo',
  Visit == 3 ~ '12 mo'
)) %>% 
  group_by(custom_subpop) %>% pairwise_t_test(freq ~ timepoint_group, paired = T)

# 4b
ggplot(d %>% filter(custom_subpop == 'Naive CD4'), aes(x=timepoint, y=freq, group=SubjectID)) +
  geom_point() + geom_path() + 
  theme_bw() +
  scale_x_continuous(breaks=c(0,3, 6, 9, 12)) +
  ggtitle('Naive CD4')

# 4f
ggplot(d %>% filter(custom_subpop == 'Treg'), aes(x=timepoint, y=freq, group=SubjectID)) +
  geom_point() + geom_path() + 
  theme_bw() +
  scale_x_continuous(breaks=c(0,3, 6, 9, 12)) +
  ggtitle('Treg')

# 4c
d <- data_freq %>% filter(flowSOM_level1 == 'CD8Tcells') %>%
  mutate(custom_subpop = case_when(
    flowSOM_level2 %in% c('CD8Tcells_Naive', 'CD8Tcells_Naive_CD24') ~ 'Naive CD8',
    flowSOM_level2 %in% c('CD8Tcells_EM', 'CD8Tcells_TCM', 'CD8Tcells_TCM_CD24', 'CD8Tcells_TEMRA', 'CD8Tcells_TEMRA_CD24', 'CD8Tcells_TEMRA_CD57') ~ 'Memory CD8',
    flowSOM_level2 == 'CD8Tcells_MAIT' ~ 'MAIT CD8'
  )) %>% 
  group_by(Subject, custom_subpop) %>% 
  summarise(freq = sum(freq)) %>%
  ungroup() %>% 
  tidyr::complete(Subject, custom_subpop, fill = list(freq=0)) %>%
  group_by(Subject) %>%
  mutate(freq_of_prev = freq/sum(freq)) %>%
  left_join(info, by='Subject')

tb <- d %>% mutate(timepoint_group = case_when(
  Visit == 1 ~ 'baseline',
  Visit == 2 ~ '3 mo',
  Visit == 3 ~ '12 mo'
)) %>% 
  group_by(custom_subpop) %>% pairwise_t_test(freq ~ timepoint_group, paired = T)

ggplot(d %>% filter(custom_subpop == 'Naive CD8'), aes(x=timepoint, y=freq, group=SubjectID)) +
  geom_point() + geom_path() + 
  theme_bw() +
  scale_x_continuous(breaks=c(0,3, 6, 9, 12)) +
  ggtitle('Naive CD8')

