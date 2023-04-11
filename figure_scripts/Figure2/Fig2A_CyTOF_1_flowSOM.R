# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## CyTOF flowSOM, for FIGURES 2A, 2B, 2C ##


# 1_flowSOM
library(FlowSOM)
library(flowCore)
library(tidyverse)
library(magrittr)
library(ggrepel)
library(cowplot)
library(data.table)
options(ggrepel.max.overlaps = Inf)

# Analysis of CyTOF immune cell frequencies with all 4 batches

# Batches:
# EXP20CL8796
# EXP20DG3622
# EXP21DG3628
# EXP22DG3665

# Load flowframe
setwd("Please insert working directory here to point at directory '../data/Figure2/'")

ff = read.FCS('220811_CyTOF_1234batches_ComBat_corrected.fcs', transformation = FALSE, truncate_max_range = FALSE)
combat_edata_B_label <- read.csv("220811_CyTOF_1234batches_ComBat_corrected_labels.csv", row.names = 1)

# 1. FlowSOM 30 clusters ----------
common_marker = colnames(ff@exprs) %>% unname()

# Run FlowSOM
fSOM_30 = FlowSOM(ff,
                  compensate = F,
                  transform = F,
                  scale = F,
                  colsToUse = common_marker,
                  nClus = 10,
                  xdim = 5, ydim = 6,
                  seed = 824)

# Median marker expression of FlowSOM cell clusters
MFI_30 = GetClusterMFIs(fSOM_30, colsUsed = T) #gets MFI for each one of the 30 clusters
clusters_30 = GetClusters(fSOM_30) # vector of cluster numbers
freqClusters_30 = data.frame(clusters = clusters_30) %>% #calculates freq of these 30 clusters
  dplyr::count(.data$clusters) %>%
  dplyr::mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  as.data.frame()
res_30 = cbind(MFI_30, freqClusters_30) #merges MFI and frequency of 30 clusters

# 2. Separate Neutrophil clusters ----------
# Identify Neutrophil clusters according to heatmap and pdf output (where metacluster 7 is neutrophils)
nrow_neutrophils = which(clusters_30 %in% c(12,13,16,17,18,19,21,23,26,27,28,29,30))

# Separate Neutrophil data
label_neutrophils = combat_edata_B_label[nrow_neutrophils,] #labels of just neutr clusters
clusters_neutrophils = clusters_30[nrow_neutrophils] #vector of clusters for neutrophils
label_neutrophils$lineage = 'Neutrophils'
label_neutrophils$subtype = 'Neutrophils'
label_neutrophils$clusters = clusters_neutrophils
data_neutrophils = ff@exprs[nrow_neutrophils,] #get data for neutrohils only
res_neutrophils = res_30 %>% dplyr::filter(clusters %in% c(12,13,16,17,18,19,21,23,26,27,28,29,30)) #results just for neutrophils (mfi and %)
res_neutrophils$lineage = 'Neutrophils'
res_neutrophils$subtype = 'Neutrophils'

# Remaining cells after removing Neutrophils
data_remaining = ff@exprs[-nrow_neutrophils,] #data from all except for neutrophils
label_remaining = combat_edata_B_label[-nrow_neutrophils,] # label from all except for neutrophils

# 3. FlowSOM 100 clusters on remaining cells ----------

ff_remaining = ff
ff_remaining@exprs = data_remaining

# Run FlowSOM in the remaining cells
fSOM_100 = FlowSOM(ff_remaining,
                   compensate = F,
                   transform = F,
                   scale = F,
                   colsToUse = common_marker,
                   nClus = 10,
                   xdim = 10, ydim = 10,
                   seed = 824)

# Median marker expression of FlowSOM cell clusters
MFI_100 = GetClusterMFIs(fSOM_100, colsUsed = T)
clusters_100 = GetClusters(fSOM_100)
freqClusters_100 = data.frame(clusters = clusters_100) %>%
  dplyr::count(.data$clusters) %>%
  dplyr::mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  as.data.frame()
res_100 = cbind(MFI_100, freqClusters_100)

# Heatmap of FlowSOM clusters
q99 = quantile(as.matrix(MFI_100), 0.99)
myheatmap <- pheatmap::pheatmap(as.matrix(MFI_100),
                                scale = 'none',
                                labels_row = paste(rownames(MFI_100),' (', round(freqClusters_30$percentage,1), '%', ')',sep = ''),
                                display_numbers = TRUE,
                                angle_col = 45,
                                breaks = seq(0, q99, q99/90),
                                main = "100 clusters: Median marker expression per cluster")


# Manual annotation
mypopulations <- data.frame(
  clusters=seq(1:100),
  metaclusters=fSOM_100$metaclustering)

mypopulations <- cbind(mypopulations, MFI_100)

mypopulations <- mypopulations %>% mutate(major_clusters=case_when(
  metaclusters %in% 1 ~ "Tcells", #CD3, CD4 CD8 CD5 etc; many subclusters
  metaclusters %in% 2 ~ "Mono_DCs", #CD11c, CD14
  metaclusters %in% 3 ~ "NKs", #CD56 
  metaclusters %in% 4 ~ "Tcells", #CD3, CD4, CD8
  metaclusters %in% 5 ~ "dead", #positive for many things; only 1 cluster for this metcluster
  metaclusters %in% 6 ~ "Basophils", #CD123+ HLADR-; only 1 subcluster
  metaclusters %in% 7 ~ "CD34_CD9",
  metaclusters %in% 8 ~ "Activated", #HLADR
  metaclusters %in% 9 ~ "Eos_Tcells",
  metaclusters %in% 10 ~ "Bcells")) #HLADR

# LEVEL 2

mypopulations <- mypopulations %>% mutate(flowSOM_level2=case_when(
  # within metacluster 1 of Tcells
  clusters %in% 1 ~ "CD4Tcells_TCM_CD127", 
  clusters %in% 2 ~ "CD4Tcells_TCM", 
  clusters %in% 3 ~ "CD4Tcells_TCM", 
  clusters %in% 4 ~ "CD4Tcells_TCM", 
  clusters %in% 5 ~ "CD4Tcells_TCM", 
  clusters %in% 11 ~ "CD4Tcells_TCM", 
  clusters %in% 12 ~ "CD4Tcells_TCM", 
  clusters %in% 13 ~ "CD4Tcells_Treg", 
  clusters %in% 14 ~ "CD4Tcells_TCM", 
  clusters %in% 15 ~ "CD4Tcells_EM", 
  clusters %in% 21 ~ "CD8Tcells_MAIT", 
  clusters %in% 22 ~ "gdTcells", 
  clusters %in% 23 ~ "CD4Tcells_TCM", 
  clusters %in% 24 ~ "CD8Tcells_TCM", 
  clusters %in% 25 ~ "CD4Tcells_TCM", 
  clusters %in% 31 ~ "CD8Tcells_Naive", 
  clusters %in% 32 ~ "CD8Tcells_EM", 
  clusters %in% 33 ~ "gdTcells", 
  clusters %in% 34 ~ "CD8Tcells_MAIT", 
  clusters %in% 35 ~ "unknown_35_CD99hi", 
  clusters %in% 41 ~ "CD8Tcells_EM", 
  clusters %in% 42 ~ "CD8Tcells_TCM", 
  clusters %in% 45 ~ "CD4Tcells_Naive", 
  clusters %in% 46 ~ "CD4Tcells_TCM", 
  clusters %in% 51 ~ "CD8Tcells_TEMRA", 
  clusters %in% 61 ~ "CD8Tcells_TEMRA", 
  clusters %in% 62 ~ "CD8Tcells_TEMRA", 
  clusters %in% 63 ~ "CD8Tcells_MAIT", 
  clusters %in% 71 ~ "CD4Tcells_Naive", 
  clusters %in% 72 ~ "CD4Tcells_Naive", 
  clusters %in% 81 ~ "CD4Tcells_Naive", 
  clusters %in% 82 ~ "CD4Tcells_Naive", 
  clusters %in% 83 ~ "CD4Tcells_Naive", 
  clusters %in% 85 ~ "CD8Tcells_Naive", 
  clusters %in% 86 ~ "CD8Tcells_Naive", 
  clusters %in% 91 ~ "CD4Tcells_Naive", 
  clusters %in% 92 ~ "CD4Tcells_Naive", 
  clusters %in% 93 ~ "CD4Tcells_Naive", 
  clusters %in% 94 ~ "CD8Tcells_Naive", 
  clusters %in% 95 ~ "CD8Tcells_Naive", 
  clusters %in% 96 ~ "CD8Tcells_Naive", 
  
  # within metacluster 2 of Mono_DCs
  clusters %in% 6 ~ "Monocytes_CM", 
  clusters %in% 7 ~ "Monocytes_CM", 
  clusters %in% 8 ~ "Monocytes_CM", 
  clusters %in% 9 ~ "Monocytes_CM", 
  clusters %in% 10 ~ "Monocytes_CM",
  clusters %in% 16 ~ "Monocytes_CM",
  clusters %in% 17 ~ "Monocytes_CM",
  clusters %in% 18 ~ "DC_CD141DC",
  clusters %in% 19 ~ "Monocytes_CM",
  clusters %in% 20 ~ "Monocytes_CM",
  clusters %in% 26 ~ "DC_CD1cDC",
  clusters %in% 27 ~ "Monocytes_CM",
  clusters %in% 28 ~ "Monocytes_CM",
  clusters %in% 29 ~ "Monocytes_CM",
  clusters %in% 30 ~ "Monocytes_NCM",
  clusters %in% 36 ~ "MDSC_M",
  clusters %in% 37 ~ "Monocytes_CM",
  clusters %in% 38 ~ "Monocytes_CM",
  clusters %in% 39 ~ "DC_CD141DC",
  clusters %in% 40 ~ "DC_CD141DC",
  clusters %in% 47 ~ "DPTcell",
  clusters %in% 48 ~ "DC_CD141DC",
  clusters %in% 49 ~ "DC_CD141DC",
  clusters %in% 56 ~ "CD8Tcells_TEMRA_CD57",
  clusters %in% 57 ~ "CD8Tcells_TEMRA_CD57",
  clusters %in% 58 ~ "CD4Tcells_TCM",
  clusters %in% 65 ~ "Neutrophils_HLADR",
  clusters %in% 66 ~ "DPTcell",
  clusters %in% 67 ~ "DC_preDC",
  clusters %in% 75 ~ "CD4Tcells_TEMRA",
  clusters %in% 76 ~ "CD8Tcells_TEMRA",
  
  # within metacluster 3 of NKs
  clusters %in% 43 ~ "NKcells_CD56bright", 
  clusters %in% 44 ~ "NKcells_CD56dim",
  clusters %in% 52 ~ "NKcells_CD56dim",
  clusters %in% 53 ~ "NKcells_CD56bright",
  clusters %in% 54 ~ "NKcells_CD56dim",
  clusters %in% 55 ~ "NKcells_CD56dim",
  clusters %in% 64 ~ "pDC",
  
  # within metacluster 4 of Tcells
  clusters %in% 50 ~ "Monocytes_IM", #CD14+HLADR+CD16+
  clusters %in% 60 ~ "CD8Tcells_TCM_CD24",
  clusters %in% 70 ~ "CD8Tcells_TEMRA_CD24",
  clusters %in% 74 ~ "Monocytes_NCM",
  clusters %in% 77 ~ "CD8Tcells_TEMRA_CD24",
  clusters %in% 84 ~ "CD8Tcells_Naive_CD24",
  clusters %in% 99 ~ "CD8Tcells_Naive",
  
  # within metacluster 5 of dead
  clusters %in% 59 ~ "dead",
  
  # within metacluster 6 of Basophils
  clusters %in% 68 ~ "Basophils", 
  
  # within metacluster 7 of CD34_CD9; CD20, IgD, CD34, CD9
  clusters %in% 69 ~ "Bcells_progenitor", 
  clusters %in% 78 ~ "Bcells_progenitor", 
  clusters %in% 89 ~ "Bcells_progenitor", 
  
  # within metacluster 8 of Activated
  clusters %in% 73 ~ "DC_CD141DC", 
  
  # within metacluster 9 of Eos_Tcells
  clusters %in% 79 ~ "Eosinophils", 
  clusters %in% 80 ~ "Eosinophils",
  clusters %in% 90 ~ "Eosinophils",
  clusters %in% 100 ~ "Eosinophils",
  
  # within metacluster 10 of Bcells
  clusters %in% 87 ~ "Bcells_memory", #CD20 IgD-
  clusters %in% 88 ~ "Bcells_naive", #CD20 IgD+
  clusters %in% 97 ~ "Bcells_naive", #CD20 IgD+
  clusters %in% 98 ~ "Bcells_memoryCD27" #CD20 IgD- CD27+
))


# CHECK and adjust annotations
colnames(mypopulations)
q99 = quantile(as.matrix(mypopulations[,3:50]), 0.99)
myplotc <- pheatmap::pheatmap(as.matrix(mypopulations[,3:50]),
                              scale = 'none',
                              labels_row = paste(mypopulations$level2,' (', round(freqClusters_30$percentage,1), '%', ')',sep = ''),
                              display_numbers = TRUE,
                              angle_col = 45,
                              breaks = seq(0, q99, q99/90),
                              main = "100 clusters: Median marker expression per cluster")

# update a proper level1
mypopulations <- mypopulations %>% mutate(flowSOM_level1=sapply(str_split(mypopulations$flowSOM_level2, "_"), "[", 1)) 

## save annotation ##
mypopulations <- mypopulations[,c(1:2,52:53,3:50)]
write.csv(mypopulations, "221003_flowsom_annotations.csv")

# levels 1 and 2
res_100$flowSOM_level1 <- mypopulations$flowSOM_level1[match(res_100$clusters, mypopulations$clusters)]
res_100$flowSOM_level2 <- mypopulations$flowSOM_level2[match(res_100$clusters, mypopulations$clusters)]

label_remaining$clusters = clusters_100
label_remaining$flowSOM_level1 <- mypopulations$flowSOM_level1[match(label_remaining$clusters, mypopulations$clusters)]
label_remaining$flowSOM_level2 <- mypopulations$flowSOM_level2[match(label_remaining$clusters, mypopulations$clusters)]

# 4. Bind remaining cells with neutrophils ----------
data_bind = rbind(data_neutrophils, data_remaining)
label_neutrophils <- label_neutrophils[,c(1:8,11)]
label_neutrophils$flowSOM_level1 <- "Neutrophils"
label_neutrophils$flowSOM_level2 <- "Neutrophils"
identical(colnames(label_neutrophils), colnames(label_remaining))
label_bind = rbind(label_neutrophils, label_remaining)

all_bind = cbind(data_bind, label_bind) 

res_neutrophils <- res_neutrophils[,c(1:51)]
res_neutrophils$flowSOM_level1 <- "Neutrophils"
res_neutrophils$flowSOM_level2 <- "Neutrophils"
res_bind = rbind(res_neutrophils, res_100)
res_bind %<>% mutate(frequency = (n/sum(res_bind$n))*100)

# Heatmap of all cell clusters
q99 = quantile(as.matrix(res_bind[,1:length(common_marker)]), 0.99)
myplotd <- pheatmap::pheatmap(as.matrix(res_bind[,1:length(common_marker)]),
                              scale = 'none',
                              labels_row = paste(res_bind$flowSOM_level2,' (', round(res_bind$percentage,1), '%', ')',sep = ''),
                              display_numbers = TRUE,
                              angle_col = 45,
                              breaks = seq(0, q99, q99/90),
                              main = "Median marker expression per cluster")

# Calculate relative frequency of each celltype in each sample

# flowSOM level2
freq_lineage = label_bind %>%
  group_by(FCSfile) %>%
  count(.data$flowSOM_level2) %>%
  mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  dplyr::select(-n) %>%
  spread(flowSOM_level2, percentage, fill=0)

# flowSOM level1
freq_lev1 = label_bind %>%
  group_by(FCSfile) %>%
  count(.data$flowSOM_level1) %>%
  mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  dplyr::select(-n) %>%
  spread(flowSOM_level1, percentage, fill=0)

# 5. Save results ----------
write.csv(freq_lineage, file = '221003_flowsom_frequency_level2_populations.csv', row.names = F)
write.csv(freq_lev1, file = '221003_flowsom_frequency_level1_populations.csv', row.names = F)
write.csv(all_bind, file = '221003_CyTOF_1234batches_ComBat_corrected_FlowSOM.csv', row.names = F)

# Save median marker expression data for Network plot -- these have FtM and MtF combined
write.table(res_bind, file = '221003_flowsom_clustered.txt', sep="\t", row.names=F, col.names = T, quote = F)

# Save sampleinfo
sampleinfo <- combat_edata_B_label[,c(2,3,4,5)] %>% unique()
rownames(sampleinfo) <- NULL
write.csv(sampleinfo, "221003_sampleinfo.csv")

### Select df with GATT only for further analysis ###
setwd("/Users/camilaconsiglio/Documents/BrodinLab/SRT/Data/Metadata/")
meta <- read.csv("220829_Metadata.csv",row.names = 1)
FtM <- meta$Subject[meta$Transition %in% "FtM"]

combat_edata_B_label$Subject <- paste(combat_edata_B_label$SubjectID, combat_edata_B_label$Visit, sep = "_")
label_FtM <- combat_edata_B_label$label[combat_edata_B_label$Subject %in% FtM]
combat_edata_B_label$rownames <- as.numeric(rownames(combat_edata_B_label))
which_FtM <- combat_edata_B_label$rownames[combat_edata_B_label$Subject %in% FtM]

# filter FtM from flowsom30
fSOM_30_FtM <- FlowSOMSubset(fSOM_30, which_FtM)
MFI_30_FtM = GetClusterMFIs(fSOM_30_FtM, colsUsed = T) #gets MFI for each one of the 30 clusters
clusters_30FtM = GetClusters(fSOM_30_FtM) # vector of cluster numbers
freqClusters_30FtM = data.frame(clusters = clusters_30FtM) %>% #calculates freq of these 30 clusters
  dplyr::count(.data$clusters) %>%
  dplyr::mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  as.data.frame()
res_30_FtM = cbind(MFI_30_FtM, freqClusters_30FtM) #merges MFI and frequency of 30 clusters

nrow_neutrophils = which(clusters_30 %in% c(12,13,16,17,18,19,21,23,26,27,28,29,30))

# Separate Neutrophil 
res_neutrophils_FtM = res_30_FtM %>% dplyr::filter(clusters %in% c(12,13,16,17,18,19,21,23,26,27,28,29,30)) #results just for neutrophils (mfi and %)

# filter FtM from flowsom100
rownames(label_remaining) <- NULL
label_remaining$rownames <- as.numeric(rownames(label_remaining))
unique(label_remaining$Subject)
label_remaining$Subject <- paste(label_remaining$SubjectID, label_remaining$Visit, sep = "_")
colnames(label_remaining)
which_FtM_100c <- label_remaining %>% filter(Subject %in% FtM) %>% pull(rownames)

fSOM_100_FtM <- FlowSOMSubset(fsom = fSOM_100, ids = which_FtM_100c)

MFI_100_FtM = GetClusterMFIs(fSOM_100_FtM, colsUsed = T) #gets MFI for each one of the 100 clusters
clusters_100FtM = GetClusters(fSOM_100_FtM) # vector of cluster numbers
freqClusters_100FtM = data.frame(clusters = clusters_100FtM) %>% #calculates freq of these 100 clusters
  dplyr::count(.data$clusters) %>%
  dplyr::mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  as.data.frame()
res_100_FtM = cbind(MFI_100_FtM, freqClusters_100FtM) #merges MFI and frequency of 100 clusters

# merge res FtM
colnames(res_neutrophils_FtM)
res_neutrophils_FtM$flowSOM_level1 <- "Neutrophils"
res_neutrophils_FtM$flowSOM_level2 <- "Neutrophils"
colnames(res_100_FtM)
res_100_FtM$flowSOM_level1 <- mypopulations$flowSOM_level1[match(res_100_FtM$clusters, mypopulations$clusters)]
res_100_FtM$flowSOM_level2 <- mypopulations$flowSOM_level2[match(res_100_FtM$clusters, mypopulations$clusters)]
res_bind_FtM = rbind(res_neutrophils_FtM, res_100_FtM)
res_bind_FtM %<>% mutate(frequency = (n/sum(res_bind_FtM$n))*100)

# relabel debris
res_bind$flowSOM_level2[res_bind$flowSOM_level2 %in% "Monocytes_NCM"] <- "debris"

# relabel of neutrophil cluster numbers so they dont overlab with nonneutrophil cluster numbers
length(res_bind$flowSOM_level2[res_bind$flowSOM_level2 == "Neutrophils"])
length(101:113)
res_bind$clusters[res_bind$flowSOM_level2 == "Neutrophils"] <- 101:113

# save sampleinfo_FtM
sampleinfo_FtM <- sampleinfo[paste(sampleinfo$SubjectID, sampleinfo$Visit, sep = "_") %in% FtM,]
write.csv(sampleinfo_FtM, "221003_sampleinfo_FtM.csv")

# save frequency of levels 1 and 2 FtM
freq_lineage_FtM <- freq_lineage %>% filter(FCSfile %in% sampleinfo_FtM$FCSfile)
freq_lev1_FtM <- freq_lev1 %>% filter(FCSfile %in% sampleinfo_FtM$FCSfile)
write.csv(freq_lineage_FtM, '221003_flowsom_frequency_level2_populations_FtM.csv')
write.csv(freq_lev1_FtM, '221003_flowsom_frequency_level1_populations_FtM.csv')

# save all_bind_FtM
all_bind_FtM <- all_bind %>% filter(FCSfile %in% sampleinfo_FtM$FCSfile)
dim(all_bind_FtM) #12377068 x 59
write.csv(all_bind_FtM, '221003_CyTOF_1234batches_ComBat_corrected_FlowSOM_FtM.csv')

# save
write.table(res_bind, file = '221013_flowsom_clustered_FtM_updatedclusternumbers.txt', sep="\t", row.names=F, col.names = T, quote = F)







