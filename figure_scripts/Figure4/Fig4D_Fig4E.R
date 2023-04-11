# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## FIGURES 4D & 4E ##
library(Seurat)
library(stringr)
library(tidyverse)
library(ggrepel)
library(ggthemes)
library(cowplot)
library(qusage)
options(ggrepel.max.overlaps = Inf)

# scRNAseq 10X, first run/Camila
setwd("Please insert working directory here to point at directory '../data/Figure4/'")
# load data
# GE
ge.data <- Read10X('4D_4E/2203_STO11/GE/')
# HTO
hto.data <- Read10X('4D_4E/2203_STO11/HTO/raw/')

# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
joint.bcs <- intersect(colnames(ge.data), colnames(hto.data$`Antibody Capture`))

# Subset RNA and HTO counts by joint cell barcodes
ex.umis <- ge.data[, joint.bcs]
ex.htos <- as.matrix(hto.data$`Antibody Capture`[, joint.bcs])
# Confirm that the HTO have the correct names
#rownames(ex.htos)

# Setup Seurat object
ex.hashtag <- CreateSeuratObject(counts = ex.umis)

# Add HTO data as a new assay independent from RNA
ex.hashtag[["HTO"]] <- CreateAssayObject(counts = ex.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
ex.hashtag <- NormalizeData(ex.hashtag, assay = "HTO", normalization.method ="CLR")

# QC and selecting cells for further analysis
ex.hashtag[["percent.mt"]] <- PercentageFeatureSet(ex.hashtag, pattern = "^MT-")
ex.hashtag <- subset(ex.hashtag, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10) #8895 samples
ex.hashtag <- HTODemux(ex.hashtag, assay = "HTO", positive.quantile = 0.99) # , verbose=T
Idents(ex.hashtag) <- "HTO_classification.global"

# Filtering
ex.hashtag.subset <- subset(ex.hashtag, idents ="Negative", invert =TRUE)
ex.singlets <- subset(ex.hashtag, idents = "Singlet")
ex.singlets$hash.ID <- droplevels(ex.singlets$hash.ID)
ex.singlets$Subject <- case_when(
  ex.singlets$hash.ID %in% "TS-B-0251" ~ "STO11",
  ex.singlets$hash.ID %in% "TS-B-0252" ~ "STO11",
  ex.singlets$hash.ID %in% "TS-B-0253" ~ "STO11",
  ex.singlets$hash.ID %in% "TS-B-0254" ~ "STO11",
  ex.singlets$hash.ID %in% "TS-B-0255" ~ "STO11")

ex.singlets$VisitMonths <- case_when(
  ex.singlets$hash.ID %in% "TS-B-0251" ~ 0,
  ex.singlets$hash.ID %in% "TS-B-0252" ~ 0,
  ex.singlets$hash.ID %in% "TS-B-0253" ~ 0,
  ex.singlets$hash.ID %in% "TS-B-0254" ~ 3,
  ex.singlets$hash.ID %in% "TS-B-0255" ~ 3)

ex.singlets$Stim <- case_when(
  ex.singlets$hash.ID %in% "TS-B-0251" ~ "NTC",
  ex.singlets$hash.ID %in% "TS-B-0252" ~ "LPS",
  ex.singlets$hash.ID %in% "TS-B-0253" ~ "R848",
  ex.singlets$hash.ID %in% "TS-B-0254" ~ "NTC",
  ex.singlets$hash.ID %in% "TS-B-0255" ~ "LPS")

ex.singlets$Stim_Visit <- paste(ex.singlets$Stim, ex.singlets$VisitMonths, sep = "_")

ex.singlets <- SCTransform(ex.singlets, method = "glmGamPoi", verbose = FALSE) 
ex.singlets <- RunPCA(ex.singlets, features =VariableFeatures(ex.singlets)) #npcs = 30, or features =VariableFeatures(ex.singlets), or  features = myfeats

ElbowPlot(ex.singlets)
mydims <- 1:15

# UMAP and Clustering
ex.singlets <- RunUMAP(ex.singlets, reduction = "pca", dims = mydims, verbose=FALSE)
ex.singlets <- FindNeighbors(ex.singlets, reduction = "pca", dims = mydims)
ex.singlets <- FindClusters(ex.singlets, resolution = 1.2) #0.8 is default - 0.9 separates some NKs from CD8Ts
DimPlot(ex.singlets, label = TRUE) 
#DimPlot(ex.singlets, group.by ="HTO_classification") + facet_wrap(~HTO_classification)
#DimPlot(ex.singlets, group.by ="VisitMonths") + facet_wrap(~VisitMonths)
#DimPlot(ex.singlets, group.by ="Stim") + facet_wrap(~Stim)

# cluster annotation
pbmc.markers <- FindAllMarkers(ex.singlets, assay = 'SCT', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3)

ex.singlets$new_clusters <- case_when(
  ex.singlets$seurat_clusters %in% c(0,1,2,3,4,6,10,11,16,18,21) ~ "CD4_Tcells",
  ex.singlets$seurat_clusters %in% c(7,13,17,19) ~ "NKcells",
  ex.singlets$seurat_clusters %in% c(5,12,14) ~ "Bcells",
  ex.singlets$seurat_clusters %in% c(8,15) ~ "CD8_Tcells",
  ex.singlets$seurat_clusters %in% c(20,22) ~ "Monocytes",
  ex.singlets$seurat_clusters %in% c(9) ~ "CCR7_Tcells")

DimPlot(ex.singlets, group.by = "new_clusters", label = TRUE)
#DimPlot(ex.singlets, group.by = 'new_clusters', split.by = "Stim_Visit")

Idents(ex.singlets) <- ex.singlets$new_clusters

# B_cell
B_cell <- subset(ex.singlets, idents = 'Bcells') 
DefaultAssay(B_cell) <- "RNA"
B_cell <- ScaleData(B_cell)
B_cell <- SCTransform(B_cell, method = "glmGamPoi", verbose = FALSE) 
B_cell <- RunPCA(B_cell, features = VariableFeatures(B_cell))
dims<-1:15
B_cell <- FindNeighbors(B_cell, dims = dims)
B_cell <- FindClusters(B_cell, resolution = 0.3)
B_cell <- RunUMAP(B_cell, dims = dims, verbose=FALSE)
Idents(B_cell) <- B_cell$Stim_Visit
B_cell_NTC_V2V1 <- B_cell %>% FindMarkers(ident.1 = "NTC_3", ident.2 = "NTC_0", verbose = FALSE, logfc.threshold = 0, min.pct = 0) %>% mutate(Cell="B_cell", Stim="NTC", Comparison="V2vsV1") %>% rownames_to_column()
B_cell_LPS_V2V1 <- B_cell %>% FindMarkers(ident.1 = "LPS_3", ident.2 = "LPS_0", verbose = FALSE, logfc.threshold = 0, min.pct = 0) %>% mutate(Cell="B_cell", Stim="LPS", Comparison="V2vsV1") %>% rownames_to_column()
B_cell_summary <- rbind(B_cell_NTC_V2V1, B_cell_LPS_V2V1)
B_cell_summary$stars = cut(B_cell_summary$p_val_adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " "))

# CD4_Tcells
CD4_Tcells <- subset(ex.singlets, idents = 'CD4_Tcells') 
DefaultAssay(CD4_Tcells) <- "RNA"
CD4_Tcells <- ScaleData(CD4_Tcells)
CD4_Tcells <- SCTransform(CD4_Tcells, method = "glmGamPoi", verbose = FALSE) 
CD4_Tcells <- RunPCA(CD4_Tcells, features = VariableFeatures(CD4_Tcells))
dims<-1:15
CD4_Tcells <- FindNeighbors(CD4_Tcells, dims = dims)
CD4_Tcells <- FindClusters(CD4_Tcells, resolution = 0.3)
CD4_Tcells <- RunUMAP(CD4_Tcells, dims = dims, verbose=FALSE)
Idents(CD4_Tcells) <- CD4_Tcells$Stim_Visit
CD4_Tcells_NTC_V2V1 <- CD4_Tcells %>% FindMarkers(ident.1 = "NTC_3", ident.2 = "NTC_0", verbose = FALSE, logfc.threshold = 0, min.pct = 0) %>% mutate(Cell="CD4_Tcells", Stim="NTC", Comparison="V2vsV1") %>% rownames_to_column()
CD4_Tcells_LPS_V2V1 <- CD4_Tcells %>% FindMarkers(ident.1 = "LPS_3", ident.2 = "LPS_0", verbose = FALSE, logfc.threshold = 0, min.pct = 0) %>% mutate(Cell="CD4_Tcells", Stim="LPS", Comparison="V2vsV1") %>% rownames_to_column()
CD4_Tcells_summary <- rbind(CD4_Tcells_NTC_V2V1, CD4_Tcells_LPS_V2V1)
CD4_Tcells_summary$stars = cut(CD4_Tcells_summary$p_val_adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " "))

# NKcells
NKcells <- subset(ex.singlets, idents = 'NKcells') 
DefaultAssay(NKcells) <- "RNA"
NKcells <- ScaleData(NKcells)
NKcells <- SCTransform(NKcells, method = "glmGamPoi", verbose = FALSE) 
NKcells <- RunPCA(NKcells, features = VariableFeatures(NKcells))
dims<-1:15
NKcells <- FindNeighbors(NKcells, dims = dims)
NKcells <- FindClusters(NKcells, resolution = 0.3)
NKcells <- RunUMAP(NKcells, dims = dims, verbose=FALSE)
Idents(NKcells) <- NKcells$Stim_Visit
NKcells_NTC_V2V1 <- NKcells %>% FindMarkers(ident.1 = "NTC_3", ident.2 = "NTC_0", verbose = FALSE, logfc.threshold = 0, min.pct = 0) %>% mutate(Cell="NKcells", Stim="NTC", Comparison="V2vsV1") %>% rownames_to_column()
NKcells_LPS_V2V1 <- NKcells %>% FindMarkers(ident.1 = "LPS_3", ident.2 = "LPS_0", verbose = FALSE, logfc.threshold = 0, min.pct = 0) %>% mutate(Cell="NKcells", Stim="LPS", Comparison="V2vsV1") %>% rownames_to_column()
NKcells_summary <- rbind(NKcells_NTC_V2V1, NKcells_LPS_V2V1)
NKcells_summary$stars = cut(NKcells_summary$p_val_adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " "))

# Monocytes
Monocytes <- subset(ex.singlets, idents = 'Monocytes') 
DefaultAssay(Monocytes) <- "RNA"
Monocytes <- ScaleData(Monocytes)
Monocytes <- SCTransform(Monocytes, method = "glmGamPoi", verbose = FALSE) 
Monocytes <- RunPCA(Monocytes, features = VariableFeatures(Monocytes))
dims<-1:15
Monocytes <- FindNeighbors(Monocytes, dims = dims)
Monocytes <- FindClusters(Monocytes, resolution = 0.3)
Monocytes <- RunUMAP(Monocytes, dims = dims, verbose=FALSE)
Idents(Monocytes) <- Monocytes$Stim_Visit
Monocytes_NTC_V2V1 <- Monocytes %>% FindMarkers(ident.1 = "NTC_3", ident.2 = "NTC_0", verbose = FALSE, logfc.threshold = 0, min.pct = 0) %>% mutate(Cell="Monocytes", Stim="NTC", Comparison="V2vsV1") %>% rownames_to_column()
Monocytes_LPS_V2V1 <- Monocytes %>% FindMarkers(ident.1 = "LPS_3", ident.2 = "LPS_0", verbose = FALSE, logfc.threshold = 0, min.pct = 0) %>% mutate(Cell="Monocytes", Stim="LPS", Comparison="V2vsV1") %>% rownames_to_column()
Monocytes_summary <- rbind(Monocytes_NTC_V2V1, Monocytes_LPS_V2V1) #no Monocytes_NTC_V2V1
Monocytes_summary <- Monocytes_LPS_V2V1
Monocytes_summary$stars = cut(Monocytes_summary$p_val_adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " "))

# CD8_Tcells
CD8_Tcells <- subset(ex.singlets, idents = 'CD8_Tcells') 
DefaultAssay(CD8_Tcells) <- "RNA"
CD8_Tcells <- ScaleData(CD8_Tcells)
CD8_Tcells <- SCTransform(CD8_Tcells, method = "glmGamPoi", verbose = FALSE) 
CD8_Tcells <- RunPCA(CD8_Tcells, features = VariableFeatures(CD8_Tcells))
dims<-1:15
CD8_Tcells <- FindNeighbors(CD8_Tcells, dims = dims)
CD8_Tcells <- FindClusters(CD8_Tcells, resolution = 0.3)
CD8_Tcells <- RunUMAP(CD8_Tcells, dims = dims, verbose=FALSE)
Idents(CD8_Tcells) <- CD8_Tcells$Stim_Visit
CD8_Tcells_NTC_V2V1 <- CD8_Tcells %>% FindMarkers(ident.1 = "NTC_3", ident.2 = "NTC_0", verbose = FALSE, logfc.threshold = 0, min.pct = 0) %>% mutate(Cell="CD8_Tcells", Stim="NTC", Comparison="V2vsV1") %>% rownames_to_column()
CD8_Tcells_LPS_V2V1 <- CD8_Tcells %>% FindMarkers(ident.1 = "LPS_3", ident.2 = "LPS_0", verbose = FALSE, logfc.threshold = 0, min.pct = 0) %>% mutate(Cell="CD8_Tcells", Stim="LPS", Comparison="V2vsV1") %>% rownames_to_column()
CD8_Tcells_summary <- rbind(CD8_Tcells_NTC_V2V1, CD8_Tcells_LPS_V2V1)
CD8_Tcells_summary$stars = cut(CD8_Tcells_summary$p_val_adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " "))

summary_all <- rbind(B_cell_summary, Monocytes_summary,
                     CD4_Tcells_summary, CD8_Tcells_summary,
                     NKcells_summary)

summary_all$Stim <- factor(summary_all$Stim, levels = c("NTC", "LPS"))

## AR response elements ##
TFBS <- read.gmt("c3.tft.v2022.1.Hs.symbols.gmt")
AR <- unique(c(TFBS$AR_01, TFBS$AR_02, TFBS$AR_03, TFBS$AR_Q2, TFBS$AR_Q6))

# filter AR responsive genes in CD4T and CD8T
summary_all <- summary_all %>% 
  filter(rowname %in% AR) %>% filter(Stim %in% "NTC") %>% 
  filter(Cell %in% c("CD4_Tcells","CD8_Tcells"))

summary_all %>%
  ggplot(aes(x=avg_log2FC, y=rowname)) + geom_point(aes(alpha=p_val_adj< 0.05)) + 
  facet_wrap(~Cell, ncol=2) + geom_vline(xintercept = 0, linetype="dashed") +
  theme(axis.text.y = element_blank()) + 
  labs(y=NULL, alpha="significance", color=NULL, x="average Log2FC (3m/baseline)") + 
  geom_text_repel(data=subset(summary_all, p_val_adj< 0.05), aes(label=rowname), size=2) +
  scale_color_manual(values = "black") + scale_alpha_manual(values=c(0.1,1)) + 
  ggtitle("scRNAseq 10X expression of predicted AR targets (unstim)") +
  scale_x_continuous(limits = c(-2,2))





