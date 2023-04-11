# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## FIGURE 4I ##
install.packages("Seurat", version=4, force=TRUE)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86) 
library(ggplot2)
library(patchwork)
set.seed(1234)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(paletteer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)
options(ggrepel.max.overlaps = Inf)
library(ggthemes)
library(cowplot)
library(scales)
library(stringr)

## scATACseq analysis ##

## Preprocessing
# Load annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'

# Load data from cellranger aggr output
setwd("Please insert working directory here to point at directory '../data/Figure4/'")

counts <- Read10X_h5(filename = "4i/1_aggr_out/filtered_peak_bc_matrix.h5")
metadata <- read.csv(file = "4i/1_aggr_out/singlecell.csv", header = T, row.names = 1)
metadata$sampleID <- sapply(str_split(rownames(metadata), "-"), "[", 2)
# fragment file: full list of unique fragments across all single cells, and the fragment index file (fragments.tsv.gz.tbi) needs to be in the same folder for this to work
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38', 
  fragments = "4i/1_aggr_out/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200)

# Create seurat object using peak/cell matrix;
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata)
Annotation(pbmc) <- annotations 

### QC ###
# QC: TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = TRUE)  
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')

# QC: nucleosome signal
pbmc <- NucleosomeSignal(object = pbmc)
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

# QC: Fraction of fragments in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

# QC: Blacklist ratio
pbmc$blacklist_fraction <- FractionCountsInRegion(
  object = pbmc, 
  assay = 'peaks',
  regions = blacklist_hg38)

pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_fraction < 0.05 & 
    nucleosome_signal < 4 &
    TSS.enrichment > 2)

#pbmc #143624 features across 12773 samples within 1 assay 

### Normalization and dimensionality reduction ###
pbmc <- RunTFIDF(pbmc)
# Feature selection
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
# Dimension reduction
pbmc <- RunSVD(pbmc)
# The first LSI component captures technical variation rather than biological, so it will be removed from downstream analysis. 
#DepthCor(pbmc)
# Nonlinear dim red
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

DimPlot(object = pbmc, label = TRUE) + NoLegend()

pbmc$sampleID <- sapply(str_split(colnames(pbmc), "-"), "[", 2)
pbmc$Subject <- case_when( 
  pbmc$sampleID %in% c("1", "2", "3") ~ "STO10",
  pbmc$sampleID %in% c("4", "5") ~ "STO14",
  pbmc$sampleID %in% c("6", "7", "8") ~ "STO11")

pbmc$VisitMonths <- case_when( 
  pbmc$sampleID %in% c("1", "4", "6") ~ 0,
  pbmc$sampleID %in% c("2", "5", "7") ~ 3,
  pbmc$sampleID %in% c("3", "8") ~ 12)

### Determine gene activities
gene.activities <- GeneActivity(pbmc)
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA))

rm(gene.activities)

DefaultAssay(pbmc) <- 'RNA'

# annotate
pbmc$new_clusters <- case_when(
  pbmc$seurat_clusters %in% c(0, 1, 5, 8, 17) ~ "CD4_Tcells",
  pbmc$seurat_clusters %in% c(2, 9) ~ "Bcells",
  pbmc$seurat_clusters %in% c(3, 11) ~ "CD8_Tcells",
  pbmc$seurat_clusters %in% c(4, 7, 13, 15) ~ "CD8_Tcells",
  pbmc$seurat_clusters %in% 10 ~ "CD8_Tcells",
  pbmc$seurat_clusters %in% c(6, 12, 16) ~ "NKcells",
  pbmc$seurat_clusters %in% c(14, 19, 20) ~ "monocytes_DCs",
  pbmc$seurat_clusters %in% 18 ~ "Bcells")

Idents(pbmc) <- pbmc$new_clusters
Idents(pbmc) <- factor(Idents(pbmc), levels = c(
  "Bcells", "CD4_Tcells", "CD8_Tcells", 
  "NKcells", "monocytes_DCs"))

pbmc$new_clusters_visit <- paste(pbmc$new_clusters, pbmc$VisitMonths, sep="_")
Idents(pbmc) <- pbmc$new_clusters
DimPlot(pbmc, label = TRUE, cols = paletteer::paletteer_d("ggthemes::excel_Vapor_Trail")[2:6]) + NoLegend()

## TF activity analysis
DefaultAssay(pbmc) <- 'peaks'
Idents(pbmc) <- pbmc$new_clusters_visit

#Finding overrepresented motifs in a set of differentially accessible peaks, 

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE, species = 9606) ) #species = 9606 is for human
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- as.logical(seqnames(granges(pbmc)) %in% main.chroms)
human_pbmc <- pbmc[keep.peaks, ]

# add motif information, now stored in a new object called human_pbmc
human_pbmc <- AddMotifs(object = human_pbmc, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
Idents(human_pbmc) <- human_pbmc$new_clusters_visit

# Function to find overrepresented motifs in a set of differentially accessible peaks
meta.feature <- GetAssayData(human_pbmc, assay = "peaks", slot = "meta.features")
DefaultAssay(human_pbmc) <- 'peaks'

motifs_analysis <- function(new_clusters_visitPost, new_clusters_visitPre) {
  #Goal: to find overrepresented motifs in a set of differentially accessible peaks
  # step1: differentially accessible peaks in a cluster between pre and post SRT
  da_peaks <- FindMarkers(
    object = human_pbmc,
    ident.1 = new_clusters_visitPost,
    ident.2 = new_clusters_visitPre,
    only.pos = TRUE, 
    test.use = 'LR',
    min.pct = 0.05,
    latent.vars = 'nCount_peaks')
  top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ]) #top differentially accessible peaks
  
  # step2: determine all open peaks in this population (background)
  open.peaks <- AccessiblePeaks(human_pbmc, idents = c(new_clusters_visitPost, new_clusters_visitPre))
  # match the overall GC content in the peak set
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[top.da.peak, ],
    n = 50000)
  
  # step3: find which motifs are overrepresented in Post vs Pre using the differntially accessible peaks as features, and the open peaks in this cell population as background
  enriched.motifs <- FindMotifs(
    object = subset(human_pbmc, idents= c(new_clusters_visitPost, new_clusters_visitPre)),
    features = top.da.peak,
    background=peaks.matched,
    p.adjust.method = "BH")
  # annotate which cell and visit comparison
  enriched.motifs$comparison <- paste(new_clusters_visitPost, "vs", new_clusters_visitPre, sep = "__" )
  # return motifs df
  return(enriched.motifs)
}

# create dfs
CD4Tcells_3x0_ma <- motifs_analysis(new_clusters_visitPost = "CD4_Tcells_3", new_clusters_visitPre = "CD4_Tcells_0")
CD4Tcells_12x0_ma <- motifs_analysis(new_clusters_visitPost = "CD4_Tcells_12", new_clusters_visitPre = "CD4_Tcells_0")
Bcells_3x0_ma <- motifs_analysis(new_clusters_visitPost = "Bcells_3", new_clusters_visitPre = "Bcells_0")
Bcells_12x0_ma <- motifs_analysis(new_clusters_visitPost = "Bcells_12", new_clusters_visitPre = "Bcells_0")
NKcells_3x0_ma <- motifs_analysis(new_clusters_visitPost = "NKcells_3", new_clusters_visitPre = "NKcells_0")
NKcells_12x0_ma <- motifs_analysis(new_clusters_visitPost = "NKcells_12", new_clusters_visitPre = "NKcells_0")
CD8Tcells_3x0_ma <- motifs_analysis(new_clusters_visitPost = "CD8_Tcells_3", new_clusters_visitPre = "CD8_Tcells_0")
CD8Tcells_12x0_ma <- motifs_analysis(new_clusters_visitPost = "CD8_Tcells_12", new_clusters_visitPre = "CD8_Tcells_0")
Myeloidcells_3x0_ma <- motifs_analysis(new_clusters_visitPost = "monocytes_DCs_3", new_clusters_visitPre = "monocytes_DCs_0")
Myeloidcells_12x0_ma <- motifs_analysis(new_clusters_visitPost = "monocytes_DCs_12", new_clusters_visitPre = "monocytes_DCs_0")

ma_all <- rbind(CD4Tcells_3x0_ma, CD4Tcells_12x0_ma,
                Bcells_3x0_ma, Bcells_12x0_ma,
                NKcells_3x0_ma, NKcells_12x0_ma,
                CD8Tcells_3x0_ma, CD8Tcells_12x0_ma,
                Myeloidcells_3x0_ma, Myeloidcells_12x0_ma)

ma_all$significance <- ifelse(ma_all$p.adjust < 0.05, "Sig", "ns")
ma_all$Cell <- str_remove(sapply(str_split(ma_all$comparison, "\\d__"), "[", 1), "1")
ma_all$visit_comparison <- ifelse(str_detect(ma_all$comparison, "12"),  ma_all$visit_comparison <- "12m_vs_0m", ma_all$visit_comparison <- "3m_vs_0m")
ma_all$Cell_comparison <- paste(ma_all$Cell, ma_all$visit_comparison, sep = "_")
ma_all$visit_comparison <- factor(ma_all$visit_comparison, levels = c("3m_vs_0m","12m_vs_0m"))
ma_all$log_fold.enrichment <- log(ma_all$fold.enrichment)

# TFs of interest
NFKB <- unique(c(unique(ma_all$motif.name)[str_detect(unique(ma_all$motif.name), "REL")],
                 unique(ma_all$motif.name)[str_detect(unique(ma_all$motif.name), "NFK")]))

STATs <- c("STAT1","STAT1::STAT2","STAT3")

more <- c("CREB1", "JUN", "FOS")

ma_all %>% 
  filter(motif.name %in% c(NFKB, STATs, more)) %>% 
  mutate(motif.name=factor(motif.name, levels = c(NFKB, STATs, more))) %>%
  filter(visit_comparison %in% "12m_vs_0m") %>%
  mutate(stars=cut(p.adjust, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " "))) %>%
  ggplot(aes(x=Cell, y=motif.name)) + 
  geom_point(aes(fill = log_fold.enrichment), shape = 21, alpha = 0.8, size=10) +
  scale_fill_gradientn(colors = paletteer::paletteer_c("grDevices::Blue-Red", n=20, direction = -1), oob=scales::oob_squish_any, limits=c(-1,1)) +
  theme_bw() + geom_text(aes(label=stars), color="black", size=5, nudge_x = 0, nudge_y = -0.05) + ggtitle("12m vs baseline") +
  labs(x=NULL, y=NULL, size=NULL, fill="log(fold enrichment)")  + theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))



