# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## FIGURE 4B ##
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

counts <- Read10X_h5(filename = "4B/1_aggr_out/filtered_peak_bc_matrix.h5")
metadata <- read.csv(file = "4B/1_aggr_out/singlecell.csv", header = T, row.names = 1)
metadata$sampleID <- sapply(str_split(rownames(metadata), "-"), "[", 2)
# fragment file: full list of unique fragments across all single cells, and the fragment index file (fragments.tsv.gz.tbi) needs to be in the same folder for this to work
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38', 
  fragments = "4B/1_aggr_out/fragments.tsv.gz",
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
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] scales_1.2.1                      cowplot_1.1.1                    
# [3] ggthemes_4.2.4                    TFBSTools_1.36.0                 
# [5] JASPAR2020_0.99.10                BSgenome.Hsapiens.UCSC.hg38_1.4.5
# [7] BSgenome_1.66.3                   rtracklayer_1.58.0               
# [9] Biostrings_2.66.0                 XVector_0.38.0                   
# [11] paletteer_1.5.0                   ggpubr_0.6.0                     
# [13] ggrepel_0.9.3                     lubridate_1.9.2                  
# [15] forcats_1.0.0                     stringr_1.5.0                    
# [17] dplyr_1.1.1                       purrr_1.0.1                      
# [19] readr_2.1.4                       tidyr_1.3.0                      
# [21] tibble_3.2.1                      tidyverse_2.0.0                  
# [23] patchwork_1.1.2                   ggplot2_3.4.2                    
# [25] EnsDb.Hsapiens.v86_2.99.0         ensembldb_2.22.0                 
# [27] AnnotationFilter_1.22.0           GenomicFeatures_1.50.4           
# [29] AnnotationDbi_1.60.2              Biobase_2.58.0                   
# [31] GenomicRanges_1.50.2              GenomeInfoDb_1.34.9              
# [33] IRanges_2.32.0                    S4Vectors_0.36.2                 
# [35] BiocGenerics_0.44.0               SeuratObject_4.1.3               
# [37] Seurat_4.3.0                      Signac_1.9.0                     
# 
# loaded via a namespace (and not attached):
# [1] rappdirs_0.3.3              scattermore_0.8             R.methodsS3_1.8.2          
# [4] bit64_4.0.5                 knitr_1.42                  irlba_2.3.5.1              
# [7] DelayedArray_0.24.0         R.utils_2.12.2              data.table_1.14.8          
# [10] rpart_4.1.19                flowCore_2.10.0             KEGGREST_1.38.0            
# [13] RCurl_1.98-1.12             generics_0.1.3              RSQLite_2.3.1              
# [16] RANN_2.6.1                  future_1.32.0               bit_4.0.5                  
# [19] tzdb_0.4.0                  spatstat.data_3.0-1         xml2_1.3.3                 
# [22] httpuv_1.6.9                SummarizedExperiment_1.28.0 DirichletMultinomial_1.40.0
# [25] xfun_0.38                   hms_1.1.3                   evaluate_0.20              
# [28] promises_1.2.0.1            fansi_1.0.4                 restfulr_0.0.15            
# [31] progress_1.2.2              caTools_1.18.2              dbplyr_2.3.2               
# [34] igraph_1.4.1                DBI_1.1.3                   htmlwidgets_1.6.2          
# [37] spatstat.geom_3.1-0         ellipsis_0.3.2              backports_1.4.1            
# [40] cytolib_2.10.1              prismatic_1.1.1             annotate_1.76.0            
# [43] biomaRt_2.54.1              deldir_1.0-6                MatrixGenerics_1.10.0      
# [46] vctrs_0.6.1                 ROCR_1.0-11                 abind_1.4-5                
# [49] cachem_1.0.7                withr_2.5.0                 progressr_0.13.0           
# [52] checkmate_2.1.0             sctransform_0.3.5           GenomicAlignments_1.34.1   
# [55] prettyunits_1.1.1           goftest_1.2-3               cluster_2.1.4              
# [58] lazyeval_0.2.2              seqLogo_1.64.0              crayon_1.5.2               
# [61] hdf5r_1.3.8                 spatstat.explore_3.1-0      labeling_0.4.2             
# [64] pkgconfig_2.0.3             nlme_3.1-162                ProtGenerics_1.30.0        
# [67] nnet_7.3-18                 rlang_1.1.0                 globals_0.16.2             
# [70] lifecycle_1.0.3             miniUI_0.1.1.1              filelock_1.0.2             
# [73] BiocFileCache_2.6.1         dichromat_2.0-0.1           polyclip_1.10-4            
# [76] matrixStats_0.63.0          lmtest_0.9-40               Matrix_1.5-4               
# [79] carData_3.0-5               zoo_1.8-11                  base64enc_0.1-3            
# [82] pheatmap_1.0.12             ggridges_0.5.4              png_0.1-8                  
# [85] viridisLite_0.4.1           rjson_0.2.21                bitops_1.0-7               
# [88] R.oo_1.25.0                 KernSmooth_2.23-20          blob_1.2.4                 
# [91] parallelly_1.35.0           spatstat.random_3.1-4       rstatix_0.7.2              
# [94] ggsignif_0.6.4              CNEr_1.34.0                 memoise_2.0.1              
# [97] magrittr_2.0.3              plyr_1.8.8                  ica_1.0-3                  
# [100] zlibbioc_1.44.0             compiler_4.2.1              BiocIO_1.8.0               
# [103] RColorBrewer_1.1-3          fitdistrplus_1.1-8          Rsamtools_2.14.0           
# [106] cli_3.6.1                   listenv_0.9.0               pbapply_1.7-0              
# [109] htmlTable_2.4.1             Formula_1.2-5               MASS_7.3-58.3              
# [112] tidyselect_1.2.0            RProtoBufLib_2.10.0         stringi_1.7.12             
# [115] yaml_2.3.7                  grid_4.2.1                  VariantAnnotation_1.44.1   
# [118] fastmatch_1.1-3             tools_4.2.1                 timechange_0.2.0           
# [121] future.apply_1.10.0         parallel_4.2.1              rstudioapi_0.14            
# [124] TFMPvalue_0.0.9             foreign_0.8-85              gridExtra_2.3              
# [127] farver_2.1.1                Rtsne_0.16                  digest_0.6.31              
# [130] shiny_1.7.4                 pracma_2.4.2                motifmatchr_1.20.0         
# [133] Rcpp_1.0.10                 car_3.1-2                   broom_1.0.4                
# [136] later_1.3.0                 RcppAnnoy_0.0.20            httr_1.4.5                 
# [139] biovizBase_1.46.0           colorspace_2.1-0            XML_3.99-0.14              
# [142] tensor_1.5                  reticulate_1.28             splines_4.2.1              
# [145] uwot_0.1.14                 RcppRoll_0.3.0              rematch2_2.1.2             
# [148] spatstat.utils_3.0-2        sp_1.6-0                    plotly_4.10.1              
# [151] xtable_1.8-4                jsonlite_1.8.4              poweRlaw_0.70.6            
# [154] R6_2.5.1                    Hmisc_5.0-1                 pillar_1.9.0               
# [157] htmltools_0.5.5             mime_0.12                   glue_1.6.2                 
# [160] fastmap_1.1.1               BiocParallel_1.32.6         codetools_0.2-19           
# [163] utf8_1.2.3                  lattice_0.21-8              spatstat.sparse_3.0-1      
# [166] curl_5.0.0                  leiden_0.4.3                gtools_3.9.4               
# [169] GO.db_3.16.0                survival_3.5-5              rmarkdown_2.21             
# [172] munsell_0.5.0               GenomeInfoDbData_1.2.9      reshape2_1.4.4             
# [175] gtable_0.3.3 

