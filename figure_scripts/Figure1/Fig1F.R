# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## FIGURES 1D & 1E ##
library(tidyverse)
library(stringr)
library(GenomicFeatures)
library(tximport)
library(DESeq2)
library(limma)
library(biomaRt)
library(fgsea)
library(cowplot)
library(ggrepel)
options(ggrepel.max.overlaps = Inf)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggthemes)
library(msigdbr)
library(ggraph)
library(plyr)

setwd("Please insert working directory here to point at directory './data/Figure1/'")
## DESeq2 ##

# load run info
meta <- read.csv("Fig1_metadata.csv", row.names = 1) 

# NOTE: Please go to folder and uncompress the files in "Figure1/Figure1F-folder".
# kallisto output files
dir <- "./Figure1/Figure1F/"
run <- list.files() 
files <- file.path(dir, run, "abundance.tsv") 
#59 files

#files
#file.exists(files)
#all(file.exists(files))
#length(run) #59
names(files) <- substr(run, 1, 11)

# load GTF file, can be downloaded from http://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/
gtffile <- file.path(getwd(), 'Homo_sapiens.GRCh38.90.gtf')
txdb <- makeTxDbFromGFF(gtffile, format = 'gtf')
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# prep meta for design
meta$Visitb <- meta$Visit
meta$Visitb <- case_when(
  meta$Visitb == 1 ~ "V1",
  meta$Visitb == 2 ~ "V2",
  meta$Visitb == 3 ~ "V3")
meta$Visitb <- factor(meta$Visitb)

meta <- meta %>% filter(!is.na(Age)) #59

# tximport kallisto output files
txi <- tximport(files, type='kallisto', tx2gene=tx2gene, ignoreTxVersion = TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ SubjectID + Age + Visitb) #

# filtering
rs_treat_int <- rowSums(counts(dds))
minrs = 100
#dim(dds) #pre-filtering:  34865    59
dds <- dds[ rs_treat_int >= minrs, ]
#dim(dds) #post-filtering: 26622    59

dds <- estimateSizeFactors(dds)

nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= ncol(dds)/4 #keep normalized count of at least 10 in at least one fourth of samples
dds <- dds[filter,]
#dim(dds) #post-filtering: 18053    59

#Run DESeq
dds <- DESeq(dds) 

# Results
resV2 <- results(dds, contrast = list("Visitb_V2_vs_V1"))
resV3 <- results(dds, contrast = list("Visitb_V3_vs_V1"))

#External gene IDs
ensembl_ids <- rownames(resV2)
mart_r <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl") 
genemap <- getBM(mart = mart_r,values = ensembl_ids,filter = c("ensembl_gene_id"),attributes = c("ensembl_gene_id", "external_gene_name") )
idx_rV2 <- match(rownames(resV2), genemap$ensembl_gene_id)
resV2$gene <- genemap$external_gene_name[idx_rV2]

idx_rV3 <- match(rownames(resV3), genemap$ensembl_gene_id)
resV3$gene <- genemap$external_gene_name[idx_rV3]

resV2$Comparison <- "Visit_V2_vs_V1"
resV3$Comparison <- "Visit_V3_vs_V1"

resV2 <- resV2 %>% as.data.frame() %>% filter(!is.na(gene))  %>% filter(!gene %in% "") 
resV3 <- resV3 %>% as.data.frame() %>% filter(!is.na(gene))  %>% filter(!gene %in% "") 
resV2 <- resV2 %>% filter(!duplicated(gene))
resV3 <- resV3 %>% filter(!duplicated(gene))

## GSEA with Hallmarks ##
m_df <- msigdbr(species = "Homo sapiens")
C3_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)

# V3: cnetplot using clusterprofiler
resV3
Hids <- bitr(rownames(resV3), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
dedup_Hids = Hids[!duplicated(Hids[c("ENSEMBL")]),] # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
Hstat.test = resV3[rownames(resV3) %in% dedup_Hids$ENSEMBL,] # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
Hstat.test$Y = dedup_Hids$ENTREZID # Create a new column in df2 with the corresponding ENTREZ IDs
Hstat.test$ENSEMBL <- rownames(Hstat.test)
rownames(Hstat.test) <- NULL
H3kegg_gene_list <- Hstat.test$log2FoldChange # Create a vector of the gene unuiverse
names(H3kegg_gene_list) <- Hstat.test$Y # Name vector with ENTREZ ids
H3kegg_gene_list <- na.omit(H3kegg_gene_list) # omit any NA values 
H3kegg_gene_list = sort(H3kegg_gene_list, decreasing = TRUE) # sort the list in decreasing order (required for clusterProfiler)
set.seed(2244)
em3 <- GSEA(H3kegg_gene_list, TERM2GENE = C3_t2g, pvalueCutoff = 1)
em3_g <- setReadable(em3, 'org.Hs.eg.db', 'ENTREZID')

mycats <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE")

pother <- cnetplot(em3_g, showCategory = mycats,  color_category="black", foldChange = H3kegg_gene_list) 

ggraph(pother$data) +
  geom_edge_link(alpha=.4) + 
  geom_node_point(data = pother$data, aes(size= abs(color), fill = color), shape=21,  alpha = 0.9) +
  scale_size_continuous(range = c(4, 14)) +
  scale_fill_gradientn(colours = c("#EE5A45", "white", "#1F8F89"), values = scales::rescale(c(-3.5, 0, 3.5)), limits=c(-3.5, 3.5)) +
  geom_node_text(data = 
                   filter(pother$data, name %in% c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE")), 
                 aes(label = name), size = 5, repel = T, max.overlaps = Inf) +
  geom_node_text(data = filter(pother$data, !name %in% c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE")),
                 aes(label = name), size = 3, repel = T, max.overlaps = Inf) +
  theme_graph(base_family = 'Helvetica') + labs(title = "FtM 12m", size="abs(NES)", fill="NES")
setwd("~/Library/CloudStorage/OneDrive-LundUniversity/Systems_Immunology_Lab/Projects/GAHT_MHT/Nature2024/Figures")
ggsave(filename = 'Nature2024_Fig1F.pdf', width = 15, height = 10)


sessionInfo()
#R version 4.2.1 (2022-06-23)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Monterey 12.5

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

#locale:
#[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] plyr_1.8.8                  ggraph_2.1.0                msigdbr_7.5.1              
#[4] ggthemes_4.2.4              org.Hs.eg.db_3.16.0         clusterProfiler_4.6.2      
#[7] ggrepel_0.9.3               cowplot_1.1.1               fgsea_1.24.0               
#[10] biomaRt_2.54.1              limma_3.54.2                DESeq2_1.38.3              
#[13] SummarizedExperiment_1.28.0 MatrixGenerics_1.10.0       matrixStats_0.63.0         
#[16] tximport_1.26.1             GenomicFeatures_1.50.4      AnnotationDbi_1.60.2       
#[19] Biobase_2.58.0              GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
#[22] IRanges_2.32.0              S4Vectors_0.36.2            BiocGenerics_0.44.0        
#[25] lubridate_1.9.2             forcats_1.0.0               stringr_1.5.0              
#[28] dplyr_1.1.1                 purrr_1.0.1                 readr_2.1.4                
#[31] tidyr_1.3.0                 tibble_3.2.1                ggplot2_3.4.2              
#[34] tidyverse_2.0.0            

#loaded via a namespace (and not attached):
# [1] utf8_1.2.3               OlinkAnalyze_3.6.2       tidyselect_1.2.0        
#[4] lme4_1.1-32              RSQLite_2.3.1            htmlwidgets_1.6.2       
#[7] grid_4.2.1               BiocParallel_1.32.6      scatterpie_0.1.8        
#[10] devtools_2.4.5           munsell_0.5.0            codetools_0.2-19        
#[13] miniUI_0.1.1.1           withr_2.5.0              colorspace_2.1-0        
#[16] GOSemSim_2.24.0          filelock_1.0.2           rstudioapi_0.14         
#[19] DOSE_3.24.2              labeling_0.4.2           emmeans_1.8.5           
#[22] GenomeInfoDbData_1.2.9   polyclip_1.10-4          farver_2.1.1            
#[25] bit64_4.0.5              datawizard_0.7.1         downloader_0.4          
#[28] treeio_1.22.0            vctrs_0.6.1              generics_0.1.3          
#[31] gson_0.1.0               timechange_0.2.0         BiocFileCache_2.6.1     
#[34] randomForest_4.7-1.1     R6_2.5.1                 graphlayouts_0.8.4      
#[37] locfit_1.5-9.7           gridGraphics_0.5-1       bitops_1.0-7            
#[40] cachem_1.0.7             DelayedArray_0.24.0      vroom_1.6.1             
#[43] promises_1.2.0.1         BiocIO_1.8.0             scales_1.2.1            
#[46] enrichplot_1.18.4        gtable_0.3.3             processx_3.8.0          
#[49] tidygraph_1.2.3          rlang_1.1.0              splines_4.2.1           
#[52] lazyeval_0.2.2           rtracklayer_1.58.0       rstatix_0.7.2           
#[55] prismatic_1.1.1          broom_1.0.4              BiocManager_1.30.20     
#[58] yaml_2.3.7               reshape2_1.4.4           abind_1.4-5             
#[61] backports_1.4.1          httpuv_1.6.9             qvalue_2.30.0           
#[64] tree_1.0-43              tools_4.2.1              usethis_2.1.6           
#[67] ggplotify_0.1.0          ellipsis_0.3.2           RColorBrewer_1.1-3      
#[70] sessioninfo_1.2.2        Rcpp_1.0.10              progress_1.2.2          
#[73] zlibbioc_1.44.0          RCurl_1.98-1.12          ps_1.7.4                
#[76] prettyunits_1.1.1        viridis_0.6.2            urlchecker_1.0.1        
#[79] fs_1.6.1                 magrittr_2.0.3           data.table_1.14.8       
#[82] lmerTest_3.1-3           mvtnorm_1.1-3            ggnewscale_0.4.8        
#[85] pkgload_1.3.2            patchwork_1.1.2          hms_1.1.3               
#[88] mime_0.12                xtable_1.8-4             HDO.db_0.99.1           
#[91] XML_3.99-0.14            readxl_1.4.2             gridExtra_2.3           
#[94] compiler_4.2.1           shadowtext_0.1.2         crayon_1.5.2            
#[97] minqa_1.2.5              htmltools_0.5.5          ggfun_0.0.9             
#[100] later_1.3.0              tzdb_0.4.0               aplot_0.1.10            
#[103] geneplotter_1.76.0       DBI_1.1.3                tweenr_2.0.2            
#[106] dbplyr_2.3.2             MASS_7.3-58.3            rappdirs_0.3.3          
#[109] boot_1.3-28.1            babelgene_22.9           Matrix_1.5-4            
#[112] car_3.1-2                cli_3.6.1                parallel_4.2.1          
#[115] insight_0.19.1           igraph_1.4.1             pkgconfig_2.0.3         
#[118] GenomicAlignments_1.34.1 numDeriv_2016.8-1.1      xml2_1.3.3              
#[121] paletteer_1.5.0          ggtree_3.6.2             annotate_1.76.0         
#[124] XVector_0.38.0           estimability_1.4.1       yulab.utils_0.0.6       
#[127] callr_3.7.3              digest_0.6.31            parameters_0.20.3       
#[130] Biostrings_2.66.0        cellranger_1.1.0         fastmatch_1.1-3         
#[133] tidytree_0.4.2           restfulr_0.0.15          curl_5.0.0              
#[136] shiny_1.7.4              Rsamtools_2.14.0         rjson_0.2.21            
#[139] nloptr_2.0.3             jsonlite_1.8.4           lifecycle_1.0.3         
#[142] nlme_3.1-162             reprtree_0.6             carData_3.0-5           
#[145] viridisLite_0.4.1        fansi_1.0.4              pillar_1.9.0            
#[148] lattice_0.21-8           KEGGREST_1.38.0          fastmap_1.1.1           
#[151] httr_1.4.5               pkgbuild_1.4.0           GO.db_3.16.0            
#[154] glue_1.6.2               remotes_2.4.2            bayestestR_0.13.1       
#[157] zip_2.2.2                png_0.1-8                bit_4.0.5               
#[160] ggforce_0.4.1            stringi_1.7.12           profvis_0.3.7           
#[163] rematch2_2.1.2           blob_1.2.4               memoise_2.0.1           
#[166] ape_5.7-1  

