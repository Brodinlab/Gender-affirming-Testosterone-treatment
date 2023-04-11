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

## DESeq2 ##

# load run info
meta <- read.csv("Fig1_metadata.csv", row.names = 1) 

# kallisto output files
dir <- "./"
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
  scale_fill_gradientn(colours = paletteer::paletteer_c("grDevices::Blue-Red 2", n=20, direction = -1), limits=c(-3.5,3.5)) +
  geom_node_text(data = 
                   filter(pother$data, name %in% c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE")), 
                 aes(label = name), size = 5, repel = T) +
  geom_node_text(data = filter(pother$data, !name %in% c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE")),
                 aes(label = name), size = 3, repel = T) +
  theme_graph(base_family = 'Helvetica') + labs(title = "FtM 12m", size="abs(NES)", fill="NES")

