# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## EXT FIG 1B &  EXT FIG 1C ##
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

# Design including all cells for EXT FIG 1B
setwd("Please insert working directory here to point at directory '../data/Figure2/' - also please note that you first need to run Fig2A_CyTOF_1_flowSOM.R to generate files needed in this script")

# Download all data from Gender-affirming-Testosterone-treatment/data/Figure1
# load run info
meta <- read.csv("220829_Metadata_FtM.csv", row.names = 1) #59 samples

# load major population info
# Load frequency data generated with flowSOM 
imm.freq <-  read.csv('221003_flowsom_frequency_level1_populations_FtM.csv', row.names = 1)
sampleinfo <- read.csv("221003_sampleinfo.csv")
m_imm.freq <- merge(sampleinfo, imm.freq, by="FCSfile")
m_imm.freq$Subject <- paste(m_imm.freq$SubjectID, m_imm.freq$Visit, sep = "_") #60 samples

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

# load GTF file
gtffile <- file.path(getwd(), 'Homo_sapiens.GRCh38.90.gtf')
txdb <- makeTxDbFromGFF(gtffile, format = 'gtf')
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# tximport kallisto output files
txi <- tximport(files, type='kallisto', tx2gene=tx2gene, ignoreTxVersion = TRUE)

# merge meta with freq cell type df
meta <- merge(meta, m_imm.freq, by="Subject") #59 samples

# prep meta for design
meta$Visitb <- meta$Visit.x
meta$Visitb <- case_when(
  meta$Visitb == 1 ~ "V1",
  meta$Visitb == 2 ~ "V2",
  meta$Visitb == 3 ~ "V3")
meta$Visitb <- factor(meta$Visitb)

meta <- meta %>% filter(!is.na(Age)) #59
meta[,c(12,22:36)] <- sapply(meta[,c(12,22:36)], scale)

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ SubjectID.x + Age + Bcells + CD4Tcells + CD8Tcells + DC + Eosinophils + Monocytes + Neutrophils + NKcells + pDC + Visitb) #

# filtering steps
rs_treat_int <- rowSums(counts(dds))
minrs = 100
dds <- dds[ rs_treat_int >= minrs, ]

dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= ncol(dds)/4 #we want a normalized count of at least 10 in at least one fourth of samples
dds <- dds[filter,]

#Run DESeq
dds <- DESeq(dds) 


# Results
resV2 <- results(dds, contrast = list("Visitb_V2_vs_V1"))
resV3 <- results(dds, contrast = list("Visitb_V3_vs_V1"))
resAge <- results(dds, contrast = list("Age"))
resBcells <- results(dds, contrast = list("Bcells"))
resCD4Tcells <- results(dds, contrast = list("CD4Tcells"))
resCD8Tcells <- results(dds, contrast = list("CD8Tcells"))
resDC <- results(dds, contrast = list("DC"))
resEosinophils <- results(dds, contrast = list("Eosinophils"))
resMonocytes <- results(dds, contrast = list("Monocytes"))
resNeutrophils <- results(dds, contrast = list("Neutrophils"))
resNKcells <- results(dds, contrast = list("NKcells"))
respDC<- results(dds, contrast = list("pDC"))

#External gene IDs
ensembl_ids <- rownames(resV2)
mart_r <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl") 
genemap <- getBM(mart = mart_r,values = ensembl_ids,filter = c("ensembl_gene_id"),attributes = c("ensembl_gene_id", "external_gene_name") )

add_extIDs <- function(myresults, mycomparison) {
  idx <- match(rownames(myresults), genemap$ensembl_gene_id)
  myresults$gene <- genemap$external_gene_name[idx]
  myresults$Comparison <- paste(mycomparison)
  myresults <- myresults  %>%
    as.data.frame() 
}

resV2 <- add_extIDs(myresults = resV2, mycomparison = "V2_vs_V1")
resV3 <- add_extIDs(myresults = resV3, mycomparison = "V3_vs_V1")
resAge <- add_extIDs(myresults = resAge, mycomparison = "Age")
resBcells <- add_extIDs(myresults = resBcells, mycomparison = "Bcells")
resCD4Tcells <- add_extIDs(myresults = resCD4Tcells, mycomparison = "CD4Tcells")
resCD8Tcells <- add_extIDs(myresults = resCD8Tcells, mycomparison = "CD8Tcells")
resDC <- add_extIDs(myresults = resDC, mycomparison = "DC")
resEosinophils <- add_extIDs(myresults = resEosinophils, mycomparison = "Eosinophils")
resMonocytes <- add_extIDs(myresults = resMonocytes, mycomparison = "Monocytes")
resNeutrophils <- add_extIDs(myresults = resNeutrophils, mycomparison = "Neutrophils")
resNKcells <- add_extIDs(myresults = resNKcells, mycomparison = "NKcells")
respDC <- add_extIDs(myresults = respDC, mycomparison = "pDC")

res_all <- rbind(
  resV2, resV3, resAge,
  resBcells, resCD4Tcells, resCD8Tcells,
  resDC, resEosinophils, resMonocytes,
  resNeutrophils, resNKcells, respDC
)

### GSEA with Hallmarks 

m_df <- msigdbr(species = "Homo sapiens")
C3_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%  # H for hallmarks
  dplyr::select(gs_name, entrez_gene)
head(C3_t2g)

run_gsea <- function(myresults, mycomparison) {
  Hids <- bitr(rownames(myresults), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
  dedup_Hids = Hids[!duplicated(Hids[c("ENSEMBL")]),] 
  Hstat.test = myresults[rownames(myresults) %in% dedup_Hids$ENSEMBL,] 
  Hstat.test$Y = dedup_Hids$ENTREZID 
  Hstat.test$ENSEMBL <- rownames(Hstat.test)
  rownames(Hstat.test) <- NULL
  H3kegg_gene_list <- Hstat.test$log2FoldChange 
  names(H3kegg_gene_list) <- Hstat.test$Y 
  H3kegg_gene_list <- na.omit(H3kegg_gene_list) 
  H3kegg_gene_list = sort(H3kegg_gene_list, decreasing = TRUE)
  set.seed(2244)
  em3 <- GSEA(H3kegg_gene_list, TERM2GENE = C3_t2g, pvalueCutoff = 1)
  em3_g <- setReadable(em3, 'org.Hs.eg.db', 'ENTREZID')
  gsea_res <- em3_g@result[1:11] 
  gsea_res$Comparison <- mycomparison
  gsea_res
}

resV3_gsea <- run_gsea(myresults = resV3, mycomparison = "V3vsV1")
resV2_gsea <- run_gsea(myresults = resV2, mycomparison = "V2vsV1")
resAge_gsea <- run_gsea(myresults = resAge, mycomparison = "Age")
resBcells_gsea <- run_gsea(myresults = resBcells, mycomparison = "Bcells")
resCD4Tcells_gsea <- run_gsea(myresults = resCD4Tcells, mycomparison = "CD4Tcells")
resCD8Tcells_gsea <- run_gsea(myresults = resCD8Tcells, mycomparison = "CD8Tcells")
resDC_gsea <- run_gsea(myresults = resDC, mycomparison = "DC")
resEosinophils_gsea <- run_gsea(myresults = resEosinophils, mycomparison = "Eosinophils")
resMonocytes_gsea <- run_gsea(myresults = resMonocytes, mycomparison = "Monocytes")
resNeutrophils_gsea <- run_gsea(myresults = resNeutrophils, mycomparison = "Neutrophils")
resNKcells_gsea <- run_gsea(myresults = resNKcells, mycomparison = "NKcells")
respDC_gsea <- run_gsea(myresults = respDC, mycomparison = "pDC")

res_gsea_all <- rbind(
  resV2_gsea, resV3_gsea, resAge_gsea,
  resBcells_gsea, resCD4Tcells_gsea, resCD8Tcells_gsea,
  resDC_gsea, resEosinophils_gsea, resMonocytes_gsea,
  resNeutrophils_gsea, resNKcells_gsea, respDC_gsea
)

ids_interest <- c("HALLMARK_INFLAMMATORY_RESPONSE",
                  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                  "HALLMARK_FATTY_ACID_METABOLISM",
                  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                  "HALLMARK_GLYCOLYSIS",
                  "HALLMARK_CHOLESTEROL_HOMEOSTASIS")

# EXT FIG 1C
P1 <- res_gsea_all %>%
  filter(ID %in% ids_interest) %>%
  mutate(stars=cut(p.adjust, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " "))) %>%
  mutate(Comparison=factor(Comparison, levels = c(unique(res_gsea_all$Comparison)))) %>%
  mutate(ID=factor(ID, levels = rev(ids_interest) )) %>%
  ggplot(aes(x=Comparison, y=ID, fill=NES)) +
  geom_point(shape=21, size=6) +
  geom_text(aes(label=stars), color="black", size=5) +
  scale_fill_gradientn(colours = paletteer::paletteer_c("grDevices::Blue-Red 2", n=20, direction = -1), limits=c(-3.5,3.5)) +
  labs(x=NULL, y=NULL) +  theme_classic() + ggtitle("full model pDC")
P1


## Design including all cells execpt pDC for EXT FIG 1C ##
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ SubjectID.x + Age + Bcells + CD4Tcells + CD8Tcells + DC + Eosinophils + Monocytes + Neutrophils + NKcells + Visitb) #

# filtering steps
rs_treat_int <- rowSums(counts(dds))
minrs = 100
dds <- dds[ rs_treat_int >= minrs, ]
dds <- estimateSizeFactors(dds)

# more filtering
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= ncol(dds)/4 #we want a normalized count of at least 10 in at least one fourth of samples
dds <- dds[filter,]

#Run DESeq
dds <- DESeq(dds) 

# Results
resV2 <- results(dds, contrast = list("Visitb_V2_vs_V1"))
resV3 <- results(dds, contrast = list("Visitb_V3_vs_V1"))
resAge <- results(dds, contrast = list("Age"))
resBcells <- results(dds, contrast = list("Bcells"))
resCD4Tcells <- results(dds, contrast = list("CD4Tcells"))
resCD8Tcells <- results(dds, contrast = list("CD8Tcells"))
resDC <- results(dds, contrast = list("DC"))
resEosinophils <- results(dds, contrast = list("Eosinophils"))
resMonocytes <- results(dds, contrast = list("Monocytes"))
resNeutrophils <- results(dds, contrast = list("Neutrophils"))
resNKcells <- results(dds, contrast = list("NKcells"))

#External gene IDs
ensembl_ids <- rownames(resV2)
mart_r <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl") 
genemap <- getBM(mart = mart_r,values = ensembl_ids,filter = c("ensembl_gene_id"),attributes = c("ensembl_gene_id", "external_gene_name") )

add_extIDs <- function(myresults, mycomparison) {
  idx <- match(rownames(myresults), genemap$ensembl_gene_id)
  myresults$gene <- genemap$external_gene_name[idx]
  myresults$Comparison <- paste(mycomparison)
  myresults <- myresults  %>%
    as.data.frame()
}

resV2 <- add_extIDs(myresults = resV2, mycomparison = "V2_vs_V1")
resV3 <- add_extIDs(myresults = resV3, mycomparison = "V3_vs_V1")
resAge <- add_extIDs(myresults = resAge, mycomparison = "Age")
resBcells <- add_extIDs(myresults = resBcells, mycomparison = "Bcells")
resCD4Tcells <- add_extIDs(myresults = resCD4Tcells, mycomparison = "CD4Tcells")
resCD8Tcells <- add_extIDs(myresults = resCD8Tcells, mycomparison = "CD8Tcells")
resDC <- add_extIDs(myresults = resDC, mycomparison = "DC")
resEosinophils <- add_extIDs(myresults = resEosinophils, mycomparison = "Eosinophils")
resMonocytes <- add_extIDs(myresults = resMonocytes, mycomparison = "Monocytes")
resNeutrophils <- add_extIDs(myresults = resNeutrophils, mycomparison = "Neutrophils")
resNKcells <- add_extIDs(myresults = resNKcells, mycomparison = "NKcells")

res_all <- rbind(
  resV2, resV3, resAge,
  resBcells, resCD4Tcells, resCD8Tcells,
  resDC, resEosinophils, resMonocytes,
  resNeutrophils, resNKcells
)

### GSEA with Hallmarks 

m_df <- msigdbr(species = "Homo sapiens")
C3_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%  # H for hallmarks
  dplyr::select(gs_name, entrez_gene)
head(C3_t2g)

run_gsea <- function(myresults, mycomparison) {
  Hids <- bitr(rownames(myresults), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
  dedup_Hids = Hids[!duplicated(Hids[c("ENSEMBL")]),] 
  Hstat.test = myresults[rownames(myresults) %in% dedup_Hids$ENSEMBL,] 
  Hstat.test$Y = dedup_Hids$ENTREZID 
  Hstat.test$ENSEMBL <- rownames(Hstat.test)
  rownames(Hstat.test) <- NULL
  H3kegg_gene_list <- Hstat.test$log2FoldChange 
  names(H3kegg_gene_list) <- Hstat.test$Y 
  H3kegg_gene_list <- na.omit(H3kegg_gene_list) 
  H3kegg_gene_list = sort(H3kegg_gene_list, decreasing = TRUE) 
  set.seed(2244)
  em3 <- GSEA(H3kegg_gene_list, TERM2GENE = C3_t2g, pvalueCutoff = 1)
  em3_g <- setReadable(em3, 'org.Hs.eg.db', 'ENTREZID')
  gsea_res <- em3_g@result[1:11] 
  gsea_res$Comparison <- mycomparison
  gsea_res
}

resV3_gsea <- run_gsea(myresults = resV3, mycomparison = "V3vsV1")
resV2_gsea <- run_gsea(myresults = resV2, mycomparison = "V2vsV1")
resAge_gsea <- run_gsea(myresults = resAge, mycomparison = "Age")
resBcells_gsea <- run_gsea(myresults = resBcells, mycomparison = "Bcells")
resCD4Tcells_gsea <- run_gsea(myresults = resCD4Tcells, mycomparison = "CD4Tcells")
resCD8Tcells_gsea <- run_gsea(myresults = resCD8Tcells, mycomparison = "CD8Tcells")
resDC_gsea <- run_gsea(myresults = resDC, mycomparison = "DC")
resEosinophils_gsea <- run_gsea(myresults = resEosinophils, mycomparison = "Eosinophils")
resMonocytes_gsea <- run_gsea(myresults = resMonocytes, mycomparison = "Monocytes")
resNeutrophils_gsea <- run_gsea(myresults = resNeutrophils, mycomparison = "Neutrophils")
resNKcells_gsea <- run_gsea(myresults = resNKcells, mycomparison = "NKcells")

res_gsea_all <- rbind(
  resV2_gsea, resV3_gsea, resAge_gsea,
  resBcells_gsea, resCD4Tcells_gsea, resCD8Tcells_gsea,
  resDC_gsea, resEosinophils_gsea, resMonocytes_gsea,
  resNeutrophils_gsea, resNKcells_gsea
)

# EXT FIG 1B
P2 <- res_gsea_all %>%
  filter(ID %in% ids_interest) %>%
  mutate(stars=cut(p.adjust, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " "))) %>%
  mutate(Comparison=factor(Comparison, levels = c(unique(res_gsea_all$Comparison)))) %>%
  mutate(ID=factor(ID, levels = rev(ids_interest) )) %>%
  ggplot(aes(x=Comparison, y=ID, fill=NES)) +
  geom_point(shape=21, size=6) +
  geom_text(aes(label=stars), color="black", size=5) +
  scale_fill_gradientn(colours = paletteer::paletteer_c("grDevices::Blue-Red 2", n=20, direction = -1), limits=c(-3.5,3.5)) +
  labs(x=NULL, y=NULL) +  theme_classic() + ggtitle("excluding pDC")


plot_grid(P2,P1,ncol=1)

