#Script to generate Figure4A in Lakshmikanth, Consiglio et al - Immune system adaptation during gender affirming testosterone therapy
#Author of script Rikard Forlin - rikard.forlin@ki.se
library(nichenetr)
library(tidyverse)
library(Seurat)
library(anndata)
library(ggplot2)
library(circlize)
library(colormap)
library(SeuratDisk)


#### GENELISTS ####
#Th1-response from https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/BOSCO_TH1_CYTOTOXIC_MODULE.html
th1_response <- c(
  "ABCB1", "ADAM19", "ADAM3A", "ADORA3", "ANKRD22", "APOL3", "ARNT2", "ATF5", "BATF2", "BCL2",
  "BCL2L14", "TMEM229B", "ANXA2R", "KNL1", "CCL7", "CCL8", "CCR5", "CD163L1", "CD28", "CD300E",
  "CD38", "CD5", "CDC7", "CDK6", "CH25H", "CMKLR1", "CMPK2", "CRABP1", "CTLA4", "CXCL10", "CXCL11",
  "CXCL9", "DNA2", "DPCD", "ENPP2", "F13A1", "F2R", "FAM20A", "CALHM6", "STRIP2", "FBXO39", "FCGR2B",
  "FFAR3", "FUT2", "GAS6", "GBP6", "GIMAP4", "GIMAP5", "GIMAP7", "GIMAP8", "GNLY", "GPR171", "GZMB",
  "GZMK", "HAMP", "HAPLN3", "HERC6", "HESX1", "HS3ST3B1", "IFNB1", "IFNG", "IL10", "IL12RB2", "IL15",
  "IL15RA", "IL21", "IL2RA", "IL4I1", "IDO1", "IRF4", "ACOD1", "ISG15", "ITGA9", "CEMIP", "KLHDC1",
  "KLRD1", "CERS4", "LILRB5", "MERTK", "MKI67", "MS4A6E", "MYBL1", "NKG7", "OR2A5", "OR5D14", "OR6K6",
  "P2RX5", "PLEKHO1", "PRF1", "PRKCA", "PRKCQ", "PTGER2", "RGL1", "RHBDF2", "RNASE2", "RRM2", "RTP4",
  "RUNX3", "SAMD3", "SDS", "SH2D1A", "SLA2", "SLAMF1", "SMOX", "SOCS1", "SSTR2", "STAB1", "STAT4",
  "TNFSF18", "TRAT1", "TSHZ3", "TXK", "TYMS", "USP18"
)


th1_drivers <- c('TNF','IL15', 'IL6', 'IL12B') #IL12A too sparsely expressed

######## Run NicheNet-analysis ########
options(warn=-1)
# First load the ligand-target matrices from NicheNet and the weighted networks.
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

#Read in adata made for NicheNet analysis (found in jupyter script for Figure4_scRNA-seq, "preparing for NicheNet")
count_data_str <- 'insert pathway generated from jupyter script "Figure4Aprep_4F.ipynb" '
meta_data_str <- 'insert pathway generated from jupyter script "Figure4Aprep_4F.ipynb" '

expression_data <- read.csv(count_data_str, row.names = 1)
expression_data <- t(expression_data)

metadata <- read.csv(meta_data_str, row.names = 1)

seuratObj <- CreateSeuratObject(counts = expression_data)
seuratObj@meta.data <- metadata


#If you want to convert into SeuratObject to not read the .csv files every time:
#Convert("/Users/rikardforlin/Forskning/Data/GenderAffirmingTestosteroneTreatment/NicheNet-data/LPS/ForNicheNetLPS_GeneCountOver3kNoLayers.h5ad", dest = "h5seurat", overwrite = FALSE)
#seuratObj <- LoadH5Seurat("/Users/rikardforlin/Forskning/Data/GenderAffirmingTestosteroneTreatment/NicheNet-data/LPS/ForNicheNetLPS_GeneCountOver3kNoLayers.h5seurat", meta.data = FALSE, misc = FALSE, assays = "RNA")
#seuratObj@meta.data <- metadata

#Normalize the data
seuratObj <- NormalizeData(seuratObj)

#####################################################################

all_cells <- c('CD8T', 'NK', 'Monocyte')
inds = rownames(seuratObj@meta.data[seuratObj@meta.data$CellType %in% all_cells, ])
seuratObj <- seuratObj[, rownames(seuratObj@meta.data) %in% inds]
#Set identity to cell type
Idents(object = seuratObj) <- seuratObj@meta.data$CellType

# List of all cell types
#all_receiver <- c('CD8T', 'CD4T', 'NK')
all_receiver <- c('CD8T', 'NK')

# Get all all expressed genes into a single vector
all_expressed_genes_receiver <- rownames(seuratObj@assays$RNA$data)

############################################################################################
# target genes in rows, ligands in columns
background_expressed_genes = all_expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

#Define what cell type is the sender - and look for expressed genes for that cell type
ct_sender = 'Monocyte'
sender_celltypes = ct_sender

# target genes in rows, ligands in columns
expressed_genes_sender <- all_expressed_genes_receiver %>% .[. %in% colnames(ligand_target_matrix)]

seurat_obj_receiver= subset(seuratObj, idents = all_receiver)

Idents(seurat_obj_receiver) <- "stimulation"
condition_oi = "lps-V2"
condition_reference = "lps-V1" 

#Define a gene set of interest: here we take differentially expressed genes in pre and post testosterone-treatment
DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.1,logfc.threshold = 0.1) %>% rownames_to_column("gene") #change_here?
DE_table_receiver <- DE_table_receiver %>% drop_na()
DE_table_receiver <- DE_table_receiver[DE_table_receiver$gene %in% th1_response, ]

#Filter the DEG-geneset
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_log2FC >= 0.15) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#Define gene set of interest for the sender cells (Monocytes)

seurat_obj_sender = subset(seuratObj, idents = ct_sender)
#seurat_obj_sender  = SetIdent(seurat_obj_sender , value = seurat_obj_sender [["stimulation"]])
Idents(seurat_obj_sender) <- "stimulation"

DE_table_sender = FindMarkers(object = seurat_obj_sender , ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.1, logfc.threshold = 0.1) %>% rownames_to_column("gene") # logfc.threshold = 0.25
DE_table_sender <- DE_table_sender %>% drop_na()
DE_table_sender <- DE_table_sender[DE_table_sender$gene %in% th1_drivers, ]

geneset_oi_ligands = DE_table_sender %>% pull(gene)
geneset_oi_ligands = geneset_oi_ligands %>% .[. %in% colnames(ligand_target_matrix)]


#Define a set of potential ligands expressed by the sender populations and bind a receptor expressed on the receiver/target population
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,geneset_oi_ligands)
expressed_receptors = intersect(receptors,all_expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()


#Perform NicheNet ligand activity analysis: rank the potential ligands based on the presence of their target genes in the gene set of intereset (compared to the background set of genes)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))

# Auroc, Aupr, pearson are measure for how well a ligand can predict the observed DEGs compared to the background set of expressed genes
# In their validation study, they showed that pearson correlation coefficient between a ligands target predictions and the observed transcriptional response was the most informative measure to define ligand activity.
#Take the top 100 of pearson correlation coefficient - this is very high but we'll filter later on
best_upstream_ligands = ligand_activities %>% top_n(100, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

########################################################################################
##Circos plot visualization to show active ligand-target links between interacting cells
##https://github.com/saeyslab/nichenetr/blob/2cb19de12d7ff8d3d8491302a43e5b3a6d3c0ae9/vignettes/circos.md

#We only show ligands with a weight higher than a predefined cutoff: links belonging to the 66% of lowest scores were removed. 
#Find and set cutoff for the highest weighted ligands.
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 20000) %>% bind_rows() %>% drop_na()
active_ligand_target_links_df <- active_ligand_target_links_df[order(active_ligand_target_links_df$weight, decreasing = TRUE), ]


active_ligand_target_links_df <- head(active_ligand_target_links_df, 100)
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.66)

active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

#Filter ligands with the cutoff
ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())
circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

circos_links$ligand_type = ct_sender
circos_links$target_type = 'all'

#Prepare the circos visualization: give each segment of ligands and targets a specific color and order
grid_col_tbl_ligand = tibble(ligand_type = ct_sender)
grid_col_tbl_target = tibble(target_type = "all")

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand)
ccolors <- colormap(colormap = colormaps$cubehelix, length(ligand_color[[1]]))
ligand_color$color_ligand_type <- ccolors
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
####################################################################################################
###Create specific colors according to which celltype has the highest mean gene count of target gene
target_color = circos_links %>% distinct(target)

target_color$color_target_type = colormap(colormap = colormaps$viridis, length(target_color[[1]]))

#### Inserted color according to mean expression of gene (calculated in python)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)
grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

#Prepare the circos visualization: order ligands and targets
target_order = circos_links$target %>% unique()
ligand_order = circos_links$ligand
order = c(ligand_order,target_order)
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == ct_sender) %>% distinct(ligand) %>% nrow() -1)),
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "all") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
)


circos.clear()
circos.par(gap.degree = 1)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #

unique(circos_links$target)
