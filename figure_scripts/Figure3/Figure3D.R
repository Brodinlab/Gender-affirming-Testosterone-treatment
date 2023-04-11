#Script to generate Figure3D in Lakshmikanth, Consiglio et al - Immune system adaptation during Gender affirming Testosterone treatment
#Author Rikard Forlin - rikard.forlin@ki.se
library(nichenetr)
library(tidyverse)
library(SeuratData)
library(SeuratDisk)
library(Seurat)
library(anndata)
library(ggplot2)
library(circlize)
library(colormap)

options(warn=-1)
# First load the ligand-target matrices from NicheNet and the weighted networks.
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#ligand_target_matrix[1:5,1:5]
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))


#Read in adata made for NicheNet analysis (found in jupyter script, "Prepare for NicheNet")
# and make it into a seurat object
adata <- read_h5ad("nichenet_adata.h5ad")#Insert pathway to annotated dataframe here
rownames(adata$obs) <- make.names(rownames(adata$obs), unique = TRUE)

X <- as.matrix(adata$X)
seuratObj <- CreateSeuratObject(counts = t(X), meta.data = adata$obs)
seuratObj@meta.data$celltype = adata$obs$celltype
seuratObj@meta.data$celltype_sub = adata$obs$celltype_sub
seuratObj@meta.data$V = adata$obs$V


#Remove cell types that we are making up a small percentage of the population, X is unknown cell types
seuratObj@meta.data$celltype = as.character(seuratObj@meta.data$celltype)
inds = rownames(seuratObj@meta.data[seuratObj@meta.data$celltype %in% c('MAIT', 'Megakaryocyte', 'gdT', 'X'), ])
seuratObj <- seuratObj[, ! rownames(seuratObj@meta.data) %in% inds]

#Set identity to cell type
Idents(object = seuratObj) <- seuratObj@meta.data$celltype

#This is a bit messy and I'm sorry for that - my R skills are under improvement...
#But here we get expressed genes from each cell type, and look into the ligand_target_matrix from NicheNet to define a background gene set
#Only genes that is a target for a ligand is included in the background_expressed_genes for downstream analysis
receiver = 'Monocyte'
expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10)
receiver2 = 'CD8T'
expressed_genes_receiver2 = get_expressed_genes(receiver2, seuratObj, pct = 0.10)
receiver3 = 'CD4T'
expressed_genes_receiver3 = get_expressed_genes(receiver3, seuratObj, pct = 0.10)
receiver4 = 'NK'
expressed_genes_receiver4 = get_expressed_genes(receiver4, seuratObj, pct = 0.10)
receiver5 = 'pDC'
expressed_genes_receiver5 = get_expressed_genes(receiver5, seuratObj, pct = 0.10)
receiver6 = 'DC'
expressed_genes_receiver6= get_expressed_genes(receiver6, seuratObj, pct = 0.10)
receiver7 = 'B'
expressed_genes_receiver7= get_expressed_genes(receiver7, seuratObj, pct = 0.10)
all_receiver = c('CD8T', 'CD4T', 'NK', 'pDC','DC','B', 'Monocyte')



expressed_genes_receiver = c(expressed_genes_receiver, expressed_genes_receiver2,expressed_genes_receiver3,expressed_genes_receiver4, expressed_genes_receiver5, expressed_genes_receiver6, expressed_genes_receiver7)
expressed_genes_receiver = unique(expressed_genes_receiver)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

#Define what cell type is the sender - and look for expressed genes for that cell type
#This is the same as expressed_genes_receiver, but for clarification I left it here
sender_celltypes = 'Monocyte'
expressed_genes_sender = get_expressed_genes(sender_celltypes, seuratObj, pct = 0.10)
expressed_genes_sender = unique(expressed_genes_sender)

#Define a gene set of interest: here we take differentially expressed genes in pre and post testosterone-treatment
seurat_obj_receiver= subset(seuratObj, idents = all_receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["V"]])
condition_oi = "V2"
condition_reference = "V1" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")
DE_table_receiver <- DE_table_receiver %>% drop_na()

#Filter the DEG-geneset
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#Define a set of potential ligands expressed by the sender populations and bind a receptor expressed on the receiver/target population
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
#head(potential_ligands)

#Perform NicheNet ligand activity analysis: rank the potential ligands based on the presence of their target genes in the gene set of intereset (compared to the background set of genes)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))

# Auroc, Aupr, pearson are measure for how well a ligand can predict the observed DEGs compared to the background set of expressed genes
# In their validation study, they showed that pearson correlation coefficient between a ligands target predictions and the observed transcriptional response was the most informative measure to define ligand activity.
#ligand_activities[1:10,]
#Take the top 20 of pearson correlation coefficient
best_upstream_ligands = ligand_activities %>% top_n(15, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

########################################################################################
##Circos plot visualization to show active ligand-target links between interacting cells
##https://github.com/saeyslab/nichenetr/blob/2cb19de12d7ff8d3d8491302a43e5b3a6d3c0ae9/vignettes/circos.md

#We only show ligands with a weight higher than a predefined cutoff: links belonging to the 63% of lowest scores were removed. 
#Find and set cutoff for the highest weighted ligands.
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 20000) %>% bind_rows() %>% drop_na()
active_ligand_target_links_df <- active_ligand_target_links_df[order(active_ligand_target_links_df$weight, decreasing = TRUE), ]
active_ligand_target_links_df <- head(active_ligand_target_links_df, 100)
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.66)

#Filter ligands with the cutoff
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)
ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())
circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

circos_links$ligand_type = 'Monocyte'
circos_links$target_type = 'all'


#Prepare the circos visualization: give each segment of ligands and targets a specific color and order
grid_col_tbl_ligand = tibble(ligand_type = "Monocyte")
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
#monocyte = purple
#NK = #FF3300 (red)
#pDC = #3333FF (blue)
#CD4T = #339966 (green)
#DC = orange
#CD8T = yellow
#B = black

target_color$color_target_type <- c("#3333FF", "#FF3300", "purple", "purple", "purple", "purple", "purple", "purple","yellow", 
                                    "purple", "purple", "#FF3300", "purple", "#3333FF", "#3333FF", "purple", 
                                    "orange", "orange", "purple","purple", "purple","purple", "purple", "#FF3300", "purple")


#target_color$color_target_type = colormap(colormap = colormaps$viridis, length(target_color[[1]]))
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
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Monocyte") %>% distinct(ligand) %>% nrow() -1)),
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




