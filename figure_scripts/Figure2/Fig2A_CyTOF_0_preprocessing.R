# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## CyTOF preprocessing, for FIGURES 2A, 2B, 2C ##
library(dplyr)
library(stringr)
library(flowCore)
library(ggplot2)
library(tibble)
library(sva)
library(data.table)
library(BiocParallel)
library(cowplot)
library(ggrastr)

# 0 - Preprocessing of all SRT CyTOF data

# Batches:
# EXP20CL8796
# EXP20DG3622
# EXP21DG3628
# EXP22DG3665

setwd("Please insert working directory here to point at directory '../data/Figure2/'")

#download CyTOF data from https://flowrepository.org/id/RvFrXOT8W7LbhnWqDpab7sf10EIFURBQBJp1AtW99dFK7jjAObTKsHvWl3DjLMwN
# info on which samples are on which batch is found in the metadata file
meta_all <- read.csv("SRT_allmetadata.csv", row.names = 1)

# save fcs files under their batch directory, so batches can first be processed separately to then be combined

## EXP-20-CL8796 ##
parent.folder <- "/EXP-20-CL8796/"
EXP20CL8796_sub.folders <- list.dirs(parent.folder, recursive =TRUE)[-1] #36 subfolders
EXP20CL8796 <- list()

for(i in EXP20CL8796_sub.folders) {
  setwd(i)
  # read FCS files
  myfcsfile1 <- list.files(pattern = ".fcs")
  myfcsfile <- flowCore::read.FCS(myfcsfile1, truncate_max_range = FALSE) 
  # add markers as colnames
  colnames(myfcsfile) <- flowCore::parameters(myfcsfile)$desc
  # exclude unnecessary markers
  toexclude <- c("Time", "Event_length", 
                 "102Pd","104Pd","105Pd","106Pd","108Pd","116Cd","131Xe","133Cs","140Ce",
                 "DNAIr191","DNAIr193",
                 "Center", "Offset", "Width", "Residual")
  myfcsfile <- myfcsfile[,!(colnames(myfcsfile) %in% toexclude)]
  # rename markers so they match with other batches
  colnames(myfcsfile)[colnames(myfcsfile) %in% "4-1BB"] <- "CD137" #this is the case for EXP20CL8796
  colnames(myfcsfile)[colnames(myfcsfile) %in% "CD3"] <- "CD3e" #not the case for EXP20CL8796
  colnames(myfcsfile)[colnames(myfcsfile) %in% "TCRgd"] <- "gdTCR" #not the case for EXP20CL8796
  colnames(myfcsfile) <- str_remove(colnames(myfcsfile), "-")
  # transform data 
  asinhTransNEWcof <- arcsinhTransform(transformationId="defaultArcsinhTransform", a=0, b=0.2, c=0)
  translistNEWcof <- transformList(colnames(myfcsfile), asinhTransNEWcof) 
  dataTransformNEWcof <- transform(myfcsfile, translistNEWcof)
  # get expression
  dataset <- exprs(dataTransformNEWcof) %>% as.data.frame()
  # read gridslim annotation
  csvfile <- list.files(pattern = ".csv")
  mycsvfile <- read.csv(csvfile, sep=',', row.names = 1)
  finalframe <- data.frame(mycsvfile$level0, mycsvfile$level1, mycsvfile$level2, dataset)
  # change rownames
  rownames(finalframe) <- paste(
    # batch
    "EXP20CL8796", 
    # subject 
    sapply(str_split(sapply(str_split(unique(sapply(str_split(list.files(), "\\."), "[", 1)), "_d"), "[", 1), " "), "[", 1),
    #visit
    sapply(str_split(sapply(str_split(unique(sapply(str_split(list.files(), "\\."), "[", 1)), "_d"), "[", 1), " "), "[", 2) %>% str_remove("\\(") %>% str_remove("\\)"),
    # fcs name
    str_remove(paste(myfcsfile1), ".fcs"),
    # grid slim annotation
    paste(finalframe$mycsvfile.level1, finalframe$mycsvfile.level2, sep="_"),
    # rownames
    rownames(dataset), 
    sep = "__")
  #remove beads and dead cells
  finalframe <- finalframe %>% dplyr::filter(mycsvfile.level0 %in% "cells") %>% dplyr::filter(!(mycsvfile.level1 %in% " ")) 
  finalframe <- finalframe %>% select(!c("mycsvfile.level0", "mycsvfile.level1", "mycsvfile.level2"))
  EXP20CL8796[[i]] <- finalframe 
  print(i)
}

final_EXP20CL8796 <- bind_rows(EXP20CL8796)
final_EXP20CL8796[1:5,1:5]
colnames(final_EXP20CL8796)
ncol(final_EXP20CL8796)
unique(sapply(str_split(rownames(final_EXP20CL8796), "__"), "[", 1)) #batch
unique(sapply(str_split(rownames(final_EXP20CL8796), "__"), "[", 2)) #subject
unique(sapply(str_split(rownames(final_EXP20CL8796), "__"), "[", 3)) #visit
unique(sapply(str_split(rownames(final_EXP20CL8796), "__"), "[", 4)) #fcsfile
unique(sapply(str_split(rownames(final_EXP20CL8796), "__"), "[", 5)) #cell
rm(dataset, dataTransformNEWcof, finalframe, mycsvfile, myfcsfile,translistNEWcof, asinhTransNEWcof, csvfile, i, parent.folder, toexclude)

## EXP-20-DG3622 ##
parent.folder <- "/EXP-20-DG3622/"
EXP20DG3622sub.folders <- list.dirs(parent.folder, recursive =TRUE)[-1] #31 subfolders
EXP20DG3622 <- list()

for(i in EXP20DG3622sub.folders) {
  setwd(i)
  # read FCS files
  myfcsfile1 <- list.files(pattern = ".fcs")
  myfcsfile <- flowCore::read.FCS(myfcsfile1, truncate_max_range = FALSE) 
  # add markers as colnames
  colnames(myfcsfile) <- flowCore::parameters(myfcsfile)$desc
  # exclude unnecessary markers
  toexclude <- c("Time", "Event_length", 
                 "102Pd","104Pd","105Pd","106Pd","108Pd","116Cd","131Xe","133Cs","140Ce",
                 "DNAIr191","DNAIr193",
                 "Center", "Offset", "Width", "Residual")
  myfcsfile <- myfcsfile[,!(colnames(myfcsfile) %in% toexclude)]
  # rename markers so they match with other batches
  colnames(myfcsfile)[colnames(myfcsfile) %in% "4-1BB"] <- "CD137" 
  colnames(myfcsfile)[colnames(myfcsfile) %in% "CD3"] <- "CD3e" 
  colnames(myfcsfile)[colnames(myfcsfile) %in% "TCRgd"] <- "gdTCR" 
  colnames(myfcsfile) <- str_remove(colnames(myfcsfile), "-")
  # transform data 
  asinhTransNEWcof <- arcsinhTransform(transformationId="defaultArcsinhTransform", a=0, b=0.2, c=0)
  translistNEWcof <- transformList(colnames(myfcsfile), asinhTransNEWcof) 
  dataTransformNEWcof <- transform(myfcsfile, translistNEWcof)
  # get expression
  dataset <- exprs(dataTransformNEWcof) %>% as.data.frame()
  # read gridslim annotation
  csvfile <- list.files(pattern = ".csv")
  mycsvfile <- read.csv(csvfile, sep=',', row.names = 1)
  finalframe <- data.frame(mycsvfile$level0, mycsvfile$level1, mycsvfile$level2, dataset)
  # change rownames
  rownames(finalframe) <- paste(
    # batch
    "EXP20DG3622", 
    # subject 
    str_remove(sapply(str_split(unique(sapply(str_split(list.files(), "\\."), "[", 1)), " "), "[", 1),"-"),
    #visit
    sapply(str_split(sapply(str_split(unique(sapply(str_split(list.files(), "\\."), "[", 1)), " "), "[", 2), "_"), "[", 2),
    # fcs name
    str_remove(paste(myfcsfile1), ".fcs"),
    # grid slim annotation
    paste(finalframe$mycsvfile.level1, finalframe$mycsvfile.level2, sep="_"),
    # rownames
    rownames(dataset), 
    sep = "__")
  #remove beads and dead cells
  finalframe <- finalframe %>% dplyr::filter(mycsvfile.level0 %in% "cells") %>% dplyr::filter(!(mycsvfile.level1 %in% " ")) 
  finalframe <- finalframe %>% select(!c("mycsvfile.level0", "mycsvfile.level1", "mycsvfile.level2"))
  EXP20DG3622[[i]] <- finalframe 
  print(i)
}

final_EXP20DG3622 <- bind_rows(EXP20DG3622)
final_EXP20DG3622[1:5,1:5]
colnames(final_EXP20DG3622)
ncol(final_EXP20DG3622)
unique(sapply(str_split(rownames(final_EXP20DG3622), "__"), "[", 1)) #batch
unique(sapply(str_split(rownames(final_EXP20DG3622), "__"), "[", 2)) #subject -- corrected name of files sto15 and sto18, which didnt follow same naming strategy as other files in this batch; other files follow STO-10 FtM_2; and these followed STO-15_1
unique(sapply(str_split(rownames(final_EXP20DG3622), "__"), "[", 3)) #visit
unique(sapply(str_split(rownames(final_EXP20DG3622), "__"), "[", 4)) #fcs file
unique(sapply(str_split(rownames(final_EXP20DG3622), "__"), "[", 5)) #cell
setdiff(colnames(final_EXP20DG3622), colnames(final_EXP20CL8796))
rm(dataset, dataTransformNEWcof, finalframe, mycsvfile, myfcsfile,translistNEWcof, asinhTransNEWcof, csvfile, i, parent.folder, toexclude)


## EXP-21-DG3628 ##
parent.folder <- "/EXP-21-DG3628/"
EXP21DG3628sub.folders <- list.dirs(parent.folder, recursive =TRUE)[-1] #16 subfolders
EXP21DG3628 <- list()

meta_all <- read.csv("meta_EXP21DG3628.csv", row.names = 1)

for(i in EXP21DG3628sub.folders) {
  setwd(i)
  # read FCS files
  myfcsfile1 <- list.files(pattern = ".fcs")
  myfcsfile <- flowCore::read.FCS(myfcsfile1, truncate_max_range = FALSE) 
  # add markers as colnames
  colnames(myfcsfile) <- flowCore::parameters(myfcsfile)$desc
  # exclude unnecessary markers
  toexclude <- c("Time", "Event_length", 
                 "102Pd","104Pd","105Pd","106Pd","108Pd","116Cd","131Xe","133Cs","140Ce",
                 "DNAIr191","DNAIr193",
                 "Center", "Offset", "Width", "Residual")
  myfcsfile <- myfcsfile[,!(colnames(myfcsfile) %in% toexclude)]
  # rename markers so they match with other batches
  colnames(myfcsfile)[colnames(myfcsfile) %in% "4-1BB"] <- "CD137" 
  colnames(myfcsfile)[colnames(myfcsfile) %in% "CD3"] <- "CD3e" 
  colnames(myfcsfile)[colnames(myfcsfile) %in% "TCRgd"] <- "gdTCR" 
  colnames(myfcsfile) <- str_remove(colnames(myfcsfile), "-")
  # transform data 
  asinhTransNEWcof <- arcsinhTransform(transformationId="defaultArcsinhTransform", a=0, b=0.2, c=0)
  translistNEWcof <- transformList(colnames(myfcsfile), asinhTransNEWcof) 
  dataTransformNEWcof <- transform(myfcsfile, translistNEWcof)
  # get expression
  dataset <- exprs(dataTransformNEWcof) %>% as.data.frame()
  # read gridslim annotation
  csvfile <- list.files(pattern = ".csv")
  mycsvfile <- read.csv(csvfile, sep=',', row.names = 1)
  finalframe <- data.frame(mycsvfile$level0, mycsvfile$level1, mycsvfile$level2, dataset)
  # change rownames
  rownames(finalframe) <- paste(
    # batch
    "EXP21DG3628", 
    # subject 
    meta_all$SubjectID[match(unique(sapply(str_split(list.files(), "\\."), "[", 1)), meta_all$GRID_ID)],
    #visit
    meta_all$Visit[match(unique(sapply(str_split(list.files(), "\\."), "[", 1)), meta_all$GRID_ID)],
    # fcs name
    str_remove(paste(myfcsfile1), ".fcs"),
    # grid slim annotation
    paste(finalframe$mycsvfile.level1, finalframe$mycsvfile.level2, sep="_"),
    # rownames
    rownames(dataset), 
    sep = "__")
  #remove beads and dead cells
  finalframe <- finalframe %>% dplyr::filter(mycsvfile.level0 %in% "cells") %>% dplyr::filter(!(mycsvfile.level1 %in% " ")) 
  finalframe <- finalframe %>% select(!c("mycsvfile.level0", "mycsvfile.level1", "mycsvfile.level2"))
  EXP21DG3628[[i]] <- finalframe 
  print(i)
}

final_EXP21DG3628 <- bind_rows(EXP21DG3628)
final_EXP21DG3628[1:5,1:5]
colnames(final_EXP21DG3628)
ncol(final_EXP21DG3628)
unique(sapply(str_split(rownames(final_EXP21DG3628), "__"), "[", 1)) #batch
unique(sapply(str_split(rownames(final_EXP21DG3628), "__"), "[", 2)) #subject
unique(sapply(str_split(rownames(final_EXP21DG3628), "__"), "[", 3)) #visit
unique(sapply(str_split(rownames(final_EXP21DG3628), "__"), "[", 4)) #fcs file
unique(sapply(str_split(rownames(final_EXP21DG3628), "__"), "[", 5)) #cells
setdiff(colnames(final_EXP21DG3628), colnames(final_EXP20CL8796))
setdiff(colnames(final_EXP20DG3622), colnames(final_EXP21DG3628))
rm(dataset, dataTransformNEWcof, finalframe, mycsvfile, myfcsfile,translistNEWcof, asinhTransNEWcof, csvfile, i, parent.folder, toexclude)


## EXP-22-DG3665 ##
setwd("./EXP22DG3665/")
mydf <- list()
EXP22DG3665_files <- list.files(pattern = ".fcs")

for(i in EXP22DG3665_files) {
  FCSfile <- flowCore::read.FCS(i, truncate_max_range = FALSE) 
  # add marker names to columns
  colnames(FCSfile) <- flowCore::parameters(FCSfile)$desc
  # exclude markers
  toexclude <- c("Time", "Event_length", "Center", "Offset", "Width", "Residual", 
                 "Amplitude", '88Sr', '190BCKG',
                 "102Pd", "103Rh", "104Pd", "105Pd", "106Pd", "108Pd", # these are now excluded here 
                 '120Sn', '127I', '138Ba', '208Pb', #'113In'  is IgD in this experiment
                 "116Cd", "131Xe", "X133Cs") #not excluding "140Ce_EQ" and "191Ir", "193Ir_Intercalator-Ir" to later get rid of beads
  FCSfile <- FCSfile[,!(colnames(FCSfile) %in% toexclude)] 
  # arcsin h transform to stabilize variance
  asinhTransNEWcof <- arcsinhTransform(transformationId="defaultArcsinhTransform", a=0, b=0.2, c=0)
  translistNEWcof <- transformList(colnames(FCSfile), asinhTransNEWcof) 
  dataTransformNEWcof <- transform(FCSfile, translistNEWcof)
  # get expression
  dataset <- exprs(dataTransformNEWcof) %>% as.data.frame() %>% mutate(Sample=str_remove(paste(i), ".fcs"))
  mydf[[i]] <- dataset 
  print(i)
}

myfinaldf <- bind_rows(mydf)

protocol <- read.csv("protocol.csv", sep = ";")
protocol <- protocol[1:19,] %>% dplyr::filter(!Subject.ID %in% "Empty")
protocol$Subject.ID <- str_remove(protocol$Subject.ID, "\n")

# add step here to name rows based on metadata of file
# format for the first 3 batches: # batch __ subject __ visit __ fcsname __ level1 __ level2 __ rownames
# format for the fourth batch: # batch __ subject __ visit __ fcsname __ rownames
rownames(myfinaldf) <- paste("EXP22DG3665", 
                             protocol$Subject.ID[match(myfinaldf$Sample, protocol$Sample.ID)],
                             str_remove(protocol$Sample.info[match(myfinaldf$Sample, protocol$Sample.ID)], "Visit "),
                             myfinaldf$Sample,
                             rownames(myfinaldf), sep = "__")

# filter out beads
myfinaldf <- myfinaldf %>% dplyr::filter(`140Ce_EQ Bead` < 2.5)

# filter out dead cells, with no 191 or 193
myfinaldf <- myfinaldf %>% dplyr::filter(`193Ir_Intercalator-Ir` > 1 & `191Ir` > 0.5)

setdiff(colnames(myfinaldf), colnames(final_EXP20CL8796)) 
# "113In" "HLA-DR" "140Ce_EQ Bead" "CD3" "TCRgd" "Siglec-8" "191Ir" "193Ir_Intercalator-Ir" "Sample"
setdiff(colnames(final_EXP20CL8796), colnames(myfinaldf)) 
# "IgD"     "HLADR"   "CD3e"    "gdTCR"   "Siglec8"

# rename colnames(myfinaldf)
colnames(myfinaldf)[colnames(myfinaldf) %in% "113In"] <- "IgD"
colnames(myfinaldf)[colnames(myfinaldf) %in% "HLA-DR"] <- "HLADR"
colnames(myfinaldf)[colnames(myfinaldf) %in% "CD3"] <- "CD3e"
colnames(myfinaldf)[colnames(myfinaldf) %in% "TCRgd"] <- "gdTCR"
colnames(myfinaldf)[colnames(myfinaldf) %in% "Siglec-8"] <- "Siglec8"

setdiff(colnames(myfinaldf), colnames(final_EXP20CL8796)) 
# "140Ce_EQ Bead"         "191Ir"                 "193Ir_Intercalator-Ir" "Sample"    
setdiff(colnames(final_EXP20CL8796), colnames(myfinaldf)) 
# none     

myfinaldf <- myfinaldf %>% select(!c("140Ce_EQ Bead","191Ir","193Ir_Intercalator-Ir","Sample"))
colnames(myfinaldf)
ncol(myfinaldf)
rownames(myfinaldf)[1:5]

unique(sapply(str_split(rownames(myfinaldf), "__"), "[", 1)) #batch
unique(sapply(str_split(rownames(myfinaldf), "__"), "[", 2)) #subject
unique(sapply(str_split(rownames(myfinaldf), "__"), "[", 3)) #visit
unique(sapply(str_split(rownames(myfinaldf), "__"), "[", 4)) #fcs file


## merge the batches ##
colnames(final_EXP20CL8796)
colnames(final_EXP20DG3622)
colnames(final_EXP21DG3628)
colnames(myfinaldf)

all_batches <- rbind(final_EXP20CL8796,
                     final_EXP20DG3622,
                     final_EXP21DG3628,
                     myfinaldf)


rm(dataset, dataTransformNEWcof, EXP20CL8796, EXP20DG3622, EXP21DG3628, FCSfile, meta_all, mydf, protocol, translistNEWcof, asinhTransNEWcof,
   i, EXP20CL8796_sub.folders, EXP20DG3622sub.folders, EXP21DG3628sub.folders, myfcsfile1, toexclude,
   EXP22DG3665_files)

## save cell by marker matrix - without batch correction ##
write.csv(all_batches, "220811_CyTOF_1234batches_uncorrected.csv")
# 19 516 904 cells by 48 markers

### Batch correct with Combat ###

# batch correction using sva's ComBat
batch=sapply(str_split(rownames(all_batches), "__"), "[", 1)
#unique(batch)
combat_edata_B = ComBat(dat=t(all_batches), batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
#combat_edata_B[1:5,1:5]
combat_edata_B <- as.data.frame(t(combat_edata_B))

## Save files following batch correction ##
# flowframe
ff = flowFrame(data.matrix(combat_edata_B)) 
write.FCS(ff, filename = '220811_CyTOF_1234batches_ComBat_corrected.fcs')

# cell by marker - batch corrected matrix
write.csv(combat_edata_B, "220811_CyTOF_1234batches_ComBat_corrected.csv") 

# labels of the cell by marker - batch corrected matrix
combat_edata_B_label <- data.frame(
  label=rownames(combat_edata_B),
  batch=sapply(str_split(rownames(combat_edata_B), "__"), "[", 1),
  SubjectID=sapply(str_split(rownames(combat_edata_B), "__"), "[", 2),
  Visit=sapply(str_split(rownames(combat_edata_B), "__"), "[", 3),
  FCSfile=sapply(str_split(rownames(combat_edata_B), "__"), "[", 4)
)

colnames(combat_edata_B_label)
unique(combat_edata_B_label$batch)

combat_edata_B_label$level1 <- ifelse(
  combat_edata_B_label$batch %in% c("EXP20CL8796","EXP20DG3622","EXP21DG3628"), 
  combat_edata_B_label$level1 <- sapply(str_split(combat_edata_B_label$label, "__"), "[", 5),
  combat_edata_B_label$level1 <- NA)

unique(combat_edata_B_label$level1)

write.csv(combat_edata_B_label, '220811_CyTOF_1234batches_ComBat_corrected_labels.csv')


