# Gender affirming Testosterone treatment
In this repository you will find all code to reproduce figures in Lakshmikanth, Consiglio et al - Immune system adaptation during Gender affirming Testosterone treatment (*insert link to article here*).

For all the data in the article, please visit this link and download the zip-file *https://drive.google.com/drive/u/0/folders/10BXCUDLt7sUAl4KxeqyTPaeAnzFoigqT*. The data folder will have multiple folders, some that are Figure specific data and some that are more general for several figures (eg scRNA-seq data that we reuse over multiple figures).

To download all the reproducible code, please clone this repository using *git clone git@github.com:Brodinlab/Gender-affirming-Testosterone-treatment.git*.

For jupyter scripts: For the file-pathways to work directly, you have to open Jupyter Lab/Notebook in the specific figure-folder (eg figure-scripts/Figure2) and put the data-directory at the top of the folder (so the path is Gender-affirming-Testosterone-treatment/data), otherwise you can specify your own data directory.

For the R-scripts, please specify your working directory where the data-folder is located before loading in the data file.

Please note that for some of the figures to be correctly reproduced, first run the jupyter script ending in \\_ preparation and then the R-file.

For any questions regarding data analyses, please contact *petter.brodin@ki.se*


