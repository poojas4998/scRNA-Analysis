# scRNA-Analysis
scRNA sequencing analysis

Here I will share functions I made to get certain type of plots through Seurat. 

Here are a few changes i made to the dataset.
I use sample not "samples" or "treatments" as a column.

After annotation, i make a column in the metadata 
seurat_object$cell_annotation <- Idents(seurat_object)

The functions are the following:
save_plot_as_png : saves the plots as a png 
save_plot_as_png_wide : saves the plot as a png that has a larger x axis.

# Load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr)) #for separate function to work in find_combined_markers function
suppressPackageStartupMessages(library(plyr)) #dont delete, this library is very important for mapvalues
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(patchwork)) # violin plot printing all at once in one image
suppressPackageStartupMessages(library(devtools)) #upgraded FindMarkers function
#devtools::install_github('immunogenomics/presto')
suppressPackageStartupMessages(library(scCustomize))

current_date <- format(Sys.Date(), "%y%m%d")
