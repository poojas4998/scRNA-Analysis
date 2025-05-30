# scRNA-Analysis
scRNA sequencing analysis

Here I will share functions I made to get certain type of plots through Seurat. 

Here are a few changes i made to the dataset.
I use sample not "samples" or "treatments" as a column. Rename columns as this to make the following functions to run.

After annotation, i make a column in the metadata 
seurat_object$cell_annotation <- Idents(seurat_object)

The functions are the following:

save_plot_as_png : saves the plots as a png 

save_plot_as_png_wide : saves the plot as a png that has a larger x axis.

find_combined_markers_singletreatment: Runs Find Markers. Ensure you define all_combinations as for example 
                                      all_combinations <- unique(oIL9R_main$samples)[!unique(oIL9R_main$samples) %in%  c("9R")] 
                                      all_combinations <- unique(seurat_object$sample)[!unique(seurat_object$sample) %in%  c("PBS")]
                                      
subset_genes: finds Significant differentially expressed genes, upregulated significant genes and downregulated significant genes and saves them in the folder. Also this result is later divided into the samples and saved in the folder.

create_dot_plot: Creates dotplots which help in annotation. By default if marker genes are not presented it shows the genes found in the immunodictionary paper presented by the cui lab: https://www.nature.com/articles/s41586-023-06816-9#data-availability

plot_dotplot_multimodal_gd_gex: fancy way to plot a dotpot based on https://www.nature.com/articles/s41590-023-01710-y  (Figure 2 A)

plot_feature_plots_to_pdf: All featureplots of interest in a pdf file

create_violin_plot : All violin plots of interest in a pdf file

plot_volcano_scrna: Volcano plots rendered

perform_go_enrichment: Gene Ontology Enrichment (provide seurat object and marker data) - Outdated

perform_single_go_enrichment: Gene Ontology Enrichment (provide dataframe) (Updated)

plot_deg_distribution_by_condition: DEG Distribution ggplot

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
