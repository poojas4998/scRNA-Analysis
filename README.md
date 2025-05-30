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
