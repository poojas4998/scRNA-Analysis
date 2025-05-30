############ NASTIEST OF NASTIEST CODE TO MAKE EPIC PLOTS ##############


### HEATMAP for cell annotation vs gene for each treatment vs pbs (inspiration gut project 2025)
library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
degs_list <- data.frame()
# Function to find DEGs for a given cytokine
find_degs <- function(seurat_slice, cytokine, control = "box", min_cells = 3) {
  if (length(which(seurat_slice$samples == cytokine)) >= min_cells &&
      length(which(seurat_slice$samples == control)) >= min_cells) {
    
    markers <- FindMarkers(
      seurat_slice, ident.1 = cytokine, ident.2 = control,
      min.pct = 0.2, logfc.threshold = 0.1, test.use = "wilcox"
    )
    
    if (nrow(markers) > 0) {
      markers$fdr <- p.adjust(markers$p_val, "fdr")
      markers$sample <- cytokine
      markers$gene <- rownames(markers)
      return(markers)
    }
  }
  return(NULL)
}


# # Initialize empty list to store results
# degs_list <- list()
# seurat_object <- oIL9R_main
# # Iterate over cell types
# for (cell_type in unique(seurat_object$seurat_clusters)) {
#   print(cell_type)
#   
#   # Subset Seurat object for the current cell type
#   seurat_slice <- subset(seurat_object, seurat_clusters == cell_type)
#   
#   # Apply DEG analysis to each cytokine
#   for (cytokine in all_combinations) {
#     deg_results <- find_degs(seurat_slice, cytokine)
#     if (!is.null(deg_results)) {
#       deg_results$cell_type <- cell_type
#       degs_list <- append(degs_list, list(deg_results))
#     }
#   }
# }
# # Combine all results into a single data frame
# all_markers <- do.call(rbind, degs_list)

# 
# # Set sample identity
# Idents(oIL9R_main) <- "samples"
# 
# seurat_object <- oIL9R_main
# 
# # Function to save plots with increased figure height for better readability
# save_plot_as_png <- function(plot, filename, width = 8, height = 15, dpi = 300) {
#   ggsave(filename, plot = plot, width = width, height = height, dpi = dpi, units = "in")
# }
# 
# # Load necessary libraries
# library(Seurat)
# library(ggplot2)
# library(dplyr)
# library(reshape2)
# 
# # sample_ii <- "21R" 
# # Function to compute and plot expression differences
# for (sample_ii in all_combinations) {
#   
#   # Subset for selected cytokine condition
#   seurat_sub <- subset(seurat_object, subset = samples == sample_ii)
#   seurat_pbs <- subset(seurat_object, subset = samples == "box") # Control condition
#   
#   seurat_sub = subset(seurat_object, subset = samples == sample_ii) 
#   all_markers_sub <- subset(all_markers, sample == sample_ii)
#   all_markers_sub <- all_markers_sub[order(all_markers_sub$avg_log2FC, decreasing = TRUE), ]
#   
#   all_markers_sub <- all_markers_sub %>% group_by(cell_type) %>% slice_head(n = 10)
#   all_markers_sub_unique <- all_markers_sub[!duplicated(all_markers_sub$gene), ]
#   top_genes_unique <- all_markers_sub_unique$gene
#   
#   # Ensure `cell_type` is a factor with correct levels
#   all_markers_sub$cell_type <- factor(all_markers_sub$cell_type, levels = unique(all_markers_sub$cell_type))
#   
#   # Get all unique cell annotations
#   all_cell_annotations <- unique(seurat_object$seurat_clusters)
#   
#   # Extract only genes present in RNA assay
#   valid_genes <- all_markers_sub$gene[all_markers_sub$gene %in% rownames(GetAssayData(seurat_object, slot = "counts"))]
#   
#   # If no valid genes are found, skip iteration
#   if (length(valid_genes) == 0) next
#   
#   # Compute average expression per cell type
#   avg_exp_cytokine <- AverageExpression(seurat_sub, features = valid_genes, group.by = "seurat_clusters")
#   avg_exp_pbs <- AverageExpression(seurat_pbs, features = valid_genes, group.by = "seurat_clusters")
#   
#   # Ensure matrices have the same genes and ordering
#   genes_intersect <- intersect(rownames(avg_exp_cytokine$RNA), rownames(avg_exp_pbs$RNA))
#   avg_exp_cytokine_matched <- avg_exp_cytokine$RNA[genes_intersect, , drop = FALSE]
#   avg_exp_pbs_matched <- avg_exp_pbs$RNA[genes_intersect, , drop = FALSE]
#   
#   # Ensure both matrices have the same cell types
#   common_cell_types <- intersect(colnames(avg_exp_cytokine_matched), colnames(avg_exp_pbs_matched))
#   avg_exp_cytokine_matched <- avg_exp_cytokine_matched[, common_cell_types, drop = FALSE]
#   avg_exp_pbs_matched <- avg_exp_pbs_matched[, common_cell_types, drop = FALSE]
#   
#   # Compute expression differences using Log2 Fold Change (LFC)
#   exp_diff <- log2((avg_exp_cytokine_matched + 1) / (avg_exp_pbs_matched + 1))
#   
#   # Replace NA values with 0
#   exp_diff[is.na(exp_diff)] <- 0
#   # Modify expression differences to cap at -2 and 2
#   exp_diff[exp_diff > 2] <- 2
#   exp_diff[exp_diff < -2] <- -2
#   
#   # Convert to long format for plotting
#   plot_df <- melt(as.matrix(exp_diff))
#   colnames(plot_df) <- c("Gene", "Cell_Type", "Log2_Fold_Change")
#   
#   # Plot heatmap of Log2 Fold Change with capped values
#   p1 <- ggplot(plot_df, aes(x = Cell_Type, y = Gene, fill = Log2_Fold_Change)) + 
#     geom_tile() +
#     theme_classic() +
#     scale_fill_gradientn(
#       name = "Log2 Fold Change", 
#       colors = colorRampPalette(c("navy", "white", "firebrick3"))(60), 
#       limits = c(-2, 2), breaks = c(-2, -1, 0, 1, 2),
#       na.value = "white"  # Avoids grey for missing values
#     ) +
#     theme(
#       panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
#       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
#       axis.text.y = element_text(face = "italic", size = 12),
#       legend.title = element_text(size = 12), 
#       legend.position = "right"
#     )  +
#     scale_y_discrete(position = "right") +
#     ggtitle(paste("Log2 Fold Change:", sample_ii, "vs box")) + 
#     coord_fixed()
#   
#   # Save heatmap
#   save_plot_as_png(p1, paste0("roughwork/heatmap_exp_diff_capped_", sample_ii, ".png"))
#   
#   # Display plot
#   print(p1)
# }
# 
# ###### end of heatmap code
# 
