# Immune dictionary like fig 2b

plot_deg_distribution_by_condition <- function(seurat_object, marker_list, output_dir = "roughwork", palette = barbie_colors) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  combined_deg_summary <- data.frame()
  individual_plot_list <- list()
  
  conditions <- names(marker_list)
  first_condition <- conditions[1]
  other_conditions <- conditions[-1]
  
  # Internal reusable function
  process_condition <- function(condition_name, show_y_axis = TRUE) {
    cat("🔍 Processing:", condition_name, "\n")
    
    genes_of_interest <- marker_list[[condition_name]]$gene_name
    cells_of_interest <- WhichCells(seurat_object, expression = samples == condition_name)
    
    if (length(genes_of_interest) == 0 || length(cells_of_interest) == 0) {
      cat("⚠️ Skipping", condition_name, "- no genes or cells found.\n")
      return(NULL)
    }
    
    meta_subset <- seurat_object@meta.data[cells_of_interest, ]
    meta_subset$seurat_clusters <- factor(meta_subset$seurat_clusters)
    
    expr_matrix <- GetAssayData(seurat_object, assay = "RNA", slot = "data")[genes_of_interest, cells_of_interest]
    binary_expr <- expr_matrix > 0
    
    cluster_vec <- meta_subset$seurat_clusters
    names(cluster_vec) <- rownames(meta_subset)
    cluster_per_cell <- cluster_vec[colnames(binary_expr)]
    
    gene_expr_per_cluster <- t(binary_expr) %>%
      as.data.frame() %>%
      dplyr::mutate(cluster = cluster_per_cell) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(across(everything(), ~ any(. > 0))) %>%
      dplyr::mutate(Num_DEGs = rowSums(across(where(is.logical)))) %>%
      dplyr::select(cluster, Num_DEGs) %>%
      dplyr::mutate(condition = condition_name)
    
    # Append global summary
    combined_deg_summary <<- bind_rows(combined_deg_summary, gene_expr_per_cluster)
    
    # Create plot
    p <- ggplot(gene_expr_per_cluster, aes(x = cluster, y = Num_DEGs, fill = cluster)) +
      geom_bar(stat = "identity", color = "black", size = 0.5) +
      coord_flip() +
      scale_fill_manual(values = palette) +
      labs(title = paste("DEGs per Cluster -", condition_name), x = NULL, y = "# DEGs") +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "none",
        axis.text.y = if (show_y_axis) element_text(size = 10) else element_blank(),
        axis.text.x = element_text(size = 10)
      )
    
    # Save PNG
    safe_name <- gsub("[^a-zA-Z0-9]", "_", condition_name)
    png_name <- file.path(output_dir, paste0("DEGs_cluster_", safe_name, ".png"))
    ggsave(png_name, p, width = 4, height = 6)
    
    individual_plot_list[[condition_name]] <<- p
  }
  
  # Process first with y-axis
  process_condition(first_condition, TRUE)
  
  # Process rest without y-axis
  for (cond in other_conditions) {
    process_condition(cond, FALSE)
  }
  
  # Combine plots
  combined_plot <- wrap_plots(individual_plot_list, nrow = 1)
  
  # Save combined
  ggsave(file.path(output_dir, "facet_deg_clusters_by_condition.png"), combined_plot, width = 6 * length(individual_plot_list), height = 6)
  
  return(combined_plot)
}
