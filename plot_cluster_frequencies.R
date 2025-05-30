plot_cluster_frequencies <- function(datasets, output_dir_base = "preliminary") {
  library(ggplot2)
  library(dplyr)
  library(scales)
  
  # Helper to save plots
  save_plot_as_png <- function(plot, path, width = 10, height = 6, dpi = 300) {
    ggsave(filename = path, plot = plot, width = width, height = height, dpi = dpi)
  }
  
  dir.create(output_dir_base, showWarnings = FALSE, recursive = TRUE)
  
  for (dataset_name in names(datasets)) {
    obj <- datasets[[dataset_name]]
    
    # Compute cluster frequencies
    freq_df <- obj@meta.data %>%
      dplyr::count(sample, cell_annotation) %>%
      group_by(sample) %>%
      mutate(freq = n / sum(n)) %>%
      ungroup() %>%
      mutate(dataset = dataset_name)
    
    # Sort sample so PBS comes first
    sample_levels <- unique(freq_df$sample)
    sample_levels <- c("PBS", setdiff(sort(sample_levels), "PBS"))
    freq_df$sample <- factor(freq_df$sample, levels = sample_levels)
    
    # Generate pastel color palette
    n_clusters <- length(unique(freq_df$cell_annotation))
    pastel_colors <- hue_pal(l = 70)(n_clusters)
    
    # Create bar plot
    p <- ggplot(freq_df, aes(x = sample, y = freq, fill = cell_annotation)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = pastel_colors) +
      theme_minimal(base_size = 14) +
      ylab("Cluster Frequency") +
      xlab("Cytokine Treatment (sample)") +
      ggtitle(paste("Cluster Frequencies -", dataset_name)) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 16),
        legend.title = element_blank()
      )
    
    print(p)
    
    # Save plot
    output_dir <- file.path(output_dir_base, dataset_name)
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    save_plot_as_png(p, file.path(output_dir, "cluster_frequency_barplot.png"))
  }
  
  cat("✅ Cluster frequency bar plots saved in:", output_dir_base, "\n")
}
