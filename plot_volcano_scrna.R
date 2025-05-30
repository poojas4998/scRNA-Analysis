plot_volcano_scrna <- function(deg_results, 
                               label_threshold_log2FC = 1.3, 
                               label_threshold = 1.3, 
                               label_n = 2, 
                               title = "Volcano Plot: Differential Expression") {
  
  library(ggplot2)
  library(dplyr)
  
  # Ensure the data has the expected columns
  if (!all(c("avg_log2FC", "p_val_adj", "gene_name") %in% colnames(deg_results))) {
    stop("The data must have columns: 'avg_log2FC', 'p_val_adj', and 'gene_name'")
  }
  
  
  # Cap values for display
  deg_results <- deg_results %>%
    mutate(
      capped_log2FC = pmax(pmin(avg_log2FC, 10), -10),
      capped_log10p = pmin(-log10(p_val_adj), 600),
      volc_plot_status = case_when(
        avg_log2FC > label_threshold_log2FC & p_val_adj < 10^(-label_threshold) ~ "UP",
        avg_log2FC < -label_threshold_log2FC & p_val_adj < 10^(-label_threshold) ~ "DOWN",
        TRUE ~ "NS"
      )
    )
  
  # Filter points to be labeled
  labeled_points <- deg_results %>%
    filter(abs(avg_log2FC) > label_threshold_log2FC & -log10(p_val_adj) > label_threshold)
  
  # Label top genes
  up_labeled_points <- labeled_points %>%
    filter(avg_log2FC > 0) %>%
    slice_max(order_by = abs(avg_log2FC) * -log10(p_val_adj), n = label_n, with_ties = FALSE)
  
  down_labeled_points <- labeled_points %>%
    filter(avg_log2FC < 0) %>%
    slice_max(order_by = abs(avg_log2FC) * -log10(p_val_adj), n = label_n, with_ties = FALSE)
  
  labeled_points <- rbind(up_labeled_points, down_labeled_points) %>%
    mutate(
      capped_log2FC = pmax(pmin(avg_log2FC, 10), -10),
      capped_log10p = pmin(-log10(p_val_adj), 600)
    )
  
  # Generate volcano plot
  volcano_plot <- ggplot(deg_results, aes(x = capped_log2FC, y = capped_log10p)) +
    geom_point(aes(color = volc_plot_status), size = 1.5, alpha = 0.7) +
    geom_text(data = labeled_points, aes(label = gene_name), 
              nudge_x = 0.2, nudge_y = 0.2, size = 3, color = "black") +
    xlab(expression("Average Log2 Fold Change (Capped at ±10)")) + 
    ylab(expression("-Log10 Adjusted P-Value (Capped at 600)")) +
    ggtitle(title) +
    geom_hline(yintercept = label_threshold, linetype = 'dashed', color = 'gray') +
    geom_vline(xintercept = c(-label_threshold_log2FC, label_threshold_log2FC), 
               linetype = 'dashed', color = 'gray') +
    scale_color_manual(values = c("UP" = "#D53E4F", "DOWN" = "#3288BD", "NS" = "gray80")) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      legend.position = "right",
      legend.title = element_blank()
    )
  
  return(volcano_plot)
}
