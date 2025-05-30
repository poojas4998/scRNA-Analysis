# Define a function to create a violin plot for selected genes
create_violin_plot <- function(seurat_obj, genes_of_interest, group_by_column,output_pdf, point_size = 0.1) {
  pdf(output_pdf, width = 8, height = 6)
  for (gene in genes_of_interest) {
    p <- VlnPlot_scCustom(
      seurat_obj,
      features = gene,
      group.by = group_by_column,
      plot_boxplot = TRUE
    ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)
      ) +
      labs(
        title = paste("Gene Expression for", gene,""),
        x = group_by_column,
        y = "Expression Level"
      ) + NoLegend()
    scale_fill_brewer(palette = "Pastel2")  # Use a visually appealing color palette
    print(p)
  }
  dev.off()
  cat("Violin plots have been saved to", output_pdf, "\n")
}
