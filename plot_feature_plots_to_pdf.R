# Function to generate feature plots for the top markers and save to a single PDF
plot_feature_plots_to_pdf <- function(seurat_obj, selected_genes, output_pdf) {
  # Start a PDF device
  pdf(output_pdf, width = 8, height = 6) # Adjust the width and height as needed
  
  # Generate and save feature plots
  for (gene in selected_genes) {
    # Generate FeaturePlot for the gene
    p <- FeaturePlot(seurat_obj, features = gene) + ggtitle(paste("FeaturePlot of", gene))
    
    # Print the plot to the PDF
    print(p)
  }
  
  # Close the PDF device
  dev.off()
  
  # Notify the user that the PDF has been created
  cat("Feature plots have been saved to", output_pdf, "\n")
}
