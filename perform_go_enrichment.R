library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(enrichplot)
library(DOSE)

# Function to perform GO enrichment for each treatment
perform_go_enrichment <- function(seurat_object, marker_data, pval_cutoff = 0.05, qval_cutoff = 0.2, save_results = TRUE) {
  
  # Extract unique treatment conditions
  treatments <- unique(seurat_object@meta.data$updated_condition)
  go_results_list <- list()
  
  # Loop through each treatment condition
  for (treatment in treatments) {
    cat("Processing treatment:", treatment, "\n")
    
    # Subset markers for the current treatment
    treatment_markers <- marker_data[marker_data$updated_condition == treatment, ]
    
    # Filter for significant marker genes
    treatment_markers_filtered <- subset(treatment_markers, p_val_adj < pval_cutoff)
    
    # Check if there are significant markers
    if (nrow(treatment_markers_filtered) == 0) {
      cat("No significant markers for treatment:", treatment, "\n")
      next
    }
    
    # Sort by fold change
    treatment_markers_sorted <- treatment_markers_filtered[order(-treatment_markers_filtered$avg_log2FC), ]
    
    # Extract gene symbols
    treatment_genes <- treatment_markers_sorted$gene
    
    # Convert gene symbols to Entrez IDs
    gene_entrez <- bitr(treatment_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    
    # Check if gene conversion was successful
    if (nrow(gene_entrez) == 0) {
      cat("No Entrez IDs found for treatment:", treatment, "\n")
      next
    }
    
    # Perform GO enrichment analysis
    go_enrichment <- enrichGO(
      gene = gene_entrez$ENTREZID,
      OrgDb = org.Mm.eg.db,
      ont = "BP",
      pvalueCutoff = pval_cutoff,
      qvalueCutoff = qval_cutoff,
      readable = TRUE
    )
    
    # Store results
    go_results_list[[paste0("Treatment_", treatment)]] <- go_enrichment
    
    # Save results if enabled
    if (save_results) {
      output_prefix <- paste0("GO_Treatment_", treatment)
      
      # Save CSV
      write.csv(as.data.frame(go_enrichment), paste0(output_prefix, ".csv"), row.names = FALSE)
      
      # Generate and save barplot
      barplot_filename <- paste0(output_prefix, "_Barplot.png")
      png(barplot_filename, width = 800, height = 600)
      print(barplot(go_enrichment, showCategory = 15, title = paste("GO Terms (Treatment:", treatment, ")")))
      dev.off()
      
      # Generate and save dotplot
      dotplot_filename <- paste0(output_prefix, "_Dotplot.png")
      png(dotplot_filename, width = 800, height = 600)
      print(dotplot(go_enrichment, showCategory = 15, title = paste("GO Terms (Treatment:", treatment, ")")))
      dev.off()
      
      cat("Saved results for treatment:", treatment, "\n")
    }
  }
  
  cat("GO enrichment analysis complete.\n")
  return(go_results_list)
}
