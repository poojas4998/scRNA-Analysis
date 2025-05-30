# Function to Rename Clusters and Save Seurat Object with Same Variable Name
annotate_and_save_seurat <- function(seurat_object_name, new_cluster_ids, file_name, umap_output_path) {
  # Use get to dynamically access the object
  seurat_object <- get(seurat_object_name)
  
  # Set cluster identities and rename clusters based on the provided list
  Idents(seurat_object) <- "seurat_clusters"
  seurat_object <- RenameIdents(seurat_object, new_cluster_ids)
  seurat_object$cell_annotation <- Idents(seurat_object)
  # Assign the modified object back to the same variable name
  assign(seurat_object_name, seurat_object, envir = .GlobalEnv)
  
  # Save the Seurat object with the same variable name
  save(list = seurat_object_name, file = file_name)
  message("Seurat object saved at: ", file_name)
  
  # Generate UMAP plot with annotations
  umap_plot <- DimPlot(seurat_object, reduction = "umap", label = TRUE, label.size = 2.5, repel = TRUE)
  print(umap_plot)
  
  # Save UMAP plot as PNG
  save_plot_as_png(umap_plot,umap_output_path)
  message("UMAP plot saved at: ", umap_output_path)
  
  # Return a summary of the new cluster annotations
  return(table(seurat_object$samples))
}
