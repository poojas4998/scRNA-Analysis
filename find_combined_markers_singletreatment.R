find_combined_markers_singletreatment <- function(seurat_object, all_combinations,group_column = "sample",control_name = "PBS") {
  # seurat_nk_cells need to be a seurat object
  DefaultAssay(seurat_object) <- "RNA"
  cluster_results <- lapply(all_combinations, function(cluster_id) {
    markers <- FindMarkers(seurat_object,
                           ident.1 = cluster_id,
                           ident.2 = control_name,
                           group.by = group_column,
                           verbose = FALSE) %>%
      rownames_to_column(var = "gene_name") %>%
      mutate(
        Comparison = paste(cluster_id, "vs." , paste0(control_name)),
        updated_condition = cluster_id
      ) 
    return(as.data.frame(markers))
  })
  combined_markers <- do.call(rbind, cluster_results)
  return(combined_markers)
}
