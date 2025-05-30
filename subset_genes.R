library(dplyr)
library(purrr)

subset_genes <- function(df, p_val_adj_threshold = 0.05, project_folder = "DEG_exports") {
  # Step 1: Subset genes
  upregulated <- df %>% filter(avg_log2FC > 0.5 & p_val_adj < p_val_adj_threshold)
  downregulated <- df %>% filter(avg_log2FC < -0.5 & p_val_adj < p_val_adj_threshold)
  significant <- df %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < p_val_adj_threshold)
  
  # Step 2: Prepare folder
  padj_name <- sprintf("%03d", as.numeric(gsub("\\.", "", as.character(p_val_adj_threshold * 100))))
  # output_folder <- file.path(project_folder, paste0("padj_", padj_name))
  output_folder <- file.path(paste0("padj_", padj_name))
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Step 3: Save combined CSVs
  write.csv(upregulated, file = file.path(output_folder, "upsigDEGs.csv"), row.names = FALSE)
  write.csv(downregulated, file = file.path(output_folder, "downsigDEGs.csv"), row.names = FALSE)
  write.csv(significant, file = file.path(output_folder, "allsigDEGs.csv"), row.names = FALSE)
  
  # Step 4: Helper to split by updated_condition and export
  split_sort_save <- function(df, category_suffix) {
    df_list <- df %>%
      split(.$updated_condition) %>%
      map(~ arrange(.x, desc(avg_log2FC)))

    # Save each condition's file
    walk2(df_list, names(df_list), function(data, condition) {
      safe_name <- paste0(gsub("[^a-zA-Z0-9]", "_", condition), "_", category_suffix)
      write.csv(data, file = file.path(output_folder, paste0(safe_name, ".csv")), row.names = FALSE)
    })
    
    return(df_list)
  }
  
  # Step 5: Split, save, and return all in a structured list
  result <- list(
    combined = list(
      upregulated = upregulated,
      downregulated = downregulated,
      significant = significant
    ),
    per_condition = list(
      upsig = split_sort_save(upregulated, "upsig"),
      downsig = split_sort_save(downregulated, "downsig"),
      sig = split_sort_save(significant, "sig")
    )
  )
  
  return(result)
}
