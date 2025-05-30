perform_single_go_enrichment <- function(df, 
                                         name = "Sample", 
                                         output_dir = "plots", 
                                         pval_cutoff = 0.05, 
                                         top_n = 10, 
                                         species = "mouse") {
  
  # Step 1: Check for essential columns
  if (!"gene_name" %in% colnames(df)) {
    stop("⚠️ The input data frame must contain a 'gene_name' column.")
  }
  
  if (!"p_val_adj" %in% colnames(df)) {
    stop("⚠️ The input data frame must contain a 'p_val_adj' column.")
  }
  
  # Step 2: Filter significant genes
  df_filtered <- df %>%
    filter(p_val_adj < pval_cutoff) %>%
    arrange(desc(avg_log2FC))
  
  if (nrow(df_filtered) == 0) {
    cat("⚠️ No significant genes in:", name, "\n")
    return(NULL)
  }
  
  # Step 3: Gene ID conversion (SYMBOL -> ENTREZ)
  OrgDb <- if (species == "mouse") org.Mm.eg.db else org.Hs.eg.db
  
  gene_entrez <- bitr(
    df_filtered$gene_name, 
    fromType = "SYMBOL", 
    toType = "ENTREZID", 
    OrgDb = OrgDb
  )
  
  if (nrow(gene_entrez) == 0) {
    cat("⚠️ No Entrez IDs mapped for:", name, "\n")
    return(NULL)
  }
  
  # Step 4: GO enrichment analysis
  go_enrichment <- enrichGO(
    gene = gene_entrez$ENTREZID,
    OrgDb = OrgDb,
    ont = "BP",  # Biological Process
    pvalueCutoff = pval_cutoff,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  if (is.null(go_enrichment) || nrow(go_enrichment) == 0) {
    cat("⚠️ GO enrichment returned no terms for:", name, "\n")
    return(NULL)
  }
  
  # Step 5: Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Step 6: Create GO plots (Barplot & Dotplot)
  bp <- barplot(go_enrichment, showCategory = top_n, title = paste("GO Barplot -", name)) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12))
  
  dp <- dotplot(go_enrichment, showCategory = top_n, title = paste("GO Dotplot -", name)) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12))
  
  # Step 7: Save plots as PNG
  ggsave(filename = file.path(output_dir, paste0("GO_Barplot_", name, ".png")), plot = bp, width = 8, height = 6)
  ggsave(filename = file.path(output_dir, paste0("GO_Dotplot_", name, ".png")), plot = dp, width = 8, height = 6)
  
  # Save enriched GO results as CSV
  result_file <- file.path(output_dir, paste0("GO_Enrichment_", name, ".csv"))
  write.csv(as.data.frame(go_enrichment), file = result_file, row.names = FALSE)
  
  cat("✅ GO enrichment done and results saved for:", name, "\n")
  
  return(go_enrichment)
}
