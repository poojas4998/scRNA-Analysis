# Load necessary libraries
# if (!requireNamespace("gtools", quietly = TRUE)) install.packages("gtools")
library(gtools)

# Define function to compute all set intersections and optionally save them
generate_venn_intersections <- function(gene_sets, output_file = NULL) {
  venn_intersections <- list()
  
  # Loop over all combination lengths
  for (r in 1:length(gene_sets)) {
    combos <- combinations(n = length(gene_sets), r = r, v = names(gene_sets), set = TRUE)
    for (i in seq_len(nrow(combos))) {
      combo <- combos[i, ]
      group_name <- paste(combo, collapse = " & ")
      intersection_genes <- Reduce(intersect, gene_sets[combo])
      venn_intersections[[group_name]] <- intersection_genes
    }
  }
  
  # Save to file if requested
  if (!is.null(output_file)) {
    writeLines(
      unlist(lapply(names(venn_intersections), function(x) {
        paste0(x, ":\n", paste(venn_intersections[[x]], collapse = ", "), "\n")
      })),
      output_file
    )
    cat("✅ Intersections written to:", output_file, "\n")
  }
  
  return(venn_intersections)
}
