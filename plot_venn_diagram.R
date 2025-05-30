plot_venn_diagram <- function(
    gene_sets,
    set_names = NULL,
    fill_colors = NULL,
    stroke_size = 0.5,
    set_name_size = 6,
    seed = 20190708
) {
  set.seed(seed)
  
  # Check and install ggvenn if needed
  if (!requireNamespace("ggvenn", quietly = TRUE)) {
    devtools::install_github("yanlinlin82/ggvenn")
  }
  library(ggvenn)
  
  n_sets <- length(gene_sets)
  if (is.null(set_names)) {
    set_names <- names(gene_sets)
  } else {
    names(gene_sets) <- set_names
  }
  
  if (is.null(fill_colors)) {
    default_palette <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#00A087FF", "#DC0000FF")
    fill_colors <- default_palette[1:n_sets]
  }
  
  # Plot with ggvenn
  ggvenn(
    gene_sets,
    fill_color = fill_colors,
    stroke_size = stroke_size,
    set_name_size = set_name_size
  )
}

