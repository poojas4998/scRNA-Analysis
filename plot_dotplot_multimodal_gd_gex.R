plot_dotplot_multimodal_gd_gex <- function(seurat_object, 
                                           features,title_text,
                                           color_palette = "RdYlBu", ident_name = "seurat_clusters") {
  
  # Step 1: Set cluster identity and order
  Idents(seurat_object) <- ident_name
  seurat_object$seurat_clusters <- factor(
    seurat_object$seurat_clusters
  )
  
  # Step 2: Plot DotPlot
  dotplot <- DotPlot(seurat_object, features = features,cols = c("RdYlBu")) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + 
    guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
    RotatedAxis() +
    labs(title = title_text) +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA),
      text = element_text(size = 10),
      panel.grid.major.x = element_line(color = "grey80"),
      panel.grid.major.y = element_line(color = "grey80") 
    )
  
  return(dotplot)
}
