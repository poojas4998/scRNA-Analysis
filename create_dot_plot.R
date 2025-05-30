create_dot_plot <- function(seurat_object, 
                            marker_genes = NULL, 
                            scale_plot = TRUE, 
                            angle = 90, 
                            vjust_value = 0.5,
                            color_palette = c("blue", "white", "red"), 
                            size_factor = 4, title_name) {
  
  # Default marker genes if none provided
  if (is.null(marker_genes)) {
    marker_genes <- c(
      "Ptprc", "Cd19", "Ms4a1", "Cd79a", "Cd79b", "H2-Ab1", "H2-Eb1", "Cd3g", "Cd3e",
      "Cd8a", "Cd8b1", "Nkg7", "Tcrg-C1", "Torg-C2", "Itgae", "Foxp3", "Ctla4", "112ra",
      "Tnfrsf4", "Pdcd1", "ICOS", "Cd28", "Mki67", "Gzma", "Gzmb", "Ncr1", "Cd4", "Ccr6",
      "Ccdc184", "KIrb1b", "Fam184b", "Cd68", "Ly6c2", "Lag3", "Havcr1", "Sell", "Siglech",
      "Irf8", "Itgax", "Btla", "Havcr2", "Cd207", "Cd86", "Cd40", "Ly75", "Cd24a", "Xcr1",
      "Clec9a", "Vwf", "Lyz2", "Mg/2", "Sirpa", "Adora2a", "Ccr7", "Apol7c", "Cd63",
      "Tnfrsf9", "H2-M2", "Fut4", "Cd80", "Aire", "C1qa", "Lyz1", "Marco", "Folr2",
      "Siglec1", "Ms4a7", "Mertk", "Itgam", "Csf1r", "Adgre1", "Ccr2", "Cd14", "Fegr3",
      "Ifih1", "Isg15", "Cd33", "Ly6g", "Stfa211", "Mmp9", "Mcpt4", "Fcerla", "Cd69",
      "Ly6c1", "Nrp1", "Pecam1", "Emcn", "Thod", "Lyve1", "Prox1", "Ecscr", "Col1a2",
      "C4b", "Dcn"
    )
  }
  
  # Create the dot plot using Seurat's DotPlot function
  dot_plot <- DotPlot(seurat_object, features = marker_genes, scale = scale_plot) +
    scale_colour_gradientn(colours = color_palette) +
    ggtitle(title_name)
    theme(axis.text.x = element_text(angle = angle, vjust = vjust_value),
          axis.text.y = element_text(size = 10, face = "bold"),
          plot.title = element_text(size = 14, face = "bold")) 
  
  return(dot_plot)
}
