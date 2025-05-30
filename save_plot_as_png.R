# Function to save plot as PNG
save_plot_as_png <- function(plot, filename) {
  png(filename, res = 250, width = 4000, height = 2000)
  print(plot)
  dev.off()
}
