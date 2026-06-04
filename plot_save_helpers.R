# Shared helper for saving AIM 1 figures.
# All plots are written as 600 dpi TIFF images in the shared figures directory.
aim1_figures_dir <- "/phi/sbi/sbi_blake/figures_aim_1"

save_aim1_plot <- function(plot, filename, width = 8, height = 6) {
  if (is.null(plot)) {
    return(invisible(NULL))
  }

  safe_filename <- gsub("[\\\\/]+", "_", filename)
  if (!grepl("\\.tiff?$", safe_filename, ignore.case = TRUE)) {
    safe_filename <- paste0(safe_filename, ".tiff")
  }

  dir.create(aim1_figures_dir, recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    filename = file.path(aim1_figures_dir, safe_filename),
    plot = plot,
    device = "tiff",
    dpi = 600,
    width = width,
    height = height,
    units = "in",
    compression = "lzw"
  )

  invisible(plot)
}
