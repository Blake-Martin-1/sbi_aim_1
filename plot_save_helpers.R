# Shared helper for saving AIM 1 figures.
# All plots are written as 600 dpi TIFF images in the configured figures directory.
if (!exists("load_path_config")) {
  source("path_config_helpers.R")
}
load_path_config()

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
