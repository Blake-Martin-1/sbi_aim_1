# Shared helper for loading file path configuration from a CSV.
# The CSV must include path_key and path_value columns.
load_path_config <- function(config_path = NULL, envir = parent.frame()) {
  cached_paths <- getOption("sbi.path_config_paths")
  if (!is.null(cached_paths)) {
    list2env(cached_paths, envir = envir)
    return(invisible(cached_paths))
  }

  if (is.null(config_path) || !nzchar(config_path)) {
    config_path <- readline(prompt = "Enter the path to the file path configuration CSV: ")
  }

  config_path <- trimws(config_path)
  if (!nzchar(config_path)) {
    stop("A file path configuration CSV is required.")
  }
  if (!file.exists(config_path)) {
    stop("File path configuration CSV does not exist: ", config_path)
  }

  path_config <- read.csv(config_path, stringsAsFactors = FALSE)
  required_cols <- c("path_key", "path_value")
  missing_cols <- setdiff(required_cols, names(path_config))
  if (length(missing_cols) > 0) {
    stop("File path configuration CSV is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  path_config$path_key <- trimws(path_config$path_key)
  path_config$path_value <- trimws(path_config$path_value)
  path_config <- path_config[nzchar(path_config$path_key) & nzchar(path_config$path_value), , drop = FALSE]
  if (anyDuplicated(path_config$path_key)) {
    duplicated_keys <- unique(path_config$path_key[duplicated(path_config$path_key)])
    stop("File path configuration CSV contains duplicate path_key value(s): ", paste(duplicated_keys, collapse = ", "))
  }

  invalid_keys <- path_config$path_key[make.names(path_config$path_key) != path_config$path_key]
  if (length(invalid_keys) > 0) {
    stop("File path configuration CSV contains invalid R variable name(s): ", paste(invalid_keys, collapse = ", "))
  }

  paths <- as.list(stats::setNames(path_config$path_value, path_config$path_key))
  options(sbi.path_config_csv = config_path, sbi.path_config_paths = paths)
  list2env(paths, envir = envir)
  invisible(paths)
}
