# ---------------------------------------------------------------------------- #
# file:    generate_IER.R
# purpose: Generate the NIH Inclusion Enrollment Report (IER) participant-level
#          CSV from the prospective AIM 1 dataset.
# ---------------------------------------------------------------------------- #

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

# Generate dataframe to use for IER table generation when the prospective
# subgroup dataframe is already loaded. Each row should represent a unique
# participant/MRN and include: mrn, age_years, race, ethnicity, is_female.
if (exists("pros_combined_for_subgroups")) {
  ier_df <- pros_combined_for_subgroups %>%
    dplyr::select(mrn, age_years, race, ethnicity, is_female) %>%
    distinct() # 2,694 unique patients
}

# The NIH CSV upload expects exactly these column names.
IER_COLUMNS <- c("Race", "Ethnicity", "Sex", "Age", "Age Unit")

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || !nzchar(as.character(x[[1]]))) y else x
}

get_arg_value <- function(args, flag, default = NULL) {
  hit <- match(flag, args)
  if (!is.na(hit) && length(args) >= hit + 1) {
    return(args[[hit + 1]])
  }
  default
}

get_output_path <- function(args = commandArgs(trailingOnly = TRUE)) {
  output_path <- get_arg_value(args, "--output") %||%
    Sys.getenv("IER_OUTPUT_FILE_PATH", unset = NA_character_)

  if (is.null(output_path) || is.na(output_path) || !nzchar(output_path)) {
    if (!interactive()) {
      stop("Provide an output CSV path with --output or IER_OUTPUT_FILE_PATH.", call. = FALSE)
    }
    output_path <- readline(prompt = "Enter the filepath where the NIH IER CSV should be saved: ")
  }

  output_path <- trimws(output_path)
  if (!nzchar(output_path)) {
    stop("An output CSV filepath is required.", call. = FALSE)
  }

  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist and could not be created: ", output_dir, call. = FALSE)
  }

  output_path
}

standardize_race <- function(x) {
  y <- str_squish(as.character(x))
  y_lower <- str_to_lower(y)

  dplyr::case_when(
    is.na(y) | y == "" | y_lower %in% c("unknown", "not reported", "unknown/not reported", "declined", "patient refused", "other_or_unknown", "other or unknown", "other") ~ "Unknown",
    str_detect(y_lower, "more than one|multiple|multi") ~ "More than one race",
    str_detect(y_lower, "american indian|alaska native|native american") ~ "American Indian",
    str_detect(y_lower, "asian") ~ "Asian",
    str_detect(y_lower, "black|african american") ~ "Black",
    str_detect(y_lower, "hawaiian|pacific islander") ~ "Hawaiian",
    str_detect(y_lower, "white|caucasian") ~ "White",
    TRUE ~ "Unknown"
  )
}

standardize_ethnicity <- function(x) {
  y <- str_squish(as.character(x))
  y_lower <- str_to_lower(y)

  dplyr::case_when(
    is.na(y) | y == "" | y_lower %in% c("unknown", "not reported", "unknown/not reported", "declined", "patient refused", "other_or_unknown", "other or unknown") ~ "Unknown",
    str_detect(y_lower, "not hispanic|non[- ]?hispanic") ~ "Not Hispanic or Latino",
    str_detect(y_lower, "hispanic|latino|latina|latinx") ~ "Hispanic or Latino",
    TRUE ~ "Unknown"
  )
}

standardize_sex_from_is_female <- function(is_female) {
  y <- suppressWarnings(as.numeric(is_female))

  dplyr::case_when(
    is.na(y) ~ "Unknown",
    y == 1 ~ "Female",
    y == 0 ~ "Male",
    TRUE ~ "Unknown"
  )
}

age_to_nih_fields <- function(age_years) {
  age_years <- suppressWarnings(as.numeric(age_years))

  age_unit <- dplyr::case_when(
    is.na(age_years) ~ "Unknown",
    age_years >= 90 ~ "Ninety Plus",
    age_years >= 1 ~ "Years",
    age_years >= 1 / 12 ~ "Months",
    age_years >= 0 ~ "Days",
    TRUE ~ "Unknown"
  )

  age_value <- dplyr::case_when(
    is.na(age_years) | age_unit %in% c("Unknown", "Ninety Plus") ~ NA_real_,
    age_unit == "Years" ~ floor(age_years),
    age_unit == "Months" ~ floor(age_years * 12),
    age_unit == "Days" ~ floor(age_years * 365.25),
    TRUE ~ NA_real_
  )

  tibble(Age = as.integer(age_value), `Age Unit` = age_unit)
}

make_ier <- function(ier_df) {
  required_cols <- c("mrn", "age_years", "race", "ethnicity", "is_female")
  missing_cols <- setdiff(required_cols, names(ier_df))
  if (length(missing_cols) > 0) {
    stop("ier_df is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  ier_df %>%
    arrange(.data$mrn) %>%
    distinct(.data$mrn, .keep_all = TRUE) %>%
    transmute(
      Race = standardize_race(.data$race),
      Ethnicity = standardize_ethnicity(.data$ethnicity),
      Sex = standardize_sex_from_is_female(.data$is_female),
      age_years_for_ier = suppressWarnings(as.numeric(.data$age_years))
    ) %>%
    bind_cols(age_to_nih_fields(.$age_years_for_ier)) %>%
    select(all_of(IER_COLUMNS))
}

write_ier_csv <- function(ier_df, output_path) {
  ier <- make_ier(ier_df)
  readr::write_csv(ier, output_path, na = "")
  message("Wrote ", nrow(ier), " participant rows to ", output_path)
  invisible(ier)
}

main <- function() {
  if (!exists("ier_df", envir = .GlobalEnv)) {
    stop("ier_df must exist before running generate_IER.R.", call. = FALSE)
  }

  output_path <- get_output_path()
  write_ier_csv(get("ier_df", envir = .GlobalEnv), output_path)
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
