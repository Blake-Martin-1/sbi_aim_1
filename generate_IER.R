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

# Generate dataframe to use for IER table generation
ier_df <- pros_combined_for_subgroups %>% dplyr::select(mrn, age_years, race, ethnicity, is_female) %>% distinct() #2,694 unique patients

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

find_first_column <- function(data, candidates, required_for = NULL) {
  hit <- intersect(candidates, names(data))
  if (length(hit) > 0) {
    return(hit[[1]])
  }
  if (!is.null(required_for)) {
    stop(
      "Could not find a ", required_for, " column. Tried: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  NA_character_
}

standardize_race <- function(x) {
  y <- str_squish(as.character(x))
  y_lower <- str_to_lower(y)

  dplyr::case_when(
    is.na(y) | y == "" | y_lower %in% c("unknown", "not reported", "unknown/not reported", "declined", "patient refused", "other_or_unknown", "other or unknown") ~ "Unknown",
    str_detect(y_lower, "more than one|multiple|multi") ~ "More than one race",
    str_detect(y_lower, "american indian|alaska native|native american") ~ "American Indian/Alaska Native",
    str_detect(y_lower, "asian") ~ "Asian",
    str_detect(y_lower, "black|african american") ~ "Black or African American",
    str_detect(y_lower, "hawaiian|pacific islander") ~ "Native Hawaiian or Other Pacific Islander",
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

standardize_sex <- function(x) {
  y <- str_squish(as.character(x))
  y_lower <- str_to_lower(y)

  dplyr::case_when(
    is.na(y) | y == "" | y_lower %in% c("unknown", "not reported", "unknown/not reported", "u") ~ "Unknown",
    y_lower %in% c("female", "f", "woman", "girl") ~ "Female",
    y_lower %in% c("male", "m", "man", "boy") ~ "Male",
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

make_ier <- function(prospective_data) {
  if ("epoch" %in% names(prospective_data)) {
    prospective_data <- prospective_data %>%
      filter(str_to_lower(as.character(.data$epoch)) %in% c("pros", "prospective"))
  }

  id_col <- find_first_column(
    prospective_data,
    c("study_id", "pat_mrn_id", "mrn", "MRN", "PAT_MRN_ID", "PAT_ENC_CSN_ID", "csn"),
    required_for = "participant identifier"
  )
  race_col <- find_first_column(prospective_data, c("race", "Race", "RACE", "old_race", "OLD_RACE"), required_for = "race")
  ethnicity_col <- find_first_column(prospective_data, c("ethnicity", "Ethnicity", "ETHNICITY"), required_for = "ethnicity")
  sex_col <- find_first_column(prospective_data, c("sex", "Sex", "SEX", "gender", "Gender", "is_female"), required_for = "sex")
  age_col <- find_first_column(prospective_data, c("age", "age_years", "Age", "AGE", "age_at_admission", "AGE_YEARS"), required_for = "age in years")

  prospective_data %>%
    arrange(.data[[id_col]]) %>%
    distinct(.data[[id_col]], .keep_all = TRUE) %>%
    transmute(
      Race = standardize_race(.data[[race_col]]),
      Ethnicity = standardize_ethnicity(.data[[ethnicity_col]]),
      Sex = if (sex_col == "is_female") {
        case_when(is.na(.data[[sex_col]]) ~ "Unknown", as.logical(.data[[sex_col]]) ~ "Female", TRUE ~ "Male")
      } else {
        standardize_sex(.data[[sex_col]])
      },
      age_years_for_ier = suppressWarnings(as.numeric(.data[[age_col]]))
    ) %>%
    bind_cols(age_to_nih_fields(.$age_years_for_ier)) %>%
    select(all_of(IER_COLUMNS))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  input_path <- get_arg_value(args, "--input") %||%
    Sys.getenv("PROSPECTIVE_MODEL_DATA_FILE_PATH", unset = NA_character_) %||%
    getOption("prospective_model_data_file_path")
  output_path <- get_arg_value(args, "--output", "IER_participant_level.csv")

  if (is.null(input_path) || is.na(input_path) || !nzchar(input_path)) {
    stop("Provide an input CSV with --input or PROSPECTIVE_MODEL_DATA_FILE_PATH.", call. = FALSE)
  }
  if (!file.exists(input_path)) {
    stop("Input CSV does not exist: ", input_path, call. = FALSE)
  }

  prospective_data <- readr::read_csv(input_path, show_col_types = FALSE)
  ier <- make_ier(prospective_data)
  readr::write_csv(ier, output_path, na = "")
  message("Wrote ", nrow(ier), " prospective participant rows to ", output_path)
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
