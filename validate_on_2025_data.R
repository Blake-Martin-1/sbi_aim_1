###### Code to test the model developed on prospective 2023-2024 data on the first 4-5 months of 2025 ###
library(dplyr)
library(tidyr)
library(lubridate)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)


epic_5_25 <- readr::read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_May2025.csv",
  na = c("", "NULL"),
  col_types = readr::cols(.default = readr::col_character()),
  progress = FALSE,
  lazy = FALSE
)

## Drop columns ending in _RAW
epic_5_25_noraw <- epic_5_25[, !grepl("_RAW$", names(epic_5_25)), drop = FALSE]
epic_5_25_noraw <- epic_5_25_noraw[, !grepl("_TS_LEN$", names(epic_5_25_noraw)), drop = FALSE]
epic_5_25_noraw <- epic_5_25_noraw[, !grepl("_IS_COMPLETE$", names(epic_5_25_noraw)), drop = FALSE]

# Remove remainder of NA columns
all_null_cols <- names(epic_5_25_noraw)[
  sapply(epic_5_25_noraw, function(x) all(is.na(x)))
]

all_null_cols
length(all_null_cols)

## Create a mapping table so you can inspect old -> new names
name_map <- data.frame(
  old_name = names(epic_5_25_noraw),
  new_name = tolower(names(epic_5_25_noraw)),
  in_predictors = tolower(names(epic_5_25_noraw)) %in% predictors,
  stringsAsFactors = FALSE
)

## Apply the rename
epic_5_25_renamed <- epic_5_25_noraw
names(epic_5_25_renamed) <- name_map$new_name

## Optional: move predictor columns to the front, keep all others after them
epic_5_25_renamed <- epic_5_25_renamed %>%
  dplyr::select(dplyr::any_of(predictors), dplyr::everything())

## Optional checks
missing_predictors <- setdiff(predictors, names(epic_5_25_renamed))
extra_columns <- setdiff(names(epic_5_25_renamed), predictors)

## View results
name_map %>% dplyr::arrange(dplyr::desc(in_predictors), new_name)
missing_predictors
extra_columns

## Ok now process data to make sure format matches the dataframe used by the models ##
pros_25 <- epic_5_25_renamed %>% rename(score_time = score_time_utc, model_score = total_score, model_type = acuity_system_id)

pros_25$model_type[pros_25$model_type == "100133"] <- "LR_no_abx"
pros_25$model_type[pros_25$model_type == "100134"] <- "LR_yes_abx"
pros_25$model_type[pros_25$model_type == "100148"] <- "RF_no_abx"
pros_25$model_type[pros_25$model_type == "100149"] <- "RF_yes_abx"

# Lower case of column names
colnames(pros_25) <- str_to_lower(colnames(pros_25))

# Fix relevant date/times
pros_25 <- pros_25 %>% mutate(score_time = as.POSIXct(score_time, tz = "UTC", format = "%Y-%m-%d %H:%M:%S"))
pros_25$score_time <- with_tz(time = pros_25$score_time, tzone = "America/Denver")

# Rename ICU admission time and fix the timezone
pros_25 <- pros_25 %>% mutate(picu_adm_date_time = as.POSIXct(icu_start_instant, format = "%Y-%m-%d %H:%M:%S"))
pros_25 <- pros_25 %>% mutate(picu_adm_date_time = force_tz(picu_adm_date_time, tzone = "America/Denver"))


# Filter for correct study start time
pros_val_start_time <- as.POSIXct("2023-09-01 00:00:01", tz = "America/Denver")
pros_25 <- pros_25 %>% filter(score_time >= pros_val_start_time)

nrow(pros_25) # 163790 predictions
n_distinct(pros_25$pat_enc_csn_id) #1,239 distinct PICU encounters



## Helper: assign greedy 5-minute timepoint clusters within each encounter
add_timepoint_clusters <- function(df, window_minutes = 5) {

  df <- df %>%
    dplyr::arrange(score_time)

  n <- nrow(df)

  if (n == 0) {
    df$timepoint_id <- integer(0)
    df$timepoint_anchor_time <- df$score_time
    return(df)
  }

  timepoint_id <- integer(n)
  timepoint_anchor_time <- df$score_time[rep(NA_integer_, n)]

  i <- 1
  k <- 0

  while (i <= n) {
    k <- k + 1
    anchor_time <- df$score_time[i]

    j <- i

    while (
      j <= n &&
      as.numeric(
        difftime(df$score_time[j], anchor_time, units = "mins")
      ) <= window_minutes
    ) {
      j <- j + 1
    }

    idx <- i:(j - 1)
    timepoint_id[idx] <- k
    timepoint_anchor_time[idx] <- anchor_time

    i <- j
  }

  df$timepoint_id <- timepoint_id
  df$timepoint_anchor_time <- timepoint_anchor_time

  df
}

pros_25_timepointed <- pros_25 %>%
  dplyr::filter(
    model_type %in% c("LR_no_abx", "LR_yes_abx", "RF_no_abx", "RF_yes_abx"),
    !is.na(pat_enc_csn_id),
    !is.na(score_time)
  ) %>%
  dplyr::arrange(pat_enc_csn_id, score_time, model_type) %>%
  dplyr::group_by(pat_enc_csn_id) %>%
  dplyr::group_modify(~ add_timepoint_clusters(.x, window_minutes = 5)) %>%
  dplyr::ungroup()

# Design collapse predictor function
collapse_predictor <- function(x) {

  if (is.character(x)) {
    x <- stringr::str_trim(x)
    x <- dplyr::na_if(x, "")
  }

  nonmissing_vals <- x[!is.na(x)]

  if (length(nonmissing_vals) == 0) {
    return(x[NA_integer_][1])
  }

  unique(nonmissing_vals)[1]
}

#Run the function
pros_25_collapsed_predictors <- pros_25_timepointed %>%
  dplyr::group_by(pat_enc_csn_id, timepoint_id, timepoint_anchor_time) %>%
  dplyr::summarise(
    n_model_rows_collapsed = dplyr::n(),
    n_models_collapsed = dplyr::n_distinct(model_type),
    first_score_time = min(score_time, na.rm = TRUE),
    last_score_time = max(score_time, na.rm = TRUE),
    score_time_span_minutes = as.numeric(
      difftime(last_score_time, first_score_time, units = "mins")
    ),
    models_included = paste(sort(unique(model_type)), collapse = ", "),
    dplyr::across(
      dplyr::all_of(predictors),
      collapse_predictor
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(pat_enc_csn_id, timepoint_anchor_time)
