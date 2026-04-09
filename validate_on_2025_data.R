###### Code to test the model developed on prospective 2023-2024 data on the first 4-5 months of 2025 ###

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
pros_25 <- epic_5_25_renamed %>% rename(score_time = instant_utc_dttm, model_score = total_score, model_type = acuity_system_id)

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
pos_val_start_time <- as.POSIXct("2023-09-01 00:00:01", tz = "America/Denver")
pros_25 <- pros_25 %>% filter(score_time >= study_start_time)



