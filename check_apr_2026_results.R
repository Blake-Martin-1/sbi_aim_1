#### Apply new nmodel to prospective results from Apr 2026 ####

apr_25 <- readr::read_csv(
  file = "/phi/sbi/sbi_blake/model_events_abxrf_parsed_4_9_26.csv",
  na = c("", "NULL"),
  col_types = readr::cols(.default = readr::col_character()),
  progress = FALSE,
  lazy = FALSE
)

## Drop columns ending in _RAW
apr_25_noraw <- apr_25[, !grepl("_RAW$", names(apr_25)), drop = FALSE]
apr_25_noraw <- apr_25_noraw[, !grepl("_TS_LEN$", names(apr_25_noraw)), drop = FALSE]
apr_25_noraw <- apr_25_noraw[, !grepl("_IS_COMPLETE$", names(apr_25_noraw)), drop = FALSE]

# Remove remainder of NA columns
all_null_cols <- names(apr_25_noraw)[
  sapply(apr_25_noraw, function(x) all(is.na(x)))
]

all_null_cols
length(all_null_cols)

## Create a mapping table so you can inspect old -> new names
name_map <- data.frame(
  old_name = names(apr_25_noraw),
  new_name = tolower(names(apr_25_noraw)),
  in_predictors = tolower(names(apr_25_noraw)) %in% predictors,
  stringsAsFactors = FALSE
)

## Apply the rename
apr_25_renamed <- apr_25_noraw
names(apr_25_renamed) <- name_map$new_name

## Optional: move predictor columns to the front, keep all others after them
apr_25_renamed <- apr_25_renamed %>%
  dplyr::select(dplyr::any_of(predictors), dplyr::everything())

## Optional checks
missing_predictors <- setdiff(predictors, names(apr_25_renamed))
extra_columns <- setdiff(names(apr_25_renamed), predictors)

## View results
name_map %>% dplyr::arrange(dplyr::desc(in_predictors), new_name)
missing_predictors
extra_columns

## Ok now process data to make sure format matches the dataframe used by the models ##
pros_25 <- apr_25_renamed %>% rename( , score_time = instant_utc_dttm, model_score = total_score, model_type = acuity_system_id)

pros_25$model_type[pros_25$model_type == "100133"] <- "LR_no_abx"
pros_25$model_type[pros_25$model_type == "100134"] <- "LR_yes_abx"
pros_25$model_type[pros_25$model_type == "100148"] <- "RF_no_abx"
pros_25$model_type[pros_25$model_type == "100149"] <- "RF_yes_abx"

# Lower case of column names
colnames(pros_25) <- str_to_lower(colnames(pros_25))