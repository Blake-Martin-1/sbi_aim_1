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

# Fix score_time column
apr_25_renamed <- apr_25_renamed %>%
  dplyr::mutate(
    picu_adm_date_time = lubridate::ymd_hms(
      icu_start_instant_str,
      tz = "America/Denver"
    )
  )

# Create score time via picu_adm_date_time and time since PICU admission
apr_25_renamed <- apr_25_renamed %>%
  dplyr::mutate(
    hours_since_icu_num = as.numeric(hours_since_icu),
    score_time = picu_adm_date_time + lubridate::dhours(hours_since_icu_num)
  )

## Ok now process data to make sure format matches the dataframe used by the models ##
pros_25 <- apr_25_renamed

# Lower case of column names
colnames(pros_25) <- str_to_lower(colnames(pros_25))

# Rearrange important cols to front
pros_25 <- pros_25 %>% relocate(csn, picu_adm_date_time, score_time, score_value)


