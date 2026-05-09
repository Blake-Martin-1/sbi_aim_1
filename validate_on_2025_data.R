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

# Add in the base excess present identifier
pros_25 <- pros_25 %>%
  dplyr::mutate(
    base_excess_present =
      !is.na(base_excess_arterial_ts) |
      !is.na(base_excess_capillary_ts) |
      !is.na(base_excess_venous_ts)
  )


## Helper: assign greedy 5-minute timepoint clusters within each encounter
add_timepoint_clusters <- function(df, window_minutes = 15) {

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
  dplyr::group_modify(~ add_timepoint_clusters(.x, window_minutes = 15)) %>%
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

## Ok now filter to just those predictions in the frst 24 hours
pros_25_collapsed_predictors_24h <- pros_25_collapsed_predictors %>%
  dplyr::group_by(pat_enc_csn_id) %>%
  dplyr::mutate(
    initial_timepoint_anchor_time = timepoint_anchor_time[timepoint_id == 1][1],
    hours_since_initial_timepoint = as.numeric(
      difftime(
        timepoint_anchor_time,
        initial_timepoint_anchor_time,
        units = "hours"
      )
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(hours_since_initial_timepoint <= 24)

# Filter to only those rows where 4 models were collapsed
pros_25_4_model <- pros_25_collapsed_predictors_24h %>% filter(n_models_collapsed == 4)

# Convert classes to match those needed by the model
## Confirm all required predictors exist before attempting conversion
missing_from_pros_25 <- setdiff(predictors, names(pros_25_4_model))
missing_from_test_df <- setdiff(predictors, names(test_df))

if (length(missing_from_pros_25) > 0) {
  stop(
    "These predictors are missing from pros_25_4_model: ",
    paste(missing_from_pros_25, collapse = ", ")
  )
}

if (length(missing_from_test_df) > 0) {
  stop(
    "These predictors are missing from test_df: ",
    paste(missing_from_test_df, collapse = ", ")
  )
}

## Helper function: coerce a new-data column to match the class of the corresponding test_df column
coerce_to_test_class <- function(x, template_col) {

  if (inherits(template_col, "POSIXct")) {
    tz_use <- attr(template_col, "tzone")
    if (is.null(tz_use) || length(tz_use) == 0 || is.na(tz_use)) {
      tz_use <- "UTC"
    }
    return(as.POSIXct(x, tz = tz_use))
  }

  if (is.factor(template_col)) {
    return(factor(as.character(x), levels = levels(template_col)))
  }

  if (is.integer(template_col)) {
    if (is.logical(x)) {
      return(as.integer(x))
    } else {
      return(as.integer(readr::parse_number(as.character(x))))
    }
  }

  if (is.numeric(template_col)) {
    if (is.logical(x)) {
      return(as.numeric(x))
    } else {
      return(readr::parse_number(as.character(x)))
    }
  }

  if (is.logical(template_col)) {
    x_chr <- tolower(trimws(as.character(x)))
    return(
      dplyr::case_when(
        x_chr %in% c("true", "t", "1", "yes", "y") ~ TRUE,
        x_chr %in% c("false", "f", "0", "no", "n") ~ FALSE,
        TRUE ~ NA
      )
    )
  }

  if (is.character(template_col)) {
    return(as.character(x))
  }

  ## Fallback: return unchanged if class is not covered above
  x
}

## Create prediction-ready dataframe with predictor classes matched to test_df
new_test_data_25 <- pros_25_4_model %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(predictors),
      ~ coerce_to_test_class(
        x = .x,
        template_col = test_df[[dplyr::cur_column()]]
      )
    )
  )

## Optional QA: verify classes now match test_df for all predictors
new_test_data_25_class_check <- tibble::tibble(
  predictor = predictors,
  test_df_class = vapply(
    predictors,
    function(x) paste(class(test_df[[x]]), collapse = "/"),
    character(1)
  ),
  new_test_data_25_class = vapply(
    predictors,
    function(x) paste(class(new_test_data_25[[x]]), collapse = "/"),
    character(1)
  ),
  class_matches = test_df_class == new_test_data_25_class
)

new_test_data_25_class_check %>%
  dplyr::filter(!class_matches)


rf_pred_prob_future <- predict(rf_tune, new_test_data_25 %>% dplyr::select(all_of(predictors)), type = "prob")[, "pos"]

##### QC Overlay SBI probabilities between the prior prospective test set and the future test set ####
library(dplyr)
library(tibble)
library(ggplot2)
library(scales)

pred_dist_df <- dplyr::bind_rows(
  tibble::tibble(
    predicted_probability = rf_pred_prob_test,
    dataset = "Test"
  ),
  tibble::tibble(
    predicted_probability = rf_pred_prob_future,
    dataset = "Future"
  )
)

p_pred_density <- ggplot(
  pred_dist_df,
  aes(
    x = predicted_probability,
    color = dataset,
    fill = dataset
  )
) +
  geom_density(
    alpha = 0.25,
    linewidth = 1.1,
    adjust = 1
  ) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title = "Distribution of SBI Model-Predicted Probabilities",
    x = "Predicted probability of SBI",
    y = "Density",
    color = "Dataset",
    fill = "Dataset"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "top"
  )

p_pred_density


## Join the test data from 2025 with the predicted probabilities
final_2025_df <- new_test_data_25 %>% bind_cols(data.frame("pred_sbi" = rf_pred_prob_future))

# Start loading in and cleaning extracted data (to get PICU admit time and SBI components)
adt_raw <- read_csv(file = "/phi/sbi/prospective_data/Prospective/data/adt_export_pros_100125.csv")

adt_df <- adt_raw
adt_df$intime <- as.POSIXct(adt_df$intime, format = "%m/%d/%y %H:%M")
adt_df$outtime <- as.POSIXct(adt_df$outtime, format = "%m/%d/%y %H:%M")

adt_df <- adt_df %>% mutate(study_id = paste0(pat_enc_csn_id, "_", cicu_event_count))
adt_df$intime <- force_tz(adt_df$intime, tzone =  "America/Denver")

first_adt <- adt_df %>% filter(cicu_event_count == 1)
first_adt <- first_adt %>% rename(picu_adm_date_time = intime)

### Limit adt dataframe to only those encounters in the test set
first_adt <- first_adt %>% mutate(pat_enc_csn_id = as.character(pat_enc_csn_id))
adt_future <- first_adt %>% filter(pat_enc_csn_id %in% final_2025_df$pat_enc_csn_id) # 1201 csn's remaining, of note final-2025_df has 1231 so some are missing

# Limit the prediction dataset to only those encounters that are in adt_future
future_interm <- final_2025_df %>% filter(pat_enc_csn_id %in% adt_future$pat_enc_csn_id) # now only same 1201 CSNs left, 25,265 rows

# join in admit info
future_w_admit <- future_interm %>% left_join(adt_future, by = "pat_enc_csn_id")

# Filter out fows where prediction was > 24 hours after pcu_adm_date_time, drops from 1201 to 1186 csn's

future_w_admit_24h <- future_w_admit %>%
  dplyr::mutate(
    hours_from_picu_adm_to_anchor = as.numeric(
      difftime(
        timepoint_anchor_time,
        picu_adm_date_time,
        units = "hours"
      )
    )
  ) %>%
  dplyr::filter(
    hours_from_picu_adm_to_anchor <= 24
  )

### Add in labs to determine lactate and culture negative sepsis outcome, will need to combine with blood culture data loaded from micro file

lab_raw <- read.csv(file = "/phi/sbi/prospective_data/Prospective/data/lab_export_pros_100125.csv")
lab_fut <- lab_raw
colnames(lab_fut) <- str_to_lower(colnames(lab_fut))

lact_df <- lab_fut %>% filter(component_name %in% c(
  "LACTATE (POCT)",
  "LACTATE WHOLE BLOOD",
  "LACTATE WB CONVERSION"
  ))

# Need to load in the demographics dataframe to be able to join by har instead of csn (or at least link the)
demo_raw <- read.csv(file = "/phi/sbi/prospective_data/Prospective/data/demog_export_pros_100125.csv")

just_ids <- demo_raw %>% dplyr::select(pat_mrn_id, hsp_account_id, pat_enc_csn_id) %>% distinct()
just_ids$pat_enc_csn_id <- as.character(just_ids$pat_enc_csn_id)
slim_ids <- just_ids %>% filter(pat_enc_csn_id %in% future_w_admit_24h$pat_enc_csn_id)

# Add on har
future_w_har <- future_w_admit_24h %>% left_join(slim_ids, by = "pat_enc_csn_id")

# Ok now will identify encounters with lactate > 2 within 24 hours of picu admission

## 1. Create an encounter-level table from future_w_har
##    This preserves one row per pat_enc_csn_id / PICU admission.
future_encounters <- future_w_har %>%
  dplyr::distinct(
    pat_enc_csn_id,
    hsp_account_id,
    picu_adm_date_time
  )

## 2. Clean lactate dataframe:
##    - Convert result_value to numeric
##    - Convert specimn_taken_time from character to POSIXct
##    - Assume times are America/Denver wall-clock times
lact_df_clean <- lact_df %>%
  dplyr::mutate(
    lactate_value = readr::parse_number(result_value),
    specimn_taken_time_dt = lubridate::ymd_hms(
      specimn_taken_time,
      tz = "America/Denver"
    )
  ) %>%
  dplyr::filter(
    !is.na(lactate_value),
    !is.na(specimn_taken_time_dt)
  )

## 3. For each encounter, determine whether any lactate >2 occurred
##    within +/- 24 hours of PICU admission
lact_over_two_by_encounter <- future_encounters %>%
  dplyr::left_join(
    lact_df_clean,
    by = "hsp_account_id"
  ) %>%
  dplyr::mutate(
    hours_from_picu_adm_to_lactate = as.numeric(
      difftime(
        specimn_taken_time_dt,
        picu_adm_date_time,
        units = "hours"
      )
    ),
    lactate_over_two_in_window =
      lactate_value > 2 &
      abs(hours_from_picu_adm_to_lactate) <= 24
  ) %>%
  dplyr::group_by(pat_enc_csn_id) %>%
  dplyr::summarise(
    lact_over_two = as.integer(any(lactate_over_two_in_window, na.rm = TRUE)),
    .groups = "drop"
  )

## 4. Join the encounter-level lact_over_two flag back onto future_w_har
future_w_har <- future_w_har %>%
  dplyr::left_join(
    lact_over_two_by_encounter,
    by = "pat_enc_csn_id"
  ) %>%
  dplyr::mutate(
    lact_over_two = dplyr::coalesce(lact_over_two, 0L)
  )

# Now load in VPS data and determine pneumonia on admission
vps_pna_raw_1 <- read_excel("/phi/sbi/sbi_blake/2025_Q1_pna.xlsx")
vps_pna_raw_2 <- read_excel("/phi/sbi/sbi_blake/2025_Q2_pna.xlsx")

vps_pna_raw_full <- bind_rows(vps_pna_raw_1, vps_pna_raw_2)

# Fix names and classes of columns, and then bind all vps df's together
vps_pna_fixed <- vps_pna_raw_full %>% dplyr::select(MRN, ICUAdmDateTime, AccountNum, `Present On Admission`)
colnames(vps_pna_fixed) <- c("mrn", "picu_adm_date_time", "hsp_account_id", "pna_on_admit")

# Fix all of the pna timezone info
pna_all <- vps_pna_fixed %>% mutate(picu_adm_date_time = force_tz(picu_adm_date_time, tzone = "America/Denver"))
pna_all$hsp_account_id <- as.integer(pna_all$hsp_account_id)

# Filter for pneumonias present on admission
pna_admit <- pna_all %>% filter(pna_on_admit == "Yes") %>% distinct() # from 577 pna rows to 420 pneumonias

# Add pneumonia label to current dataframe
future_w_pna <- future_w_har %>% mutate(pna_1_0 = ifelse(hsp_account_id %in% pna_admit$hsp_account_id, yes = 1, no = 0))

# Now call micro script to determine micro sbi presence
source(file = "/phi/sbi/sbi_blake/micro_tidying_2025.R")

## Identify study id's with a micro sbi
pros_all <- future_w_micro %>% mutate(micro_sbi_1_0 = ifelse(study_id %in% micro_df$study_id, yes = 1, no = 0))

# Clarify that cnss has to occur in setting of drawing a blood culture, drops cnss rate from 7% to 3.7%
pros_all <- pros_all %>% group_by(study_id) %>% mutate(cx_neg_sepsis = ifelse(min(sbp_min) < 5 & lact_over_two == 1, yes = 1, no = 0)) %>% ungroup()
pros_all <- pros_all %>% mutate(cnss_true = ifelse(cx_neg_sepsis == 1 & study_id %in% bcx_mrn_hosp_adm$study_id, yes = 1, no = 0))

# Now identify sbi overall with combo of pneumonia, cnss, and micro sbi's
pros_all <- pros_all %>% mutate(sbi_present = ifelse((pna_1_0 + cnss_true + micro_sbi_1_0) > 0, yes = 1, no = 0)) # 27.5% of all predictions have an SBI present in 24hr period around PICU admit

# Identify study_id level rate of SBI
# BY prediction
sbi_rate_by_prediction <- round(sum(pros_all$sbi_present) / nrow(pros_all) * 100, digits = 2)

# By study_id
n_with_sbi <- pros_all %>%
  group_by(study_id) %>%                 # one group per PICU stay
  summarise(has_sbi = any(sbi_present == 1), .groups = "drop") %>%
  filter(has_sbi) %>%                    # keep stays with ≥1 “1”
  nrow() # 309 with sbi out of 1186

sbi_rate_by_study_id <- round(n_with_sbi / n_distinct(pros_all$study_id) * 100, digits = 2) # 26% of all PICU admissions

# Now filter to appropriate time period based on available data
study_end <- as.POSIXct("2025-05-31 23:59:59", tz = "America/Denver")
study_start_time <- as.POSIXct("2025-01-01 00:00:01", tz = "America/Denver")

# Change anchor to score time
pros_all <- pros_all %>% rename(score_time = timepoint_anchor_time)

pros_full_data <- pros_all # store full dataset in case needed
pros_all <- pros_all %>% filter(score_time <= study_end) %>% filter (score_time >= study_start_time)# filter for end of study period

# Load in antibiotic data
source(file = "/phi/sbi/sbi_blake/abx_pros_2025.R")

