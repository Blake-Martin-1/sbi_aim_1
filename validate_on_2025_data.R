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

#Set seed for reproducibility
set.seed(2025)

# Generate predictions
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

# Filter out predictions that occur >24 hours after picu admission
pros_all<- pros_all %>% mutate(t_diff = difftime(score_time, picu_adm_date_time, units = "hours"))
pros_all <- pros_all %>% filter(t_diff <= 24) # n of prediction rows unchanged

# Recalculate number of encounters
n_enc <- n_distinct(pros_all$study_id) #1186 study encounters
n_pts <- nrow(pros_all %>% dplyr::select(pat_mrn_id) %>% distinct()) #1131 unique hospitalizations
n_pred_pre <- nrow(pros_all) #24,512 predictions

# Fix it so that sbi_present is 1 if cx_neg_sepsis is 1 at any point during the PICU encounter
# Step 1: Identify study_ids where cx_neg_sepsis is 1 in any row
study_ids_with_cx_neg_sepsis <- pros_all %>%
  group_by(study_id) %>%
  summarise(max_cx_neg_sepsis = max(cnss_true)) %>%
  filter(max_cx_neg_sepsis == 1) %>%
  pull(study_id)

# Create another column to show if cx_neg_sepsis was every present (instead of just becoming 1 once)
pros_all <- pros_all %>%
  mutate("ever_cx_neg_sepsis" = if_else(study_id %in% study_ids_with_cx_neg_sepsis, 1, 0))

pros_all <- pros_all %>% mutate(sbi_present = ifelse(sbi_present == 1 | ever_cx_neg_sepsis == 1, yes = 1, no = sbi_present))

# Load in demographic info to help match csn and har
demo_ref <- read_csv(file = "/phi/sbi/prospective_data/Prospective/data/demog_export_pros_100125.csv")
id_df <- demo_ref %>% dplyr::select(hsp_account_id, pat_enc_csn_id, pat_mrn_id) %>% distinct()
id_df$pat_enc_csn_id <- as.character(id_df$pat_enc_csn_id)

# Add in race, ethnicity, language, insurance data from demo_slim after adding hsp_account_id via id_df,
pros_all <- pros_all %>% left_join(demo_slim %>% dplyr::select(hsp_account_id, race, ethnicity, language, insurance_type) %>% distinct(), by = "hsp_account_id")

# Fix ethnicity and other SDOH categories
pros_all$ethnicity[pros_all$ethnicity == "Decline to Answer"] <- "Other_Or_Unknown"
pros_all$ethnicity[pros_all$ethnicity == "Not Reported"] <- "Other_Or_Unknown"
pros_all$ethnicity[pros_all$ethnicity == "Unknown"] <- "Other_Or_Unknown"
pros_all$ethnicity[is.na(pros_all$ethnicity)] <- "Other_Or_Unknown"

pros_all$race[is.na(pros_all$race)] <- "Other_Or_Unknown"
pros_all$race[pros_all$race == "Not Reported"] <- "Other_Or_Unknown"
pros_all$race[pros_all$race == "Decline to Answer"] <- "Other_Or_Unknown"
pros_all$race[pros_all$race == "Other"] <- "Other_Or_Unknown"
pros_all$race[pros_all$race == "Unknown"] <- "Other_Or_Unknown"

pros_all$language[!(pros_all$language %in% c("English", "Spanish"))] <- "Other_Or_Unknown"
pros_all$insurance_type[pros_all$insurance_type == "Other" | is.na(pros_all$insurance_type)] <- "Other_Or_Unknown"

# Establish slim dataframe with unique study id and other identifiers
adm_and_csn <- pros_all %>% dplyr::select(study_id, pat_enc_csn_id, hsp_account_id, pat_mrn_id, picu_adm_date_time) %>% distinct()

# Now calculate number of encounters with SBI and without SBI
n_sbi_with <- pros_all %>%
  group_by(study_id) %>%
  summarise(has_sbi = any(sbi_present == 1), .groups = "drop") %>%
  filter(has_sbi) %>%
  nrow()

n_sbi_without <- pros_all %>%
  group_by(study_id) %>%
  summarise(has_sbi = any(sbi_present == 1), .groups = "drop") %>%
  filter(!has_sbi) %>%
  nrow()


# QC: Output the results
n_sbi_with #309
n_sbi_without #877
n_sbi_with + n_sbi_without  # This should now equal 1186 which is the number of study_ids that are unique
n_enc <- nrow(pros_all %>% dplyr::select(study_id) %>% distinct()); n_enc # 1186 total encounters

# Ensure that if a patient has a pneumonia, micro, or cnss SBI during one prediction that that study_id has all rows with pna present, etc
pros_all <- pros_all %>%
  group_by(study_id) %>%
  mutate(pna_1_0 = if_else(any(pna_1_0 == 1), 1, 0)) %>%
  ungroup()


pros_all <- pros_all %>%
  group_by(study_id) %>%
  mutate(micro_sbi_1_0 = if_else(any(micro_sbi_1_0 == 1), 1, 0)) %>%
  ungroup()

pros_all <- pros_all %>%
  group_by(study_id) %>%
  mutate(ever_cx_neg_sepsis = if_else(any(ever_cx_neg_sepsis == 1), 1, 0)) %>%
  ungroup()

pros_all <- pros_all %>%
  group_by(study_id) %>%
  mutate(sbi_present = if_else(any(sbi_present == 1), 1, 0)) %>%
  ungroup()

# Ensure the culture negative sepsis flag is only positive if pneumonia and micro are == 0
pros_all$ever_cx_neg_sepsis[pros_all$pna_1_0 == 1] <- 0
pros_all$ever_cx_neg_sepsis[pros_all$micro_sbi_1_0 == 1] <- 0

# # Write / read fst of pros_all before modeling
# write.csv(x = pros_all, file = "/phi/sbi/sbi_blake/pros_2025_validation_pros_all_just_b4_modeling_5_8_26.csv")
# pros_all <- read.csv(file =  "/phi/sbi/sbi_blake/pros_2025_validation_pros_all_just_b4_modeling_5_8_26.csv")

# Evaluate predictions
colAUC(pros_all$pred_sbi, pros_all$sbi_present, plotROC = TRUE) # determine AUROC in test set = 0.83

# Rename df and column to enable cuse of prior code to evaluate future set performance
rf_future_df <- pros_all %>% rename(rf_prob = pred_sbi)

# Create similar SBI column
rf_future_df <- rf_future_df %>%
  dplyr::mutate(
    SBI = dplyr::case_when(
      sbi_present == 1 ~ "pos",
      sbi_present == 0 ~ "neg",
      TRUE ~ NA_character_
    ),
    SBI = factor(
      SBI,
      levels = levels(rf_test_df$SBI)
    )
  )

# Change study_id to factor to match rf_test_df structure
rf_future_df$study_id <- as.factor(rf_future_df$study_id)

# Create hours_since_picu_adm column
rf_future_df <- rf_future_df %>%
  dplyr::mutate(
    hours_since_picu_adm = as.numeric(
      round(
        as.numeric(
          difftime(
            score_time,
            picu_adm_date_time,
            units = "hours"
          )
        ),
        digits = 0
      )
    )
  )

# Filter out any predictions made before PICU admission
rf_future_df <- rf_future_df %>% filter(hours_since_picu_adm >= 0)




######### Future test set performance analyses #############
# These analyses mirror the rf_test_df evaluations from create_new_pros_model.R,
# but use rf_future_df from the 2025 prospective validation data.

# Preserve the original prospective test-set performance objects so combined plots
# and tables can compare the create_new_pros_model.R test set against this future
# test set.
rf_test_hourly_metrics_boot <- rf_hourly_metrics_boot
split_performance_table_original <- split_performance_table

# -----------------------------
# Calibration plot on future test set
# Predicted probability = Pr(SBI present)
# Observed outcome = SBI present
# -----------------------------
calib_future_df <- rf_future_df %>%
  dplyr::filter(
    !is.na(rf_prob),
    !is.na(sbi_present)
  ) %>%
  dplyr::mutate(
    sbi_present_num = to_sbi_numeric(sbi_present)
  )

calib_future_plot_df <- calib_future_df %>%
  dplyr::mutate(
    pred_bin = dplyr::ntile(rf_prob, 10)
  ) %>%
  dplyr::group_by(pred_bin) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_predicted_risk = mean(rf_prob, na.rm = TRUE),
    observed_sbi_rate = mean(sbi_present_num == 1, na.rm = TRUE),
    .groups = "drop"
  )

p_calibration_future <- ggplot2::ggplot(
  calib_future_plot_df,
  ggplot2::aes(x = mean_predicted_risk, y = observed_sbi_rate)
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.9
  ) +
  ggplot2::geom_line(
    color = "blue3",
    linewidth = 0.9
  ) +
  ggplot2::geom_point(
    color = "blue3",
    size = 3
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "Calibration Plot for New Prospective Random Forest Model",
    subtitle = "Future test set; points represent deciles of predicted SBI risk",
    x = "Mean predicted probability of SBI",
    y = "Observed SBI rate"
  ) +
  ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    axis.title = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank()
  )

p_calibration_future

# -----------------------------
# Future test set timepoint-level summary at the same 0.12 threshold
# -----------------------------
threshold <- 0.12
n_boot <- 1000
set.seed(123)

rf_future_hour_df <- rf_future_df %>%
  dplyr::filter(hours_since_picu_adm %in% 1:24)

overall_sbi_neg_timepoint_summary_future <- rf_future_hour_df %>%
  dplyr::filter(
    !is.na(sbi_present),
    !is.na(rf_prob)
  ) %>%
  dplyr::mutate(
    sbi_present_num = to_sbi_numeric(sbi_present),
    pred_sbi_negative = rf_prob <= threshold,
    actual_sbi_negative = sbi_present_num == 0,
    correctly_identified_sbi_negative = pred_sbi_negative & actual_sbi_negative
  ) %>%
  dplyr::summarise(
    n_patient_timepoint_pairs = dplyr::n(),
    n_actual_sbi_negative_timepoint_pairs = sum(actual_sbi_negative, na.rm = TRUE),
    n_predicted_sbi_negative_timepoint_pairs = sum(pred_sbi_negative, na.rm = TRUE),
    n_correctly_identified_sbi_negative_timepoint_pairs = sum(correctly_identified_sbi_negative, na.rm = TRUE),
    pct_actual_sbi_negative_timepoints_identified =
      n_correctly_identified_sbi_negative_timepoint_pairs / n_actual_sbi_negative_timepoint_pairs,
    npv_among_predicted_sbi_negative_timepoints =
      n_correctly_identified_sbi_negative_timepoint_pairs / n_predicted_sbi_negative_timepoint_pairs
  )

overall_sbi_neg_timepoint_summary_future # NPV 0.94 and 48% of SBI- pt-timepoint pairs identified

# Recreate overall_sbi_neg_timepoint_summary using rf_future_df.
overall_sbi_neg_timepoint_summary <- overall_sbi_neg_timepoint_summary_future

future_dup_check <- rf_future_hour_df %>%
  dplyr::count(study_id, hours_since_picu_adm, name = "n_rows") %>%
  dplyr::filter(n_rows > 1)

if (nrow(future_dup_check) > 0) {
  warning("Some future study_id + hours_since_picu_adm combinations had >1 row. Keeping the first row per patient-hour.")

  rf_future_hour_df <- rf_future_hour_df %>%
    dplyr::group_by(study_id, hours_since_picu_adm) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
}

# -----------------------------
# Split-level performance, keeping the original test-set row and adding the
# future test set row. AUROC/AUPRC use SBI-negative as the case/event class via
# calc_metrics_once() from create_new_pros_model.R.
# -----------------------------
future_split_performance_row <- calc_split_metrics_boot(
  rf_future_df,
  threshold = threshold,
  n_boot = n_boot
) %>%
  dplyr::mutate(
    split = "future test set",
    `AUROC (95% CI)` = format_metric_ci(auroc, auroc_lower, auroc_upper),
    `AUPRC (95% CI)` = format_metric_ci(auprc, auprc_lower, auprc_upper),
    `NPV (95% CI)` = format_metric_ci(npv, npv_lower, npv_upper),
    .before = 1
  ) %>%
  dplyr::select(split, `AUROC (95% CI)`, `AUPRC (95% CI)`, `NPV (95% CI)`)

split_performance_table <- dplyr::bind_rows(
  split_performance_table_original %>%
    dplyr::filter(split == "test set") %>%
    dplyr::select(split, `AUROC (95% CI)`, `AUPRC (95% CI)`, `NPV (95% CI)`),
  future_split_performance_row
)

split_performance_table_future <- split_performance_table
split_performance_table

# -----------------------------
# Future hourly metrics and SBI-negative prevalence
# -----------------------------
rf_future_hourly_metrics_boot <- purrr::map_dfr(
  1:24,
  function(hr) {
    dat_hr <- rf_future_hour_df %>%
      dplyr::filter(hours_since_picu_adm == hr)

    calc_hour_metrics_boot(
      dat = dat_hr,
      threshold = threshold,
      n_boot = n_boot
    ) %>%
      dplyr::mutate(hours_since_picu_adm = hr, .before = 1)
  }
) %>%
  dplyr::mutate(
    auroc = round(auroc, 2),
    auroc_lower = round(auroc_lower, 2),
    auroc_upper = round(auroc_upper, 2),
    auprc = round(auprc, 2),
    auprc_lower = round(auprc_lower, 2),
    auprc_upper = round(auprc_upper, 2),
    npv_label = dplyr::if_else(
      !is.na(n_negative_pred) & n_negative_pred > 0,
      paste0(tn, "/", n_negative_pred),
      NA_character_
    ),
    pct_sbi_neg_identified_label = dplyr::if_else(
      !is.na(pct_sbi_neg_identified),
      paste0(round(100 * pct_sbi_neg_identified), "%"),
      NA_character_
    ),
    auroc_label = dplyr::if_else(
      !is.na(auroc),
      sprintf("%.2f", auroc),
      NA_character_
    ),
    auprc_label = dplyr::if_else(
      !is.na(auprc),
      sprintf("%.2f", auprc),
      NA_character_
    )
  )

# Recreate rf_hourly_metrics_boot using rf_future_df, while the original test-set
# metrics remain available as rf_test_hourly_metrics_boot for combined plots.
rf_hourly_metrics_boot <- rf_future_hourly_metrics_boot
rf_hourly_metrics_boot

future_sbi_consistency_check <- rf_future_hour_df %>%
  dplyr::mutate(
    sbi_present_num = to_sbi_numeric(sbi_present)
  ) %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(
    n_outcome_values = dplyr::n_distinct(sbi_present_num[!is.na(sbi_present_num)]),
    .groups = "drop"
  )

if (any(future_sbi_consistency_check$n_outcome_values > 1)) {
  warning("Some future study_id values have more than one sbi_present value across hours.")
}

future_patient_level_outcome <- rf_future_hour_df %>%
  dplyr::mutate(
    sbi_present_num = to_sbi_numeric(sbi_present)
  ) %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(
    sbi_present_num = dplyr::first(sbi_present_num[!is.na(sbi_present_num)]),
    .groups = "drop"
  )

sbi_neg_prevalence_df_future <- future_patient_level_outcome %>%
  dplyr::summarise(
    n_encounters = dplyr::n(),
    n_sbi_negative = sum(sbi_present_num == 0, na.rm = TRUE),
    n_sbi_positive = sum(sbi_present_num == 1, na.rm = TRUE),
    sbi_neg_prevalence = n_sbi_negative / n_encounters,
    sbi_pos_prevalence = n_sbi_positive / n_encounters
  )

# Recreate sbi_neg_prevalence_df using rf_future_df.
sbi_neg_prevalence_df <- sbi_neg_prevalence_df_future
sbi_neg_prevalence_future <- sbi_neg_prevalence_df$sbi_neg_prevalence
sbi_neg_prevalence_label_future <- paste0(
  "Future SBI-negative prevalence = ",
  scales::percent(sbi_neg_prevalence_future, accuracy = 1)
)

sbi_neg_prevalence_df

# -----------------------------
# Plot helpers
# -----------------------------
plot_npv_by_hour <- function(plot_df, plot_title, show_legend = FALSE) {
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = hours_since_picu_adm, color = dataset, fill = dataset)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = npv_lower, ymax = npv_upper),
      alpha = 0.25,
      color = NA,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = npv),
      linewidth = 0.9,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = npv),
      size = 2.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = pct_sbi_neg_identified, linetype = dataset),
      linewidth = 0.9,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = pct_sbi_neg_identified),
      size = 2.3,
      shape = 17,
      na.rm = TRUE
    ) +
    ggplot2::scale_x_continuous(breaks = 1:24) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1.05),
      breaks = seq(0, 1, by = 0.1),
      expand = ggplot2::expansion(mult = c(0.01, 0.10)),
      name = "NPV",
      sec.axis = ggplot2::sec_axis(
        transform = ~ . * 100,
        name = "% of SBI-negative patients correctly identified"
      )
    ) +
    ggplot2::scale_color_manual(values = c("Test set" = "gray35", "Future test set" = "blue3")) +
    ggplot2::scale_fill_manual(values = c("Test set" = "gray70", "Future test set" = "lightblue")) +
    ggplot2::labs(
      title = plot_title,
      x = "Hours Since PICU Admission",
      color = "Dataset",
      fill = "Dataset",
      linetype = "Dataset"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.title.x = ggplot2::element_text(face = "bold"),
      axis.title.y.left = ggplot2::element_text(face = "bold"),
      axis.title.y.right = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = if (show_legend) "bottom" else "none"
    )
}

plot_discrimination_by_hour <- function(plot_df, metric, lower, upper, label, plot_title, y_axis_label, show_legend = FALSE) {
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = hours_since_picu_adm, y = .data[[metric]], color = dataset, fill = dataset)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[[lower]], ymax = .data[[upper]]),
      alpha = 0.25,
      color = NA,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      linewidth = 0.9,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      size = 2.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data[[label]]),
      vjust = -0.8,
      size = 3.5,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    ggplot2::scale_x_continuous(breaks = 1:24) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1.05),
      breaks = seq(0, 1, by = 0.1),
      expand = ggplot2::expansion(mult = c(0.01, 0.08))
    ) +
    ggplot2::scale_color_manual(values = c("Test set" = "gray35", "Future test set" = "blue3")) +
    ggplot2::scale_fill_manual(values = c("Test set" = "gray70", "Future test set" = "lightblue")) +
    ggplot2::labs(
      title = plot_title,
      x = "Hours Since PICU Admission",
      y = y_axis_label,
      color = "Dataset",
      fill = "Dataset"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = if (show_legend) "bottom" else "none"
    )
}

rf_future_plot_df <- rf_future_hourly_metrics_boot %>%
  dplyr::mutate(dataset = "Future test set")

rf_test_future_plot_df <- dplyr::bind_rows(
  rf_test_hourly_metrics_boot %>% dplyr::mutate(dataset = "Test set"),
  rf_future_plot_df
) %>%
  dplyr::mutate(
    dataset = factor(dataset, levels = c("Test set", "Future test set"))
  )

# NPV plots: future-only and test-set comparison
p_npv_future <- plot_npv_by_hour(
  rf_future_plot_df,
  "Negative Predictive Value by PICU Hour: Future Test Set"
)
p_npv_future

p_npv_test_vs_future <- plot_npv_by_hour(
  rf_test_future_plot_df,
  "Negative Predictive Value by PICU Hour: Test vs Future Test Set",
  show_legend = TRUE
)
p_npv_test_vs_future

# AUROC plots: future-only and test-set comparison. SBI-negative patients are the
# case/event class because calc_hour_metrics_boot() calls calc_metrics_once().
p_auroc_future <- plot_discrimination_by_hour(
  rf_future_plot_df,
  metric = "auroc",
  lower = "auroc_lower",
  upper = "auroc_upper",
  label = "auroc_label",
  plot_title = "AUROC by PICU Hour: Future Test Set",
  y_axis_label = "AUROC"
)
p_auroc_future

p_auroc_test_vs_future <- plot_discrimination_by_hour(
  rf_test_future_plot_df,
  metric = "auroc",
  lower = "auroc_lower",
  upper = "auroc_upper",
  label = "auroc_label",
  plot_title = "AUROC by PICU Hour: Test vs Future Test Set",
  y_axis_label = "AUROC",
  show_legend = TRUE
)
p_auroc_test_vs_future

# AUPRC plots: future-only and test-set comparison. SBI-negative patients are the
# case/event class, so prevalence reference lines use SBI-negative prevalence.
p_auprc_future <- plot_discrimination_by_hour(
  rf_future_plot_df,
  metric = "auprc",
  lower = "auprc_lower",
  upper = "auprc_upper",
  label = "auprc_label",
  plot_title = "AUPRC by PICU Hour: Future Test Set",
  y_axis_label = "AUPRC"
) +
  ggplot2::geom_hline(
    yintercept = sbi_neg_prevalence_future,
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +
  ggplot2::annotate(
    "text",
    x = 24,
    y = sbi_neg_prevalence_future,
    label = sbi_neg_prevalence_label_future,
    hjust = 1,
    vjust = -0.6,
    size = 3.8,
    color = "black",
    fontface = "bold"
  )
p_auprc_future

sbi_neg_prevalence_test <- patient_level_outcome %>%
  dplyr::summarise(
    sbi_neg_prevalence = mean(sbi_present_num == 0, na.rm = TRUE)
  ) %>%
  dplyr::pull(sbi_neg_prevalence)

p_auprc_test_vs_future <- plot_discrimination_by_hour(
  rf_test_future_plot_df,
  metric = "auprc",
  lower = "auprc_lower",
  upper = "auprc_upper",
  label = "auprc_label",
  plot_title = "AUPRC by PICU Hour: Test vs Future Test Set",
  y_axis_label = "AUPRC",
  show_legend = TRUE
) +
  ggplot2::geom_hline(
    yintercept = sbi_neg_prevalence_test,
    linetype = "dashed",
    color = "gray35",
    linewidth = 0.8
  ) +
  ggplot2::geom_hline(
    yintercept = sbi_neg_prevalence_future,
    linetype = "dashed",
    color = "blue3",
    linewidth = 0.8
  )
p_auprc_test_vs_future


###### Section to apply decision policy to new future data ######
######### Future decision-policy evaluation #############
# This section applies the single best policy selected in
# sbi_decision_policy_explore.R to the 2025 future test set. It assumes that the
# decision-policy script has already been run in the current R session so that
# best_policy and the policy helper functions are available.
required_policy_objects <- c(
  "best_policy",
  "prep_seq_data",
  "build_patient_list",
  "apply_policy_dataset",
  "summarise_policy"
)
missing_policy_objects <- required_policy_objects[
  !vapply(required_policy_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_policy_objects) > 0) {
  stop(
    "Run sbi_decision_policy_explore.R before this section. Missing objects: ",
    paste(missing_policy_objects, collapse = ", ")
  )
}

# Rebuild the complete patient-hour interval dataframe using the same helper
# used for rf_test_df in sbi_decision_policy_explore.R. Missing hours remain NA
# and therefore break multi-hour rules.
rf_future_seq <- prep_seq_data(rf_future_df)
future_patients <- build_patient_list(rf_future_seq)

future_decisions_best <- apply_policy_dataset(
  patient_list = future_patients,
  params = best_policy
)

future_best_summary <- summarise_policy(
  decisions_df = future_decisions_best,
  params = best_policy
) %>%
  dplyr::mutate(
    dataset = "future_test",
    meets_npv_target = !is.na(npv) & npv >= 0.95,
    coverage_score = dplyr::if_else(
      is.na(prop_true_neg_ruled_out),
      -Inf,
      prop_true_neg_ruled_out
    ),
    speed_score = dplyr::if_else(
      is.na(median_hour_ruleout_true_neg),
      Inf,
      median_hour_ruleout_true_neg
    )
  )

# Table analogous to best_policy, but evaluated on rf_future_df rather than the
# validation set used to identify the policy.
best_policy_future_performance <- future_best_summary %>%
  dplyr::select(
    policy_id,
    min_ruleout_hour,
    low_threshold,
    low_rule,
    low_k,
    low_m,
    low_n,
    high_threshold,
    high_rule,
    high_k,
    high_m,
    high_n,
    n_patients,
    n_true_neg,
    n_true_pos,
    n_ruled_out,
    n_ruled_out_true_neg,
    n_ruled_out_false_neg,
    npv,
    npv_lower,
    npv_upper,
    prop_all_patients_ruled_out,
    prop_true_neg_ruled_out,
    prop_true_pos_ruleout_error,
    median_hour_ruleout_all,
    median_hour_ruleout_true_neg,
    n_not_eligible,
    prop_true_pos_not_eligible,
    prop_true_neg_not_eligible,
    ppv_not_eligible,
    prop_indeterminate,
    prop_true_neg_indeterminate,
    prop_true_pos_indeterminate,
    meets_npv_target,
    coverage_score,
    speed_score,
    dataset
  )

best_policy_future_performance

future_trajectories_with_decision <- rf_future_seq %>%
  dplyr::left_join(
    future_decisions_best,
    by = c("study_id", "sbi_present")
  )

## -----------------------------
## Prep labels for final states
## -----------------------------
future_plot_decisions <- future_decisions_best %>%
  dplyr::mutate(
    truth_label = dplyr::case_when(
      sbi_present == 0 ~ "SBI-negative",
      sbi_present == 1 ~ "SBI-positive",
      TRUE ~ NA_character_
    ),
    final_state_label = dplyr::case_when(
      final_state == "ruled_out" ~ "Ruled out",
      final_state == "not_eligible" ~ "High-risk / not eligible",
      final_state == "indeterminate" ~ "Indeterminate",
      TRUE ~ final_state
    ),
    truth_label = factor(truth_label, levels = c("SBI-negative", "SBI-positive")),
    final_state_label = factor(
      final_state_label,
      levels = c("Ruled out", "High-risk / not eligible", "Indeterminate")
    )
  )

## -----------------------------
## Panel A: disposition by truth
## -----------------------------
future_panel_a_df <- future_plot_decisions %>%
  dplyr::count(truth_label, final_state_label, name = "n") %>%
  dplyr::group_by(truth_label) %>%
  dplyr::mutate(
    prop = n / sum(n),
    label = paste0(scales::percent(prop, accuracy = 1), "\n(n=", n, ")")
  ) %>%
  dplyr::ungroup()

p_policy_future_a <- ggplot2::ggplot(
  future_panel_a_df,
  ggplot2::aes(x = truth_label, y = prop, fill = final_state_label)
) +
  ggplot2::geom_col(width = 0.7, color = "white", linewidth = 0.6) +
  ggplot2::geom_text(
    ggplot2::aes(label = label),
    position = ggplot2::position_stack(vjust = 0.5),
    size = 4.2,
    lineheight = 0.95,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Ruled out" = "forestgreen",
      "High-risk / not eligible" = "orange",
      "Indeterminate" = "grey70"
    )
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Patients within true SBI group",
    fill = NULL,
    title = "Future test set policy disposition by true SBI status"
  ) +
  ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 15),
    axis.title = ggplot2::element_text(face = "bold", size = 14),
    axis.text = ggplot2::element_text(size = 12),
    legend.text = ggplot2::element_text(size = 12),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

## ----------------------------------------
## Panel B: cumulative rule-out over time
## ----------------------------------------
future_cum_df <- future_plot_decisions %>%
  dplyr::filter(!is.na(truth_label)) %>%
  dplyr::select(study_id, truth_label, final_state, decision_hour)

future_hours_df <- tidyr::expand_grid(
  study_id = unique(future_cum_df$study_id),
  hour = 0:24
) %>%
  dplyr::left_join(
    future_cum_df,
    by = "study_id"
  ) %>%
  dplyr::mutate(
    cum_ruled_out = dplyr::case_when(
      final_state == "ruled_out" & !is.na(decision_hour) & decision_hour <= hour ~ 1,
      TRUE ~ 0
    )
  )

future_panel_b_df <- future_hours_df %>%
  dplyr::group_by(truth_label, hour) %>%
  dplyr::summarise(
    prop_cum_ruled_out = mean(cum_ruled_out, na.rm = TRUE),
    .groups = "drop"
  )

p_policy_future_b <- ggplot2::ggplot(
  future_panel_b_df,
  ggplot2::aes(x = hour, y = prop_cum_ruled_out)
) +
  ggplot2::geom_area(
    data = future_panel_b_df %>% dplyr::filter(truth_label == "SBI-negative"),
    fill = "blue",
    alpha = 0.18
  ) +
  ggplot2::geom_area(
    data = future_panel_b_df %>% dplyr::filter(truth_label == "SBI-positive"),
    fill = "red",
    alpha = 0.22
  ) +
  ggplot2::geom_line(
    ggplot2::aes(color = truth_label),
    linewidth = 1.2
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = truth_label),
    size = 2.2
  ) +
  ggplot2::scale_color_manual(
    values = c("SBI-negative" = "blue", "SBI-positive" = "red")
  ) +
  ggplot2::scale_x_continuous(breaks = 0:24) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    x = "Hours since PICU admission",
    y = "Cumulative proportion ruled out",
    color = NULL,
    title = "Future test set cumulative rule-out over time"
  ) +
  ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 15),
    axis.title = ggplot2::element_text(face = "bold", size = 14),
    axis.text = ggplot2::element_text(size = 11),
    legend.text = ggplot2::element_text(size = 12),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

policy_future_figure <- p_policy_future_a + p_policy_future_b +
  patchwork::plot_layout(widths = c(1, 1.15))

policy_future_figure


