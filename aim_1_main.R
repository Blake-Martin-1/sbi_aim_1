# Script to evaluate changes over time in input availability, input values, and patient mix. To be used for prospective validation paper looking at reasons why the model had poor
# performance when evaluated prospectively.


# Call setup script to import needed filepaths
setwd(dir = "/phi/sbi/sbi_blake/aim_1_paper_materials/")

source("setup_aim_1.R")

# Load additional scripts
library(scales)
library(gbm)
library(lubridate)
library(data.table)
library(pedbp)
library(DescTools)
library(stringi)
library(stringr)
library(readxl)
library(devtools)
library(rriskDistributions)
library(mltools)
library(tidyr)
library(rsample)
library(parsnip)
library(ranger)
library(pROC)      # ROC + AUC
library(caret)
library(iml)
library(purrr)
library(glmnet)
library(tidyverse)
library(survey)
library(forcats)
library(gtsummary)
library(flextable)
library(officer)
library(ggplot2)
library(dplyr)
library(tibble)
library(yardstick)
library(patchwork)
library(PRROC)     # PR curve + AUPRC
library(WeightIt)
library(cobalt)
library(emmeans)
library(ggrepel)
library(WeightedROC)
library(gtsummary)
library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)

# Run script to generate models and rf_df and rf_df_abx
source(file = "retro_models_to_2026.R")

# Load in similar dataframe for the prospective model data
prospective_model_data_file_path <- getOption(
  "prospective_model_data_file_path",
  "/phi/sbi/sbi_blake/pros_all_just_b4_modeling_1_15_26_all_models.csv"
)
pros_all <- read_csv(file = prospective_model_data_file_path)

# Filter to just the antibiotic unexposed random forest model
pros_rf_ue <- pros_all %>% filter(model_type == "RF_no_abx")

# Filter to just the 2hr mark
pros_2hr_rf_ue <- pros_rf_ue %>% filter(hours_since_picu_adm >= 2) %>% filter(hours_since_picu_adm < 4)

pros_2hr_rf_ue <- pros_2hr_rf_ue %>%
  group_by(study_id) %>%
  slice_min(order_by = hours_since_picu_adm, n = 1, with_ties = FALSE) %>%
  ungroup()


#Now select only the same input and output columns as present in the retrospective dataset, make sure names and format match, scales match, etc
thresh_rf_no_abx <- 0.05

# Fix name of case_id
rf_df <- rf_df %>% rename(study_id = case_id)
retro_rf_same <- rf_df

pros_rf_same <- pros_2hr_rf_ue %>%
  dplyr::select(all_of(intersect(colnames(pros_2hr_rf_ue),
                                 colnames(retro_rf_same))))

retro_rf_same <- retro_rf_same %>%
  dplyr::select(all_of(sort(colnames(retro_rf_same))))


retro_rf_same <- retro_rf_same %>% relocate(sbi_present, model_score)
pros_rf_same <- pros_rf_same %>% relocate(sbi_present, model_score)

retro_rf_same$model_score <- retro_rf_same$model_score * 100
retro_rf_same <- retro_rf_same %>% mutate(sbi_present = ifelse(sbi_present == "yes", yes = 1, no = 0))

#Add epoch column
retro_rf_same$epoch <- "retro"
pros_rf_same$epoch <- "prospective"

# In addition we'll need to know how many labs, fio2 values etc were available vs. imputed and if available how many measurements were present
# will start with the retrospective dataset

# First need to get back the case_id values (whoops)
# load model_data_df
full_retro_df <- read.fst(path = "/home/martinbl/sbi_blake/pre_micro_vps_vitals_11_21_23.fst")

# Change names of columns I'll match on
full_retro_name_fix <- full_retro_df %>% rename(age = age_years, hr_min = min_hr)

# str(full_retro_name_fix %>% dplyr::select(age, los_before_icu_days, hr_min, min_sbp, min_dbp))
# str(retro_rf_same %>% dplyr::select(age, los_before_icu_days, hr_min, sbp_min, dbp_min))

#
retro_rf_with_case <- retro_rf_same

names(retro_rf_with_case) |> grep("hr", x = _, value = TRUE)



retro_rf_with_case <- retro_rf_with_case %>% relocate(study_id)

# Read pt_sex into df pt_sex, then correct column classes,
df_demographics <- read_csv(demographics_file_path)
pt_sex <- df_demographics %>%
  dplyr::transmute(
    mrn = as.character(MRN),
    sex = as.factor(SEX),
    ethnicity = as.factor(ETHNICITY)
  )

encounter_df <- read_csv(file =  adm_dx_file_path)
colnames(encounter_df) <- c("mrn", "hosp_adm_date_time", "hosp_dc_date_time", "picu_adm_date_time", "admission_dx", "billing_dx", "proc_dx")
encounter_df <- encounter_df %>% mutate(mrn = as.character(encounter_df$mrn),
                                        hosp_adm_date_time = as.POSIXct(encounter_df$hosp_adm_date_time, tz = "GMT",
                                                                        format = c("%m/%d/%Y %H:%M")),
                                        hosp_dc_date_time = as.POSIXct(encounter_df$hosp_dc_date_time, tz = "GMT", format = c("%m/%d/%Y %H:%M")),
                                        picu_adm_date_time = as.POSIXct(encounter_df$picu_adm_date_time, tz = "GMT", format = c("%m/%d/%Y %H:%M")))

encounter_df <- encounter_df %>% mutate(hosp_adm_date_time = force_tz(hosp_adm_date_time, tz = "America/Denver"),
                                        hosp_dc_date_time = force_tz(hosp_dc_date_time, tz = "America/Denver"),
                                        picu_adm_date_time = force_tz(picu_adm_date_time, tz = "America/Denver"))


vitals_full_file_path <- getOption("vitals_full_file_path", vitals_full_file_path)
vitals_full <- read_csv(file = vitals_full_file_path)

# # #Correct class type for the vitals_full df
vitals_full <- vitals_full %>% mutate(mrn = as.character(vitals_full$MRN),
                                      result_value = as.character(vitals_full$RESULT_VALUE),
                                      hosp_admit_time = as.POSIXct(vitals_full$HOSP_ADMSN_TIME, tz = "GMT"),
                                      entry_time = as.POSIXct(vitals_full$ENTRY_TIME, tz = "GMT"),
                                      vital_time_stamp = as.POSIXct(vitals_full$RECORDED_TIME, tz = "GMT"),
                                      entry_time = as.POSIXct(vitals_full$ENTRY_TIME, tz = "GMT"),
                                      vital_type = as.character(vitals_full$VITAL_TYPE)) %>%
  dplyr::select(-c(PREVIOUS_VALUE, PREVIOUS_ENTRY_TIME, MRN, HOSP_ADMSN_TIME, ENTRY_TIME, RECORDED_TIME, ENTRY_TIME, VITAL_TYPE))

vitals_full <- vitals_full %>% mutate(hosp_admit_time = force_tz(hosp_admit_time, tz = "America/Denver"), entry_time = force_tz(entry_time, tz = "America/Denver"),
                                      vital_time_stamp = force_tz(vital_time_stamp, tz = "America/Denver"), entry_time = force_tz(entry_time, tz = "America/Denver"))

# Set vitals_for_death equal to df_demographics for later use in determining mortality
vitals_for_death <- encounter_df

# Filter out caseID's from the VPS patient list that don't have an MRN match in the vital sign data frame, and filter out mrn's in the vitals df not present in vps list
mrns_epic <- vitals_full %>% dplyr::select(mrn) %>% distinct()
mrns_vps <- vps_pt_list_full %>% dplyr::select(mrn) %>% distinct()

vps_pt_list <- vps_pt_list_full %>% filter(mrn %in% mrns_epic$mrn)

vps_pt_list <- vps_pt_list %>%
  mutate(
    picu_adm_date_time = force_tz(picu_adm_date_time, tz = "America/Denver")
  )

vitals_full <- vitals_full %>% filter(mrn %in% mrns_vps$mrn)

# Rename hosp admit time column
vitals_full <- vitals_full %>% mutate(hosp_adm_date_time = vitals_full$hosp_admit_time) %>% dplyr::select(-c(hosp_admit_time))

# Now will join the encounter info to the vitals_full info such that for each vital sign there is an associated ICU admission time (will aid with matching to VPS)
vitals_full <- left_join(vitals_full, encounter_df %>% dplyr::select(mrn, hosp_adm_date_time, hosp_dc_date_time, picu_adm_date_time), by = c("mrn", "hosp_adm_date_time"))

# Filter out vital sign data from before the VPS database in 2011 starts
first_vps_admit_time <- as.POSIXct("2011-01-01 00:01:00", tz = "America/Denver")
vitals_full <- vitals_full %>% filter(picu_adm_date_time >= first_vps_admit_time)

# Need to identify the combo of mrn and picu-adm_date_time in the vitals_full df that matches up (i.e. diff < 1 hour) with the
#unique mrn and hosp_admit pairs in the epic set
uniq_epic_mrn_hosp_admit <- vitals_full %>% dplyr::select(mrn, hosp_adm_date_time, picu_adm_date_time) %>% distinct()
uniq_vps_mrn_picu_adm_case_id <- vps_pt_list %>% dplyr::select(case_id, mrn, picu_adm_date_time) %>% distinct()

uniq_epic_mrn_hosp_admit <- arrange(.data = uniq_epic_mrn_hosp_admit, picu_adm_date_time)
uniq_vps_mrn_picu_adm_case_id <- uniq_vps_mrn_picu_adm_case_id %>%
  mutate(
    picu_adm_date_time = force_tz(picu_adm_date_time, tz = "America/Denver")
  )
uniq_epic_mrn_hosp_admit <- uniq_epic_mrn_hosp_admit %>% mutate(case_id = NA_character_)

for(i in c(1:nrow(uniq_epic_mrn_hosp_admit))) {
  temp_vps <- uniq_vps_mrn_picu_adm_case_id %>% filter(mrn == uniq_epic_mrn_hosp_admit$mrn[i])
  right_case_id <- temp_vps %>% filter(abs(difftime(picu_adm_date_time, uniq_epic_mrn_hosp_admit$picu_adm_date_time[i], tz = "America/Denver", units = "min")) < 240)
  if(nrow(right_case_id) > 1){print(">1 case_id")}
  uniq_epic_mrn_hosp_admit$case_id[i] <- right_case_id$case_id[1]
}

uniq_epic_mrn_hosp_admit <- uniq_epic_mrn_hosp_admit %>% filter(!is.na(case_id))

vitals_full <- vitals_full %>% left_join(uniq_epic_mrn_hosp_admit, by = c("mrn", "picu_adm_date_time"))

#Filter out vitals documented >2 hours after ICU admit time
new_vitals_df <- vitals_full
new_vitals_df <- new_vitals_df %>% filter(entry_time <= (picu_adm_date_time + 2*60*60))


# Remove extra hosp_adm_time_date column, rename RESULT_VALUE column
new_vitals_df <- new_vitals_df %>% dplyr::select(-c("hosp_adm_date_time.y")) %>% mutate(hosp_adm_date_time = new_vitals_df$hosp_adm_date_time.x,
                                                                                        vital_value = new_vitals_df$RESULT_VALUE) %>% dplyr::select(-c("hosp_adm_date_time.x", "RESULT_VALUE"))

# new_vitals_df <- new_vitals_df %>% dplyr::select(-case_id.y) %>% rename(case_id = "case_id.x")

# Initialize columns for BP measurements
new_vitals_df$SBP <- 0
new_vitals_df$DBP <- 0


# Create SBP from bp string
# First remove bp vitals from new_vitals_df
bp_df <- new_vitals_df %>% filter(vital_type == "bp")
new_vitals_df <- new_vitals_df %>% filter(vital_type != "bp")

# Now assign SBP and DBP values
bp_df <- bp_df %>% mutate(SBP = as.integer(sub("\\/.*", "", result_value)))
bp_df <- bp_df %>% mutate(DBP = as.integer(sub(".*\\/", "", result_value)))

# Re-add bp_df to new_vitals_df
new_vitals_df <- bind_rows(new_vitals_df, bp_df)

# Calculate max, min, mean, median, and number of observations for a given vital sign for each case_id via a
# data_frame for each vital sign type. Note that the first of these vital sign df's, hr_stats, includes ICU arrival time, LOS before ICU simply because it's the leftmost
# dataframe when they get joined. Weight and height were originally given in ounces and inches so these are converted to kg and cm

hr_df <- filter(new_vitals_df, vital_type == "hr"); hr_df$vital_value <- as.integer(hr_df$vital_value)
fio2_df <- filter(new_vitals_df, vital_type == "fio2"); fio2_df$vital_value <- as.integer(fio2_df$vital_value)
spo2_df <- filter(new_vitals_df, vital_type == "spo2"); spo2_df$vital_value <- as.integer(spo2_df$vital_value)
bp_df <- filter(new_vitals_df, vital_type == "bp")
temp_df <- filter(new_vitals_df, vital_type == "temp"); temp_df$vital_value <- as.numeric(temp_df$vital_value)
rr_df <- filter(new_vitals_df, vital_type == "RR"); rr_df$vital_value <- as.numeric(rr_df$vital_value)
weight_df <- filter(new_vitals_df, vital_type == "weight"); weight_df$vital_value <- round(as.numeric(weight_df$vital_value) * 28.3495 / 1000, digits = 2)
height_df <- filter(new_vitals_df, vital_type == "height"); height_df$vital_value <- round(as.numeric(height_df$vital_value) * 2.54, digits = 1)

########### Now will filter out obviously incorrect vital sign entries using clinical guidance
hr_df <- hr_df %>% filter(vital_value > 10 & vital_value < 300)
fio2_df <- fio2_df %>% filter(vital_value >=20 & vital_value <=100)
spo2_df <- spo2_df %>% filter(vital_value >=0 & vital_value <=100)
bp_df <- bp_df %>% filter(SBP > 20 & SBP < 300) %>% filter(DBP > 10 & DBP < 250)
temp_df <- temp_df %>% filter(vital_value >70 & vital_value < 110)
rr_df <- rr_df %>% filter(vital_value > 0 & vital_value < 130)
weight_df <- weight_df %>% filter(vital_value >=1 & vital_value < 300)
height_df <- height_df %>% filter(vital_value < 210 & vital_value >=40)

# Now get number of vital sign values for each case_id
hr_counts <- hr_df %>%
  filter(!is.na(case_id)) %>%
  count(case_id, name = "n_hr_rows")

fio2_counts <- fio2_df %>%
  filter(!is.na(case_id)) %>%
  count(case_id, name = "n_fio2_rows")

o2sat_counts <- spo2_df %>%
  filter(!is.na(case_id)) %>%
  count(case_id, name = "n_o2sat_rows")

bp_counts <- bp_df %>%
  filter(!is.na(case_id)) %>%
  count(case_id, name = "n_bp_rows")

temp_counts <- temp_df %>%
  filter(!is.na(case_id)) %>%
  count(case_id, name = "n_temp_rows")

rr_counts <- rr_df %>%
  filter(!is.na(case_id)) %>%
  count(case_id, name = "n_rr_rows")

# Back to creating the retrospective dataset
retro_rf_with_case <- retro_rf_with_case %>% rename(case_id = study_id)
retro_rf_vs <- retro_rf_with_case %>% left_join(hr_counts, by = "case_id") %>% left_join(fio2_counts, by = "case_id") %>% left_join(o2sat_counts, by = "case_id")%>% left_join(bp_counts, by = "case_id") %>%
                left_join(temp_counts, by = "case_id") %>% left_join(rr_counts, by = "case_id")

# Create fio2 flag using fio2_counts column
retro_rf_vs <- retro_rf_vs %>% mutate(fio2_present = ifelse(is.na(n_fio2_rows), yes = 0, no = 1))
retro_rf_vs$n_fio2_rows[is.na(retro_rf_vs$n_fio2_rows)] <- 0

# Now add in lab data
retro_labs <- read_csv(file = "/home/martinbl/sbi_blake/labs_for_comp.csv")
# retro_fio2 <- read_csv(file = "/home/martinbl/sbi_blake/fio2_data.csv")

# Now nead to add in the number of pco2 case_id values that have values and how many
co2_df <- retro_labs %>% filter(lab_time_taken <= picu_adm_date_time + 60 * 60 * 2) %>%
  group_by(case_id) %>%
  summarise(
    pco2_gas_blood_present = as.integer(any(lab_type == "pco2_gas_blood", na.rm = TRUE)),
    n_pco2_gas_blood       = sum(lab_type == "pco2_gas_blood", na.rm = TRUE),
    .groups = "drop"
  )

# Fix case_id type
co2_df$case_id <- as.character(co2_df$case_id)

# Join vs dataset with lab dataset and then fill in 0's for patients with no labs
retro_full <- left_join(retro_rf_vs, co2_df, by = "case_id")
retro_full$pco2_gas_blood_present[is.na(retro_full$pco2_gas_blood_present)] <- 0
retro_full$n_pco2_gas_blood[is.na(retro_full$n_pco2_gas_blood)] <- 0


# Ensure we have correct variables for abx exposed version of RF model

#Fix classes

rf_df_abx$sbi_present <- as.factor(rf_df_abx$sbi_present)


### Now do same changes to be able to have full retro abx exposed model variables, flags, and counts

# Filter to just the antibiotic unexposed random forest model
pros_rf_ae <- pros_all %>% filter(model_type == "RF_yes_abx")

# Filter to just the 2hr mark
pros_2hr_rf_ae <- pros_rf_ae %>% filter(hours_since_picu_adm >= 2) %>% filter(hours_since_picu_adm < 4)

pros_2hr_rf_ae <- pros_2hr_rf_ae %>%
  group_by(study_id) %>%
  slice_min(order_by = hours_since_picu_adm, n = 1, with_ties = FALSE) %>%
  ungroup()


#Now select only the same input and output columns as present in the retrospective dataset, make sure names and format match, scales match, etc
thresh_rf_yes_abx <- 0.1 # need to ensure this is correct

retro_rf_same_abx <- rf_df_abx

pros_rf_same_abx <- pros_2hr_rf_ae %>%
  dplyr::select(all_of(intersect(colnames(pros_2hr_rf_ae),
                                 colnames(retro_rf_same_abx))))

retro_rf_same_abx <- retro_rf_same_abx %>%
  dplyr::select(all_of(sort(colnames(retro_rf_same_abx))))


retro_rf_same_abx <- retro_rf_same_abx %>% relocate(sbi_present, model_score)
pros_rf_same_abx <- pros_rf_same_abx %>% relocate(sbi_present, model_score)

#Add epoch column
retro_rf_same_abx$epoch <- "retro"
pros_rf_same_abx$epoch <- "prospective"

# In addition we'll need to know how many labs, fio2 values etc were available vs. imputed and if available how many measurements were present
# will start with the retrospective dataset

# First need to get back the case_id values (whoops)
# load model_data_df
full_retro_df <- read.fst(path = "/home/martinbl/sbi_blake/pre_micro_vps_vitals_11_21_23.fst")

# Change names of columns I'll match on
full_retro_name_fix <- full_retro_df %>% rename(age = age_years, hr_min = min_hr)

# str(full_retro_name_fix %>% dplyr::select(age, los_before_icu_days, hr_min, min_sbp, min_dbp))
# str(retro_rf_same %>% dplyr::select(age, los_before_icu_days, hr_min, sbp_min, dbp_min))



retro_rf_with_case_abx <- retro_rf_same_abx

# Move case_id to front remove key columnns
retro_rf_with_case_abx <- retro_rf_with_case_abx %>% relocate(case_id)

# Join with vs counts
retro_rf_vs_abx <- retro_rf_with_case_abx %>% left_join(hr_counts, by = "case_id") %>% left_join(fio2_counts, by = "case_id") %>% left_join(o2sat_counts, by = "case_id")%>% left_join(bp_counts, by = "case_id") %>%
  left_join(temp_counts, by = "case_id") %>% left_join(rr_counts, by = "case_id")

# Create fio2 flag using fio2_counts column
retro_rf_vs_abx <- retro_rf_vs_abx %>% mutate(fio2_present = ifelse(is.na(n_fio2_rows), yes = 0, no = 1))
retro_rf_vs_abx$n_fio2_rows[is.na(retro_rf_vs_abx$n_fio2_rows)] <- 0

# Now add in lab data
retro_labs <- read_csv(file = "/home/martinbl/sbi_blake/labs_for_comp.csv")
# retro_fio2 <- read_csv(file = "/home/martinbl/sbi_blake/fio2_data.csv")

# Now nead to add in the number of lactate case_id values that have values and how many
lactate_df <- retro_labs %>% filter(lab_time_taken <= picu_adm_date_time + 2 * 3600) %>%
  group_by(case_id) %>%
  summarise(
    lactate_present = as.integer(any(lab_type == "lactate_blood", na.rm = TRUE)),
    n_lactate       = sum(lab_type == "lactate_blood", na.rm = TRUE),
    .groups = "drop"
  )

# Fix case_id type
lactate_df$case_id <- as.character(lactate_df$case_id)

# Join vs dataset with lab dataset and then fill in 0's for patients with no labs
retro_full_abx <- left_join(retro_rf_vs_abx, lactate_df, by = "case_id")
retro_full_abx$lactate_present[is.na(retro_full_abx$lactate_present)] <- 0
retro_full_abx$n_lactate[is.na(retro_full_abx$n_lactate)] <- 0

# Now nead to add in the number of lactate case_id values that have values and how many
hematocrit_df <- retro_labs %>% filter(lab_time_taken <= picu_adm_date_time + 2 * 3600) %>%
  group_by(case_id) %>%
  summarise(
    hematocrit_present = as.integer(any(lab_type == "hematocrit_blood", na.rm = TRUE)),
    n_hematocrit       = sum(lab_type == "hematocrit_blood", na.rm = TRUE),
    .groups = "drop"
  )

# Fix case_id type
hematocrit_df$case_id <- as.character(hematocrit_df$case_id)

# Join vs dataset with lab dataset and then fill in 0's for patients with no labs
retro_full_abx <- left_join(retro_rf_vs_abx, hematocrit_df, by = "case_id")
retro_full_abx$hematocrit_present[is.na(retro_full_abx$hematocrit_present)] <- 0
retro_full_abx$n_hematocrit[is.na(retro_full_abx$n_hematocrit)] <- 0


### Now need to attach all of the vital sign information and lab information to the prospective data ###
# Read in prospective datasets from Epic
# Get csn's and time of score to quickly filter available labs and vitals

csn_time <- pros_2hr_rf_ue %>% dplyr::select(study_id, pat_enc_csn_id, hsp_account_id, score_time, picu_adm_date_time) %>% distinct()

# Load needed pros data
vitals_pros <- read_csv("/phi/sbi/prospective_data/Prospective/data/vitals_export_pros_100125.csv")


# Copy to avoid modifying originals by reference
csn_dt <- as.data.table(csn_time)
vit_dt <- as.data.table(vitals_pros)

# Drop vitals with missing type (optional but recommended)
vit_dt <- vit_dt[!is.na(vital_type)]

# Define window per study_id
csn_dt[, window_start := picu_adm_date_time - as.difftime(24, units = "hours")]
csn_dt[, window_end   := score_time]

# Optional: prefilter vitals by relevant accounts
vit_dt <- vit_dt[hsp_account_id %in% csn_dt$hsp_account_id]

# Ensure types
stopifnot(inherits(csn_dt$window_start, "POSIXct"),
          inherits(csn_dt$window_end,   "POSIXct"),
          inherits(vit_dt$recorded_time,"POSIXct"))

# Identify rows with a valid window (we'll match only these)
csn_dt[, valid_window := !is.na(window_start) & !is.na(window_end) & (window_end >= window_start)]

# Index to speed joins without sorting (often faster than setkey on huge tables)
setindex(vit_dt, hsp_account_id, recorded_time)

# Non-equi join only for rows with valid windows
matched <- vit_dt[
  csn_dt[valid_window == TRUE],
  on = .(hsp_account_id,
         recorded_time >= window_start,
         recorded_time <= window_end),
  nomatch = 0L,
  allow.cartesian = TRUE
]

# Count vitals per study_id x vital_type
counts_long <- matched[, .N, by = .(study_id, vital_type)]

# Cast to wide with 0 fill
counts_wide <- dcast(counts_long, study_id ~ vital_type, value.var = "N", fill = 0L)

# Rename to n_<vital_type>_rows
vital_cols <- setdiff(names(counts_wide), "study_id")
setnames(counts_wide, vital_cols, paste0("n_", vital_cols, "_rows"))

# Merge back onto csn_time and fill missing with 0
csn_time_vitals <- merge(csn_dt, counts_wide, by = "study_id", all.x = TRUE)

new_count_cols <- grep("^n_.*_rows$", names(csn_time_vitals), value = TRUE)
for (cc in new_count_cols) {
  set(csn_time_vitals, which(is.na(csn_time_vitals[[cc]])), cc, 0L)
}

# For rows with invalid windows, ensure counts exist and are 0
# (If no vital types existed at all in matched data, new_count_cols could be empty; guard it.)
if (length(new_count_cols) > 0) {
  for (cc in new_count_cols) {
    set(csn_time_vitals, which(csn_time_vitals$valid_window == FALSE), cc, 0L)
  }
}

# Clean up helper cols
csn_time_vitals[, c("window_start", "window_end", "valid_window") := NULL]

csn_time_vitals <- as_tibble(csn_time_vitals)


# Join together by study_id
pros_full_no_abx <- pros_2hr_rf_ue %>% left_join(csn_time_vitals %>% dplyr::select(-pat_enc_csn_id,  -hsp_account_id, -score_time, -picu_adm_date_time), by = "study_id")

# Make fio2_present column
pros_full_no_abx <- pros_full_no_abx %>% mutate(fio2_present = ifelse(n_fio2_rows == 0, yes = 0, no = 1))

# Now need to identify the RF labs, which were present, and how many were used for the score

labs_pros <- read_csv(file = "/phi/sbi/prospective_data/Prospective/data/lab_export_pros_100125.csv")


# Load the prospective data:
epic_7_23  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_Sep2023.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_9_23 <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_Oct2023.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_10_23 <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_Nov2023.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_11_23 <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_Dec2023.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_12_23  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_Jan2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_1_24  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_Feb2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_2_24  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_actual_feb2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_3_24  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_March2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_4_24  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_April2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_5_24  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_May2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_6_24  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_June2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_7_24  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_July2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_8_24  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_August2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

# Fix a few issues with datasets to ensure can bind the rows (i.e. data class)

epic_pros_raw <- bind_rows(epic_9_23, epic_10_23, epic_11_23, epic_12_23, epic_1_24, epic_2_24, epic_3_24, epic_4_24, epic_5_24, epic_6_24,
                           epic_7_24, epic_8_24)

epic_pros_data <- epic_pros_raw

# Rename certain columns
pros_df <- epic_pros_data %>% rename(score_time = INSTANT_UTC_DTTM, model_score = TOTAL_SCORE, model_type = ACUITY_SYSTEM_ID)

pros_df$model_type[pros_df$model_type == "100133"] <- "LR_no_abx"
pros_df$model_type[pros_df$model_type == "100134"] <- "LR_yes_abx"
pros_df$model_type[pros_df$model_type == "100148"] <- "RF_no_abx"
pros_df$model_type[pros_df$model_type == "100149"] <- "RF_yes_abx"

# Lower case of column names
colnames(pros_df) <- str_to_lower(colnames(pros_df))

# Fix relevant date/times
pros_df <- pros_df %>% mutate(score_time = as.POSIXct(score_time, tz = "UTC", format = "%Y-%m-%d %H:%M:%S"))
pros_df$score_time <- with_tz(time = pros_df$score_time, tzone = "America/Denver")

# Rename ICU admission time and fix the timezone
pros_df <- pros_df %>% mutate(picu_adm_date_time = as.POSIXct(icu_start_instant, format = "%Y-%m-%d %H:%M:%S"))
pros_df <- pros_df %>% mutate(picu_adm_date_time = force_tz(picu_adm_date_time, tzone = "America/Denver"))


# Filter for correct study start time
study_start_time <- as.POSIXct("2023-09-01 00:00:01", tz = "America/Denver")
pros_df <- pros_df %>% filter(score_time >= study_start_time)

# Obtain pco2 from BG columns
pco2_pros_raw <- pros_df %>% dplyr::select(pat_enc_csn_id,  model_type, icu_start_instant, score_time, pco2_arterial_tc_ts, pco2_capillary_tc_ts, pco2_venous_tc_ts, pco2_venous_nonpoct_tc_ts)

pco2_pros_raw <- pco2_pros_raw %>% rename(picu_adm_date_time = icu_start_instant)


# Retain only matching rows
# 1. Fix picu_adm_date_time class in pco2_pros_raw and pros_full_no_abx
#    (timezone assumed identical, so we don't force tz conversion)
pco2_pros_raw_fixed <- pco2_pros_raw %>%
  mutate(
    picu_adm_date_time = as.POSIXct(
      picu_adm_date_time,
      format = "%Y-%m-%d %H:%M:%OS"
    )
  )


pros_full_no_abx_fixed <- pros_full_no_abx %>%
  mutate(score_time = force_tz(score_time, tzone = "America/Denver"))

pco2_slim_no_abx <- pco2_pros_raw_fixed %>%
  filter(model_type == "RF_no_abx") %>%
  semi_join(
    pros_full_no_abx_fixed %>% dplyr::select(pat_enc_csn_id, score_time),
    by = c("pat_enc_csn_id", "score_time")
  )

# Helper function to count PCO2 measurements in a packed string
count_pco2 <- function(x) {
  dplyr::case_when(
    is.na(x)    ~ 0L,
    x == "NULL" ~ 0L,
    x == ""     ~ 0L,
    TRUE        ~ stringr::str_count(x, stringr::fixed("~,~")) + 1L
  )
}

# Create new dataframe with counts
count_pco2_no_abx <- pco2_slim_no_abx %>%
  dplyr::mutate(
    n_pco2_arterial  = count_pco2(pco2_arterial_tc_ts),
    n_pco2_capillary = count_pco2(pco2_capillary_tc_ts),

    # Count both venous sources and sum
    n_pco2_venous = count_pco2(pco2_venous_tc_ts) +
      count_pco2(pco2_venous_nonpoct_tc_ts),

    n_pco2_total = n_pco2_arterial +
      n_pco2_capillary +
      n_pco2_venous
  )

count_pco2_no_abx <- count_pco2_no_abx %>% mutate(pco2_gas_blood_present = ifelse(n_pco2_total == 0, yes = 0, no = 1))
count_pco2_no_abx <- count_pco2_no_abx %>% rename(n_pco2_gas_blood = n_pco2_total)

## Ok now combine the co2 data with the main pros no abx df
pros_full_no_abx_fixed <- pros_full_no_abx_fixed %>% left_join(count_pco2_no_abx %>% dplyr::select(pat_enc_csn_id, score_time, n_pco2_gas_blood, pco2_gas_blood_present), by = c("pat_enc_csn_id", "score_time"))




############## Now prospective abx exposed
# Get csn's and time of score to quickly filter available labs and vitals

csn_time_abx <- pros_2hr_rf_ae %>% dplyr::select(study_id, pat_enc_csn_id, hsp_account_id, score_time, picu_adm_date_time) %>% distinct()

csn_dt_abx <- as.data.table(csn_time_abx)
vit_dt_abx <- as.data.table(vitals_pros)

# Drop vitals with missing type (optional but recommended)
vit_dt_abx <- vit_dt_abx[!is.na(vital_type)]

# Define window per study_id
csn_dt_abx[, window_start := picu_adm_date_time - as.difftime(24, units = "hours")]
csn_dt_abx[, window_end   := score_time]

# Optional: prefilter vitals by relevant accounts
vit_dt_abx <- vit_dt_abx[hsp_account_id %in% csn_dt$hsp_account_id]

# Ensure types
stopifnot(inherits(csn_dt_abx$window_start, "POSIXct"),
          inherits(csn_dt_abx$window_end,   "POSIXct"),
          inherits(vit_dt_abx$recorded_time,"POSIXct"))

# Identify rows with a valid window (we'll match only these)
csn_dt_abx[, valid_window := !is.na(window_start) & !is.na(window_end) & (window_end >= window_start)]

# Index to speed joins without sorting (often faster than setkey on huge tables)
setindex(vit_dt_abx, hsp_account_id, recorded_time)

# Non-equi join only for rows with valid windows
matched_abx <- vit_dt_abx[
  csn_dt_abx[valid_window == TRUE],
  on = .(hsp_account_id,
         recorded_time >= window_start,
         recorded_time <= window_end),
  nomatch = 0L,
  allow.cartesian = TRUE
]

# Count vitals per study_id x vital_type
counts_long_abx <- matched_abx[, .N, by = .(study_id, vital_type)]

# Cast to wide with 0 fill
counts_wide_abx <- dcast(counts_long_abx, study_id ~ vital_type, value.var = "N", fill = 0L)

# Rename to n_<vital_type>_rows
vital_cols_abx <- setdiff(names(counts_wide_abx), "study_id")
setnames(counts_wide_abx, vital_cols_abx, paste0("n_", vital_cols_abx, "_rows"))

# Merge back onto csn_time and fill missing with 0
csn_time_vitals_abx <- merge(csn_dt_abx, counts_wide_abx, by = "study_id", all.x = TRUE)

new_count_cols_abx <- grep("^n_.*_rows$", names(csn_time_vitals_abx), value = TRUE)
for (cc in new_count_cols_abx) {
  set(csn_time_vitals_abx, which(is.na(csn_time_vitals_abx[[cc]])), cc, 0L)
}

# For rows with invalid windows, ensure counts exist and are 0
# (If no vital types existed at all in matched data, new_count_cols could be empty; guard it.)
if (length(new_count_cols_abx) > 0) {
  for (cc in new_count_cols_abx) {
    set(csn_time_vitals_abx, which(csn_time_vitals_abx$valid_window == FALSE), cc, 0L)
  }
}

# Clean up helper cols
csn_time_vitals_abx[, c("window_start", "window_end", "valid_window") := NULL]

csn_time_vitals_abx <- as_tibble(csn_time_vitals_abx)


# Join together by study_id
pros_full_yes_abx <- pros_2hr_rf_ae %>% left_join(csn_time_vitals_abx %>% dplyr::select(-pat_enc_csn_id,  -hsp_account_id, -score_time, -picu_adm_date_time), by = "study_id")

# Make fio2_present column
pros_full_yes_abx <- pros_full_yes_abx %>% mutate(fio2_present = ifelse(n_fio2_rows == 0, yes = 0, no = 1))

# Now need to identify the RF labs, which were present, and how many were used for the score
# Obtain hct from _TS column
hct_pros_raw <- pros_df %>% dplyr::select(pat_enc_csn_id,  model_type, picu_adm_date_time, score_time, hematocrit_poc_ts, hematocrit_ts)


# Retain only matching rows
# 1. Fix picu_adm_date_time class in pco2_pros_raw and pros_full_no_abx
#    (timezone assumed identical, so we don't force tz conversion)
hct_pros_raw_fixed <- hct_pros_raw %>%
  mutate(
    picu_adm_date_time = as.POSIXct(
      picu_adm_date_time,
      format = "%Y-%m-%d %H:%M:%OS"
    )
  )


pros_full_yes_abx_fixed <- pros_full_yes_abx %>%
  mutate(score_time = force_tz(score_time, tzone = "America/Denver"))

hct_slim_yes_abx <- hct_pros_raw_fixed %>%
  filter(model_type == "RF_yes_abx") %>%
  semi_join(
    pros_full_yes_abx_fixed %>% dplyr::select(pat_enc_csn_id, score_time),
    by = c("pat_enc_csn_id", "score_time")
  )

# Helper function to count PCO2 measurements in a packed string
count_hct <- function(x) {
  dplyr::case_when(
    is.na(x)    ~ 0L,
    x == "NULL" ~ 0L,
    x == ""     ~ 0L,
    TRUE        ~ stringr::str_count(x, stringr::fixed("~,~")) + 1L
  )
}

# Create new dataframe with counts
count_hct_yes_abx <- hct_slim_yes_abx %>%
  dplyr::mutate(
    n_hct_poc  = count_hct(hematocrit_poc_ts),
    n_hct_reg = count_hct(hematocrit_ts),

    n_hct_total = n_hct_poc +
      n_hct_reg
  )

count_hct_yes_abx <- count_hct_yes_abx %>% mutate(hematocrit_present = ifelse(n_hct_total == 0, yes = 0, no = 1))
count_hct_yes_abx <- count_hct_yes_abx %>% rename(n_hematocrit_blood = n_hct_total)

# Now repeat with lactate data
# Obtain hct from _TS column
lactate_pros_raw <- pros_df %>% dplyr::select(pat_enc_csn_id,  model_type, picu_adm_date_time, score_time, lactate_ts)


# Retain only matching rows
#    (timezone assumed identical, so we don't force tz conversion)
lactate_pros_raw_fixed <- lactate_pros_raw %>%
  mutate(
    picu_adm_date_time = as.POSIXct(
      picu_adm_date_time,
      format = "%Y-%m-%d %H:%M:%OS"
    )
  )

lactate_slim_yes_abx <- lactate_pros_raw_fixed %>%
  filter(model_type == "RF_yes_abx") %>%
  semi_join(
    pros_full_yes_abx_fixed %>% dplyr::select(pat_enc_csn_id, score_time),
    by = c("pat_enc_csn_id", "score_time")
  )

# Helper function to count PCO2 measurements in a packed string
count_lactate <- function(x) {
  dplyr::case_when(
    is.na(x)    ~ 0L,
    x == "NULL" ~ 0L,
    x == ""     ~ 0L,
    TRUE        ~ stringr::str_count(x, stringr::fixed("~,~")) + 1L
  )
}

# Create new dataframe with counts
count_lactate_yes_abx <- lactate_slim_yes_abx %>%
  dplyr::mutate(
    n_lactate  = count_lactate(lactate_ts)
  )

count_lactate_yes_abx <- count_lactate_yes_abx %>% rename(n_lactate_blood = n_lactate)

## Ok now combine the hct and lactate data with the main pros no abx df
pros_full_yes_abx_fixed <- pros_full_yes_abx_fixed %>% left_join(count_hct_yes_abx %>% dplyr::select(pat_enc_csn_id, score_time, n_hematocrit_blood, hematocrit_present), by = c("pat_enc_csn_id", "score_time"))
pros_full_yes_abx_fixed <- pros_full_yes_abx_fixed %>% left_join(count_lactate_yes_abx %>% dplyr::select(pat_enc_csn_id, score_time, n_lactate_blood), by = c("pat_enc_csn_id", "score_time"))

# Add epoch column
pros_full_no_abx_fixed$epoch <- "prospective"
pros_full_yes_abx_fixed$epoch <- "prospective"

# Fix case_id to study_id in retro datasets, fix other column names
retro_full <- retro_full %>% rename(study_id = case_id)
pros_full_no_abx_fixed <- pros_full_no_abx_fixed %>% rename(n_o2sat_rows = n_spo2_rows, n_rr_rows = n_RR_rows)

colnames(retro_full_abx)
retro_full_abx <- retro_full_abx %>% rename(study_id = case_id)
setdiff(names(retro_full_abx), names(pros_full_yes_abx_fixed))

# To do: add picu_adm_date_time back in, add in number of total vital sign rows, PCCC, broad diagnosisi categories, then at 2H the following: IMV, vasoactives in first 2h (calculate # of vasoactives)
# infection labels (SBI type, viral infections, fungal infections), number of cultures drawn (how many before prediction time), timing of antibiotic initiation,

# First get relevant retrospective data
add_retro <- read.fst(path = "/phi/sbi/sbi_blake/model_data_w_ids_9_13_24.fst")
add_slim_retro <- add_retro %>% dplyr::select(case_id, picu_adm_date_time, picu_dc_date_time, icu_los_days, race, ethnicity, death_date,
                                              prism_3_score, pim_2_score, post_op, ccc, cancer, sbi_pneumonia, blood, cns, urine, virus, virus_24,
                                              sbi_cx_neg_sepsis, death, inv_vent)



retro_full_w_vps <- retro_full %>% left_join(add_slim_retro %>% rename(study_id = case_id), by = "study_id")
retro_full_w_vps$picu_adm_date_time <- force_tz(retro_full_w_vps$picu_adm_date_time, tz = "America/Denver")
retro_full_w_vps$picu_dc_date_time <- force_tz(retro_full_w_vps$picu_dc_date_time, tz = "America/Denver")

retro_full_w_vps_abx <- retro_full_abx %>% left_join(add_slim_retro %>% rename(study_id = case_id), by = "study_id")
retro_full_w_vps_abx$picu_adm_date_time <- force_tz(retro_full_w_vps_abx$picu_adm_date_time, tz = "America/Denver")
retro_full_w_vps_abx$picu_dc_date_time <- force_tz(retro_full_w_vps_abx$picu_dc_date_time, tz = "America/Denver")


# Ok now we also want the VPS primary problem data including diagnosis, category, and subcategory for problems present at admission
#Load in the vps data for diagnosis codes and procedures
vps_diagnosis_11_17 <- read_excel(vps_file_path, sheet = "Diagnoses")
vps_diagnosis_18_20 <- readxl::read_xlsx(
  path  = vps_file_path_10_yr,
  sheet = "Diagnoses"
)

vps_diagnosis <- bind_rows(vps_diagnosis_11_17, vps_diagnosis_18_20)

vps_primary_dx_admission <- vps_diagnosis %>% dplyr::select(case_id, STAR_code_descr, category, subcategory, primary_diag, present_on_admit)
vps_primary_dx_admission <- vps_primary_dx_admission %>% filter(primary_diag == "Yes") %>% filter(present_on_admit == "Yes") %>% distinct()
vps_primary_dx_admission$case_id <- as.character(vps_primary_dx_admission$case_id)

# Bind diagnosis data to other retro data
retro_full_w_vps <- retro_full_w_vps %>% left_join(vps_primary_dx_admission %>% rename(study_id = case_id), by = "study_id")
retro_full_w_vps_abx <- retro_full_w_vps_abx %>% left_join(vps_primary_dx_admission %>% rename(study_id = case_id), by = "study_id")

retro_full_w_vps <- retro_full_w_vps %>% dplyr::select(-present_on_admit, -primary_diag)
retro_full_w_vps_abx <- retro_full_w_vps_abx %>% dplyr::select(-present_on_admit, -primary_diag)

# Add the year, quarter, month columns
add_calendar_cols <- function(df, dt_col = "picu_adm_date_time") {
  df %>%
    mutate(
      adm_year    = year(.data[[dt_col]]),
      adm_quarter = quarter(.data[[dt_col]]),              # 1–4
      adm_month   = month(.data[[dt_col]]),                # 1–12
      adm_ym      = format(.data[[dt_col]], "%Y-%m"),      # e.g., "2018-06" (useful for plots)
      adm_yq      = paste0(adm_year, "-Q", adm_quarter)    # e.g., "2018-Q2"
    )
}

retro_full_w_vps_abx <- add_calendar_cols(retro_full_w_vps_abx)
retro_full_w_vps <- add_calendar_cols(retro_full_w_vps)
pros_full_yes_abx_fixed <- add_calendar_cols(pros_full_yes_abx_fixed)
pros_full_no_abx_fixed <- add_calendar_cols(pros_full_no_abx_fixed)

### Placeholder for VPS addition to prospective data

pros_full_w_vps_abx <- pros_full_yes_abx_fixed
pros_full_w_vps <- pros_full_no_abx_fixed


### Add in number of cultures sent in 26 hour period as available for blood, urine, csf, start with retro ###
micro_data_retro <- read_csv("/phi/sbi/updated_retro_data/Retrospective/micro_export_retro_041624_readonly.csv")

# Filter to only those cultures I want

proc_list <- c("MENINGITIS ENCEPHALITIS PANEL (MEP)", "MENINGITIS ENCEPHALITIS PANEL CONDITIONAL (MEP)")
comp_list <- c("CULTURE, ANAEROBIC BLOOD", "CULTURE, BLOOD", "CULTURE, BLOOD ANAEROBIC", "CULTURE, CSF", "CULTURE, URINE", "CULTURE, URINE", "CULTURE, URINE QUANT",
               "MEP CONDITIONAL", "MEP RESULTS", "RESULTS FOR MEP")
micro_data_retro_cultures <- micro_data_retro %>% filter(proc_name %in% proc_list | component_name %in% comp_list) # 266k -> ~19k rows

# Now need to join to case_id via the MRN and time the culture was taken
slim_mrn_retro <- full_retro_name_fix %>% dplyr::select(case_id, mrn, picu_adm_date_time) %>% distinct()

# Add culture category
# Assign tests to culture categories
b_cat <- c("CULTURE, ANAEROBIC BLOOD", "CULTURE, BLOOD", "CULTURE, BLOOD ANAEROBIC")
u_cat <- c("CULTURE, URINE", "CULTURE, URINE", "CULTURE, URINE QUANT")
c_cat <- c( "CULTURE, CSF", "MEP CONDITIONAL", "MEP RESULTS", "RESULTS FOR MEP", "MENINGITIS ENCEPHALITIS PANEL (MEP)")

micro_data_retro_cultures <- micro_data_retro_cultures %>% mutate(cult_cat = case_when(
  component_name %in% b_cat ~ "blood_cx",
  component_name %in% u_cat ~ "urine_cx",
  (component_name %in% c_cat | proc_name %in% c_cat) ~ "csf_cx",
  TRUE ~ "none"
))

micro_slim_retro <- micro_data_retro_cultures %>% dplyr::select(pat_mrn_id, cult_cat, specimn_taken_time)

# Compare culture times to PICU admission times and identify patients who had each culture obtained
micro_slim_retro <- micro_slim_retro %>% rename(mrn = pat_mrn_id)
micro_slim_retro$mrn <- as.character(micro_slim_retro$mrn)

# Need to add in mrn to the main dfs
retro_full_w_vps_w_mrn <- retro_full_w_vps %>% left_join(slim_mrn_retro %>% dplyr::select(case_id, mrn) %>% distinct(), by = c("study_id" = "case_id"))
retro_full_w_vps_w_mrn <- retro_full_w_vps_w_mrn %>% relocate(study_id, mrn, picu_adm_date_time, sbi_present, model_score, adm_month, adm_quarter, adm_year)

# Ensure micro times correct
micro_slim_retro$specimn_taken_time <- force_tz(micro_slim_retro$specimn_taken_time, tz = "America/Denver")

# Now go through and see which retro rows have each type of culture present
# Convert to data.table
retro_dt <- as.data.table(retro_full_w_vps_w_mrn)
micro_dt <- as.data.table(micro_slim_retro)

# Define time window per PICU admission row
retro_dt[, window_start := picu_adm_date_time - as.difftime(24, units = "hours")]
retro_dt[, window_end   := picu_adm_date_time + as.difftime(2,  units = "hours")]

# Optional but helpful: keep only culture rows for MRNs in retro_dt
micro_dt <- micro_dt[mrn %in% retro_dt$mrn]

# Ensure POSIXct (defensive)
stopifnot(inherits(retro_dt$window_start, "POSIXct"),
          inherits(retro_dt$window_end,   "POSIXct"),
          inherits(micro_dt$specimn_taken_time,   "POSIXct"))

# Create a stable row id so we can aggregate back per admission row
retro_dt[, row_id := .I]

# Index micro for speed (doesn't reorder like setkey)
setindex(micro_dt, mrn, specimn_taken_time)

# Non-equi join: cultures within window for each admission row
matched <- micro_dt[
  retro_dt,
  on = .(mrn,
         specimn_taken_time >= window_start,
         specimn_taken_time <= window_end),
  nomatch = 0L,
  allow.cartesian = TRUE
]

# Create presence flags per row_id
flags_long <- matched[, .(
  blood_y_n = as.integer(any(cult_cat == "blood_cx")),
  urine_y_n = as.integer(any(cult_cat == "urine_cx")),
  csf_y_n   = as.integer(any(cult_cat == "csf_cx"))
), by = row_id]

# Merge flags back onto retro, defaulting to 0 when no matches
retro_out <- merge(retro_dt, flags_long, by = "row_id", all.x = TRUE)

for (cc in c("blood_y_n", "urine_y_n", "csf_y_n")) {
  retro_out[is.na(get(cc)), (cc) := 0L]
}

# Cleanup helper columns
retro_out[, c("row_id", "window_start", "window_end") := NULL]

# Final dataframe
retro_full_w_vps_w_mrn <- as_tibble(retro_out)
retro_no_abx_final <- retro_full_w_vps_w_mrn


######### Redo with the abx exposed retrospective dataframe #######
# Need to add in mrn to the main dfs
retro_full_w_vps_w_mrn_abx <- retro_full_w_vps_abx %>% left_join(slim_mrn_retro %>% dplyr::select(case_id, mrn) %>% distinct(), by = c("study_id" = "case_id"))
retro_full_w_vps_w_mrn_abx <- retro_full_w_vps_w_mrn_abx %>% relocate(study_id, mrn, picu_adm_date_time, sbi_present, model_score, adm_month, adm_quarter, adm_year)

# Now go through and see which retro rows have each type of culture present
# Convert to data.table
retro_dt_abx <- as.data.table(retro_full_w_vps_w_mrn_abx)
micro_dt <- as.data.table(micro_slim_retro)

# Define time window per PICU admission row
retro_dt_abx[, window_start := picu_adm_date_time - as.difftime(24, units = "hours")]
retro_dt_abx[, window_end   := picu_adm_date_time + as.difftime(2,  units = "hours")]

# Optional but helpful: keep only culture rows for MRNs in retro_dt
micro_dt <- micro_dt[mrn %in% retro_dt_abx$mrn]

# Ensure POSIXct (defensive)
stopifnot(inherits(retro_dt_abx$window_start, "POSIXct"),
          inherits(retro_dt_abx$window_end,   "POSIXct"),
          inherits(micro_dt$specimn_taken_time,   "POSIXct"))

# Create a stable row id so we can aggregate back per admission row
retro_dt_abx[, row_id := .I]

# Index micro for speed (doesn't reorder like setkey)
setindex(micro_dt, mrn, specimn_taken_time)

# Non-equi join: cultures within window for each admission row
matched_abx <- micro_dt[
  retro_dt_abx,
  on = .(mrn,
         specimn_taken_time >= window_start,
         specimn_taken_time <= window_end),
  nomatch = 0L,
  allow.cartesian = TRUE
]

# Create presence flags per row_id
flags_long_abx <- matched_abx[, .(
  blood_y_n = as.integer(any(cult_cat == "blood_cx")),
  urine_y_n = as.integer(any(cult_cat == "urine_cx")),
  csf_y_n   = as.integer(any(cult_cat == "csf_cx"))
), by = row_id]

# Merge flags back onto retro, defaulting to 0 when no matches
retro_out_abx <- merge(retro_dt_abx, flags_long_abx, by = "row_id", all.x = TRUE)

for (cc in c("blood_y_n", "urine_y_n", "csf_y_n")) {
  retro_out_abx[is.na(get(cc)), (cc) := 0L]
}

# Cleanup helper columns
retro_out_abx[, c("row_id", "window_start", "window_end") := NULL]

# Final dataframe
retro_full_w_vps_w_mrn_abx <- as_tibble(retro_out_abx)
retro_yes_abx_final <- retro_full_w_vps_w_mrn_abx


######## Now get cultures numbers for prospective cohort ########
micro_raw_pros <- read_csv(file = "/phi/sbi/prospective_data/Prospective/data/micro_export_pros_100125.csv")

View(micro_raw_pros %>% dplyr::select(proc_name, component_name) %>% distinct())
proc_abx <- c("ANAEROBIC BLOOD BACTERIAL", "ANAEROBIC BLOOD BACTERIAL CULTURE (SECOND)", "BLOOD BACTERIAL CULTURE", "BLOOD BACTERIAL CULTURE-SECOND", "MENINGITIS ENCEPHALITIS PANEL (MEP)",
              "URINE BACTERIAL CULTURE")
comp_abx <- c("CULTURE, BLOOD", "CULTURE, BLOOD ANAEROBIC", "CULTURE, URINE", "RESULTS FOR MEP")

# Slim down prospective micro dataset
micro_data_pros_cultures <- micro_raw_pros %>% filter(proc_name %in% proc_abx | component_name %in% comp_abx) # #68,902 -> 15,803 rows

# Add culture category
# Assign tests to culture categories
b_cat_pros_proc <- c("ANAEROBIC BLOOD BACTERIAL", "ANAEROBIC BLOOD BACTERIAL CULTURE (SECOND)", "BLOOD BACTERIAL CULTURE", "BLOOD BACTERIAL CULTURE-SECOND")
b_cat_pros_comp <- c("CULTURE, BLOOD", "CULTURE, BLOOD ANAEROBIC")

u_cat_pros_proc <- c("URINE BACTERIAL CULTURE")
u_cat_pros_comp <- c("CULTURE, URINE")

c_cat_pros_proc <- c("MENINGITIS ENCEPHALITIS PANEL (MEP)")
c_cat_pros_comp <- c("RESULTS FOR MEP")

micro_data_pros_cultures <- micro_data_pros_cultures %>% mutate(cult_cat = case_when(
  (component_name %in% b_cat_pros_comp | proc_name %in% b_cat_pros_proc)  ~ "blood_cx",
  (component_name %in% u_cat_pros_comp | proc_name %in% u_cat_pros_proc) ~ "urine_cx",
  (component_name %in% c_cat_pros_comp | proc_name %in% c_cat_pros_proc) ~ "csf_cx",
  TRUE ~ "none"
))

micro_slim_pros <- micro_data_pros_cultures %>% dplyr::select(pat_mrn_id, cult_cat, specimn_taken_time)

# Compare culture times to PICU admission times and identify patients who had each culture obtained
micro_slim_pros <- micro_slim_pros %>% rename(mrn = pat_mrn_id)
micro_slim_pros$mrn <- as.character(micro_slim_pros$mrn)

# slim_mrn_pros definition to add in study_id
pros_full_w_vps_abx <- pros_full_w_vps_abx %>% rename(mrn = pat_mrn_id)
pros_full_w_vps_abx$mrn <- as.character(pros_full_w_vps_abx$mrn)

pros_full_w_vps <- pros_full_w_vps %>% rename(mrn = pat_mrn_id)
pros_full_w_vps$mrn <- as.character(pros_full_w_vps$mrn)

# Update name, Relocate a few things
pros_full_w_vps_w_mrn <- pros_full_w_vps
pros_full_w_vps_w_mrn_abx <- pros_full_w_vps_abx

pros_full_w_vps_w_mrn <- pros_full_w_vps_w_mrn %>% relocate(study_id, mrn, picu_adm_date_time, sbi_present, model_score, adm_month, adm_quarter, adm_year)
pros_full_w_vps_w_mrn_abx <- pros_full_w_vps_w_mrn_abx %>% relocate(study_id, mrn, picu_adm_date_time, sbi_present, model_score, adm_month, adm_quarter, adm_year)

# Now go through and see which pros rows have each type of culture present
# Convert to data.table
pros_dt <- as.data.table(pros_full_w_vps_w_mrn)
micro_dt <- as.data.table(micro_slim_pros)

# Define time window per PICU admission row
pros_dt[, window_start := picu_adm_date_time - as.difftime(24, units = "hours")]
pros_dt[, window_end   := picu_adm_date_time + as.difftime(2,  units = "hours")]

# Optional but helpful: keep only culture rows for MRNs in pros_dt
micro_dt <- micro_dt[mrn %in% pros_dt$mrn]

# Ensure POSIXct (defensive)
stopifnot(inherits(pros_dt$window_start, "POSIXct"),
          inherits(pros_dt$window_end,   "POSIXct"),
          inherits(micro_dt$specimn_taken_time,   "POSIXct"))

# Create a stable row id so we can aggregate back per admission row
pros_dt[, row_id := .I]

# Index micro for speed (doesn't reorder like setkey)
setindex(micro_dt, mrn, specimn_taken_time)

# Non-equi join: cultures within window for each admission row
matched <- micro_dt[
  pros_dt,
  on = .(mrn,
         specimn_taken_time >= window_start,
         specimn_taken_time <= window_end),
  nomatch = 0L,
  allow.cartesian = TRUE
]

# Create presence flags per row_id
flags_long <- matched[, .(
  blood_y_n = as.integer(any(cult_cat == "blood_cx")),
  urine_y_n = as.integer(any(cult_cat == "urine_cx")),
  csf_y_n   = as.integer(any(cult_cat == "csf_cx"))
), by = row_id]

# Merge flags back onto pros, defaulting to 0 when no matches
pros_out <- merge(pros_dt, flags_long, by = "row_id", all.x = TRUE)

for (cc in c("blood_y_n", "urine_y_n", "csf_y_n")) {
  pros_out[is.na(get(cc)), (cc) := 0L]
}

# Cleanup helper columns
pros_out[, c("row_id", "window_start", "window_end") := NULL]

# Final dataframe
pros_full_w_vps_w_mrn <- as_tibble(pros_out)

# Filter to the unexposed (we already only had the rows for the unexposed model)
pros_no_abx_final <- pros_full_w_vps_w_mrn %>% filter(abx_exp == 0)

# Repeat for antibiotics exposed patients / model
pros_dt_abx <- as.data.table(pros_full_w_vps_w_mrn_abx)
micro_dt <- as.data.table(micro_slim_pros)

# Define time window per PICU admission row
pros_dt_abx[, window_start := picu_adm_date_time - as.difftime(24, units = "hours")]
pros_dt_abx[, window_end   := picu_adm_date_time + as.difftime(2,  units = "hours")]

# Optional but helpful: keep only culture rows for MRNs in pros_dt
micro_dt <- micro_dt[mrn %in% pros_dt_abx$mrn]

# Ensure POSIXct (defensive)
stopifnot(inherits(pros_dt_abx$window_start, "POSIXct"),
          inherits(pros_dt_abx$window_end,   "POSIXct"),
          inherits(micro_dt$specimn_taken_time,   "POSIXct"))

# Create a stable row id so we can aggregate back per admission row
pros_dt_abx[, row_id := .I]

# Index micro for speed (doesn't reorder like setkey)
setindex(micro_dt, mrn, specimn_taken_time)

# Non-equi join: cultures within window for each admission row
matched_abx <- micro_dt[
  pros_dt_abx,
  on = .(mrn,
         specimn_taken_time >= window_start,
         specimn_taken_time <= window_end),
  nomatch = 0L,
  allow.cartesian = TRUE
]

# Create presence flags per row_id
flags_long_abx <- matched_abx[, .(
  blood_y_n = as.integer(any(cult_cat == "blood_cx")),
  urine_y_n = as.integer(any(cult_cat == "urine_cx")),
  csf_y_n   = as.integer(any(cult_cat == "csf_cx"))
), by = row_id]

# Merge flags back onto pros, defaulting to 0 when no matches
pros_out_abx <- merge(pros_dt_abx, flags_long_abx, by = "row_id", all.x = TRUE)

for (cc in c("blood_y_n", "urine_y_n", "csf_y_n")) {
  pros_out_abx[is.na(get(cc)), (cc) := 0L]
}

# Cleanup helper columns
pros_out_abx[, c("row_id", "window_start", "window_end") := NULL]

# Final dataframe
pros_full_w_vps_w_mrn_abx <- as_tibble(pros_out_abx)

# Filter to the unexposed (we already only had the rows for the unexposed model)
pros_yes_abx_final <- pros_full_w_vps_w_mrn_abx %>% filter(abx_exp == 1)

### Ensure the SBI types identified in the retrospective dataset (cns, blood, urine) are replicated in the prospective dataset
all_pros_micro <- read_csv(file = "/phi/sbi/retro_data/all_pros_micro_w_sites_1_27_26.csv")

all_pros_micro$order_time <- force_tz(all_pros_micro$order_time, tz = "America/Denver")
all_pros_micro$result_time <- force_tz(all_pros_micro$result_time, tz = "America/Denver")
all_pros_micro$time_obtained <- force_tz(all_pros_micro$time_obtained, tz = "America/Denver")
all_pros_micro$hosp_admit_date_time <- force_tz(all_pros_micro$hosp_admit_date_time, tz = "America/Denver")
all_pros_micro$hosp_disch_time <- force_tz(all_pros_micro$hosp_disch_time, tz = "America/Denver")
all_pros_micro$picu_admit_date_time <- force_tz(all_pros_micro$picu_admit_date_time, tz = "America/Denver")


pros_micro_slim <- all_pros_micro %>% dplyr::select(mrn, hsp_account_id, pat_enc_csn_id, result, description, proc_name, component_name, organism_name, specimen_source, time_obtained,
                                                    class, infxn_site, picu_admit_date_time) %>% distinct()
pros_bact_pros <- pros_micro_slim %>% filter(class == "bacteria")

# Now fix non blood, urine, cns results
bact_sites <- c("blood", "cns", "urinary_tract")
pros_bact_pros$infxn_site[!(pros_bact_pros$infxn_site %in% bact_sites)] <- "other"
pros_bact_pros <- pros_bact_pros %>% filter(infxn_site %in% bact_sites)

# Fix time zones
pros_bact_pros$time_obtained <- force_tz(pros_bact_pros$time_obtained, tzone = "America/Denver")
pros_bact_pros$picu_admit_date_time <- force_tz(pros_bact_pros$picu_admit_date_time, tzone = "America/Denver")

## Now go back to the prospective full datasets and add in the blood, cns, urine columns that the retro ones have
pros_no_abx_final$picu_adm_date_time <- force_tz(pros_no_abx_final$picu_adm_date_time, tzone = "America/Denver")
pros_no_abx_final$score_time <- force_tz(pros_no_abx_final$score_time, tzone = "America/Denver")
pros_no_abx_final$intime <- force_tz(pros_no_abx_final$intime, tzone = "America/Denver")
pros_no_abx_final$outtime <- force_tz(pros_no_abx_final$outtime, tzone = "America/Denver")
pros_no_abx_final$icu_start_instant <- force_tz(pros_no_abx_final$icu_start_instant, tzone = "America/Denver")

pros_yes_abx_final$picu_adm_date_time <- force_tz(pros_yes_abx_final$picu_adm_date_time, tzone = "America/Denver")
pros_yes_abx_final$score_time <- force_tz(pros_yes_abx_final$score_time, tzone = "America/Denver")
pros_yes_abx_final$intime <- force_tz(pros_yes_abx_final$intime, tzone = "America/Denver")
pros_yes_abx_final$outtime <- force_tz(pros_yes_abx_final$outtime, tzone = "America/Denver")
pros_yes_abx_final$icu_start_instant <- force_tz(pros_yes_abx_final$icu_start_instant, tzone = "America/Denver")

# Now add the blood, urine, and cns columns to the no abx pros df
# Convert to data.table
pros_dt <- as.data.table(pros_no_abx_final)
bact_dt <- as.data.table(pros_bact_pros)

# Standardize join key types (MRN)
pros_dt[, mrn := as.character(mrn)]
bact_dt[, mrn := as.character(mrn)]

# Define the window per PICU admission row
pros_dt[, window_start := picu_adm_date_time - as.difftime(24, units = "hours")]
pros_dt[, window_end   := picu_adm_date_time + as.difftime(24,  units = "hours")]

# Optional: shrink bact table to relevant MRNs + sites (speeds things up)
bact_dt <- bact_dt[
  mrn %in% pros_dt$mrn &
    infxn_site %chin% c("blood", "urinary_tract", "cns")
]

# Defensive checks
stopifnot(
  inherits(pros_dt$picu_adm_date_time, "POSIXct"),
  inherits(bact_dt$time_obtained, "POSIXct")
)

# Row id so we can aggregate back per admission row
pros_dt[, row_id := .I]

# Index for speed (doesn't reorder data like setkey)
setindex(bact_dt, mrn, time_obtained)

# Non-equi join: match cultures within the window for each admission row
matched <- bact_dt[
  pros_dt,
  on = .(mrn,
         time_obtained >= window_start,
         time_obtained <= window_end),
  nomatch = 0L,
  allow.cartesian = TRUE
]

# Presence flags per row_id
flags <- matched[, .(
  blood = as.integer(any(infxn_site == "blood")),
  urine = as.integer(any(infxn_site == "urinary_tract")),
  cns   = as.integer(any(infxn_site == "cns"))
), by = row_id]

# Merge back; rows with no matches should be 0
pros_out <- merge(pros_dt, flags, by = "row_id", all.x = TRUE)

for (cc in c("blood", "urine", "cns")) {
  pros_out[is.na(get(cc)), (cc) := 0L]
}

# Drop helper columns
pros_out[, c("row_id", "window_start", "window_end") := NULL]

# Return to tibble if desired
pros_no_abx_final <- as_tibble(pros_out)

### Repeat for the yes_abx cohort
pros_dt <- as.data.table(pros_yes_abx_final)
bact_dt <- as.data.table(pros_bact_pros)

# Standardize join key types (MRN)
pros_dt[, mrn := as.character(mrn)]
bact_dt[, mrn := as.character(mrn)]

# Define the window per PICU admission row
pros_dt[, window_start := picu_adm_date_time - as.difftime(24, units = "hours")]
pros_dt[, window_end   := picu_adm_date_time + as.difftime(24,  units = "hours")]

# Optional: shrink bact table to relevant MRNs + sites (speeds things up)
bact_dt <- bact_dt[
  mrn %in% pros_dt$mrn &
    infxn_site %chin% c("blood", "urinary_tract", "cns")
]

# Defensive checks
stopifnot(
  inherits(pros_dt$picu_adm_date_time, "POSIXct"),
  inherits(bact_dt$time_obtained, "POSIXct")
)

# Row id so we can aggregate back per admission row
pros_dt[, row_id := .I]

# Index for speed (doesn't reorder data like setkey)
setindex(bact_dt, mrn, time_obtained)

# Non-equi join: match cultures within the window for each admission row
matched <- bact_dt[
  pros_dt,
  on = .(mrn,
         time_obtained >= window_start,
         time_obtained <= window_end),
  nomatch = 0L,
  allow.cartesian = TRUE
]

# Presence flags per row_id
flags <- matched[, .(
  blood = as.integer(any(infxn_site == "blood")),
  urine = as.integer(any(infxn_site == "urinary_tract")),
  cns   = as.integer(any(infxn_site == "cns"))
), by = row_id]

# Merge back; rows with no matches should be 0
pros_out <- merge(pros_dt, flags, by = "row_id", all.x = TRUE)

for (cc in c("blood", "urine", "cns")) {
  pros_out[is.na(get(cc)), (cc) := 0L]
}

# Drop helper columns
pros_out[, c("row_id", "window_start", "window_end") := NULL]

# Return to tibble
pros_yes_abx_final <- as_tibble(pros_out)

# Create suspected infection columns for use
retro_no_abx_final$suspected_infection <- 1
retro_yes_abx_final$suspected_infection <- 1

pros_no_abx_final$suspected_infection <- 0
pros_yes_abx_final$suspected_infection <- 0

# Filter out non suspected infection testing
micro_raw_pros

# micro table (convert once)
micro_dt <- as.data.table(micro_raw_pros %>% mutate(pat_mrn_id = as.character(pat_mrn_id)) %>% rename(mrn = pat_mrn_id) %>%
                            mutate(specimn_taken_time = force_tz(specimn_taken_time, tz = "America/Denver")))

# Defensive: drop rows with missing specimen times or mrn
micro_dt <- micro_dt[!is.na(mrn) & !is.na(specimn_taken_time)]

# Index for speed (doesn't reorder like setkey)
setindex(micro_dt, mrn, specimn_taken_time)

add_suspected_infection <- function(pros_df, micro_dt) {
  pros_dt <- as.data.table(pros_df)

  # Ensure MRN is character
  pros_dt[, mrn := as.character(mrn)]

  # Define ±24h window around PICU admission time
  pros_dt[, window_start := picu_adm_date_time - as.difftime(24, units = "hours")]
  pros_dt[, window_end   := picu_adm_date_time + as.difftime(24, units = "hours")]

  # Row id to aggregate back per admission row
  pros_dt[, row_id := .I]

  # (Optional) reduce micro rows to only MRNs present in this cohort
  micro_sub <- micro_dt[mrn %in% pros_dt$mrn]

  # Non-equi join: any specimen within window for that MRN
  matched <- micro_sub[
    pros_dt,
    on = .(mrn,
           specimn_taken_time >= window_start,
           specimn_taken_time <= window_end),
    nomatch = 0L,
    allow.cartesian = TRUE
  ]

  # Any match -> suspected infection = 1
  flags <- matched[, .(suspected_infection = 1L), by = row_id]

  # Merge back, fill missing with 0
  out <- merge(pros_dt, flags, by = "row_id", all.x = TRUE)
  out[is.na(suspected_infection), suspected_infection := 0L]

  # Cleanup helper cols
  out[, c("row_id", "window_start", "window_end") := NULL]

  as_tibble(out)
}

# Apply to both cohorts
# First remove old susp infxn cols
pros_no_abx_final <- pros_no_abx_final %>% dplyr::select(-suspected_infection)
pros_no_abx_final  <- add_suspected_infection(pros_no_abx_final,  micro_dt)
pros_yes_abx_final <- pros_yes_abx_final %>% dplyr::select(-suspected_infection)
pros_yes_abx_final <- add_suspected_infection(pros_yes_abx_final, micro_dt)


# Load in CXR information to incorporate into suspected infection
cxr_pros <- read_csv(file = "/phi/sbi/prospective_data/Prospective/data/hosp_chest_xray_export_pros_100125.csv")
cxr_pros$chest_x_ray_order <- as.POSIXct(x = cxr_pros$chest_x_ray_order, format = "%Y-%m-%d %H:%M:%OS", tz = "America/Denver")


# ---- 0) Prep: make sure timezones/classes are consistent ----
# (Do this once; change tz since timestamps are truly UTC.)
stopifnot(inherits(cxr_pros$chest_x_ray_order, "POSIXct"))
stopifnot(inherits(pros_no_abx_final$picu_adm_date_time, "POSIXct"))
stopifnot(inherits(pros_yes_abx_final$picu_adm_date_time, "POSIXct"))

# Make sure MRNs align type-wise (cxr is numeric; pros_* mrn is chr)
# We'll use character everywhere to avoid dropping leading zeros.
cxr_dt <- as.data.table(cxr_pros) %>%
  mutate(
    mrn = as.character(pat_mrn_id)
  ) %>%
  as.data.table()

# keep only actual CXR orders if needed
cxr_dt <- cxr_dt[chest_xray_yn == 1 & !is.na(chest_x_ray_order)]

# Helper to update suspected_infection in a given pros df
flag_cxr_suspected_infection <- function(pros_df, cxr_dt) {

  pros_dt <- as.data.table(pros_df)
  pros_dt[, mrn := as.character(mrn)]

  # Window bounds around PICU admission
  pros_dt[, `:=`(
    win_start = picu_adm_date_time - 24*60*60,
    win_end   = picu_adm_date_time + 24*60*60
  )]

  # Non-equi join: for each admission row, find any CXR order in the window
  # i.* columns come from pros_dt; x.* from cxr_dt
  hit_ids <- unique(
    cxr_dt[
      pros_dt,
      on = .(mrn, chest_x_ray_order >= win_start, chest_x_ray_order <= win_end),
      nomatch = 0L,
      .(study_id)
    ]$study_id
  )

  # Flip suspected_infection to 1 for any hit
  pros_dt[, suspected_infection := fifelse(study_id %in% hit_ids, 1L, suspected_infection)]

  # Drop helper cols
  pros_dt[, c("win_start", "win_end") := NULL]

  as_tibble(pros_dt)
}

# ---- 1) Apply to both datasets ----
pros_no_abx_final  <- flag_cxr_suspected_infection(pros_no_abx_final,  cxr_dt)
pros_yes_abx_final <- flag_cxr_suspected_infection(pros_yes_abx_final, cxr_dt)

# Fix timezones for the retro datasets
retro_no_abx_final$picu_adm_date_time <- force_tz(retro_no_abx_final$picu_adm_date_time, tz = "America/Denver")
retro_yes_abx_final$picu_adm_date_time <- force_tz(retro_yes_abx_final$picu_adm_date_time, tz = "America/Denver")

retro_no_abx_final$picu_dc_date_time <- force_tz(retro_no_abx_final$picu_dc_date_time, tz = "America/Denver")
retro_yes_abx_final$picu_dc_date_time <- force_tz(retro_yes_abx_final$picu_dc_date_time, tz = "America/Denver")

# Add in preicu location stuff to remainder of retrospective dataset
temp_retro <- read.fst(path = "~/sbi_blake/jan_25_23_model_data_df.fst")
temp_loc <- temp_retro %>% dplyr::select(case_id, pre_icu)
temp_loc <- data.table(temp_loc)
temp_loc <- one_hot(dt = temp_loc, cols = c("pre_icu"))
temp_loc <- temp_loc %>% distinct()

retro_no_abx_final <- left_join(retro_no_abx_final, temp_loc %>% rename(study_id = case_id), by = "study_id")

# Fix cx neg sepsis column name & pneumonia name & pccc and cancer
retro_no_abx_final <- retro_no_abx_final %>% rename(ever_cx_neg_sepsis = sbi_cx_neg_sepsis)
retro_yes_abx_final <- retro_yes_abx_final %>% rename(ever_cx_neg_sepsis = sbi_cx_neg_sepsis)

retro_no_abx_final <- retro_no_abx_final %>% rename(pna_1_0 = sbi_pneumonia)
retro_yes_abx_final <- retro_yes_abx_final %>% rename(pna_1_0 = sbi_pneumonia)

retro_no_abx_final <- retro_no_abx_final %>% rename(pccc = ccc, malignancy_pccc = cancer)
retro_yes_abx_final <- retro_yes_abx_final %>% rename(pccc = ccc, malignancy_pccc = cancer)

# Now get IMV info from flowsheets:
flowsheet <- read.csv("/phi/sbi/sbi_blake/flo_export_pros_100125.csv")
demog     <- read.csv("/phi/sbi/sbi_blake/demog_export_pros_100125.csv")

# Helper functions / variablesbvæ
tz_use <- "America/Denver"

# Robust parser for timestamps with/without fractional seconds
parse_ts <- function(x, tz = tz_use) {
  suppressWarnings(parse_date_time(x, orders = c("Ymd HMS", "Ymd HMSOS"), tz = tz))
}

# Create linker for mrn, csn, har
csn_xwalk <- demog %>%
  transmute(
    pat_enc_csn_id = as.integer(pat_enc_csn_id),
    hsp_account_id = as.integer(hsp_account_id),
    pat_mrn_id     = as.integer(pat_mrn_id),
    dob            = suppressWarnings(mdy(dob))
  ) %>%
  distinct()

# Create map for identifying imv
support_map <- tibble::tribble(
  ~meas_value, ~support_level,
  "Aerosol mask", "supplemental oxygen",
  "Cool mist", "supplemental oxygen",
  "Face tent", "supplemental oxygen",
  "High flow nasal cannula", "HFNC",
  "Manual bag ventilation - artificial airway", "IMV",
  "Manual bag ventilation - mask", "BiPAP",
  "Nasal cannula", "supplemental oxygen",
  "Noninvasive full face mask", "BiPAP",
  "Noninvasive nasal mask", "BiPAP",
  "Noninvasive prongs", "BiPAP",
  "Noninvasive total face", "BiPAP",
  "Oscillator", "IMV",
  "Other - see comment", "exclude",
  "OxiMask", "supplemental oxygen",
  "Oxyhood", "supplemental oxygen",
  "Partial rebreather mask", "supplemental oxygen",
  "Simple mask", "supplemental oxygen",
  "T-piece pneumatic resuscitator - airway", "IMV",
  "Trach direct connect", "supplemental oxygen",
  "Trach HME", "supplemental oxygen",
  "Ventilator", "IMV"
)

# Get df of imv events
imv_events <- flowsheet %>%
  mutate(
    recorded_time  = parse_ts(recorded_time, tz_use),
    hsp_account_id = as.integer(hsp_account_id)
  ) %>%
  inner_join(csn_xwalk %>% dplyr::select(hsp_account_id, pat_enc_csn_id),
             by = "hsp_account_id") %>%
  filter(row_type == "o2 delivery method") %>%
  left_join(support_map, by = "meas_value") %>%
  filter(!is.na(support_level), support_level == "IMV") %>%   # keep IMV only
  transmute(
    pat_enc_csn_id = as.integer(pat_enc_csn_id),
    t_event        = recorded_time
  ) %>%
  filter(!is.na(t_event)) %>%
  distinct()

imv_events_dt <- as.data.table(imv_events)

# Function to id IMV events between 0 to +2 hours after PICU admit
flag_imv_near_admission <- function(pros_df, imv_events_dt,
                                    hrs_before = 0, hrs_after = 2) {

  pros_dt <- as.data.table(pros_df)

  # create window around admission
  pros_dt[, `:=`(
    win_start = picu_adm_date_time - hrs_before * 3600,
    win_end   = picu_adm_date_time + hrs_after  * 3600
  )]

  # non-equi join: IMV event within [win_start, win_end] for same CSN
  hits <- unique(
    imv_events_dt[
      pros_dt,
      on = .(pat_enc_csn_id,
             t_event >= win_start,
             t_event <= win_end),
      nomatch = 0L,
      .(study_id)
    ]$study_id
  )

  pros_dt[, imv_at_picu_adm := fifelse(study_id %in% hits, 1L, 0L)]
  pros_dt[, c("win_start", "win_end") := NULL]

  as_tibble(pros_dt)
}

# Apply to both pros dfs
pros_no_abx_final  <- flag_imv_near_admission(pros_no_abx_final,  imv_events_dt,
                                              hrs_before = 0, hrs_after = 2)

pros_yes_abx_final <- flag_imv_near_admission(pros_yes_abx_final, imv_events_dt,
                                              hrs_before = 0, hrs_after = 2)


# Fix the pre_icu names
retro_no_abx_final <- retro_no_abx_final %>%
  rename(
    preicu_6thfloor          = pre_icu_6th_floor,
    preicu_8thfloor          = pre_icu_8th_floor,
    preicu_9thfloor          = pre_icu_9th_floor,
    preicu_er                = pre_icu_er,
    preicu_hem_onc_bmt       = pre_icu_hem_onc_bmt,
    preicu_nocer             = pre_icu_noc_er,
    preicu_nocinpatient      = pre_icu_noc_inpatient,
    preicu_nonchco           = `pre_icu_Non-CHCO Facility`,
    preicu_or                = pre_icu_operating_room,
    preicu_othericu          = pre_icu_other_CHCO_ICU,
    preicu_outpatient        = pre_icu_outpatient,
    preicu_procedure_center  = pre_icu_procedure_center
  ) %>%
  # These exist in retro but do NOT exist in pros_no_abx_final
  # For now will drop the other admission sources
  dplyr::select(-`pre_icu_Other or Unknown`)

# Fix inv_vent name in retro df
retro_no_abx_final <- retro_no_abx_final %>% rename(imv_at_picu_adm = inv_vent)
retro_yes_abx_final <- retro_yes_abx_final %>% rename(imv_at_picu_adm = inv_vent)

# save Columns in retro_no_abx_final but NOT in pros_no_abx_final for later
# Fix some columns
save_retro_no_abx_cols <- retro_no_abx_final %>% dplyr::select(c(study_id, setdiff(
  names(retro_no_abx_final),
  names(pros_no_abx_final)
)))

retro_no_abx_outcomes <- retro_no_abx_final
retro_no_abx_final <- retro_no_abx_outcomes %>% dplyr::select(-all_of(as.vector(setdiff(
  names(retro_no_abx_final),
  names(pros_no_abx_final)
))))

# Fix retro col names for abx exposed group
retro_yes_abx_final <- retro_yes_abx_final %>% rename(hematocrit_mean = hematocrit_blood, fio2_last = last_fio2, dbp_max = max_dbp, n_spo2_rows = n_o2sat_rows,
                                                      n_RR_rows = n_rr_rows, n_hematocrit_blood = n_hematocrit)

save_retro_yes_abx_cols <- retro_yes_abx_final %>% dplyr::select(c(study_id, setdiff(
  names(retro_yes_abx_final),
  names(pros_yes_abx_final)
)))

retro_yes_abx_outcomes <- retro_yes_abx_final
retro_yes_abx_final <- retro_yes_abx_outcomes %>% dplyr::select(-all_of(as.vector(setdiff(
  names(retro_yes_abx_final),
  names(pros_yes_abx_final)
))))

### Left off needing to remove unecessary columns from the prospective datasdets
save_pros_no_abx_cols <- pros_no_abx_final %>% dplyr::select(c(study_id, setdiff(
  names(pros_no_abx_final),
  names(retro_no_abx_final)
)))

rm_from_pros_no_abx <- c("pat_enc_csn_id", "score_time", "hours_since_picu_adm", "intime", "outtime", "error_id", "pat_id", "line", "model_type", "time_elapsed",
                          "abx_log_last_rdi_utc", "abx_rf_last_rdi_utc", "noabx_log_last_rdi_utc", "noabx_rf_last_rdi_utc", "contact_creation_to_admission_delta", "height",
                          "height_date", "old_race", "weight", "weight_date", "albumin_mean", "ast_present", "bacteria_urine_mean", "bands_perc_mean",
                          "bands_perc_present", "base_excess_present", "bun_mean", "chloride_mean",
                          "crp_mean", "dbp_last_raw", "dbp_max_raw", "dbp_mean_raw",
                          "dbp_median_raw", "dbp_min_raw", "fibrinogen_mean", "gcs_total_min",
                          "gcs_verbal_min", "hematocrit_mean", "hemoglobin_present", "hr_last_raw",
                          "hr_max_raw", "hr_mean_raw", "hr_median_raw", "hr_min_raw",
                          "icu_start_instant", "lactate_max", "lactate_present", "ldh_total_mean",
                          "leukocytes_urine_mean", "leukocytes_urine_present", "lipase_mean", "mcv_mean",
                          "monocytes_perc_mean", "nitrate_urine_mean", "o2sat_count", "o2sat_max",
                          "ph_gas_blood_mean", "ph_urine_mean", "platelet_count_mean", "po2_arterial_present",
                          "po2_venous_mean", "po2_venous_present", "prior_unit", "respiratory_support_any_positive",
                          "rr_mean_raw", "rr_mean", "rr_median_raw", "rr_median",
                          "rr_min_raw", "salicylates_present", "sbp_last_raw", "sbp_max_raw",
                          "sbp_mean_raw", "sbp_median_raw", "sbp_min_raw", "scheduled_admit",
                          "sodium_present", "uric_acid_present", "valid_start_instant", "valid_end_instant",
                          "wbc_csf_present", "wbc_urine_mean", "weight_present", "is_interventional_radiology",
                          "cx_neg_sepsis", "bcx_sent", "micro_sbi_1_0", "cnss_true",
                          "start", "end", "abx_exp", "abx_duration_after_score",
                          "t_diff", "hsp_account_id", "language", "insurance_type",
                          "rowid", "n_map_rows", "n_pupil_reaction_size_bilat_rows", "n_pupil_reaction_size_left_rows",
                          "n_pupil_reaction_size_right_rows", "n_resp_flo_rows", "n_spo2_hr_rows")

pros_no_abx_final <- pros_no_abx_final %>% dplyr::select(-all_of(rm_from_pros_no_abx))


# Add biologic sex to the retro dataframes to ensure we have demographic info on everyone.
get_sex_var <- read.fst(path = "/phi/sbi/sbi_blake/model_data_w_ids_9_13_24.fst") %>% dplyr::select(case_id, sex) %>% distinct()
get_sex_var <- get_sex_var %>% mutate(is_female = ifelse(sex == "Female", yes = 1, no = 0))
get_sex_var <- get_sex_var %>% dplyr::select(-sex) %>% rename(study_id = case_id)

retro_no_abx_final <- retro_no_abx_final %>% left_join(get_sex_var, by = "study_id")
retro_yes_abx_final <- retro_yes_abx_final %>% left_join(get_sex_var, by = "study_id")

# remove / rename certain columns
pros_no_abx_final <- pros_no_abx_final %>% dplyr::select(-origin)

# Now ensure the abx exposed dataframes are taken care of
save_pros_yes_abx_cols <- pros_yes_abx_final %>% dplyr::select(c(study_id, setdiff(
  names(pros_yes_abx_final),
  names(retro_yes_abx_final)
)))

pros_yes_abx_outcomes <- pros_yes_abx_final

pros_yes_abx_final <- pros_yes_abx_outcomes %>% dplyr::select(-all_of(as.vector(setdiff(
  names(pros_yes_abx_final),
  names(retro_yes_abx_final )
))))

# Now add in the reason for suspcion of infection, start with CRP
micro_data_retro$TAKEN_TIME <- force_tz(micro_data_retro$specimn_taken_time, tz = "America/Denver")
micro_data_retro <- micro_data_retro %>% rename(mrn = pat_mrn_id)
micro_data_retro$mrn <- as.character(micro_data_retro$mrn)


# 1) Create a skinny CRP event table (unique per MRN + time)
crp_events <- micro_data_retro %>%
  filter(proc_name == "C-REACTIVE PROTEIN") %>%
  transmute(
    mrn = as.character(mrn),
    taken_time = specimn_taken_time
  ) %>%
  filter(!is.na(mrn), !is.na(taken_time)) %>%
  distinct(mrn, taken_time)

# 2) Function to add crp_pres to any retro df with mrn + picu_adm_date_time
add_crp_pres_dt <- function(retro_df, crp_tbl) {
  r <- as.data.table(retro_df)
  c <- as.data.table(crp_tbl)

  r[, mrn := as.character(mrn)]
  r[, `:=`(win_start = picu_adm_date_time - 24*3600,
           win_end   = picu_adm_date_time + 24*3600)]
  r[, .row_id := .I]

  c[, mrn := as.character(mrn)]

  setkey(c, mrn, taken_time)

  # non-equi: taken_time within window for same MRN
  hits <- c[r, on = .(mrn, taken_time >= win_start, taken_time <= win_end),
            nomatch = 0L, .(.row_id)] |> unique()

  r[, crp_pres := 0L]
  r[hits, on = ".row_id", crp_pres := 1L]

  r[, c(".row_id", "win_start", "win_end") := NULL]
  as.data.frame(r)
}

retro_no_abx_final  <- add_crp_pres_dt(retro_no_abx_final,  crp_events)
retro_yes_abx_final <- add_crp_pres_dt(retro_yes_abx_final, crp_events)

pre_cxr_retro_no <- retro_no_abx_final
pre_cxr_retro_yes <- retro_yes_abx_final

### Now change to procalcitonin criteria
# 1) Create a skinny PCT event table (unique per MRN + time)
pct_events <- micro_data_retro %>%
  filter(proc_name == "PROCALCITONIN") %>%
  transmute(
    mrn = as.character(mrn),
    taken_time = specimn_taken_time
  ) %>%
  filter(!is.na(mrn), !is.na(taken_time)) %>%
  distinct(mrn, taken_time)

# 2) Function to add pct_pres to any retro df with mrn + picu_adm_date_time
add_pct_pres_dt <- function(retro_df, pct_tbl) {
  r <- as.data.table(retro_df)
  c <- as.data.table(pct_tbl)

  r[, mrn := as.character(mrn)]
  r[, `:=`(win_start = picu_adm_date_time - 24*3600,
           win_end   = picu_adm_date_time + 24*3600)]
  r[, .row_id := .I]

  c[, mrn := as.character(mrn)]

  setkey(c, mrn, taken_time)

  # non-equi: taken_time within window for same MRN
  hits <- c[r, on = .(mrn, taken_time >= win_start, taken_time <= win_end),
            nomatch = 0L, .(.row_id)] |> unique()

  r[, pct_pres := 0L]
  r[hits, on = ".row_id", pct_pres := 1L]

  r[, c(".row_id", "win_start", "win_end") := NULL]
  as.data.frame(r)
}

retro_no_abx_final  <- add_pct_pres_dt(retro_no_abx_final,  pct_events)
retro_yes_abx_final <- add_pct_pres_dt(retro_yes_abx_final, pct_events)

# Load in retro cxr data
cxr_retro <- read_csv(file = "/phi/sbi/sbi_blake/chest_xray_retro_040524.csv")
cxr_retro$chest_x_ray_order <- force_tz(cxr_retro$chest_x_ray_order, tz = "America/Denver")

cxr_slim_retro <- cxr_retro %>% filter(!is.na(chest_x_ray_order)) %>% dplyr::select(pat_mrn_id, chest_x_ray_order)

add_cxr_pres_nonequi <- function(retro_df, cxr_slim_retro, window_hours = 24) {

  # retro: one row per admission (interval)
  retro_dt <- as.data.table(retro_df)
  retro_dt[, mrn := as.character(mrn)]
  retro_dt[, picu_adm_date_time := as.POSIXct(picu_adm_date_time, tz = "America/Denver")]
  retro_dt[, .row_id := .I]

  retro_dt[, `:=`(
    win_start = picu_adm_date_time - window_hours * 3600,
    win_end   = picu_adm_date_time + window_hours * 3600
  )]

  # cxr: many rows per MRN (point-interval)
  cxr_dt <- as.data.table(cxr_slim_retro)
  cxr_dt[, mrn := as.character(pat_mrn_id)]
  cxr_dt[, chest_x_ray_order := as.POSIXct(chest_x_ray_order, tz = "America/Denver")]
  cxr_dt <- cxr_dt[!is.na(mrn) & !is.na(chest_x_ray_order),
                   .(mrn, cxr_time = chest_x_ray_order)]
  cxr_dt <- unique(cxr_dt)

  # foverlaps requires start/end columns in BOTH tables
  setnames(retro_dt, c("win_start", "win_end"), c("start", "end"))
  cxr_dt[, `:=`(start = cxr_time, end = cxr_time)]

  setkey(retro_dt, mrn, start, end)
  setkey(cxr_dt,  mrn, start, end)

  # overlaps: each matching row indicates at least one CXR in window
  hits <- foverlaps(cxr_dt, retro_dt, type = "within", nomatch = 0L)

  hit_ids <- unique(hits$.row_id)

  # attach indicator back to original retro
  retro_dt[, cxr_pres := as.integer(.row_id %in% hit_ids)]

  # cleanup + restore column names
  setnames(retro_dt, c("start", "end"), c("win_start", "win_end"))
  retro_dt[, c(".row_id", "win_start", "win_end") := NULL]

  as.data.frame(retro_dt)
}

# Apply
retro_no_abx_final  <- add_cxr_pres_nonequi(retro_no_abx_final,  cxr_slim_retro, window_hours = 24)
retro_yes_abx_final <- add_cxr_pres_nonequi(retro_yes_abx_final, cxr_slim_retro, window_hours = 24)

table(retro_no_abx_final$cxr_pres,  useNA = "ifany")
table(retro_yes_abx_final$cxr_pres, useNA = "ifany")

# Now flag the micro test suspicion of infection (i.e. non CRP labs from micro dataaset)
no_crp_micro_retro <- micro_data_retro %>% filter(!(proc_name %in% c("C-REACTIVE PROTEIN", "CHOLESTEROL, PLEURAL FLUID", "GLUCOSE, CSF", "LDH, PLEURAL FLUID",
                                                                      "MRSA CULTURE", "PERICARD FL CC&DIFF(CYTOSPIN)", "PH,PLEURAL FLUID", "PLEURAL FLD CC&DIFF(CYTOSPIN)",
                                                                      "PROCALCITONIN", "PROTEIN, PLEURAL FLUID", "SHUNT FLUID CC&DIFF(CYTOSPIN)", "SYNOVIAL FL CC&DIFF(CYTOSPIN)",
                                                                      "TRIGLYCERIDE, PLEURAL FLUID", "BAL CELL CNT & DIFF", "ALBUMIN, PLEURAL FLUID", "CSF CELL CNT & DIFF")))



add_any_micro_pres_foverlaps <- function(retro_df, micro_df, window_hours = 24) {

  # Keep input class so get back tibble if passed tibble
  is_tbl <- inherits(retro_df, "tbl_df")

  # ---- retro intervals (one per admission row)
  retro_dt <- as.data.table(copy(retro_df))
  retro_dt[, mrn := as.character(mrn)]
  retro_dt[, t0  := as.POSIXct(picu_adm_date_time, tz = "America/Denver")]
  retro_dt[, rid := .I]
  retro_dt[, `:=`(
    start = t0 - window_hours * 3600,
    end   = t0 + window_hours * 3600
  )]
  retro_int <- retro_dt[, .(rid, mrn, start, end)]

  # ---- micro "point" intervals (TAKEN_TIME treated as start=end)
  micro_dt <- as.data.table(copy(micro_df))
  micro_dt[, mrn := as.character(mrn)]
  micro_dt[, t   := as.POSIXct(specimn_taken_time, tz = "America/Denver")]
  micro_int <- unique(micro_dt[!is.na(mrn) & !is.na(t), .(mrn, start = t, end = t)])

  # ---- overlap (micro within retro window)
  setkey(retro_int, mrn, start, end)
  setkey(micro_int, mrn, start, end)

  hits <- foverlaps(micro_int, retro_int, type = "within", nomatch = 0L)
  any_hit_rid <- unique(hits[, .(rid)])

  # ---- write flag back onto full retro
  retro_dt[, any_micro_pres := as.integer(rid %in% any_hit_rid$rid)]
  retro_dt[, c("rid", "t0", "start", "end") := NULL]

  out <- as.data.frame(retro_dt)
  if (is_tbl) out <- tibble::as_tibble(out)
  out
}


# Apply to both cohorts
retro_no_abx_final  <- add_any_micro_pres_foverlaps(retro_no_abx_final,  no_crp_micro_retro, window_hours = 24)
retro_yes_abx_final <- add_any_micro_pres_foverlaps(retro_yes_abx_final, no_crp_micro_retro, window_hours = 24)

# Quick check
table(retro_no_abx_final$any_micro_pres,  useNA = "ifany")
table(retro_yes_abx_final$any_micro_pres, useNA = "ifany")
mean(retro_no_abx_final$any_micro_pres)
mean(retro_yes_abx_final$any_micro_pres)

# Ok now need to do the same for the prospective for cxr, micro, and pct/crp
nrow(labs_pros %>% filter(procedure_name == "PROCALCITONIN"))
nrow(labs_pros %>% filter(procedure_name == "C-REACTIVE PROTEIN"))
nrow(labs_pros %>% filter(procedure_name == "CRP, HIGHLY SENSITIVE, S"))

labs_pros$specimn_taken_time <- force_tz(labs_pros$specimn_taken_time, tz = "America/Denver")
labs_pros$order_time <- force_tz(labs_pros$order_time, tz = "America/Denver")
labs_pros$result_time <- force_tz(labs_pros$result_time, tz = "America/Denver")


# Helper: generic "any test within +/- window_hours" using foverlaps
# ---------------------------
add_presence_foverlaps <- function(pros_df,
                                   events_df,
                                   event_time_col,
                                   id_col_pros = "mrn",
                                   id_col_events,
                                   window_hours = 24,
                                   new_col) {

  # preserve tibble class
  is_tbl <- inherits(pros_df, "tbl_df")

  # ---- pros intervals (one per PICU admission row)
  pros_dt <- as.data.table(copy(pros_df))
  pros_dt[, (id_col_pros) := as.character(get(id_col_pros))]
  pros_dt[, t0 := as.POSIXct(picu_adm_date_time, tz = "America/Denver")]
  pros_dt[, rid := .I]
  pros_dt[, `:=`(
    start = t0 - window_hours * 3600,
    end   = t0 + window_hours * 3600
  )]
  pros_int <- pros_dt[, .(rid, mrn = get(id_col_pros), start, end)]

  # ---- event "point" intervals (time treated as start=end)
  ev_dt <- as.data.table(copy(events_df))
  ev_dt[, mrn := as.character(get(id_col_events))]
  ev_dt[, t   := as.POSIXct(get(event_time_col), tz = "America/Denver")]
  ev_int <- unique(ev_dt[!is.na(mrn) & !is.na(t), .(mrn, start = t, end = t)])

  # ---- overlap (event within pros window)
  setkey(pros_int, mrn, start, end)
  setkey(ev_int,   mrn, start, end)

  hits <- foverlaps(ev_int, pros_int, type = "within", nomatch = 0L)
  hit_rids <- unique(hits[, rid])

  # ---- write flag back
  pros_dt[, (new_col) := as.integer(rid %in% hit_rids)]
  pros_dt[, c("rid", "t0", "start", "end") := NULL]

  out <- as.data.frame(pros_dt)
  if (is_tbl) out <- tibble::as_tibble(out)
  out
}

crp_pros_df <- labs_pros %>% filter(procedure_name %in% c("C-REACTIVE PROTEIN", "CRP, HIGHLY SENSITIVE, S"))
pct_pros_df <- labs_pros %>% filter(procedure_name == "PROCALCITONIN")

acute_inf_tests_proc <- c(
  "ADENOVIRUS PCR QUAL - STOOL ONLY",
  "ADENOVIRUS PCR QUANTITATIVE",
  "AEROBIC BACTERIAL CULTURE",
  "AMOEBA ANTIBODY",
  "ANAEROBIC BLOOD BACTERIAL",
  "ANAEROBIC BLOOD BACTERIAL CULTURE (SECOND)",
  "ANAEROBIC BROTH CULTURE",
  "ANAEROBIC CULTURE",
  "BACTERIAL CULTURE",
  "BARTONELLA AB PANEL, IGG & IGM",
  "BARTONELLA PCR, B",
  "BARTONELLA PCR, TISSUE",
  "BK VIRUS PCR QUANTITATIVE",
  "BLASTOMYCES AB, EIA, S",
  "BLOOD BACTERIAL CULTURE",
  "BLOOD BACTERIAL CULTURE-SECOND",
  "BLOOD PARASITE EVALUATION PROFILE",
  "BORDETELLA PERTUSSIS AND BORDETELLA PARAPERTUSSIS MOLECULAR DETECTION, PCR - MAYO",
  "BRONCHOALVEOLAR LAVAGE (BAL) CULTURE",
  "BRUCELLA AB SCREEN IGG + IGM",
  "C. DIFFICILE PCR (TOXIGENIC)",
  "CHLAMYDIA TRACHOMATIS AND NEISSERIA GONORRHEA PCR",
  "CMV PCR QUANTITATIVE",
  "COCCIDIOIDES AB W REFLEX, S (BACKUP)",
  "COCCIDIOIDES AB, S",
  "CYSTIC FIBROSIS PATHOGEN PANEL",
  "CYTOMEGALOVIRUS IGM ANTIBODY - UCH",
  "CYTOMEGALOVIRUS PCR QUALITATIVE",
  "DENGUE VIRUS AB, IGG AND IGM,S",
  "EBV PCR QUANTITATIVE",
  "ENTEROVIRUS AND PARECHOVIRUS PCR - ARUP",
  "FUNGAL CULTURE",
  "GIP WITH CDIFF",
  "GIP WITHOUT CDIFF",
  "HANTAVIRUS AB (IGG AND IGM)",
  "HEPATITIS A ANTIBODY, IGM",
  "HEPATITIS ACUTE PANEL",
  "HEPATITIS B CORE AB, IGM",
  "HEPATITIS B CORE TOTAL AB",
  "HEPATITIS B ENVELOPE  AB",
  "HEPATITIS B VIRAL DNA, QN, S",
  "HEPATITIS C RNA DETECT/QUANT, S (BACKUP)",
  "HEPATITIS C VIRUS ANTIBODY",
  "HEPATITIS E ANTIBODY IGM",
  "HERPES SIMPLEX VIRUS PCR",
  "HHV6 PCR QUANTITATIVE",
  "HISTOPLASMA AB",
  "HIV AG/AB SCREEN",
  "HIV-1 CAP/TAQ RNA PCR QUANT",
  "HIV-1 DNA QUALITATIVE PCR",
  "HIV/HCV/HBV NAT",
  "JC VIRUS PCR QUAL",
  "KARIUS PATHOGEN TEST",
  "KINGELLA PCR",
  "LEGIONELLA SPECIES, MOLECULAR DETECTION, PCR",
  "LEPTOSPIRA IGM",
  "LYME DISEASE SEROLOGY, S",
  "MEASLES PCR QUAL",
  "MENINGITIS ENCEPHALITIS PANEL (MEP)",
  "MONO SCRN (HETEROPHILE AB)",
  "MRSA CULTURE",
  "MRSA SA SSTI PCR",
  "MUMPS PCR QUAL",
  "MYCOBACTERIAL (AFB) CULTURE",
  "MYCOBACTERIUM TUBERCULOSIS PCR - NJH",
  "MYCOPLASMA HOMINIS PCR QUAL",
  "MYCOPLASMA PNEUMONIAE PCR QUAL",
  "PARVOVIRUS B19 PCR, P",
  "PARVOVIRUS B19, IGG & IGM",
  "PNEUMOCYSTIS JIROVECII PCR",
  "Q FEVER ANTIBODIES, IGG/IGM, S",
  "RUBEOLA VIRUS ANTIBODY, IGM AND IGG, SERUM",
  "SCHISTOSOMA PARASITE EXAM, UR",
  "SPOTTED FEVER GRP AB IGG IGM",
  "STREP A ONLY CULTURE",
  "STREP A REFLEX GROUP",
  "STREPTOCOCCUS PNEUMONIAE PCR (SEMI-QUANTITATIVE)",
  "TOXO INFANT <6 MO PROFILE",
  "TOXO NONINFANT >6 MO PROFILE",
  "TOXOPLASMA PCR",
  "TRICHOMONAS PCR",
  "URINE BACTERIAL CULTURE",
  "VANCOMYCIN RESISTANT ENTEROCOCCUS CULTURE STOOL OR RECTAL SWAB SOURCE",
  "VARICELLA-ZOSTER VIRUS PCR",
  "VARICELLA-ZOSTER VIRUS,AB,IGM",
  "WEST NILE IGG AND IGM, CSF",
  "WEST NILE VIRUS ABS, IGG & IGM",
  "ZZENTEROVIRUS PCR QUAL"
)

# Filter for only acute infection tests
micro_pros_df <- micro_raw_pros %>% filter(proc_name %in% acute_inf_tests_proc)

# Fix timezone
micro_pros_df <- micro_pros_df %>%
  mutate(
    order_time = force_tz(order_time, tzone = "America/Denver"),
    specimn_taken_time = force_tz(specimn_taken_time, tzone = "America/Denver")
  )


pros_no_abx_with_susp_elements <- add_presence_foverlaps(
  pros_no_abx_final,
  crp_pros_df,
  event_time_col = "specimn_taken_time",
  id_col_pros = "mrn",
  id_col_events = "pat_mrn_id",
  window_hours = 24,
  new_col = "crp_pres"
)

pros_no_abx_with_susp_elements <- add_presence_foverlaps(
  pros_no_abx_with_susp_elements,
  pct_pros_df,
  event_time_col = "specimn_taken_time",
  id_col_pros = "mrn",
  id_col_events = "pat_mrn_id",
  window_hours = 24,
  new_col = "pct_pres"
)

pros_no_abx_with_susp_elements <- add_presence_foverlaps(
  pros_no_abx_with_susp_elements,
  micro_pros_df,
  event_time_col = "specimn_taken_time",
  id_col_pros = "mrn",
  id_col_events = "pat_mrn_id",
  window_hours = 24,
  new_col = "micro_pres"
)

# Now add CXR data
pros_no_abx_with_susp_elements <- add_presence_foverlaps(
  pros_no_abx_with_susp_elements,
  cxr_pros,
  event_time_col = "chest_x_ray_order",
  id_col_pros = "mrn",
  id_col_events = "pat_mrn_id",
  window_hours = 24,
  new_col = "cxr_pres"
)


# Obtain slim version of the dataset for plotting
pros_susp_elements_to_plot_no_abx <- pros_no_abx_with_susp_elements %>% dplyr::select(study_id, picu_adm_date_time,
                                      crp_pres, pct_pres, micro_pres, cxr_pres, suspected_infection) %>% distinct()


### Now repeat with prospective abx exposed group ###
pros_yes_abx_with_susp_elements <- add_presence_foverlaps(
  pros_yes_abx_final,
  crp_pros_df,
  event_time_col = "specimn_taken_time",
  id_col_pros = "mrn",
  id_col_events = "pat_mrn_id",
  window_hours = 24,
  new_col = "crp_pres"
)

pros_yes_abx_with_susp_elements <- add_presence_foverlaps(
  pros_yes_abx_with_susp_elements,
  pct_pros_df,
  event_time_col = "specimn_taken_time",
  id_col_pros = "mrn",
  id_col_events = "pat_mrn_id",
  window_hours = 24,
  new_col = "pct_pres"
)

pros_yes_abx_with_susp_elements <- add_presence_foverlaps(
  pros_yes_abx_with_susp_elements,
  micro_pros_df,
  event_time_col = "specimn_taken_time",
  id_col_pros = "mrn",
  id_col_events = "pat_mrn_id",
  window_hours = 24,
  new_col = "micro_pres"
)

# Now add CXR data
pros_yes_abx_with_susp_elements <- add_presence_foverlaps(
  pros_yes_abx_with_susp_elements,
  cxr_pros,
  event_time_col = "chest_x_ray_order",
  id_col_pros = "mrn",
  id_col_events = "pat_mrn_id",
  window_hours = 24,
  new_col = "cxr_pres"
)


# Obtain slim version of the dataset for plotting
pros_susp_elements_to_plot_yes_abx <- pros_yes_abx_with_susp_elements %>% dplyr::select(study_id, picu_adm_date_time,
                                                                                      crp_pres, pct_pres, micro_pres, cxr_pres, suspected_infection) %>% distinct()




# Generate first PICU stay per hospital admission
pros_no_abx_1st_infxn <- pros_no_abx_final %>%
  dplyr::filter(stringr::str_ends(study_id, "_1"))

pros_yes_abx_1st_infxn <- pros_yes_abx_final %>%
  dplyr::filter(stringr::str_ends(study_id, "_1"))

# Fix PCCC rows
pros_no_abx_1st_infxn$pccc <- as.factor(pros_no_abx_1st_infxn$pccc)
pros_no_abx_1st_infxn$malignancy_pccc <- as.factor(pros_no_abx_1st_infxn$malignancy_pccc)

pros_yes_abx_1st_infxn$pccc <- as.factor(pros_yes_abx_1st_infxn$pccc)
pros_yes_abx_1st_infxn$malignancy_pccc <- as.factor(pros_yes_abx_1st_infxn$malignancy_pccc)


# Calculate score of retro model on prospective data
source(file = "/phi/sbi/sbi_blake/aim_1_paper_materials/retro_model_to_pros_data.R")
# source(file = "/phi/sbi/sbi_blake/aim_1_paper_materials/retro_model_to_pros_susp_infxn.R")

### Replace old, incorrect prospective model scores with correct scores ###
pros_no_abx_1st_infxn <- pros_no_abx_1st_infxn %>%
  mutate(model_score = rf_pred_prob_pros[, "yes"])

pros_yes_abx_1st_infxn <- pros_yes_abx_1st_infxn %>%
  mutate(model_score = rf_pred_prob_pros_abx[, "yes"])

## Create table for suplement with list of predictors ##
## Load packages
library(dplyr)
library(stringr)
library(tibble)
library(flextable)
library(officer)

## Your predictor vectors
predictors <- c(
  "age", "is_female", "malignancy_pccc", "pccc", "albumin_mean", "ast_present",
  "base_excess_present", "bun_mean", "chloride_mean", "crp_mean", "dbp_last",
  "dbp_max", "dbp_mean", "dbp_median", "dbp_min", "dbp_slope", "fio2_count",
  "fio2_last", "fio2_max", "fio2_mean", "fio2_median", "fio2_min", "fio2_present",
  "fio2_slope", "gcs_total_min", "gcs_verbal_min", "hemoglobin_present",
  "hr_last", "hr_max", "hr_mean", "hr_median", "hr_min", "hr_slope",
  "lactate_present", "leukocytes_urine_present", "o2sat_count", "o2sat_max",
  "o2sat_mean", "o2sat_median", "o2sat_min", "o2sat_slope", "ph_urine_mean",
  "po2_arterial_present", "po2_venous_present", "preicu_9thfloor", "preicu_er",
  "preicu_nocer", "preicu_nonchco", "preicu_or", "respiratory_support_any_positive",
  "rr_mean", "rr_median", "rr_min", "rr_slope", "sbp_last", "sbp_max",
  "sbp_mean", "sbp_median", "sbp_min", "sbp_slope", "scheduled_admit",
  "sodium_present", "temp_last", "temp_max", "temp_mean", "temp_median",
  "temp_min", "temp_slope", "weight_present"
)

predictors_abx <- c(
  "age", "hematocrit_blood", "lactate_max", "dbp_last", "last_fio2",
  "hr_last", "sbp_last", "temp_last", "los_before_icu_days", "max_dbp",
  "hr_max", "sbp_max", "temp_max", "dbp_mean", "fio2_mean",
  "hr_mean", "rr_mean", "o2sat_mean", "sbp_mean", "temp_mean",
  "dbp_median", "fio2_median", "hr_median", "rr_median", "sbp_median",
  "temp_median", "dbp_min", "hr_min", "rr_min", "o2sat_min",
  "sbp_min", "o2sat_count", "dbp_slope", "hr_slope", "rr_slope",
  "o2sat_slope", "sbp_slope", "temp_slope", "preicu_6thfloor", "preicu_8thfloor",
  "preicu_9thfloor", "preicu_er", "preicu_hem_onc_bmt", "preicu_nocer",
  "preicu_nocinpatient", "preicu_nonchco", "preicu_or", "preicu_othericu",
  "preicu_outpatient", "preicu_procedure_center"
)

## ---------------------------------------------------------
## 1. Standardize variable names across the two models
## ---------------------------------------------------------
standardize_predictor_name <- function(x) {
  dplyr::case_when(
    x == "last_fio2" ~ "fio2_last",
    x == "max_dbp" ~ "dbp_max",
    TRUE ~ x
  )
}

predictors_std     <- standardize_predictor_name(predictors)
predictors_abx_std <- standardize_predictor_name(predictors_abx)

## ---------------------------------------------------------
## 2. Draft definitions for variables
##    Edit these as needed for manuscript precision
## ---------------------------------------------------------
definition_lookup <- tibble::tribble(
  ~Predictor, ~Definition,
  "age", "Age at PICU admission.",
  "is_female", "Indicator for female sex.",
  "malignancy_pccc", "Indicator for malignancy-related pediatric complex chronic condition.",
  "pccc", "Indicator for any pediatric complex chronic condition.",
  "albumin_mean", "Mean serum albumin value in the predictor window.",
  "ast_present", "Indicator that aspartate aminotransferase was measured in the predictor window.",
  "base_excess_present", "Indicator that base excess was measured in the predictor window.",
  "bun_mean", "Mean blood urea nitrogen value in the predictor window.",
  "chloride_mean", "Mean serum chloride value in the predictor window.",
  "crp_mean", "Mean C-reactive protein value in the predictor window.",
  "dbp_last", "Most recent diastolic blood pressure in the predictor window.",
  "dbp_max", "Maximum diastolic blood pressure in the predictor window.",
  "dbp_mean", "Mean diastolic blood pressure in the predictor window.",
  "dbp_median", "Median diastolic blood pressure in the predictor window.",
  "dbp_min", "Minimum diastolic blood pressure in the predictor window.",
  "dbp_slope", "Slope of diastolic blood pressure over time in the predictor window.",
  "fio2_count", "Number of recorded FiO2 measurements in the predictor window.",
  "fio2_last", "Most recent FiO2 in the predictor window.",
  "fio2_max", "Maximum FiO2 in the predictor window.",
  "fio2_mean", "Mean FiO2 in the predictor window.",
  "fio2_median", "Median FiO2 in the predictor window.",
  "fio2_min", "Minimum FiO2 in the predictor window.",
  "fio2_present", "Indicator that FiO2 was recorded in the predictor window.",
  "fio2_slope", "Slope of FiO2 over time in the predictor window.",
  "gcs_total_min", "Minimum total Glasgow Coma Scale score in the predictor window.",
  "gcs_verbal_min", "Minimum Glasgow Coma Scale verbal subscore in the predictor window.",
  "hemoglobin_present", "Indicator that hemoglobin was measured in the predictor window.",
  "hematocrit_blood", "Blood hematocrit value in the predictor window.",
  "hr_last", "Most recent heart rate in the predictor window.",
  "hr_max", "Maximum heart rate in the predictor window.",
  "hr_mean", "Mean heart rate in the predictor window.",
  "hr_median", "Median heart rate in the predictor window.",
  "hr_min", "Minimum heart rate in the predictor window.",
  "hr_slope", "Slope of heart rate over time in the predictor window.",
  "lactate_present", "Indicator that lactate was measured in the predictor window.",
  "lactate_max", "Maximum lactate value in the predictor window.",
  "leukocytes_urine_present", "Indicator that urine leukocytes were measured or present in the predictor window.",
  "los_before_icu_days", "Hospital length of stay before PICU admission, in days.",
  "o2sat_count", "Number of recorded oxygen saturation measurements in the predictor window.",
  "o2sat_max", "Maximum oxygen saturation in the predictor window.",
  "o2sat_mean", "Mean oxygen saturation in the predictor window.",
  "o2sat_median", "Median oxygen saturation in the predictor window.",
  "o2sat_min", "Minimum oxygen saturation in the predictor window.",
  "o2sat_slope", "Slope of oxygen saturation over time in the predictor window.",
  "ph_urine_mean", "Mean urine pH in the predictor window.",
  "po2_arterial_present", "Indicator that arterial PO2 was measured in the predictor window.",
  "po2_venous_present", "Indicator that venous PO2 was measured in the predictor window.",
  "preicu_6thfloor", "Indicator that patient location before PICU was the 6th floor.",
  "preicu_8thfloor", "Indicator that patient location before PICU was the 8th floor.",
  "preicu_9thfloor", "Indicator that patient location before PICU was the 9th floor.",
  "preicu_er", "Indicator that patient location before PICU was the emergency department.",
  "preicu_hem_onc_bmt", "Indicator that patient location before PICU was hematology/oncology/BMT.",
  "preicu_nocer", "Indicator that patient location before PICU was a non-CHCO emergency room.",
  "preicu_nocinpatient", "Indicator that patient location before PICU was a non-ICU inpatient unit.",
  "preicu_nonchco", "Indicator that patient location before PICU was outside the primary institution.",
  "preicu_or", "Indicator that patient location before PICU was the operating room.",
  "preicu_othericu", "Indicator that patient location before PICU was another ICU.",
  "preicu_outpatient", "Indicator that patient location before PICU was an outpatient setting.",
  "preicu_procedure_center", "Indicator that patient location before PICU was a procedure center.",
  "respiratory_support_any_positive", "Indicator for any positive-pressure respiratory support in the predictor window.",
  "rr_mean", "Mean respiratory rate in the predictor window.",
  "rr_median", "Median respiratory rate in the predictor window.",
  "rr_min", "Minimum respiratory rate in the predictor window.",
  "rr_slope", "Slope of respiratory rate over time in the predictor window.",
  "sbp_last", "Most recent systolic blood pressure in the predictor window.",
  "sbp_max", "Maximum systolic blood pressure in the predictor window.",
  "sbp_mean", "Mean systolic blood pressure in the predictor window.",
  "sbp_median", "Median systolic blood pressure in the predictor window.",
  "sbp_min", "Minimum systolic blood pressure in the predictor window.",
  "sbp_slope", "Slope of systolic blood pressure over time in the predictor window.",
  "scheduled_admit", "Indicator that PICU admission was scheduled/elective.",
  "sodium_present", "Indicator that sodium was measured in the predictor window.",
  "temp_last", "Most recent temperature in the predictor window.",
  "temp_max", "Maximum temperature in the predictor window.",
  "temp_mean", "Mean temperature in the predictor window.",
  "temp_median", "Median temperature in the predictor window.",
  "temp_min", "Minimum temperature in the predictor window.",
  "temp_slope", "Slope of temperature over time in the predictor window.",
  "weight_present", "Indicator that weight was recorded in the predictor window."
)

## ---------------------------------------------------------
## 3. Build the combined supplement table
## ---------------------------------------------------------
all_predictors <- sort(unique(c(predictors_std, predictors_abx_std)))

supp_predictor_table <- tibble(
  Predictor = all_predictors
) %>%
  dplyr::mutate(
    `Model(s) Used By` = dplyr::case_when(
      Predictor %in% predictors_std & Predictor %in% predictors_abx_std ~ "Both",
      Predictor %in% predictors_std ~ "RF Unexposed",
      Predictor %in% predictors_abx_std ~ "RF Exposed",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::left_join(definition_lookup, by = "Predictor") %>%
  dplyr::rename(`Definition` = Definition) %>%
  dplyr::arrange(Predictor)

## If any definitions are missing, fill with a placeholder for manual editing
supp_predictor_table <- supp_predictor_table %>%
  dplyr::mutate(
    Definition = dplyr::if_else(
      is.na(Definition),
      "Definition to be added.",
      Definition
    )
  )

## Write the predictors to eTable 1
write.csv(x = supp_predictor_table, file = "/phi/sbi/sbi_blake/aim_1_paper_materials/eTable_1_predictors.csv")



# =============================================================================
# 4-panel ROC + 4-panel PRC (retro/pros × abx yes/no) with bootstrap CIs
# =============================================================================
# ----------------------------
# 0) Small helpers
# ----------------------------
to_prob01 <- function(x) {
  x <- as.numeric(x)
  if (max(x, na.rm = TRUE) > 1.5) x <- x / 100  # convert 0-100 -> 0-1 if needed
  x
}

# ----------------------------
# 1) Build the 4 scenario tables
#    (THIS is the only part you should need to “wire”)
# ----------------------------

# Retro test: unexposed (already in test_pr_df: truth + p_yes)
# Retro test: unexposed (already in test_pr_df: truth + p_yes)
retro_no_abx <- test_pr_df %>%
  transmute(
    scenario  = "Retro • Abx-",
    truth     = factor(as.character(truth), levels = c("yes","no")),
    truth_num = as.integer(truth == "yes"),
    score     = to_prob01(p_yes)
  )

# Retro test: exposed (already in test_pr_df_abx: truth + p_yes)
retro_yes_abx <- test_pr_df_abx %>%
  transmute(
    scenario  = "Retro • Abx+",
    truth     = factor(as.character(truth), levels = c("yes","no")),
    truth_num = as.integer(truth == "yes"),
    score     = to_prob01(p_yes)
  )

# Prospective: unexposed
pros_no_abx <- tibble(
  scenario  = "Pros • Abx-",
  truth_num = as.integer(pros_no_abx_1st_infxn$sbi_present),
  truth     = factor(if_else(truth_num == 1L, "yes", "no"), levels = c("yes","no")),
  score     = to_prob01(rf_pred_prob_pros[, "yes"])
)

# Prospective: exposed
pros_yes_abx <- tibble(
  scenario  = "Pros • Abx+",
  truth_num = as.integer(pros_yes_temp$sbi_present),
  truth     = factor(if_else(truth_num == 1L, "yes", "no"), levels = c("yes","no")),
  score     = to_prob01(rf_pred_prob_pros_abx[, "yes"])
)

dat_all <- bind_rows(
  reto_yes_abx = retro_yes_abx,
  reto_no_abx  = retro_no_abx,
  pros_yes_abx = pros_yes_abx,
  pros_no_abx  = pros_no_abx
) %>%
  filter(!is.na(truth_num), !is.na(score), is.finite(score)) %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c("Retro • Abx+", "Retro • Abx-", "Pros • Abx+", "Pros • Abx-")
    )
  )

dat_list <- split(dat_all, dat_all$scenario)

### Create calibration summary table and plots
## Combine all four scenario dataframes into one
calib_df <- purrr::list_rbind(dat_list)

## Create 0.1 probability bins and summarize calibration within each bin
calib_summary <- calib_df %>%
  dplyr::mutate(
    prob_bin = cut(
      score,
      breaks = seq(0, 1, by = 0.1),
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  dplyr::group_by(scenario, prob_bin) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_pred = mean(score, na.rm = TRUE),
    obs_rate = mean(truth_num, na.rm = TRUE),
    events = sum(truth_num, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    bin_lower = seq(0, 0.9, by = 0.1)[as.numeric(prob_bin)],
    bin_upper = bin_lower + 0.1,
    bin_mid = bin_lower + 0.05
  )

## View the calibration summary table
calib_summary

# Make faceted plot
p_calibration <- ggplot(calib_summary, aes(x = mean_pred, y = obs_rate, color = scenario)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray50") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_text(
    aes(label = paste0("n=", n)),
    vjust = -0.8,
    size = 3.5,
    show.legend = FALSE
  ) +
  facet_wrap(~ scenario, ncol = 2) +
  scale_color_manual(
    values = c(
      "Retro • Abx+" = "red",
      "Retro • Abx-" = "red",
      "Pros • Abx+"  = "blue",
      "Pros • Abx-"  = "blue"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  labs(
    title = "Calibration plots by scenario",
    x = "Mean predicted probability within bin",
    y = "Observed event rate within bin"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

p_calibration


# ----------------------------
# 2) Curves + metrics
# ----------------------------

roc_curve_df <- function(df_one) {
  if (dplyr::n_distinct(df_one$truth_num) < 2) {
    return(tibble(scenario = df_one$scenario[1], fpr = NA_real_, tpr = NA_real_))
  }
  truth_case <- 1L - as.integer(df_one$truth_num)
  score_case <- 1 - df_one$score
  roc_obj <- pROC::roc(truth_case, score_case, quiet = TRUE, direction = "<")
  tibble(
    scenario = df_one$scenario[1],
    fpr      = 1 - roc_obj$specificities,
    tpr      = roc_obj$sensitivities
  )
}

auroc_one <- function(df_one) {
  truth_case <- 1L - as.integer(df_one$truth_num)
  score_case <- 1 - df_one$score
  roc_obj <- pROC::roc(truth_case, score_case, quiet = TRUE, direction = "<")
  as.numeric(pROC::auc(roc_obj))
}

pr_curve_df <- function(df_one) {
  if (dplyr::n_distinct(df_one$truth_num) < 2) {
    return(tibble(scenario = df_one$scenario[1], recall = NA_real_, precision = NA_real_))
  }

  df_case <- df_one %>%
    mutate(
      truth = factor(if_else(as.integer(truth_num) == 0L, "yes", "no"), levels = c("yes", "no")),
      score = 1 - score
    )
  pc <- yardstick::pr_curve(df_case, truth, score, event_level = "first")

  if (all(c("precision", "recall") %in% names(pc))) {
    out <- pc %>% transmute(recall = .data[["recall"]], precision = .data[["precision"]])
  } else if (all(c("precision", "sensitivity") %in% names(pc))) {
    out <- pc %>% transmute(recall = .data[["sensitivity"]], precision = .data[["precision"]])
  } else {
    stop("Unexpected columns from pr_curve(): ", paste(names(pc), collapse = ", "))
  }

  out %>%
    mutate(scenario = df_one$scenario[1]) %>%
    filter(is.finite(recall), is.finite(precision))
}

auprc_one <- function(df_one) {
  df_case <- df_one %>%
    mutate(
      truth = factor(if_else(as.integer(truth_num) == 0L, "yes", "no"), levels = c("yes", "no")),
      score = 1 - score
    )
  yardstick::pr_auc(df_case, truth, score, event_level = "first") %>%
    pull(.estimate) %>%
    as.numeric()
}

roc_df <- map_dfr(dat_list, roc_curve_df)
pr_df  <- map_dfr(dat_list, pr_curve_df)

metrics_df <- map_dfr(dat_list, ~ tibble(
  scenario   = .x$scenario[1],
  n          = nrow(.x),
  prevalence = mean(.x$truth_num == 0L),
  AUROC      = auroc_one(.x),
  AUPRC      = auprc_one(.x)
)) %>%
  mutate(
    scenario = factor(scenario, levels = levels(dat_all$scenario)),
    auprc_prev_ratio = AUPRC / prevalence,
    roc_label = paste0(
      "AUROC = ", sprintf("%.3f", AUROC)
    ),
    pr_label = paste0(
      "AUPRC = ", sprintf("%.3f", AUPRC), "\n",
      "AUPRC / Prev = ",
      sprintf("%.3f", AUPRC), " / ",
      sprintf("%.3f", prevalence), " = ",
      sprintf("%.1f", auprc_prev_ratio)
    )
  )

# ----------------------------
# 3) Plots
# ----------------------------

roc_df <- roc_df %>%
  mutate(
    cohort_type = case_when(
      grepl("retro", scenario, ignore.case = TRUE) ~ "Retrospective",
      grepl("pros", scenario, ignore.case = TRUE) ~ "Prospective",
      TRUE ~ "Other"
    )
  )

pr_df <- pr_df %>%
  mutate(
    cohort_type = case_when(
      grepl("retro", scenario, ignore.case = TRUE) ~ "Retrospective",
      grepl("pros", scenario, ignore.case = TRUE) ~ "Prospective",
      TRUE ~ "Other"
    )
  )

roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, color = cohort_type)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_line(linewidth = 1) +
  facet_wrap(~ scenario, ncol = 2) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_color_manual(values = c("Retrospective" = "red", "Prospective" = "blue", "Other" = "black")) +
  labs(
    title = "ROC Curves (Retrospective TEST + Prospective)",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  geom_text(
    data = metrics_df,
    aes(x = 0.35, y = 0.15, label = roc_label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 8),
    legend.position = "none"
  )

pr_plot <- ggplot(pr_df, aes(x = recall, y = precision, color = cohort_type)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ scenario, ncol = 2) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_color_manual(values = c("Retrospective" = "red", "Prospective" = "blue", "Other" = "black")) +
  labs(
    title = "Precision–Recall Curves (Retrospective TEST + Prospective)",
    x = "Recall",
    y = "Precision"
  ) +
  geom_text(
    data = metrics_df,
    aes(x = 0.10, y = 0.15, label = pr_label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 8),
    legend.position = "none"
  )

print(roc_plot)
print(pr_plot)
metrics_df

# ----------------------------
# 4) Repeat ROC/PR evaluation for prospective patients with suspected infection
# ----------------------------

pros_no_abx_susp <- tibble(
  suspected_infection = pros_no_abx_1st_infxn$suspected_infection,
  sbi_present         = pros_no_abx_1st_infxn$sbi_present,
  score_raw           = rf_pred_prob_pros[, "yes"]
) %>%
  filter(suspected_infection == 1L) %>%
  transmute(
    scenario  = "Pros • Abx-",
    truth_num = as.integer(sbi_present),
    truth     = factor(if_else(truth_num == 1L, "yes", "no"), levels = c("yes", "no")),
    score     = to_prob01(score_raw)
  )

pros_yes_abx_susp <- tibble(
  suspected_infection = pros_yes_temp$suspected_infection,
  sbi_present         = pros_yes_temp$sbi_present,
  score_raw           = rf_pred_prob_pros_abx[, "yes"]
) %>%
  filter(suspected_infection == 1L) %>%
  transmute(
    scenario  = "Pros • Abx+",
    truth_num = as.integer(sbi_present),
    truth     = factor(if_else(truth_num == 1L, "yes", "no"), levels = c("yes", "no")),
    score     = to_prob01(score_raw)
  )

dat_all_susp <- bind_rows(
  reto_yes_abx = retro_yes_abx,
  reto_no_abx  = retro_no_abx,
  pros_yes_abx = pros_yes_abx_susp,
  pros_no_abx  = pros_no_abx_susp
) %>%
  filter(!is.na(truth_num), !is.na(score), is.finite(score)) %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c("Retro • Abx+", "Retro • Abx-", "Pros • Abx+", "Pros • Abx-")
    )
  )

dat_list_susp <- split(dat_all_susp, dat_all_susp$scenario)

roc_df_susp <- map_dfr(dat_list_susp, roc_curve_df)
pr_df_susp  <- map_dfr(dat_list_susp, pr_curve_df)

metrics_df_susp <- map_dfr(dat_list_susp, ~ tibble(
  scenario   = .x$scenario[1],
  n          = nrow(.x),
  prevalence = mean(.x$truth_num == 1L),
  AUROC      = auroc_one(.x),
  AUPRC      = auprc_one(.x)
)) %>%
  mutate(
    scenario = factor(scenario, levels = levels(dat_all_susp$scenario)),
    auprc_prev_ratio = AUPRC / prevalence,
    roc_label = paste0(
      "AUROC = ", sprintf("%.3f", AUROC)
    ),
    pr_label = paste0(
      "AUPRC = ", sprintf("%.3f", AUPRC), "\n",
      "AUPRC / Prev = ",
      sprintf("%.3f", AUPRC), " / ",
      sprintf("%.3f", prevalence), " = ",
      sprintf("%.1f", auprc_prev_ratio)
    )
  )

roc_df_susp <- roc_df_susp %>%
  mutate(
    cohort_type = case_when(
      grepl("Retro", scenario, ignore.case = TRUE) ~ "Retrospective",
      grepl("Pros", scenario, ignore.case = TRUE) ~ "Prospective",
      TRUE ~ "Other"
    )
  )

pr_df_susp <- pr_df_susp %>%
  mutate(
    cohort_type = case_when(
      grepl("Retro", scenario, ignore.case = TRUE) ~ "Retrospective",
      grepl("Pros", scenario, ignore.case = TRUE) ~ "Prospective",
      TRUE ~ "Other"
    )
  )

roc_plot_susp <- ggplot(roc_df_susp, aes(x = fpr, y = tpr, color = cohort_type)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_line(linewidth = 1) +
  facet_wrap(~ scenario, ncol = 2) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_color_manual(values = c("Retrospective" = "red", "Prospective" = "blue", "Other" = "black")) +
  labs(
    title = "ROC Curves (Suspected Infection Only)",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  geom_text(
    data = metrics_df_susp,
    aes(x = 0.35, y = 0.15, label = roc_label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 8),
    legend.position = "none"
  )

pr_plot_susp <- ggplot(pr_df_susp, aes(x = recall, y = precision, color = cohort_type)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ scenario, ncol = 2) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_color_manual(values = c("Retrospective" = "red", "Prospective" = "blue", "Other" = "black")) +
  labs(
    title = "Precision-Recall Curves (Suspected Infection Only)",
    x = "Recall",
    y = "Precision"
  ) +
  geom_text(
    data = metrics_df_susp,
    aes(x = 0.10, y = 0.15, label = pr_label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 8),
    legend.position = "none"
  )

print(roc_plot_susp)
print(pr_plot_susp)
metrics_df_susp



# Ensure correct scale for model scores
rescale_to_unit <- function(x) {
  if (is.numeric(x) && max(x, na.rm = TRUE) > 1) {
    x / 100
  } else {
    x
  }
}

retro_no_abx_1st_infxn <- retro_no_abx_final
retro_yes_abx_1st_infxn <- retro_yes_abx_final

retro_no_abx_1st_infxn$model_score <-
  rescale_to_unit(retro_no_abx_1st_infxn$model_score)

retro_yes_abx_1st_infxn$model_score <-
  rescale_to_unit(retro_yes_abx_1st_infxn$model_score)

pros_no_abx_1st_infxn$model_score <-
  rescale_to_unit(pros_no_abx_1st_infxn$model_score)

pros_yes_abx_1st_infxn$model_score <-
  rescale_to_unit(pros_yes_abx_1st_infxn$model_score)

# #Write pros datasets and retro sets to file
# write_csv(x = retro_no_abx_1st_infxn, "/phi/sbi/sbi_blake/retro_1st_no_abx_2_13_26.csv")
# write_csv(x = retro_yes_abx_1st_infxn, "/phi/sbi/sbi_blake/retro_1st_yes_abx_2_13_26.csv")
# write_csv(x = pros_no_abx_1st_infxn, "/phi/sbi/sbi_blake/pros_1st_no_abx_2_13_26.csv")
# write_csv(x = pros_yes_abx_1st_infxn, "/phi/sbi/sbi_blake/pros_1st_yes_abx_2_13_26.csv")


# ------------------------------------------------------------

# ----------------------------
# Inputs: my 4 evaluation dataframes
# ----------------------------
# ------------------------------------------------------------
# Build evaluation dfs with the RIGHT score source
# ------------------------------------------------------------
dfs_eval <- list(
  "Retro • Abx-" = retro_no_abx_1st_infxn %>%
    mutate(eval_score = model_score),

  "Retro • Abx+" = retro_yes_abx_1st_infxn %>%
    mutate(eval_score = model_score),

  "Pros • Abx-" = pros_no_abx_1st_infxn %>%
    mutate(eval_score = rf_pred_prob_pros[, "yes"]),

  "Pros • Abx+" = pros_yes_temp %>%
    mutate(eval_score = rf_pred_prob_pros_abx[, "yes"])
)

score_col <- "eval_score"
y_col     <- "sbi_present"
sus_col   <- "suspected_infection"

thr_tbl <- tibble(
  abx = c("Abx-", "Abx+"),
  threshold = c(0.05, 0.074)
)

parse_scenario <- function(s) {
  epoch <- ifelse(grepl("^Retro", s), "Retro", "Pros")
  abx   <- ifelse(grepl("Abx\\+", s), "Abx+", "Abx-")
  tibble(epoch = epoch, abx = abx)
}

to_truth01 <- function(x) {
  if (is.factor(x)) x <- as.character(x)

  if (is.character(x)) {
    x2 <- tolower(trimws(x))
    return(dplyr::case_when(
      x2 %in% c("1","yes","y","true","event") ~ 1L,
      x2 %in% c("0","no","n","false","none")  ~ 0L,
      TRUE ~ NA_integer_
    ))
  }

  x_num <- suppressWarnings(as.numeric(x))
  as.integer(ifelse(is.na(x_num), NA_integer_, ifelse(x_num >= 1, 1L, 0L)))
}

ruleout_stats_one <- function(df, scenario, score_col, y_col, threshold,
                              sus_col = "suspected_infection",
                              include_all_comers = TRUE) {
  sc <- parse_scenario(scenario)

  dat <- df %>%
    transmute(
      score = as.numeric(.data[[score_col]]),
      truth = to_truth01(.data[[y_col]]),
      suspected = if (sus_col %in% names(df)) to_truth01(.data[[sus_col]]) else NA_integer_
    ) %>%
    filter(!is.na(score), !is.na(truth), truth %in% c(0L, 1L))

  if (!include_all_comers) {
    dat <- dat %>% filter(suspected == 1L)
  }

  n_total <- nrow(dat)
  n_sbi0  <- sum(dat$truth == 0L)
  n_sbi1  <- sum(dat$truth == 1L)

  low <- dat$score <= threshold

  n_low        <- sum(low)
  n_low_sbi0   <- sum(low & dat$truth == 0L)
  n_low_sbi1   <- sum(low & dat$truth == 1L)

  prop_sbineg_low <- if (n_sbi0 > 0) n_low_sbi0 / n_sbi0 else NA_real_
  npv             <- if (n_low > 0) n_low_sbi0 / n_low else NA_real_

  tibble(
    epoch = sc$epoch,
    abx   = sc$abx,
    threshold = threshold,
    n_total = n_total,
    n_sbi0  = n_sbi0,
    n_sbi1  = n_sbi1,
    n_low      = n_low,
    n_low_sbi0 = n_low_sbi0,
    n_low_sbi1 = n_low_sbi1,
    x_prop_sbineg_low = prop_sbineg_low,
    y_npv             = npv,
    has_low           = n_low > 0
  )
}

cohort_tbl <- tibble(
  cohort = c("All encounters", "Suspected infection only"),
  include_all_comers = c(TRUE, FALSE)
)

quad_df <- tidyr::crossing(
  scenario = names(dfs_eval),
  cohort_tbl
) %>%
  mutate(
    stats = purrr::map2(
      scenario,
      include_all_comers,
      ~{
        abx_lbl <- ifelse(grepl("Abx\\+", .x), "Abx+", "Abx-")
        thr <- thr_tbl$threshold[thr_tbl$abx == abx_lbl][1]
        ruleout_stats_one(
          dfs_eval[[.x]],
          scenario = .x,
          score_col = score_col,
          y_col = y_col,
          threshold = thr,
          sus_col = sus_col,
          include_all_comers = .y
        )
      }
    )
  ) %>%
  select(-include_all_comers) %>%
  tidyr::unnest(stats)

quad_df <- quad_df %>%
  mutate(
    x_plot = ifelse(is.na(x_prop_sbineg_low), 0, x_prop_sbineg_low),
    y_plot = ifelse(is.na(y_npv), 1, y_npv),
    note   = dplyr::case_when(
      has_low ~ paste0(
        "NPV=", percent(y_npv, 0.1), "\n",
        "SBI- captured=", percent(x_prop_sbineg_low, 0.1), "\n",
        "low-risk n=", n_low, " (FN=", n_low_sbi1, ")"
      ),
      TRUE ~ paste0("No patients < threshold\nlow-risk n=0")
    )
  )

seg_df <- quad_df %>%
  select(cohort, abx, epoch, x_plot, y_plot) %>%
  tidyr::pivot_wider(names_from = epoch, values_from = c(x_plot, y_plot)) %>%
  filter(!is.na(x_plot_Retro), !is.na(x_plot_Pros)) %>%
  transmute(
    cohort = cohort,
    abx = abx,
    x = x_plot_Retro, y = y_plot_Retro,
    xend = x_plot_Pros, yend = y_plot_Pros
  )

# Optional: set y-zoom based on observed values but keep full data (zoom, don't drop)
y_lo <- max(0, floor(min(quad_df$y_plot, na.rm = TRUE) * 100) / 100 - 0.02)  # e.g., 0.759 -> 0.74
y_hi <- 1.01

p_quadrant <-
  ggplot(quad_df, aes(x = x_plot, y = y_plot, color = epoch)) +

  # connecting arrows Retro -> Pros
  geom_segment(
    data = seg_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    linewidth = 1,
    alpha = 0.6,
    color = "gray50"
  ) +

  # reference lines
  geom_hline(yintercept = 0.95, linetype = 2, alpha = 0.5) +

  # points: hollow + lower alpha when n_low==0
  geom_point(
    aes(shape = has_low, alpha = has_low),
    size = 4,
    stroke = 1.1
  ) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45), guide = "none") +

  scale_color_manual(
    values = c("Retro" = "red", "Pros" = "blue")
  ) +

  ggrepel::geom_text_repel(
    aes(label = note),
    size = 4.2,
    show.legend = FALSE,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.alpha = 0.5,
    nudge_y = -0.005,
    seed = 1
  ) +

  facet_grid(cohort ~ abx) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  coord_cartesian(ylim = c(y_lo, y_hi)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.5)) +

  labs(
    title = "Clinical rule-out performance at fixed SBI risk thresholds",
    x = "Proportion of SBI-negative encounters classified as low risk",
    y = "Negative predictive value (NPV)",
    color = "Epoch"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13, face = "bold"),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_quadrant

# Now plot the distribution of scores
# ----------------------------
# Inputs
# ----------------------------
# Build evaluation dfs using the CORRECT score source for each cohort
dfs_eval <- list(
  "Retro • Abx-" = retro_no_abx_1st_infxn %>%
    mutate(eval_score = model_score),

  "Retro • Abx+" = retro_yes_abx_1st_infxn %>%
    mutate(eval_score = model_score),

  "Pros • Abx-" = pros_no_abx_1st_infxn %>%
    mutate(eval_score = rf_pred_prob_pros[, "yes"]),

  "Pros • Abx+" = pros_yes_temp %>%
    mutate(eval_score = rf_pred_prob_pros_abx[, "yes"])
)

score_col <- "eval_score"
y_col     <- "sbi_present"
sus_col   <- "suspected_infection"

thr_tbl <- tibble(
  abx = c("Abx-", "Abx+"),
  threshold = c(0.05, 0.074)
)

parse_scenario <- function(s) {
  tibble(
    epoch = ifelse(grepl("^Retro", s), "Retro", "Pros"),
    abx   = ifelse(grepl("Abx\\+", s), "Abx+", "Abx-")
  )
}

# ----------------------------
# Build long dataframe of scores
# ----------------------------
score_df <- imap_dfr(dfs_eval, ~{
  sc <- parse_scenario(.y)
  tibble(
    epoch = sc$epoch,
    abx   = sc$abx,
    score = as.numeric(.x[[score_col]])
  )
}) %>%
  filter(is.finite(score), !is.na(score)) %>%
  mutate(
    epoch = factor(epoch, levels = c("Retro", "Pros")),
    abx   = factor(abx,   levels = c("Abx-", "Abx+"))
  ) %>%
  left_join(thr_tbl, by = c("abx"))

# % below threshold per panel for annotation
annot_df <- score_df %>%
  group_by(epoch, abx, threshold) %>%
  summarise(
    n = n(),
    pct_below = mean(score <= threshold),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "N = ", n, "\n",
      "% < threshold = ", percent(pct_below, accuracy = 0.1)
    )
  )

# ----------------------------
# Plot: histogram + density + vertical threshold line
# ----------------------------
p_dist <-
  ggplot2::ggplot(score_df, ggplot2::aes(x = score, fill = epoch, color = epoch)) +
  ggplot2::geom_histogram(
    ggplot2::aes(y = after_stat(density)),
    bins = 25,
    alpha = 0.20,
    linewidth = 0.2,
    position = "identity"
  ) +
  ggplot2::geom_vline(
    ggplot2::aes(xintercept = threshold),
    linetype = "dashed",
    linewidth = 1,
    color = "black"
  ) +
  ggplot2::geom_text(
    data = annot_df,
    ggplot2::aes(x = 0.98, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1.2,
    size = 3.8
  ) +
  ggplot2::facet_grid(epoch ~ abx) +
  ggplot2::scale_fill_manual(
    values = c("Retro" = "red", "Pros" = "blue")
  ) +
  ggplot2::scale_color_manual(
    values = c("Retro" = "red", "Pros" = "blue")
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(accuracy = 0.1)
  ) +
  ggplot2::labs(
    title = "Distribution of predicted SBI probabilities by epoch and antibiotic exposure",
    subtitle = "Dashed line = fixed rule-out threshold (Abx-: 0.05, Abx+: 0.074)",
    x = "Predicted SBI probability",
    y = "Density"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 16, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 14, face = "bold"),
    axis.text = ggplot2::element_text(size = 12),
    strip.text = ggplot2::element_text(size = 13, face = "bold"),
    legend.position = "none"
  )

p_dist

to_truth01 <- function(x) {
  if (is.factor(x)) x <- as.character(x)

  if (is.character(x)) {
    x2 <- tolower(trimws(x))
    out <- dplyr::case_when(
      x2 %in% c("1","yes","y","true","event") ~ 1L,
      x2 %in% c("0","no","n","false","none")  ~ 0L,
      TRUE ~ NA_integer_
    )
    return(out)
  }

  x_num <- suppressWarnings(as.numeric(x))
  as.integer(ifelse(is.na(x_num), NA_integer_, ifelse(x_num >= 1, 1L, 0L)))
}

score_truth_df <- imap_dfr(dfs_eval, ~{
  sc <- parse_scenario(.y)
  tibble(
    epoch     = sc$epoch,
    abx       = sc$abx,
    score     = as.numeric(.x[[score_col]]),
    truth_num = to_truth01(.x[[y_col]])
  )
}) %>%
  filter(is.finite(score), !is.na(score), truth_num %in% c(0L, 1L)) %>%
  mutate(
    epoch = factor(epoch, levels = c("Retro", "Pros")),
    abx   = factor(abx,   levels = c("Abx-", "Abx+")),
    truth = factor(ifelse(truth_num == 1L, "SBI+", "SBI-"), levels = c("SBI-", "SBI+"))
  ) %>%
  left_join(thr_tbl, by = c("abx"))

score_truth_df <- score_truth_df %>%
  mutate(panel = factor(
    paste(epoch, abx, sep = " • "),
    levels = c("Retro • Abx-", "Retro • Abx+", "Pros • Abx-", "Pros • Abx+")
  ))

p_dist_by_truth <-
  ggplot2::ggplot(score_truth_df, ggplot2::aes(x = score, fill = truth)) +
  ggplot2::geom_histogram(
    position = "identity",
    alpha = 0.4,
    bins = 30
  ) +
  ggplot2::geom_vline(
    data = dplyr::distinct(score_truth_df, panel, threshold),
    ggplot2::aes(xintercept = threshold),
    linetype = "dashed",
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  ggplot2::facet_wrap(~ panel, ncol = 2, scales = "free_y") +
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  ggplot2::scale_fill_manual(values = c("SBI-" = "#6E6E6E", "SBI+" = "#E69F00")) +
  ggplot2::labs(
    title = "Distribution of predicted SBI probability by epoch, antibiotic exposure, and true SBI status",
    subtitle = "Dashed line = fixed rule-out threshold (Abx-: 0.05, Abx+: 0.074)",
    x = "Predicted SBI probability",
    y = "Count",
    fill = "True outcome"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 16, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 14, face = "bold"),
    axis.text = ggplot2::element_text(size = 12),
    strip.text = ggplot2::element_text(size = 13, face = "bold"),
    legend.title = ggplot2::element_text(size = 12, face = "bold"),
    legend.text = ggplot2::element_text(size = 11),
    legend.position = "bottom"
  )

p_dist_by_truth


# Fix all of the race categories to be the same and match prospective (more condensed)
retro_no_abx_1st_infxn$race[retro_no_abx_1st_infxn$race == "Not Reported"] <- "Other_Or_Unknown"
retro_no_abx_1st_infxn$race[retro_no_abx_1st_infxn$race == "More than one Race"] <- "Other_Or_Unknown"
retro_no_abx_1st_infxn$race[retro_no_abx_1st_infxn$race == "Other"] <- "Other_Or_Unknown"
retro_no_abx_1st_infxn$race[retro_no_abx_1st_infxn$race == "Unknown/Not Reported"] <- "Other_Or_Unknown"
retro_no_abx_1st_infxn$race[retro_no_abx_1st_infxn$race == "White"] <- "White/Caucasian"

retro_yes_abx_1st_infxn$race[retro_yes_abx_1st_infxn$race == "Not Reported"] <- "Other_Or_Unknown"
retro_yes_abx_1st_infxn$race[retro_yes_abx_1st_infxn$race == "More than one Race"] <- "Other_Or_Unknown"
retro_yes_abx_1st_infxn$race[retro_yes_abx_1st_infxn$race == "Other"] <- "Other_Or_Unknown"
retro_yes_abx_1st_infxn$race[retro_yes_abx_1st_infxn$race == "Unknown/Not Reported"] <- "Other_Or_Unknown"
retro_yes_abx_1st_infxn$race[retro_yes_abx_1st_infxn$race == "White"] <- "White/Caucasian"

retro_no_abx_1st_infxn$race <- as.factor(retro_no_abx_1st_infxn$race)
retro_yes_abx_1st_infxn$race <- as.factor(retro_yes_abx_1st_infxn$race)
pros_no_abx_1st_infxn$race <- as.factor(pros_no_abx_1st_infxn$race)
pros_yes_abx_1st_infxn$race <- as.factor(pros_yes_abx_1st_infxn$race)


retro_no_abx_1st_infxn <- retro_no_abx_1st_infxn %>%
  mutate(
    ethnicity = as.character(ethnicity),
    ethnicity = case_when(
      ethnicity %in% c("Not Reported", "Unknown") ~ "Other_Or_Unknown",
      TRUE ~ ethnicity
    ),
    ethnicity = factor(ethnicity)
  )

retro_yes_abx_1st_infxn <- retro_yes_abx_1st_infxn %>%
  mutate(
    ethnicity = as.character(ethnicity),
    ethnicity = case_when(
      ethnicity %in% c("Not Reported", "Unknown") ~ "Other_Or_Unknown",
      TRUE ~ ethnicity
    ),
    ethnicity = factor(ethnicity)
  )


# ============================================================
# PROPENSITY SCORE DECOMPOSITION: MODELS 1, 2, 3
#   Model 1 = demographics/comorbidity
#   Model 2 = Model 1 + encounter context / pre-ICU location
#   Model 3 = Model 2 + physiology / lab distribution
# ============================================================

# ============================================================
# SAFE OUTCOME STANDARDIZATION + PS ANALYSIS
# ============================================================
# ------------------------------------------------------------
# Robust conversion of outcome to 0/1
# ------------------------------------------------------------
to_truth01 <- function(x) {
  if (is.factor(x)) x <- as.character(x)

  if (is.character(x)) {
    x2 <- tolower(trimws(x))
    return(dplyr::case_when(
      x2 %in% c("1", "yes", "y", "true", "event") ~ 1L,
      x2 %in% c("0", "no", "n", "false", "none")  ~ 0L,
      TRUE ~ NA_integer_
    ))
  }

  x_num <- suppressWarnings(as.numeric(x))
  as.integer(ifelse(is.na(x_num), NA_integer_, ifelse(x_num >= 1, 1L, 0L)))
}

# ------------------------------------------------------------
# Standardize epoch dataframe before prep
# ------------------------------------------------------------
standardize_epoch_input <- function(df) {
  df %>%
    mutate(
      sbi_present = to_truth01(sbi_present),
      is_female = as.factor(is_female),
      race = as.factor(race),
      ethnicity = as.factor(ethnicity),
      pccc = as.factor(pccc),
      malignancy_pccc = as.factor(malignancy_pccc)
    )
}

# ------------------------------------------------------------
# Existing prep function, slightly cleaned
# ------------------------------------------------------------
prep_epoch_df <- function(df, epoch_label, cutoff_prob, score_vec) {

  needs_pct_to_prob <- max(score_vec, na.rm = TRUE) > 1

  df %>%
    mutate(
      epoch = epoch_label,
      A = if_else(epoch == "pros", 1L, 0L),
      model_prob = if (needs_pct_to_prob) score_vec / 100 else score_vec,
      low_risk = if_else(!is.na(model_prob) & model_prob <= cutoff_prob, 1L, 0L)
    )
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
ps_exclude_vars <- function(df) {
  c(
    "study_id", "mrn", "picu_adm_date_time",
    "epoch", "A", "A_retro",
    "sbi_present",
    "model_score", "model_prob", "low_risk"
  )
}

w_mean <- function(x, w) {
  sum(x * w, na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE)
}

ruleout_metrics <- function(df, w_col = NULL) {

  w <- if (is.null(w_col)) rep(1, nrow(df)) else df[[w_col]]

  y  <- df$sbi_present
  lr <- df$low_risk

  ok <- !is.na(y) & !is.na(lr) & !is.na(w)
  y  <- y[ok]
  lr <- lr[ok]
  w  <- w[ok]

  p_low_all <- w_mean(lr == 1, w)

  neg <- (y == 0)
  p_low_among_neg <- if (sum(w[neg]) > 0) w_mean(lr[neg] == 1, w[neg]) else NA_real_

  lr1 <- (lr == 1)
  npv <- if (sum(w[lr1]) > 0) w_mean(y[lr1] == 0, w[lr1]) else NA_real_

  fn_rate <- if (sum(w) > 0) w_mean((y == 1) & (lr == 1), w) else NA_real_

  tibble(
    n = length(y),
    p_low_all = p_low_all,
    p_low_among_sbineg = p_low_among_neg,
    npv = npv,
    fn_rate = fn_rate,
    fn_per_1000 = 1000 * fn_rate
  )
}

safe_unweighted_auroc <- function(df, truth_col = "sbi_present", score_col = "model_prob") {
  d <- df %>% filter(!is.na(.data[[truth_col]]), !is.na(.data[[score_col]]))
  if (nrow(d) == 0 || length(unique(d[[truth_col]])) < 2) return(NA_real_)

  r <- pROC::roc(
    response = 1 - d[[truth_col]],
    predictor = 1 - d[[score_col]],
    quiet = TRUE,
    direction = "<"
  )
  as.numeric(pROC::auc(r))
}

safe_weighted_auroc <- function(df,
                                truth_col = "sbi_present",
                                score_col = "model_prob",
                                w_col = "w") {
  d <- df %>%
    dplyr::filter(
      !is.na(.data[[truth_col]]),
      !is.na(.data[[score_col]]),
      !is.na(.data[[w_col]]),
      .data[[w_col]] > 0
    ) %>%
    dplyr::mutate(
      truth01 = 1 - as.integer(.data[[truth_col]]),
      score_case = 1 - as.numeric(.data[[score_col]])
    )

  if (nrow(d) == 0 || length(unique(d$truth01)) < 2) {
    return(NA_real_)
  }

  roc_obj <- WeightedROC::WeightedROC(
    guess  = d$score_case,
    label  = d$truth01,
    weight = d[[w_col]]
  )

  as.numeric(WeightedROC::WeightedAUC(roc_obj))
}

safe_unweighted_auprc <- function(df, truth_col = "sbi_present", score_col = "model_prob") {
  d <- df %>% filter(!is.na(.data[[truth_col]]), !is.na(.data[[score_col]]))
  if (nrow(d) == 0 || length(unique(d[[truth_col]])) < 2) return(NA_real_)

  scores_pos <- (1 - d[[score_col]])[d[[truth_col]] == 0]
  scores_neg <- (1 - d[[score_col]])[d[[truth_col]] == 1]

  PRROC::pr.curve(
    scores.class0 = scores_pos,
    scores.class1 = scores_neg,
    curve = FALSE
  )$auc.integral
}

safe_weighted_auprc <- function(df, truth_col = "sbi_present", score_col = "model_prob", w_col = "w") {
  d <- df %>% filter(!is.na(.data[[truth_col]]), !is.na(.data[[score_col]]), !is.na(.data[[w_col]]))
  if (nrow(d) == 0 || length(unique(d[[truth_col]])) < 2) return(NA_real_)

  scores_pos <- (1 - d[[score_col]])[d[[truth_col]] == 0]
  scores_neg <- (1 - d[[score_col]])[d[[truth_col]] == 1]
  w_pos <- d[[w_col]][d[[truth_col]] == 0]
  w_neg <- d[[w_col]][d[[truth_col]] == 1]

  PRROC::pr.curve(
    scores.class0 = scores_pos,
    scores.class1 = scores_neg,
    weights.class0 = w_pos,
    weights.class1 = w_neg,
    curve = FALSE
  )$auc.integral
}

safe_unweighted_roc_df <- function(df,
                                   scenario,
                                   truth_col = "sbi_present",
                                   score_col = "model_prob") {
  d <- df %>%
    dplyr::filter(
      !is.na(.data[[truth_col]]),
      !is.na(.data[[score_col]])
    ) %>%
    dplyr::mutate(
      truth01 = 1 - as.integer(.data[[truth_col]]),
      score_case = 1 - as.numeric(.data[[score_col]])
    )

  if (nrow(d) == 0 || length(unique(d$truth01)) < 2) {
    return(tibble::tibble(
      scenario = character(),
      fpr = numeric(),
      tpr = numeric()
    ))
  }

  roc_obj <- pROC::roc(
    response  = d$truth01,
    predictor = d$score_case,
    quiet     = TRUE,
    direction = "<"
  )

  tibble::tibble(
    scenario = scenario,
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
}

safe_weighted_roc_df <- function(df,
                                 scenario,
                                 truth_col = "sbi_present",
                                 score_col = "model_prob",
                                 w_col = "w") {
  d <- df %>%
    dplyr::filter(
      !is.na(.data[[truth_col]]),
      !is.na(.data[[score_col]]),
      !is.na(.data[[w_col]]),
      .data[[w_col]] > 0
    ) %>%
    dplyr::mutate(
      truth01 = 1 - as.integer(.data[[truth_col]]),
      score_case = 1 - as.numeric(.data[[score_col]])
    )

  if (nrow(d) == 0 || length(unique(d$truth01)) < 2) {
    return(tibble::tibble(
      scenario = character(),
      fpr = numeric(),
      tpr = numeric()
    ))
  }

  roc_obj <- WeightedROC::WeightedROC(
    guess  = d$score_case,
    label  = d$truth01,
    weight = d[[w_col]]
  )

  tibble::tibble(
    scenario = scenario,
    fpr = roc_obj$FPR,
    tpr = roc_obj$TPR
  )
}

ess <- function(w) {
  (sum(w, na.rm = TRUE)^2) / sum(w^2, na.rm = TRUE)
}

# ------------------------------------------------------------
# PS models 1-3
# ------------------------------------------------------------
get_ps_covariates_m1 <- function(df) {
  cand <- get_shared_ps_candidates(df)

  intersect(
    c("age", "is_female", "race", "ethnicity", "pccc", "malignancy_pccc"),
    cand
  )
}

get_ps_covariates_m2 <- function(df) {
  cand <- get_shared_ps_candidates(df)

  base <- intersect(
    c("age", "is_female", "race", "ethnicity", "pccc", "malignancy_pccc",
      "los_before_icu_days", "imv_at_picu_adm"),
    cand
  )

  preicu_cols <- grep("^preicu_", cand, value = TRUE)

  unique(c(base, preicu_cols))
}

get_ps_covariates_m3 <- function(df) {
  cand <- get_shared_ps_candidates(df)

  base_context <- get_ps_covariates_m2(df)

  phys_like <- intersect(
    c(
      "hr_last", "hr_max", "hr_mean", "hr_median", "hr_min", "hr_slope",
      "sbp_last", "sbp_max", "sbp_mean", "sbp_median", "sbp_min", "sbp_slope",
      "dbp_last", "dbp_max", "dbp_mean", "dbp_median", "dbp_min", "dbp_slope",
      "rr_mean", "rr_median", "rr_min", "rr_slope",
      "temp_last", "temp_max", "temp_mean", "temp_median", "temp_min", "temp_slope",
      "o2sat_count", "o2sat_mean", "o2sat_median", "o2sat_min", "o2sat_slope",
      "fio2_count", "fio2_last", "fio2_max", "fio2_mean", "fio2_median", "fio2_min", "fio2_slope",
      "lactate_max", "hematocrit_mean", "pco2_mean"
    ),
    cand
  )

  unique(c(base_context, phys_like))
}

# ------------------------------------------------------------
# Weighting
# ------------------------------------------------------------
make_weights_pros_to_retro <- function(df,
                                       covariate_fun,
                                       estimand = "ATT",
                                       ps_trim = c(0.01, 0.99),
                                       weight_cap = 20,
                                       verbose = TRUE) {

  df2 <- df %>%
    mutate(A_retro = if_else(epoch == "retro", 1L, 0L))

  covars <- covariate_fun(df2)

  if (length(covars) == 0) stop("No PS covariates found for this model.")

  ps_form <- as.formula(paste("A_retro ~", paste(covars, collapse = " + ")))

  wfit <- WeightIt::weightit(
    formula = ps_form,
    data = df2,
    method = "ps",
    estimand = estimand
  )

  df_w <- df2 %>%
    mutate(
      ps = wfit$ps,
      w_raw = wfit$weights
    ) %>%
    filter(!is.na(ps), ps >= ps_trim[1], ps <= ps_trim[2]) %>%
    mutate(
      w = if (is.null(weight_cap)) w_raw else pmin(w_raw, weight_cap)
    )

  bal_final <- cobalt::bal.tab(
    x = ps_form,
    data = df_w,
    weights = df_w$w,
    estimand = estimand,
    un = TRUE,
    m.threshold = 0.1
  )

  if (verbose) {
    cat("\n==============================\n")
    cat("Covariates in PS model:\n")
    print(covars)
    cat("\nBalance on final weighted sample:\n")
    print(bal_final)
  }

  list(
    df_w = df_w,
    weightit = wfit,
    balance = bal_final,
    covariates = covars,
    formula = ps_form
  )
}

# ------------------------------------------------------------
# Summaries
# ------------------------------------------------------------
summarize_ps_run <- function(run_obj, stratum_label, model_label) {

  df_w <- run_obj$df_w

  retro_df <- df_w %>% filter(epoch == "retro", A == 0)
  pros_df  <- df_w %>% filter(epoch == "pros",  A == 1)

  ruleout_tbl <- bind_rows(
    retro_df %>% ruleout_metrics() %>% mutate(epoch = "Retro", weighting = "Unweighted"),
    pros_df  %>% ruleout_metrics() %>% mutate(epoch = "Pros",  weighting = "Unweighted"),
    pros_df  %>% ruleout_metrics(w_col = "w") %>% mutate(epoch = "Pros", weighting = "Weighted to Retro")
  ) %>%
    mutate(stratum = stratum_label, ps_model = model_label, .before = 1)

  discrim_tbl <- bind_rows(
    tibble(
      epoch = "Retro",
      weighting = "Unweighted",
      auroc = safe_unweighted_auroc(retro_df),
      auprc = safe_unweighted_auprc(retro_df)
    ),
    tibble(
      epoch = "Pros",
      weighting = "Unweighted",
      auroc = safe_unweighted_auroc(pros_df),
      auprc = safe_unweighted_auprc(pros_df)
    ),
    tibble(
      epoch = "Pros",
      weighting = "Weighted to Retro",
      auroc = safe_weighted_auroc(pros_df, w_col = "w"),
      auprc = safe_weighted_auprc(pros_df, w_col = "w")
    )
  ) %>%
    mutate(
      stratum = stratum_label,
      ps_model = model_label,
      ess_pros_weighted = ess(pros_df$w),
      .before = 1
    )

  list(
    ruleout = ruleout_tbl,
    discrimination = discrim_tbl,
    covariates = tibble(
      stratum = stratum_label,
      ps_model = model_label,
      covariate = run_obj$covariates
    ),
    balance = run_obj$balance,
    df_w = df_w
  )
}

run_ps_model_suite <- function(df, stratum_label,
                               ps_trim = c(0.01, 0.99),
                               weight_cap = 20) {

  model_defs <- list(
    "Model 1: Demographics/comorbidity"        = get_ps_covariates_m1,
    "Model 2: + Encounter context / pre-ICU"   = get_ps_covariates_m2,
    "Model 3: + Physiology / lab distribution" = get_ps_covariates_m3
  )

  runs <- imap(model_defs, function(cov_fun, model_label) {
    run_obj <- make_weights_pros_to_retro(
      df = df,
      covariate_fun = cov_fun,
      estimand = "ATT",
      ps_trim = ps_trim,
      weight_cap = weight_cap,
      verbose = TRUE
    )

    summarize_ps_run(
      run_obj = run_obj,
      stratum_label = stratum_label,
      model_label = model_label
    )
  })

  list(
    ruleout = bind_rows(map(runs, "ruleout")),
    discrimination = bind_rows(map(runs, "discrimination")),
    covariates = bind_rows(map(runs, "covariates")),
    full_runs = runs
  )
}

# ------------------------------------------------------------
# Variables that are actually available in BOTH epochs
# within the current stratum dataframe
# ------------------------------------------------------------
get_shared_ps_candidates <- function(df) {
  exclude <- ps_exclude_vars(df)

  cand <- setdiff(names(df), exclude)

  # keep only vars that are not all NA in retro and not all NA in pros
  keep <- cand[
    vapply(cand, function(v) {
      retro_ok <- "retro" %in% df$epoch && any(!is.na(df[df$epoch == "retro", v]))
      pros_ok  <- "pros"  %in% df$epoch && any(!is.na(df[df$epoch == "pros",  v]))
      retro_ok && pros_ok
    }, logical(1))
  ]

  keep
}

# ============================================================
# BUILD STRATA SAFELY
# ============================================================

retro_no_abx_std <- standardize_epoch_input(retro_no_abx_1st_infxn)
pros_no_abx_std  <- standardize_epoch_input(pros_no_abx_1st_infxn)

retro_yes_abx_std <- standardize_epoch_input(retro_yes_abx_1st_infxn)
pros_yes_abx_std  <- standardize_epoch_input(pros_yes_temp)   # swap object if needed

df_noabx <- bind_rows(
  prep_epoch_df(
    retro_no_abx_std,
    epoch_label = "retro",
    cutoff_prob = 0.05,
    score_vec = retro_no_abx_std$model_score
  ),
  prep_epoch_df(
    pros_no_abx_std,
    epoch_label = "pros",
    cutoff_prob = 0.05,
    score_vec = rf_pred_prob_pros[, "yes"]
  )
)

df_yesabx <- bind_rows(
  prep_epoch_df(
    retro_yes_abx_std,
    epoch_label = "retro",
    cutoff_prob = 0.074,
    score_vec = retro_yes_abx_std$model_score
  ),
  prep_epoch_df(
    pros_yes_abx_std,
    epoch_label = "pros",
    cutoff_prob = 0.074,
    score_vec = rf_pred_prob_pros_abx[, "yes"]
  )
)

# Sanity checks
df_noabx %>% count(epoch, sbi_present)
df_yesabx %>% count(epoch, sbi_present)

# ============================================================
# RUN MODELS
# ============================================================

ps_noabx_results <- run_ps_model_suite(
  df = df_noabx,
  stratum_label = "Abx-",
  ps_trim = c(0.01, 0.99),
  weight_cap = 20
)

ps_yesabx_results <- run_ps_model_suite(
  df = df_yesabx,
  stratum_label = "Abx+",
  ps_trim = c(0.01, 0.99),
  weight_cap = 20
)

# ============================================================
# COMBINED OUTPUTS
# ============================================================

ps_ruleout_results <- bind_rows(
  ps_noabx_results$ruleout,
  ps_yesabx_results$ruleout
) %>%
  arrange(stratum, ps_model, factor(weighting, levels = c("Unweighted", "Weighted to Retro")), epoch)

ps_discrimination_results <- bind_rows(
  ps_noabx_results$discrimination,
  ps_yesabx_results$discrimination
) %>%
  arrange(stratum, ps_model, factor(weighting, levels = c("Unweighted", "Weighted to Retro")), epoch)

ps_covariates_used <- bind_rows(
  ps_noabx_results$covariates,
  ps_yesabx_results$covariates
)

ps_ruleout_results
ps_discrimination_results
ps_covariates_used

# Optional: write out
# readr::write_csv(ps_ruleout_results, "ps_ruleout_results_models_1_2_3.csv")
# readr::write_csv(ps_discrimination_results, "ps_discrimination_results_models_1_2_3.csv")
# readr::write_csv(ps_covariates_used, "ps_covariates_used_models_1_2_3.csv")


#### Plot weighted and unweighted pros charts
# ============================================================
# PLOT RETRO (UNWEIGHTED), PROS (UNWEIGHTED), AND
# PROS (WEIGHTED TO RETRO) FOR EACH STRATUM / MODEL
# ============================================================
# ----------------------------
# Helper: ROC points + AUROC
# Uses WeightedROC for both weighted and unweighted cases
# ----------------------------
roc_points_weighted <- function(df,
                                truth_col = "sbi_present",
                                score_col = "model_prob",
                                w_col = NULL) {
  d <- df %>%
    filter(!is.na(.data[[truth_col]]), !is.na(.data[[score_col]])) %>%
    mutate(
      y = 1 - as.integer(.data[[truth_col]]),
      s = 1 - as.numeric(.data[[score_col]]),
      w = if (is.null(w_col)) 1 else as.numeric(.data[[w_col]])
    ) %>%
    filter(is.finite(s), !is.na(w), w > 0)

  if (nrow(d) == 0 || length(unique(d$y)) < 2) {
    return(tibble(
      fpr = numeric(),
      tpr = numeric(),
      auc = numeric()
    ))
  }

  roc_obj <- WeightedROC::WeightedROC(
    guess  = d$s,
    label  = d$y,
    weight = d$w
  )

  tibble(
    fpr = roc_obj$FPR,
    tpr = roc_obj$TPR,
    auc = WeightedROC::WeightedAUC(roc_obj)
  )
}

# ----------------------------
# Helper: PR points + AUPRC
# ----------------------------
pr_points_weighted <- function(df,
                               truth_col = "sbi_present",
                               score_col = "model_prob",
                               w_col = NULL) {
  d <- df %>%
    filter(!is.na(.data[[truth_col]]), !is.na(.data[[score_col]])) %>%
    mutate(
      y = 1 - as.integer(.data[[truth_col]]),
      s = 1 - as.numeric(.data[[score_col]]),
      w = if (is.null(w_col)) 1 else as.numeric(.data[[w_col]])
    ) %>%
    filter(is.finite(s), !is.na(w), w > 0)

  if (nrow(d) == 0 || length(unique(d$y)) < 2) {
    return(tibble(
      recall = numeric(),
      precision = numeric(),
      auc = numeric()
    ))
  }

  scores_pos <- d$s[d$y == 1]
  scores_neg <- d$s[d$y == 0]
  w_pos <- d$w[d$y == 1]
  w_neg <- d$w[d$y == 0]

  pr_obj <- PRROC::pr.curve(
    scores.class0 = scores_pos,
    scores.class1 = scores_neg,
    weights.class0 = w_pos,
    weights.class1 = w_neg,
    curve = TRUE
  )

  tibble(
    recall = pr_obj$curve[, 1],
    precision = pr_obj$curve[, 2],
    auc = pr_obj$auc.integral
  )
}

curve_prevalence <- function(df,
                             truth_col = "sbi_present",
                             w_col = NULL) {
  d <- df %>%
    dplyr::filter(!is.na(.data[[truth_col]]))

  if (nrow(d) == 0) return(NA_real_)

  if (is.null(w_col)) {
    mean(as.integer(d[[truth_col]]) == 0L)
  } else {
    w <- d[[w_col]]
    y <- as.integer(d[[truth_col]]) == 0L
    sum(w * y, na.rm = TRUE) / sum(w, na.rm = TRUE)
  }
}


# ----------------------------
# One run -> ROC + PR plot
# run_obj should be something like:
# ps_noabx_results$full_runs[["Model 1: Demographics/comorbidity"]]
# ----------------------------
plot_ps_curves_one_run <- function(run_obj, title_suffix = "") {

  df_w <- run_obj$df_w

  retro_df <- df_w %>% filter(epoch == "retro", A == 0)
  pros_df  <- df_w %>% filter(epoch == "pros",  A == 1)

  # ROC data
  roc_retro <- roc_points_weighted(retro_df, w_col = NULL) %>%
    mutate(curve = "Retro (unweighted)")
  roc_pros_unw <- roc_points_weighted(pros_df, w_col = NULL) %>%
    mutate(curve = "Pros (unweighted)")
  roc_pros_w <- roc_points_weighted(pros_df, w_col = "w") %>%
    mutate(curve = "Pros (weighted to Retro)")

  roc_df <- bind_rows(roc_retro, roc_pros_unw, roc_pros_w) %>%
    mutate(
      curve = factor(
        curve,
        levels = c("Retro (unweighted)", "Pros (unweighted)", "Pros (weighted to Retro)")
      )
    )

  # PR data
  pr_retro <- pr_points_weighted(retro_df, w_col = NULL) %>%
    mutate(curve = "Retro (unweighted)")
  pr_pros_unw <- pr_points_weighted(pros_df, w_col = NULL) %>%
    mutate(curve = "Pros (unweighted)")
  pr_pros_w <- pr_points_weighted(pros_df, w_col = "w") %>%
    mutate(curve = "Pros (weighted to Retro)")

  pr_df <- bind_rows(pr_retro, pr_pros_unw, pr_pros_w) %>%
    mutate(
      curve = factor(
        curve,
        levels = c("Retro (unweighted)", "Pros (unweighted)", "Pros (weighted to Retro)")
      )
    )

  # Labels
  roc_lab <- bind_rows(
    roc_retro %>% summarise(auc = first(auc)) %>% mutate(curve = "Retro (unweighted)"),
    roc_pros_unw %>% summarise(auc = first(auc)) %>% mutate(curve = "Pros (unweighted)"),
    roc_pros_w %>% summarise(auc = first(auc)) %>% mutate(curve = "Pros (weighted to Retro)")
  ) %>%
    mutate(
      curve = factor(curve, levels = c("Retro (unweighted)", "Pros (unweighted)", "Pros (weighted to Retro)")),
      x = 0.25,
      y = c(0.16, 0.10, 0.04),
      label = sprintf("%s AUROC = %.3f", curve, auc)
    )

  pr_lab <- bind_rows(
    tibble(
      curve = "Retro (unweighted)",
      auc = dplyr::first(pr_retro$auc),
      prev = curve_prevalence(retro_df, w_col = NULL)
    ),
    tibble(
      curve = "Pros (unweighted)",
      auc = dplyr::first(pr_pros_unw$auc),
      prev = curve_prevalence(pros_df, w_col = NULL)
    ),
    tibble(
      curve = "Pros (weighted to Retro)",
      auc = dplyr::first(pr_pros_w$auc),
      prev = curve_prevalence(pros_df, w_col = "w")
    )
  ) %>%
    mutate(
      curve = factor(
        curve,
        levels = c("Retro (unweighted)", "Pros (unweighted)", "Pros (weighted to Retro)")
      ),
      x = 0.2,
      y = c(0.18, 0.12, 0.06),
      label = paste0(
        "AUPRC = ", sprintf("%.3f", auc),
        " (prev = ", sprintf("%.2f", prev), ")"
      )
  ) %>%
    mutate(
      curve = factor(
        curve,
        levels = c("Retro (unweighted)", "Pros (unweighted)", "Pros (weighted to Retro)")
      ),
      legend_label = paste0(
        as.character(curve),
        ": AUPRC = ",
        sprintf("%.3f", auc),
        " (prev = ",
        sprintf("%.2f", prev),
        ")"
      )
    )

  col_vals <- c(
    "Retro (unweighted)"      = "red",
    "Pros (unweighted)"       = "blue",
    "Pros (weighted to Retro)" = "purple"
  )

  lt_vals <- c(
    "Retro (unweighted)"      = "solid",
    "Pros (unweighted)"       = "solid",
    "Pros (weighted to Retro)" = "dashed"
  )

  p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr, color = curve, linetype = curve)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    geom_line(linewidth = 1.1, na.rm = TRUE) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_color_manual(values = col_vals, drop = FALSE) +
    scale_linetype_manual(values = lt_vals, drop = FALSE) +
    labs(
      title = paste0("ROC (", title_suffix, ")"),
      x = "False positive rate",
      y = "True positive rate",
      color = "Curve",
      linetype = "Curve"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    geom_text(
      data = roc_lab,
      aes(x = x, y = y, label = label, color = curve),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3.4,
      show.legend = FALSE
    )

  pr_breaks <- c(
    "Retro (unweighted)",
    "Pros (unweighted)",
    "Pros (weighted to Retro)"
  )

  pr_label_map <- setNames(pr_lab$legend_label, as.character(pr_lab$curve))

  p_pr <- ggplot(pr_df, aes(x = recall, y = precision, color = curve, linetype = curve)) +
    geom_line(linewidth = 1.1, na.rm = TRUE) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_color_manual(
      values = col_vals,
      breaks = c("Retro (unweighted)", "Pros (unweighted)", "Pros (weighted to Retro)"),
      labels = c("Retro", "Pros", "Pros weighted"),
      drop = FALSE
    ) +
    scale_linetype_manual(
      values = lt_vals,
      breaks = c("Retro (unweighted)", "Pros (unweighted)", "Pros (weighted to Retro)"),
      labels = c("Retro", "Pros", "Pros weighted"),
      drop = FALSE
    ) +
    labs(
      title = paste0("PR (", title_suffix, ")"),
      x = "Recall",
      y = "Precision",
      color = "Curve",
      linetype = "Curve"
    ) +
    geom_text(
      data = pr_lab,
      aes(x = x, y = y, label = label, color = curve),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3.4,
      show.legend = FALSE
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  p_roc | p_pr
}

# ============================================================
# BUILD SIX PLOTS: 2 STRATA x 3 MODELS
# ============================================================

p_noabx_m1 <- plot_ps_curves_one_run(
  ps_noabx_results$full_runs[["Model 1: Demographics/comorbidity"]],
  title_suffix = "Abx- | Model 1"
)

p_noabx_m2 <- plot_ps_curves_one_run(
  ps_noabx_results$full_runs[["Model 2: + Encounter context / pre-ICU"]],
  title_suffix = "Abx- | Model 2"
)

p_noabx_m3 <- plot_ps_curves_one_run(
  ps_noabx_results$full_runs[["Model 3: + Physiology / lab distribution"]],
  title_suffix = "Abx- | Model 3"
)

p_yesabx_m1 <- plot_ps_curves_one_run(
  ps_yesabx_results$full_runs[["Model 1: Demographics/comorbidity"]],
  title_suffix = "Abx+ | Model 1"
)

p_yesabx_m2 <- plot_ps_curves_one_run(
  ps_yesabx_results$full_runs[["Model 2: + Encounter context / pre-ICU"]],
  title_suffix = "Abx+ | Model 2"
)

p_yesabx_m3 <- plot_ps_curves_one_run(
  ps_yesabx_results$full_runs[["Model 3: + Physiology / lab distribution"]],
  title_suffix = "Abx+ | Model 3"
)

# Print individually
p_noabx_m1
p_noabx_m2
p_noabx_m3
p_yesabx_m1
p_yesabx_m2
p_yesabx_m3


####### Do some  Table 1 creation ##########
# -----------------------------
# 1) Helper: collapse preicu_* one-hot into a single factor
# -----------------------------
add_preicu_location <- function(df, prefix = "preicu_") {

  pre_cols <- grep(paste0("^", prefix), names(df), value = TRUE)

  if (length(pre_cols) == 0) {
    df$preicu_location <- NA_character_
    return(df)
  }

  # Ensure numeric 0/1
  df <- df %>%
    mutate(across(all_of(pre_cols), ~ as.numeric(.)))

  # For each row, find which preicu_* column is "1"
  loc <- apply(df[, pre_cols, drop = FALSE], 1, function(x) {
    idx <- which(x == 1)
    if (length(idx) == 0) return(NA_character_)          # none flagged
    if (length(idx) > 1) return("multiple conflicting")  # >1 flagged
    pre_cols[idx]
  })

  df %>%
    mutate(
      preicu_location = loc,
      preicu_location = str_remove(preicu_location, paste0("^", prefix)),
      preicu_location = str_replace_all(preicu_location, "_", " "),
      preicu_location = str_to_upper(preicu_location),
      preicu_location = fct_explicit_na(as.factor(preicu_location), na_level = "MISSING")
    )
}

# -----------------------------
# 2) Helper: build Table 1 for one stratum dataframe
# -----------------------------
make_table1_epoch <- function(df, title = NULL) {

  # Add preicu_location
  df <- add_preicu_location(df)

  # Make sure epoch is a factor in the correct order
  df <- df %>%
    mutate(
      epoch = factor(epoch, levels = c("retro", "pros")),
      age = as.numeric(age),

      # binary indicators coded 0/1 for one-line display
      is_female = case_when(
        as.character(is_female) %in% c("1", "TRUE", "T", "Female", "FEMALE") ~ 1,
        as.character(is_female) %in% c("0", "FALSE", "F", "Male", "MALE") ~ 0,
        TRUE ~ NA_real_
      ),

      pccc = case_when(
        as.character(pccc) %in% c("1", "TRUE", "T", "Yes", "YES") ~ 1,
        as.character(pccc) %in% c("0", "FALSE", "F", "No", "NO") ~ 0,
        TRUE ~ NA_real_
      ),

      malignancy_pccc = case_when(
        as.character(malignancy_pccc) %in% c("1", "TRUE", "T", "Yes", "YES") ~ 1,
        as.character(malignancy_pccc) %in% c("0", "FALSE", "F", "No", "NO") ~ 0,
        TRUE ~ NA_real_
      ),

      race = as.character(race),
      race = str_replace_all(race, "Other_Or_Unknown", "Other Or Unknown"),
      race = str_replace_all(race, "_", " "),
      race = str_to_upper(race),
      race = fct_explicit_na(as.factor(race), na_level = "MISSING"),

      ethnicity = as.character(ethnicity),
      ethnicity = str_replace_all(ethnicity, "Other_Or_Unknown", "Other Or Unknown"),
      ethnicity = str_replace_all(ethnicity, "_", " "),
      ethnicity = str_to_upper(ethnicity),
      ethnicity = fct_explicit_na(as.factor(ethnicity), na_level = "MISSING")
    )

  # Variables to include
  vars <- c(
    "age",
    "is_female",
    "race",
    "ethnicity",
    "pccc",
    "malignancy_pccc",
    "los_before_icu_days",
    "imv_at_picu_adm",
    "preicu_location"
  )

  # Keep only vars that exist in df
  vars <- intersect(vars, names(df))

  tab <- df %>%
    dplyr::select(all_of(c("epoch", vars))) %>%
    tbl_summary(
      by = epoch,
      type = list(
        age ~ "continuous",
        los_before_icu_days ~ "continuous",
        is_female ~ "dichotomous",
        pccc ~ "dichotomous",
        malignancy_pccc ~ "dichotomous",
        imv_at_picu_adm ~ "dichotomous",
        race ~ "categorical",
        ethnicity ~ "categorical",
        preicu_location ~ "categorical"
      ),
      value = list(
        is_female ~ 1,
        pccc ~ 1,
        malignancy_pccc ~ 1,
        imv_at_picu_adm ~ 1
      ),
      statistic = list(
        all_continuous() ~ "{median} [{p25}, {p75}]",
        all_categorical() ~ "{n} ({p}%)",
        all_dichotomous() ~ "{n} ({p}%)"
      ),
      digits = all_continuous() ~ 2,
      missing = "ifany",
      label = list(
        age ~ "AGE",
        is_female ~ "FEMALE",
        race ~ "RACE",
        ethnicity ~ "ETHNICITY",
        pccc ~ "PCCC",
        malignancy_pccc ~ "MALIGNANCY",
        los_before_icu_days ~ "LOS BEFORE ICU DAYS",
        imv_at_picu_adm ~ "IMV AT PICU ADM",
        preicu_location ~ "PREICU LOCATION"
      )
    ) %>%
    modify_header(label ~ "**CHARACTERISTIC**") %>%
    modify_spanning_header(all_stat_cols() ~ "**EPOCH**") %>%
    bold_labels()

  if (!is.null(title)) {
    tab <- tab %>% modify_caption(paste0("**", title, "**"))
  }

  tab
}

# -----------------------------
# 4) Helper: create retro vs prospective suspected infection subset table
# -----------------------------
make_table1_retro_vs_pros_si <- function(df, title = NULL) {

  # Add preicu_location
  df <- add_preicu_location(df)

  # Keep all retro, and only prospective patients with suspected infection == 1
  df <- df %>%
    mutate(
      suspected_infection = case_when(
        as.character(suspected_infection) %in% c("1", "TRUE", "T") ~ 1,
        as.character(suspected_infection) %in% c("0", "FALSE", "F") ~ 0,
        TRUE ~ NA_real_
      ),
      epoch = factor(epoch, levels = c("retro", "pros"))
    ) %>%
    filter(
      epoch == "retro" | (epoch == "pros" & suspected_infection == 1)
    ) %>%
    mutate(
      comparison_group = case_when(
        epoch == "retro" ~ "RETRO",
        epoch == "pros" & suspected_infection == 1 ~ "PROSPECTIVE, SUSPECTED INFECTION",
        TRUE ~ NA_character_
      ),
      comparison_group = factor(
        comparison_group,
        levels = c("RETRO", "PROSPECTIVE, SUSPECTED INFECTION")
      ),

      age = as.numeric(age),

      # binary indicators coded 0/1 for one-line display
      is_female = case_when(
        as.character(is_female) %in% c("1", "TRUE", "T", "Female", "FEMALE") ~ 1,
        as.character(is_female) %in% c("0", "FALSE", "F", "Male", "MALE") ~ 0,
        TRUE ~ NA_real_
      ),

      pccc = case_when(
        as.character(pccc) %in% c("1", "TRUE", "T", "Yes", "YES") ~ 1,
        as.character(pccc) %in% c("0", "FALSE", "F", "No", "NO") ~ 0,
        TRUE ~ NA_real_
      ),

      malignancy_pccc = case_when(
        as.character(malignancy_pccc) %in% c("1", "TRUE", "T", "Yes", "YES") ~ 1,
        as.character(malignancy_pccc) %in% c("0", "FALSE", "F", "No", "NO") ~ 0,
        TRUE ~ NA_real_
      ),

      imv_at_picu_adm = case_when(
        as.character(imv_at_picu_adm) %in% c("1", "TRUE", "T", "Yes", "YES") ~ 1,
        as.character(imv_at_picu_adm) %in% c("0", "FALSE", "F", "No", "NO") ~ 0,
        TRUE ~ suppressWarnings(as.numeric(as.character(imv_at_picu_adm)))
      ),

      race = as.character(race),
      race = str_replace_all(race, "Other_Or_Unknown", "Other Or Unknown"),
      race = str_replace_all(race, "_", " "),
      race = str_to_upper(race),
      race = fct_explicit_na(as.factor(race), na_level = "MISSING"),

      ethnicity = as.character(ethnicity),
      ethnicity = str_replace_all(ethnicity, "Other_Or_Unknown", "Other Or Unknown"),
      ethnicity = str_replace_all(ethnicity, "_", " "),
      ethnicity = str_to_upper(ethnicity),
      ethnicity = fct_explicit_na(as.factor(ethnicity), na_level = "MISSING")
    )

  # Variables to include
  vars <- c(
    "age",
    "is_female",
    "race",
    "ethnicity",
    "pccc",
    "malignancy_pccc",
    "los_before_icu_days",
    "imv_at_picu_adm",
    "preicu_location"
  )

  vars <- intersect(vars, names(df))

  tab <- df %>%
    dplyr::select(all_of(c("comparison_group", vars))) %>%
    tbl_summary(
      by = comparison_group,
      type = list(
        age ~ "continuous",
        los_before_icu_days ~ "continuous",
        is_female ~ "dichotomous",
        pccc ~ "dichotomous",
        malignancy_pccc ~ "dichotomous",
        imv_at_picu_adm ~ "dichotomous",
        race ~ "categorical",
        ethnicity ~ "categorical",
        preicu_location ~ "categorical"
      ),
      value = list(
        is_female ~ 1,
        pccc ~ 1,
        malignancy_pccc ~ 1,
        imv_at_picu_adm ~ 1
      ),
      statistic = list(
        all_continuous() ~ "{median} [{p25}, {p75}]",
        all_categorical() ~ "{n} ({p}%)",
        all_dichotomous() ~ "{n} ({p}%)"
      ),
      digits = all_continuous() ~ 2,
      missing = "ifany",
      label = list(
        age ~ "AGE",
        is_female ~ "FEMALE",
        race ~ "RACE",
        ethnicity ~ "ETHNICITY",
        pccc ~ "PCCC",
        malignancy_pccc ~ "MALIGNANCY",
        los_before_icu_days ~ "LOS BEFORE ICU DAYS",
        imv_at_picu_adm ~ "IMV AT PICU ADM",
        preicu_location ~ "PREICU LOCATION"
      )
    ) %>%
    modify_header(label ~ "**CHARACTERISTIC**") %>%
    modify_spanning_header(all_stat_cols() ~ "**GROUP**") %>%
    bold_labels()

  if (!is.null(title)) {
    tab <- tab %>% modify_caption(paste0("**", title, "**"))
  }

  tab
}

# -----------------------------
# 3) Create the two Table 1s
# -----------------------------
table1_noabx  <- make_table1_epoch(
  df_noabx,
  title = "Table 1. Baseline characteristics by epoch — No antibiotics stratum"
)

table1_yesabx <- make_table1_epoch(
  df_yesabx,
  title = "Table 1. Baseline characteristics by epoch — Antibiotics stratum"
)

table1_noabx
table1_yesabx

# Now make comparisons with just suspicion of infection:
table1_noabx_si <- make_table1_retro_vs_pros_si(
  df_noabx,
  title = "Table 1. Baseline characteristics: retrospective cohort vs prospective suspected infection subset — No antibiotics stratum"
)

table1_yesabx_si <- make_table1_retro_vs_pros_si(
  df_yesabx,
  title = "Table 1. Baseline characteristics: retrospective cohort vs prospective suspected infection subset — Antibiotics stratum"
)

table1_noabx_si
table1_yesabx_si

#### Ok now we look at missingness and distribution of key inputs between retro and pros.
# =============================================================================
# Dataset shift analysis (retro vs prospective), stratified by abx exposure
# NOW:
#   - DOES NOT write CSVs or save PNGs (those lines are commented out)
#   - Returns objects into the R environment
#   - Prints the shift tables and prints the plots to the active device
# =============================================================================

# ----------------------------
# USER SETTINGS
# ----------------------------

EXCLUDE_COLS <- c(
  "study_id", "mrn", "picu_adm_date_time",
  "sbi_present", "model_score", "epoch",
  "adm_month", "adm_quarter", "adm_year", "adm_ym", "adm_yq",
  "pna_1_0", "blood", "cns", "urine", "ever_cx_neg_sepsis",
  "suspected_infection"
)

# Fix some column names
retro_yes_abx_1st_infxn <- retro_yes_abx_1st_infxn %>% rename(n_rr_rows = n_RR_rows)
retro_yes_abx_1st_infxn <- retro_yes_abx_1st_infxn %>% rename(n_o2sat_rows = n_spo2_rows)

pros_yes_abx_1st_infxn <- pros_yes_abx_1st_infxn %>% rename(n_rr_rows = n_RR_rows)
pros_yes_abx_1st_infxn <- pros_yes_abx_1st_infxn %>% rename(n_o2sat_rows = n_spo2_rows)


FORCE_CATEGORICAL <- c("race", "ethnicity", "is_female")

TOP_N_SMD_PLOT        <- 30
TOP_N_AVAIL_PLOT      <- 30
TOP_N_DENSITY_FEATURE <- 6

# ----------------------------
# Helper: pretty variable labels
# ----------------------------
pretty_var_label <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\bdbp\\b", "DBP") %>%
    str_replace_all("\\bsbp\\b", "SBP") %>%
    str_replace_all("\\bhr\\b", "HR") %>%
    str_replace_all("\\brr\\b", "RR") %>%
    str_replace_all("\\bfio2\\b", "FiO2") %>%
    str_replace_all("\\bo2sat\\b", "O2 Sat") %>%
    str_replace_all("\\bpco2\\b", "PCO2") %>%
    str_replace_all("\\bpicu\\b", "PICU") %>%
    str_replace_all("\\bimv\\b", "IMV") %>%
    str_replace_all("\\bpccc\\b", "PCCC") %>%
    str_replace_all("\\bnocer\\b", "NoCER") %>%
    str_replace_all("\\bpreicu\\b", "Pre-ICU") %>%
    str_replace_all("\\by n\\b", "Y/N") %>%
    str_squish() %>%
    str_to_title()
}

# ----------------------------
# Helper: safe quantiles
# ----------------------------
.safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
}

# ----------------------------
# Helper: classify variable type
# ----------------------------
.var_type <- function(x, varname) {
  if (varname %in% FORCE_CATEGORICAL) return("categorical")
  if (str_detect(varname, "(_pres|_present)$")) return("binary")
  if (is.factor(x) || is.character(x)) return("categorical")

  if (is.numeric(x) || is.integer(x)) {
    ux <- unique(x[!is.na(x)])
    if (length(ux) > 0 && all(ux %in% c(0, 1))) return("binary")
    return("continuous")
  }

  "categorical"
}

# ----------------------------
# Availability / "not missing" rate
# ----------------------------
.availability_rate <- function(df, varname) {
  x <- df[[varname]]
  if (str_detect(varname, "(_pres|_present)$")) {
    mean(x == 1, na.rm = TRUE)
  } else {
    1 - mean(is.na(x))
  }
}

# ----------------------------
# Continuous summaries
# ----------------------------
.summarize_continuous <- function(x) {
  qs <- .safe_quantile(x, c(0.25, 0.50, 0.75))
  m  <- ifelse(all(is.na(x)), NA_real_, mean(x, na.rm = TRUE))
  s  <- ifelse(all(is.na(x)), NA_real_, stats::sd(x, na.rm = TRUE))
  list(median = qs[2], q1 = qs[1], q3 = qs[3], mean = m, sd = s)
}

# ----------------------------
# SMD + KS
# ----------------------------
.smd_continuous <- function(m0, s0, m1, s1) {
  if (any(is.na(c(m0, s0, m1, s1))) || (s0 == 0 && s1 == 0)) return(NA_real_)
  sp <- sqrt((s0^2 + s1^2) / 2)
  if (sp == 0) return(NA_real_)
  (m1 - m0) / sp
}

.smd_binary <- function(p0, p1) {
  if (any(is.na(c(p0, p1)))) return(NA_real_)
  p_pool <- (p0 + p1) / 2
  denom <- sqrt(p_pool * (1 - p_pool))
  if (denom == 0) return(NA_real_)
  (p1 - p0) / denom
}

.smd_categorical <- function(x0, x1) {
  levs <- union(unique(x0[!is.na(x0)]), unique(x1[!is.na(x1)]))
  if (length(levs) < 2) return(NA_real_)

  p0 <- map_dbl(levs, ~mean(x0 == .x, na.rm = TRUE))
  p1 <- map_dbl(levs, ~mean(x1 == .x, na.rm = TRUE))
  p_pool <- (p0 + p1) / 2

  keep <- p_pool > 0
  if (!any(keep)) return(NA_real_)

  sqrt(sum(((p1 - p0)[keep]^2) / p_pool[keep], na.rm = TRUE))
}

.ks_stat <- function(x0, x1) {
  x0 <- x0[is.finite(x0)]
  x1 <- x1[is.finite(x1)]
  if (length(x0) < 5 || length(x1) < 5) return(NA_real_)
  as.numeric(stats::ks.test(x0, x1)$statistic)
}

# ----------------------------
# Build shift table for one stratum
# ----------------------------
build_shift_table <- function(df_retro, df_pros, stratum_name) {

  pred_cols <- setdiff(intersect(names(df_retro), names(df_pros)), EXCLUDE_COLS)

  df_retro <- df_retro[, pred_cols, drop = FALSE]
  df_pros  <- df_pros[,  pred_cols, drop = FALSE]

  out <- map_dfr(pred_cols, function(v) {
    x0 <- df_retro[[v]]
    x1 <- df_pros[[v]]

    vtype <- .var_type(x0, v)

    avail0 <- .availability_rate(df_retro, v)
    avail1 <- .availability_rate(df_pros,  v)

    if (vtype == "continuous") {
      s0 <- .summarize_continuous(x0)
      s1 <- .summarize_continuous(x1)

      retro_summary <- sprintf("%.3f [%.3f, %.3f]", s0$median, s0$q1, s0$q3)
      pros_summary  <- sprintf("%.3f [%.3f, %.3f]", s1$median, s1$q1, s1$q3)

      smd <- .smd_continuous(s0$mean, s0$sd, s1$mean, s1$sd)
      ks  <- .ks_stat(x0, x1)

      tibble(
        stratum       = stratum_name,
        variable      = v,
        variable_plot = pretty_var_label(v),
        type          = vtype,
        retro_summary = retro_summary,
        pros_summary  = pros_summary,
        retro_avail   = avail0,
        pros_avail    = avail1,
        avail_diff    = (avail1 - avail0),
        smd           = smd,
        ks_stat       = ks
      )
    } else if (vtype == "binary") {
      p0 <- mean(x0 == 1, na.rm = TRUE)
      p1 <- mean(x1 == 1, na.rm = TRUE)

      n0 <- sum(!is.na(x0))
      n1 <- sum(!is.na(x1))

      retro_summary <- sprintf("%d (%.1f%%)", round(p0*n0), 100*p0)
      pros_summary  <- sprintf("%d (%.1f%%)", round(p1*n1), 100*p1)

      smd <- .smd_binary(p0, p1)

      tibble(
        stratum       = stratum_name,
        variable      = v,
        variable_plot = pretty_var_label(v),
        type          = vtype,
        retro_summary = retro_summary,
        pros_summary  = pros_summary,
        retro_avail   = avail0,
        pros_avail    = avail1,
        avail_diff    = (avail1 - avail0),
        smd           = smd,
        ks_stat       = NA_real_
      )
    } else {
      smd <- .smd_categorical(as.character(x0), as.character(x1))

      top0 <- names(sort(table(x0), decreasing = TRUE))[1]
      top1 <- names(sort(table(x1), decreasing = TRUE))[1]
      p0 <- mean(as.character(x0) == top0, na.rm = TRUE)
      p1 <- mean(as.character(x1) == top1, na.rm = TRUE)

      retro_summary <- if (!is.na(top0)) sprintf("%s: %.1f%%", top0, 100*p0) else NA_character_
      pros_summary  <- if (!is.na(top1)) sprintf("%s: %.1f%%", top1, 100*p1) else NA_character_

      tibble(
        stratum       = stratum_name,
        variable      = v,
        variable_plot = pretty_var_label(v),
        type          = vtype,
        retro_summary = retro_summary,
        pros_summary  = pros_summary,
        retro_avail   = avail0,
        pros_avail    = avail1,
        avail_diff    = (avail1 - avail0),
        smd           = smd,
        ks_stat       = NA_real_
      )
    }
  })

  out %>%
    mutate(abs_smd = abs(smd)) %>%
    arrange(desc(abs_smd))
}

# ----------------------------
# Plot helpers
# ----------------------------
# helper: identify count / row-count style variables that should be excluded
is_excluded_density_var <- function(variable, variable_plot = NA_character_) {
  variable_clean <- variable %>%
    stringr::str_trim() %>%
    stringr::str_to_lower()

  variable_plot_clean <- dplyr::coalesce(variable_plot, "") %>%
    stringr::str_trim() %>%
    stringr::str_to_lower()

  # exclude raw variable names like n_fio2_rows, number_xxx, fio2_count
  # and display labels like "N Fio2 Rows" or "Fio2 count"
  stringr::str_detect(variable_clean, "^n_|^number_|_count$|count$") |
    stringr::str_detect(variable_plot_clean, "^n\\b|^n\\s|^number\\b|count$")
}

plot_density_overlays <- function(df_retro, df_pros, shift_tbl, stratum_name) {

  shift_tbl_density <- shift_tbl %>%
    dplyr::mutate(
      exclude_from_density = is_excluded_density_var(variable, variable_plot)
    ) %>%
    dplyr::filter(
      type == "continuous",
      is.finite(abs_smd),
      !exclude_from_density
    )

  top_vars <- shift_tbl_density %>%
    dplyr::slice_head(n = TOP_N_DENSITY_FEATURE) %>%
    dplyr::pull(variable)

  if (length(top_vars) == 0) return(NULL)

  label_map <- shift_tbl_density %>%
    dplyr::filter(variable %in% top_vars) %>%
    dplyr::distinct(variable, variable_plot)

  long <- dplyr::bind_rows(
    df_retro %>% dplyr::select(dplyr::all_of(top_vars)) %>% dplyr::mutate(epoch2 = "Retro"),
    df_pros  %>% dplyr::select(dplyr::all_of(top_vars)) %>% dplyr::mutate(epoch2 = "Prospective")
  ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(top_vars),
      names_to = "variable",
      values_to = "value"
    ) %>%
    dplyr::left_join(label_map, by = "variable") %>%
    dplyr::mutate(value_plot = value) %>%
    dplyr::filter(is.finite(value_plot))

  ggplot(long, aes(x = value_plot, color = epoch2, fill = epoch2)) +
    geom_density(alpha = 0.18, linewidth = 1) +
    facet_wrap(~ variable_plot, scales = "free", ncol = 2) +
    scale_color_manual(values = c("Retro" = "#D55E00", "Prospective" = "#0072B2")) +
    scale_fill_manual(values = c("Retro" = "#D55E00", "Prospective" = "#0072B2")) +
    labs(
      title = paste0("Top shifted continuous predictors — ", stratum_name),
      subtitle = "Density overlays comparing retrospective vs prospective cohorts",
      x = NULL,
      y = "Density",
      color = "Epoch",
      fill = "Epoch"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}
plot_availability_diff <- function(shift_tbl, stratum_name, min_abs_diff = 0.01) {

  top <- shift_tbl %>%
    mutate(abs_avail_diff = abs(avail_diff)) %>%
    filter(is.finite(abs_avail_diff), abs_avail_diff >= min_abs_diff) %>%
    arrange(desc(abs_avail_diff)) %>%
    slice_head(n = TOP_N_AVAIL_PLOT) %>%
    mutate(variable_plot = fct_reorder(variable_plot, abs_avail_diff))

  if (nrow(top) == 0) {
    message("No availability diffs >= ", min_abs_diff, " in ", stratum_name)
    return(NULL)
  }

  ggplot(top, aes(x = variable_plot, y = avail_diff, fill = avail_diff > 0)) +
    geom_col(width = 0.75, alpha = 0.9) +
    coord_flip() +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "gray35") +
    scale_fill_manual(values = c(`TRUE` = "#009E73", `FALSE` = "#CC79A7"), guide = "none") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = paste0("Availability shift — ", stratum_name),
      subtitle = paste0("Prospective minus retrospective; showing |difference| ≥ ", percent(min_abs_diff, accuracy = 1)),
      x = NULL,
      y = "Difference in availability"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold")
    )
}

plot_abs_smd <- function(shift_tbl, stratum_name) {
  top <- shift_tbl %>%
    filter(
      is.finite(abs_smd),
      !variable %in% c("blood_y_n", "urine_y_n", "csf_y_n")
    ) %>%
    slice_head(n = TOP_N_SMD_PLOT) %>%
    mutate(
      variable_plot = str_replace(variable_plot, "^N\\s+", "NUMBER OF "),
      variable_plot = str_to_upper(variable_plot),
      variable_plot = fct_reorder(variable_plot, abs_smd),
      type = factor(type, levels = c("continuous", "binary", "categorical"))
    )

  ggplot(top, aes(x = variable_plot, y = abs_smd, fill = type)) +
    geom_col(width = 0.72, alpha = 0.95) +
    coord_flip() +
    geom_vline(xintercept = NULL) +
    geom_hline(yintercept = 0.1, linetype = "dashed", linewidth = 0.8, color = "gray35") +
    geom_text(
      aes(label = sprintf("%.2f", abs_smd)),
      hjust = -0.1,
      size = 3.5,
      color = "black"
    ) +
    scale_fill_manual(
      values = c(
        "continuous" = "#56B4E9",
        "binary" = "#E69F00",
        "categorical" = "#CC79A7"
      ),
      name = "Variable type"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(
      title = paste0("Most shifted predictors by |SMD| — ", stratum_name),
      subtitle = "Dashed line marks |SMD| = 0.10",
      x = NULL,
      y = "|Standardized Mean Difference|"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.major.y = element_blank()
    )
}

# ----------------------------
# Run for each abx stratum
# ----------------------------
# ----------------------------
# Create prospective suspected-infection subsets
# ----------------------------
pros_no_abx_1st_infxn_si <- pros_no_abx_1st_infxn %>%
  filter(suspected_infection == 1)

pros_yes_abx_1st_infxn_si <- pros_yes_abx_1st_infxn %>%
  filter(suspected_infection == 1)

# Optional: quick check of sample sizes
cat("No-abx prospective suspected infection subset n =", nrow(pros_no_abx_1st_infxn_si), "\n")
cat("Yes-abx prospective suspected infection subset n =", nrow(pros_yes_abx_1st_infxn_si), "\n")


run_stratum <- function(df_retro, df_pros, stratum_name) {

  shift_tbl <- build_shift_table(df_retro, df_pros, stratum_name)

  p_density <- plot_density_overlays(df_retro, df_pros, shift_tbl, stratum_name)
  p_avail   <- plot_availability_diff(shift_tbl, stratum_name)
  p_smd     <- plot_abs_smd(shift_tbl, stratum_name)

  cat("\n============================================================\n")
  cat("SHIFT TABLE — ", stratum_name, "\n", sep = "")
  cat("============================================================\n")
  print(
    shift_tbl %>%
      select(variable, variable_plot, type, retro_summary, pros_summary,
             retro_avail, pros_avail, avail_diff, smd, abs_smd, ks_stat) %>%
      slice_head(n = 50)
  )
  cat("\n(Note: only top 50 rows printed; full table is in the returned object.)\n")

  cat("\n============================================================\n")
  cat("PLOTS — ", stratum_name, "\n", sep = "")
  cat("============================================================\n")

  if (!is.null(p_density)) print(p_density)
  if (!is.null(p_avail)) print(p_avail)
  print(p_smd)

  list(
    shift_table = shift_tbl,
    p_density   = p_density,
    p_avail     = p_avail,
    p_smd       = p_smd
  )
}



# Antibiotic-unexposed stratum
no_abx_results <- run_stratum(
  df_retro = retro_no_abx_1st_infxn,
  df_pros  = pros_no_abx_1st_infxn,
  stratum_name = "No antibiotics prior to prediction"
)

# Antibiotic-exposed stratum
yes_abx_results <- run_stratum(
  df_retro = retro_yes_abx_1st_infxn,
  df_pros  = pros_yes_abx_1st_infxn,
  stratum_name = "Antibiotics prior to prediction"
)

# Convenience: expose key objects at top-level
shift_no_abx <- no_abx_results$shift_table
shift_yes_abx <- yes_abx_results$shift_table

p_no_abx_density <- no_abx_results$p_density
p_no_abx_avail   <- no_abx_results$p_avail
p_no_abx_smd     <- no_abx_results$p_smd

p_yes_abx_density <- yes_abx_results$p_density
p_yes_abx_avail   <- yes_abx_results$p_avail
p_yes_abx_smd     <- yes_abx_results$p_smd


# Now another version of shifed SMD predictor plot without abs value
plot_signed_smd <- function(shift_tbl, stratum_name) {
  top <- shift_tbl %>%
    filter(
      is.finite(smd),
      is.finite(abs_smd),
      !variable %in% c("blood_y_n", "urine_y_n", "csf_y_n")
    ) %>%
    slice_head(n = TOP_N_SMD_PLOT) %>%
    mutate(
      variable_plot = str_replace(variable_plot, "^N\\s+", "NUMBER OF "),
      variable_plot = str_to_upper(variable_plot),
      variable_plot = fct_reorder(variable_plot, smd),
      type = factor(type, levels = c("continuous", "binary", "categorical"))
    )

  ggplot(top, aes(x = variable_plot, y = smd, fill = type)) +
    geom_col(width = 0.72, alpha = 0.95) +
    coord_flip() +
    geom_hline(yintercept = 0, linewidth = 0.7, color = "black") +
    geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", linewidth = 0.8, color = "gray35") +
    geom_text(
      aes(
        label = sprintf("%.2f", smd),
        hjust = ifelse(smd >= 0, -0.1, 1.1)
      ),
      size = 3.5,
      color = "black"
    ) +
    scale_fill_manual(
      values = c(
        "continuous" = "#56B4E9",
        "binary" = "#E69F00",
        "categorical" = "#CC79A7"
      ),
      name = "Variable type"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.10, 0.10))
    ) +
    labs(
      title = paste0("Most shifted predictors by signed SMD — ", stratum_name),
      subtitle = "Bars show direction of shift; dashed lines mark SMD = ±0.10",
      x = NULL,
      y = "STANDARDIZED MEAN DIFFERENCE"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.major.y = element_blank()
    )
}

p_no_abx_smd_signed <- plot_signed_smd(
  shift_no_abx,
  "No antibiotics prior to prediction"
)

p_yes_abx_smd_signed <- plot_signed_smd(
  shift_yes_abx,
  "Antibiotics prior to prediction"
)

p_no_abx_smd_signed
p_yes_abx_smd_signed


##### Repeat all above plots with only the suspicion of infection plots #####
# ----------------------------
# Run retrospective vs prospective suspected-infection subset
# ----------------------------

# Antibiotic-unexposed stratum
no_abx_results_si <- run_stratum(
  df_retro = retro_no_abx_1st_infxn,
  df_pros  = pros_no_abx_1st_infxn_si,
  stratum_name = "No antibiotics prior to prediction: retrospective vs prospective SI+ subset"
)

# Antibiotic-exposed stratum
yes_abx_results_si <- run_stratum(
  df_retro = retro_yes_abx_1st_infxn,
  df_pros  = pros_yes_abx_1st_infxn_si,
  stratum_name = "Antibiotics prior to prediction: retrospective vs prospective SI+ subset"
)

# Convenience: expose key objects at top-level
shift_no_abx_si <- no_abx_results_si$shift_table
shift_yes_abx_si <- yes_abx_results_si$shift_table

p_no_abx_density_si <- no_abx_results_si$p_density
p_no_abx_avail_si   <- no_abx_results_si$p_avail
p_no_abx_smd_si     <- no_abx_results_si$p_smd

p_yes_abx_density_si <- yes_abx_results_si$p_density
p_yes_abx_avail_si   <- yes_abx_results_si$p_avail
p_yes_abx_smd_si     <- yes_abx_results_si$p_smd

## Signed version of SMD plots ##
p_no_abx_smd_signed_si <- plot_signed_smd(
  shift_no_abx_si,
  "No antibiotics prior to prediction: retrospective vs prospective SI+ subset"
)

p_yes_abx_smd_signed_si <- plot_signed_smd(
  shift_yes_abx_si,
  "Antibiotics prior to prediction: retrospective vs prospective SI+ subset"
)

p_no_abx_smd_signed_si
p_yes_abx_smd_signed_si


#### Now will plot components of suspicion of infection over time: CRP, PCT, CXR, micro
# Retro data exists in the below dataframes

r_plot_susp_yes_abx <- retro_yes_abx_final %>%
  dplyr::select(study_id, adm_year, adm_quarter, crp_pres, pct_pres, cxr_pres, any_micro_pres) %>%
  distinct()

r_plot_susp_no_abx <- retro_no_abx_final %>%
  dplyr::select(study_id, adm_year, adm_quarter, crp_pres, pct_pres, cxr_pres, any_micro_pres) %>%
  distinct()

joint_retro <- bind_rows(r_plot_susp_no_abx, r_plot_susp_yes_abx)

# Make the plots
make_testing_plot <- function(df, plot_title) {

  df2 <- df %>%
    distinct(
      study_id, adm_year, adm_quarter,
      crp_pres, pct_pres, cxr_pres, any_micro_pres
    ) %>%
    mutate(
      adm_quarter = as.integer(adm_quarter),
      adm_year = as.integer(adm_year),
      year_quarter = paste0(adm_year, " Q", adm_quarter)
    )

  # preserve true chronological order
  quarter_levels <- df2 %>%
    distinct(adm_year, adm_quarter, year_quarter) %>%
    arrange(adm_year, adm_quarter) %>%
    pull(year_quarter)

  long_df <- df2 %>%
    pivot_longer(
      cols = c(crp_pres, pct_pres, cxr_pres, any_micro_pres),
      names_to = "test",
      values_to = "present"
    ) %>%
    mutate(
      test = recode(
        test,
        crp_pres = "CRP",
        pct_pres = "Procalcitonin",
        cxr_pres = "CXR",
        any_micro_pres = "Micro testing"
      ),
      test = factor(test, levels = c("CRP", "Procalcitonin", "CXR", "Micro testing")),
      year_quarter = factor(year_quarter, levels = quarter_levels)
    )

  plot_df <- long_df %>%
    group_by(adm_year, adm_quarter, year_quarter, test) %>%
    summarise(
      n = n(),
      n_present = sum(present == 1, na.rm = TRUE),
      prop = n_present / n,
      .groups = "drop"
    ) %>%
    bind_cols(
      binom::binom.confint(
        x = .$n_present,
        n = .$n,
        methods = "wilson"
      ) %>%
        dplyr::select(lower, upper)
    ) %>%
    rename(
      ci_low = lower,
      ci_high = upper
    )

  test_colors <- c(
    "CRP" = "#D55E00",
    "Procalcitonin" = "#7B3294",
    "CXR" = "#0072B2",
    "Micro testing" = "#009E73"
  )

  ggplot(plot_df, aes(x = year_quarter, y = prop, color = test, fill = test, group = test)) +
    geom_ribbon(
      aes(ymin = ci_low, ymax = ci_high),
      alpha = 0.18,
      color = NA
    ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.8) +
    scale_color_manual(values = test_colors) +
    scale_fill_manual(values = test_colors) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = plot_title,
      x = "PICU admission year and quarter",
      y = "Percent of patients with test performed",
      color = "Test",
      fill = "Test"
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      axis.title.x = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold", size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 14),
      legend.title = element_text(face = "bold", size = 15),
      legend.text = element_text(size = 14),
      strip.text = element_text(face = "bold", size = 15),
      legend.position = "bottom"
    )
}

p_yes_abx <- make_testing_plot(
  r_plot_susp_yes_abx,
  "Testing over time: suspected infection, antibiotic exposed before PICU"
)

p_no_abx <- make_testing_plot(
  r_plot_susp_no_abx,
  "Testing over time: suspected infection, no antibiotic exposure before PICU"
)

p_joint <- make_testing_plot(
  joint_retro,
  "Suspicion of Infection Testing Over Time: Retrospective Epoch"
)

p_yes_abx
p_no_abx
p_joint

### Now repeat plot for the prospective group ###
# Add year of admit and fix micro column name to work with function
pros_susp_elements_to_plot_no_abx <- pros_susp_elements_to_plot_no_abx %>%
  mutate(adm_year = year(picu_adm_date_time)) %>%
  rename(any_micro_pres = micro_pres)

pros_susp_elements_to_plot_no_abx <- pros_susp_elements_to_plot_no_abx %>%
  mutate(adm_quarter = quarter(picu_adm_date_time))

p_no_abx_pros <- make_testing_plot(
  pros_susp_elements_to_plot_no_abx %>% filter(suspected_infection == 1),
  "Testing over time: —Antibiotic exposure before PICU, prospective \n Only those with suspected infection"
)

p_no_abx_pros

### Repeat for abx exposed group and then full cohort ###
# Add year of admit and fix micro column name to work with function
pros_susp_elements_to_plot_yes_abx <- pros_susp_elements_to_plot_yes_abx %>%
  mutate(adm_year = year(picu_adm_date_time)) %>%
  rename(any_micro_pres = micro_pres)

pros_susp_elements_to_plot_yes_abx <- pros_susp_elements_to_plot_yes_abx %>%
  mutate(adm_quarter = quarter(picu_adm_date_time))

p_yes_abx_pros <- make_testing_plot(
  pros_susp_elements_to_plot_yes_abx %>% filter(suspected_infection == 1),
  "Testing over time: +Antibiotic exposure before PICU, prospective \n Patients with Suspected Infection"
)

p_yes_abx_pros

### Make joint df and plot
joint_pros <- bind_rows(pros_susp_elements_to_plot_no_abx, pros_susp_elements_to_plot_yes_abx)

p_joint_pros <- make_testing_plot(
  joint_pros,
  "Suspicion of Infection Testing Over Time: Prospective Epoch"
)

p_joint_pros


##Finally will eval number of SBI neg patient who got abx, what proportion predicted to
## be SBI-negative by the model at the 2 hour mark (from all comers not just suspicion of
# infection), and the median duration of antibiotics (up to 7 days) given after the 2 hour
# mark

#-----------------------------
# Load abx dataset
#-----------------------------
abx_raw <- read_csv(file = "/phi/sbi/prospective_data/Prospective/antinfective_export_pros_091225.csv")

# Exact MEDICATION_NAME values from the new source file
# Classify each as antibiotic, antifungal, antiviral, or other

antibiotic_meds <- c(
  "AMIKACIN (5 MG/ML) IN NS (BLADDER IRRIGATION)",
  "AMOXICILLIN 125 MG PO CHEW TAB",
  "AMOXICILLIN 250 MG PO CAP",
  "AMOXICILLIN 250 MG PO CHEW TAB",
  "AMOXICILLIN 400 MG/5ML PO RECON SUSP",
  "AMOXICILLIN 500 MG PO CAP",
  "AMOXICILLIN 875 MG PO TAB",
  "AMOXICILLIN-POT CLAVULANATE 200-28.5 MG PO CHEW TAB",
  "AMOXICILLIN-POT CLAVULANATE 500-125 MG PO TAB",
  "AMOXICILLIN-POT CLAVULANATE 600-42.9 MG/5ML PO RECON SUSP",
  "AMOXICILLIN-POT CLAVULANATE 875-125 MG PO TAB",
  "AMPICILLIN IN NS INJ",
  "AMPICILLIN-SULBACTAM IN NS INJ",
  "ANTIBIOTIC LOCK THERAPY WITH HEPARIN 10 UNITS/ML",
  "AZITHROMYCIN 125 MG PO TABS - IP ONLY",
  "AZITHROMYCIN 200 MG/5ML PO RECON SUSP",
  "AZITHROMYCIN 250 MG PO TAB",
  "AZITHROMYCIN IN NS INJ",
  "CEFAZOLIN (CONC: 20MG/ML) IN NS",
  "CEFAZOLIN PF (CONC: 200 MG/ML) SUBCONJUNCTIVAL INJ (OR USE ONLY)",
  "CEFAZOLIN SODIUM 1 G INJ RECON SOLN",
  "CEFAZOLIN SODIUM-DEXTROSE 1-4 GM/50ML-% IV SOLUTION",
  "CEFDINIR 125 MG/5ML PO RECON SUSP",
  "CEFDINIR 300 MG PO CAP",
  "CEFEPIME IN NS",
  "CEFIXIME 100 MG/5ML PO RECON SUSP",
  "CEFOXITIN IN NS IV INJ",
  "CEFPODOXIME PROXETIL 100 MG PO TAB",
  "CEFPODOXIME PROXETIL 100 MG/5ML PO RECON SUSP",
  "CEFPODOXIME PROXETIL 200 MG PO TAB",
  "CEFPROZIL 250 MG/5ML PO RECON SUSP",
  "CEFTAROLINE (CONC: 12 MG/ML) IN NS",
  "CEFTAZIDIME (CONC: 2.25MG/0.1ML) INTRAVITREAL INJ",
  "CEFTAZIDIME IN NS IV INJ",
  "CEFTRIAXONE IN DEXTROSE (40MG/ML) IV SOLUTION",
  "CEFTRIAXONE IN NS (40 MG/ML) IV SOLUTION",
  "CEFTRIAXONE SODIUM 1 G VIAL",
  "CEFTRIAXONE SODIUM 500 MG INJ RECON SOLN",
  "CEFTRIAXONE SODIUM IN DEXTROSE 40 MG/ML IV SOLUTION",
  "CEFTRIAXONE WITH LIDOCAINE 1% PF (350 MG/ML) IM SOLUTION (1 G VIAL)",
  "CEPHALEXIN 250 MG PO CAP",
  "CEPHALEXIN 250 MG/5ML PO RECON SUSP",
  "CEPHALEXIN 500 MG PO CAP",
  "CIPROFLOXACIN 500 MG/5ML (10%) PO RECON SUSP",
  "CIPROFLOXACIN HCL 125 MG PO TABS - IP ONLY",
  "CIPROFLOXACIN HCL 250 MG PO TAB",
  "CIPROFLOXACIN HCL 500 MG PO TAB",
  "CIPROFLOXACIN IN D5W 200 MG/100ML IV SOLUTION",
  "CLINDAMYCIN (CONC: 18MG/ML) IN NS",
  "CLINDAMYCIN HCL 300 MG PO CAP",
  "CLINDAMYCIN HCL 75 MG PO CAP",
  "CLINDAMYCIN PALMITATE HCL 75 MG/5ML PO RECON SOLN",
  "CLINDAMYCIN PHOSPHATE IN D5W 900 MG/50ML IV SOLUTION",
  "DAPTOMYCIN IN LR IV",
  "DOXYCYCLINE (CONC: 8MG/ML) IN NS - PLEURODESIS",
  "DOXYCYCLINE HYCLATE 100 MG PO CAP",
  "DOXYCYCLINE IN NS FOR SCLEROTHERAPY",
  "DOXYCYCLINE IN NS IV INJ",
  "DOXYCYCLINE MONOHYDRATE 25 MG/5ML PO RECON SUSP",
  "ERTAPENEM INJECTABLE SOLN",
  "ERYTHROMYCIN BASE 125 MG PO TABS - IP ONLY",
  "ERYTHROMYCIN BASE 250 MG PO TAB",
  "ERYTHROMYCIN ETHYLSUCCINATE 200 MG/5ML PO RECON SUSP",
  "ERYTHROMYCIN LACTOBIONATE IN NS IV INJ",
  "FIDAXOMICIN 40 MG/ML PO RECON SUSP",
  "FIRST-METRONIDAZOLE 50 MG/ML PO RECON SUSP",
  "FIRST-VANCOMYCIN HCL 50 MG/ML PO RECON SOLN",
  "GENTAMICIN IN NS BLADDER IRRIGATION",
  "GENTAMICIN IN SALINE 2-0.9 MG/ML-% IV SOLUTION",
  "ISONIAZID 100 MG PO TAB",
  "LEVOFLOXACIN 125 MG PO TABS - IP ONLY",
  "LEVOFLOXACIN 25 MG/ML PO SOLUTION",
  "LEVOFLOXACIN 250 MG PO TAB",
  "LEVOFLOXACIN 500 MG PO TAB",
  "LEVOFLOXACIN IN D5W 500 MG/100ML IV SOLUTION",
  "LINEZOLID 100 MG/5ML PO RECON SUSP",
  "LINEZOLID 2 MG/ML IN DEXTROSE",
  "LINEZOLID 300 MG PO TABS - IP ONLY",
  "LINEZOLID 600 MG PO TAB",
  "LINEZOLID 600 MG/300ML IV SOLUTION",
  "MEROPENEM IN NS (20 MG/ML) IV INJECTION",
  "MEROPENEM IN NS (20 MG/ML) IV SOLUTION EXTENDED INFUSION",
  "MEROPENEM IV INJ",
  "METRONIDAZOLE 250 MG PO TAB",
  "METRONIDAZOLE 500 MG PO TAB",
  "METRONIDAZOLE IN NACL 5 MG/ML IV SOLUTION WRAPPER",
  "MOXIFLOXACIN HCL IN NACL 400 MG/250ML IV SOLUTION",
  "NAFCILLIN CONTINUOUS INFUSION",
  "NAFCILLIN SODIUM IN DEXTROSE 2 GM/100ML IV SOLUTION",
  "NITROFURANTOIN 25 MG/5ML PO SUSP",
  "NITROFURANTOIN MACROCRYSTAL 100 MG PO CAP",
  "NITROFURANTOIN MACROCRYSTAL 50 MG PO CAP",
  "PENICILLIN G POTASSIUM INJ CONTINUOUS INFUSION",
  "PENICILLIN V POTASSIUM 250 MG PO TAB",
  "PENICILLIN V POTASSIUM 250 MG/5ML PO RECON SOLN",
  "PIPERACILLIN - TAZOBACTAM (ZOSYN)(CONC:60 MG/ML) IN NS IV INJ",
  "PYRAZINAMIDE 100 MG/ML ORAL SUSP",
  "RIFAMPIN 25 MG/ML ORAL SUSP",
  "RIFAMPIN IN NS IV INJ",
  "RIFAXIMIN (CONC: 20MG/ML) ORAL SUSPENSION",
  "RIFAXIMIN 550 MG PO TAB",
  "SULFAMETHOXAZOLE-TRIMETHOPRIM 200-40 MG PO TABS - IP ONLY",
  "SULFAMETHOXAZOLE-TRIMETHOPRIM 200-40 MG/5ML PO SUSP",
  "SULFAMETHOXAZOLE-TRIMETHOPRIM 400-80 MG PO TAB",
  "SULFAMETHOXAZOLE-TRIMETHOPRIM 800-160 MG PO TAB",
  "SULFAMETHOXAZOLE-TRIMETHOPRIM IN D5W INFUSION (CPOE)",
  "TOBRAMYCIN 300 MG/5ML INH NEB SOLN",
  "VANCOMYCIN (CONC: 1MG/0.1ML) INTRAVITREAL INJ",
  "VANCOMYCIN 1 MG/ML IRRIGATION SOLUTION",
  "VANCOMYCIN CONTINUOUS INFUSION",
  "VANCOMYCIN HCL 1000 MG IV RECON SOLN FOR BEADS/TOPICAL USE",
  "VANCOMYCIN HCL 500 MG IV SOLR FOR BEADS/TOPICAL USE",
  "VANCOMYCIN HCL IN NACL 1-0.9 GM/200ML-% IV SOLUTION",
  "VANCOMYCIN HCL IN NACL 500-0.9 MG/100ML-% IV SOLUTION",
  "VANCOMYCIN IN NS IV INJ"
)

antifungal_meds <- c(
  "AMBISOME (AMPHOTERICIN B) IN D5W INJ (CPOE)",
  "AMPHOTERICIN B IN D5W INJ (CPOE)",
  "ANIDULAFUNGIN IN NS IV",
  "FLUCONAZOLE 100 MG PO TAB",
  "FLUCONAZOLE 150 MG PO TAB",
  "FLUCONAZOLE 40 MG/ML PO RECON SUSP",
  "FLUCONAZOLE 50 MG PO TAB",
  "FLUCONAZOLE IN SODIUM CHLORIDE 200-0.9 MG/100ML-% IV SOLUTION",
  "ISAVUCONAZONIUM SULFATE 186 MG PO CAP",
  "ISAVUCONAZONIUM SULFATE 74.5 MG PO CAP",
  "ISAVUCONAZONIUM SULFATE IN NS IV",
  "MICAFUNGIN IN NS 1.5 MG/ML (50 MG VIAL) INFUSION",
  "MICAFUNGIN SODIUM-NACL 100-0.9 MG/100ML-% IV SOLUTION",
  "MICAFUNGIN SODIUM-NACL 50-0.9 MG/50ML-% IV SOLUTION",
  "POSACONAZOLE 100 MG PO DR TAB",
  "POSACONAZOLE 40 MG/ML PO SUSP",
  "POSACONAZOLE IN NS IV",
  "TERBINAFINE HCL 125 MG PO TAB - IP ONLY",
  "TERBINAFINE HCL 250 MG PO TAB",
  "VORICONAZOLE (CONC: 100MCG/0.1ML) INTRAVITREAL INJ",
  "VORICONAZOLE 200 MG PO TAB",
  "VORICONAZOLE 40 MG/ML PO RECON SUSP",
  "VORICONAZOLE 50 MG PO TAB",
  "VORICONAZOLE/NS INFUSION"
)

antiviral_meds <- c(
  "ACYCLOVIR 200 MG PO CAP",
  "ACYCLOVIR 200 MG/5ML PO SUSP",
  "ACYCLOVIR 400 MG PO TAB",
  "ACYCLOVIR 800 MG PO TAB",
  "ACYCLOVIR IN NS INJ (CPOE)",
  "BRINCIDOFOVIR (INVESTIGATIONAL) (EIND #132354) IN D5W IV INFUSION",
  "CIDOFOVIR IN NS INJ. (CPOE)",
  "FOSCARNET SODIUM 6000 MG/250ML IV SOLUTION",
  "GANCICLOVIR IN NS INJ (CPOE)",
  "LETERMOVIR (1.2 MG/ML) IN NS IV",
  "LETERMOVIR 120 MG PO TABS - IP ONLY",
  "LETERMOVIR 240 MG IN NS 250ML IV INFUSION",
  "LETERMOVIR 240 MG PO TAB",
  "LETERMOVIR 480 MG IN NS 250ML IV INFUSION",
  "LETERMOVIR 480 MG PO TAB",
  "MARIBAVIR 200 MG PO TAB",
  "NIRMATRELVIR&RITONAVIR 300/100 20 X 150 MG & 10 X 100MG PO TAB TX PACK",
  "OSELTAMIVIR PHOSPHATE 30 MG PO CAP",
  "OSELTAMIVIR PHOSPHATE 6 MG/ML PO RECON SUSP",
  "OSELTAMIVIR PHOSPHATE 75 MG PO CAP",
  "PERAMIVIR 200 MG/20ML IV SOLUTION",
  "REMDESIVIR IN NS IV GREATER THAN 40 KG",
  "REMDESIVIR IN NS IV LESS THAN 40 KG",
  "VALACYCLOVIR (CONC: 50 MG/ML) ORAL SUSPENSION",
  "VALACYCLOVIR HCL 500 MG PO TAB",
  "VALGANCICLOVIR HCL 450 MG PO TAB",
  "VALGANCICLOVIR HCL 50 MG/ML PO RECON SOLN"
)

other_meds <- c(
  "ALBENDAZOLE 200 MG PO TAB",
  "BASE SOLN/ADDITIVES (CPOE)",
  "EXTEMPORANEOUS MIXTURE TEMPLATE",
  "HYDROXYCHLOROQUINE SULFATE 200 MG PO TAB",
  "HYDROXYCHLOROQUINE SULFATE 25 MG/ML ORAL SUSP",
  "IV ROOM MIXTURE - BLANK",
  "IVERMECTIN 3 MG PO TAB",
  "LIDOCAINE IN NS 4 ML INHALED MIXTURE",
  "PENTAMIDINE IN D5W IV INJ",
  "PENTAMIDINE ISETHIONATE 300 MG INH RECON SOLN",
  "PYRANTEL PAMOATE 144 (50 BASE) MG/ML PO SUSP",
  "TINIDAZOLE (CONC: 67 MG/ML) ORAL SUSPENSION"
)

# Build the classification dataframe
antimicrobial_class_df <- bind_rows(
  tibble(medication = antibiotic_meds, medication_class = "antibiotic"),
  tibble(medication = antifungal_meds, medication_class = "antifungal"),
  tibble(medication = antiviral_meds, medication_class = "antiviral"),
  tibble(medication = other_meds, medication_class = "other")
)

# QC checks
all_meds_in_file <- abx_raw %>%
  distinct(MEDICATION_NAME) %>%
  arrange(MEDICATION_NAME)

# Should be 178 rows in each of these
nrow(all_meds_in_file)
nrow(antimicrobial_class_df)

# These should both return 0 rows if every medication was classified exactly once
anti_join(all_meds_in_file, antimicrobial_class_df, by = c("MEDICATION_NAME" = "medication"))
anti_join(antimicrobial_class_df, all_meds_in_file, by = c("medication" = "MEDICATION_NAME"))

# Optional: make sure there are no duplicates in your classification table
antimicrobial_class_df %>%
  count(medication) %>%
  filter(n > 1)

# Filter to antibiotics only
list_abx <- antimicrobial_class_df %>%
  filter(medication_class == "antibiotic")

just_abx <- abx_raw %>%
  filter(MEDICATION_NAME %in% list_abx$medication)

#-----------------------------
# 1) Combine prospective data and derive CSN from study_id
#    Keep CSN as CHARACTER to avoid any coercion issues
#-----------------------------
pros_all_1st_infxn <- bind_rows(
  pros_no_abx_1st_infxn %>% mutate(source_df = "pros_no_abx"),
  pros_yes_abx_1st_infxn %>% mutate(source_df = "pros_yes_abx")
) %>%
  mutate(
    csn = str_remove(as.character(study_id), "_1$"),
    picu_adm_date_time = with_tz(as.POSIXct(picu_adm_date_time), tzone = "America/Denver"),
    window_start = picu_adm_date_time + lubridate::dhours(2),
    window_end   = picu_adm_date_time + lubridate::dhours(26)
  ) %>%
  distinct(study_id, .keep_all = TRUE)

#-----------------------------
# 2) Prep antibiotic administrations
#    Keep CSN as CHARACTER
#-----------------------------
abx_admin_full <- just_abx %>%
  transmute(
    csn             = as.character(PAT_ENC_CSN_ID),
    taken_time      = force_tz(as.POSIXct(taken_time), tzone = "America/Denver"),
    GENERIC_NAME    = GENERIC_NAME,
    MEDICATION_NAME = MEDICATION_NAME,
    HSP_ACCOUNT_ID  = HSP_ACCOUNT_ID,
    route_name      = route_name
  ) %>%
  filter(!is.na(csn), !is.na(taken_time)) %>%
  arrange(csn, taken_time)

abx_admin <- abx_admin_full %>% left_join(pros_all_1st_infxn %>% dplyr::select(csn, picu_adm_date_time), by = "csn")
abx_admin <- abx_admin %>% filter(taken_time > picu_adm_date_time)

# Lookup table to fill in missing/NULL generic names
generic_name_fill <- tibble::tribble(
  ~MEDICATION_NAME,                                               ~generic_fill,
  "AMPICILLIN IN NS INJ",                                         "Ampicillin",
  "AMPICILLIN-SULBACTAM IN NS INJ",                               "Ampicillin-Sulbactam",
  "ANTIBIOTIC LOCK THERAPY WITH HEPARIN 10 UNITS/ML",             "Antibiotic Lock Therapy",
  "AZITHROMYCIN IN NS INJ",                                       "Azithromycin",
  "CEFAZOLIN (CONC: 20MG/ML) IN NS",                              "Cefazolin",
  "CEFEPIME IN NS",                                               "Cefepime",
  "CEFOXITIN IN NS IV INJ",                                       "Cefoxitin",
  "CEFTAROLINE (CONC: 12 MG/ML) IN NS",                           "Ceftaroline",
  "CEFTAZIDIME IN NS IV INJ",                                     "Ceftazidime",
  "CEFTRIAXONE IN DEXTROSE (40MG/ML) IV SOLUTION",                "Ceftriaxone",
  "CEFTRIAXONE IN NS (40 MG/ML) IV SOLUTION",                     "Ceftriaxone",
  "DAPTOMYCIN IN LR IV",                                          "Daptomycin",
  "DOXYCYCLINE (CONC: 8MG/ML) IN NS - PLEURODESIS",               "Doxycycline",
  "DOXYCYCLINE IN NS FOR SCLEROTHERAPY",                          "Doxycycline",
  "DOXYCYCLINE IN NS IV INJ",                                     "Doxycycline",
  "ERTAPENEM INJECTABLE SOLN",                                    "Ertapenem",
  "GENTAMICIN IN NS BLADDER IRRIGATION",                          "Gentamicin",
  "MEROPENEM IN NS (20 MG/ML) IV INJECTION",                      "Meropenem",
  "MEROPENEM IN NS (20 MG/ML) IV SOLUTION EXTENDED INFUSION",     "Meropenem",
  "PIPERACILLIN - TAZOBACTAM (ZOSYN)(CONC:60 MG/ML) IN NS IV INJ", "Piperacillin-Tazobactam",
  "SULFAMETHOXAZOLE-TRIMETHOPRIM IN D5W INFUSION (CPOE)",         "Sulfamethoxazole-Trimethoprim",
  "VANCOMYCIN IN NS IV INJ",                                      "Vancomycin"
)

abx_admin <- abx_admin %>%
  dplyr::left_join(generic_name_fill, by = "MEDICATION_NAME") %>%
  dplyr::mutate(
    GENERIC_NAME = dplyr::if_else(
      is.na(GENERIC_NAME) | GENERIC_NAME == "NULL" | GENERIC_NAME == "",
      generic_fill,
      GENERIC_NAME
    )
  ) %>%
  dplyr::select(-generic_fill)

#-----------------------------
# 3) Join by CSN and create antibiotic episodes
#-----------------------------
all_abx_by_encounter <- pros_all_1st_infxn %>%
  dplyr::select(source_df, study_id, csn, picu_adm_date_time, window_start, window_end) %>%
  inner_join(abx_admin %>% dplyr::select(-picu_adm_date_time), by = "csn") %>%
  arrange(study_id, taken_time) %>%
  group_by(study_id) %>%
  mutate(
    gap_from_prev_hours = as.numeric(difftime(taken_time, lag(taken_time), units = "hours")),
    new_episode = case_when(
      row_number() == 1 ~ 1L,
      gap_from_prev_hours > 26 ~ 1L,
      TRUE ~ 0L
    ),
    abx_episode_id = cumsum(new_episode),
    in_target_window = taken_time >= window_start & taken_time < window_end
  ) %>%
  ungroup()

#-----------------------------
# 4) Direct QC: how many actual doses fall in the target window?
#-----------------------------
window_hits <- all_abx_by_encounter %>%
  filter(in_target_window)

window_hits %>%
  summarise(
    n_doses_in_window = n(),
    n_encounters_in_window = n_distinct(study_id)
  )

#-----------------------------
# 5) Identify the first in-window dose for each encounter
#-----------------------------
first_in_window <- window_hits %>%
  arrange(study_id, taken_time) %>%
  group_by(study_id) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    study_id,
    first_abx_time_in_window = taken_time,
    first_episode_id = abx_episode_id
  )

#-----------------------------
# 6) Build the encounter-level FLAG directly from first_in_window
#    Do NOT depend on the course summary for this
#-----------------------------
encounter_abx_flag <- first_in_window %>%
  distinct(study_id) %>%
  mutate(abx_in_24h_after_2h_flag = 1L)

#-----------------------------
# 7) Keep all rows from the antibiotic episode containing
#    the first qualifying in-window dose
#-----------------------------
course_rows <- all_abx_by_encounter %>%
  inner_join(first_in_window, by = "study_id") %>%
  filter(abx_episode_id == first_episode_id)

#-----------------------------
# 8) Summarize the first qualifying antibiotic course
#-----------------------------
abx_course_summary <- course_rows %>%
  group_by(study_id) %>%
  summarise(
    first_abx_time_in_window = first(first_abx_time_in_window),
    first_course_start_time  = min(taken_time, na.rm = TRUE),
    first_course_end_time    = max(taken_time, na.rm = TRUE),
    n_doses_first_course     = n(),
    meds_in_first_course     = paste(sort(unique(GENERIC_NAME)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    first_course_duration_hours =
      as.numeric(difftime(first_course_end_time, first_course_start_time, units = "hours")),
    first_course_duration_days =
      as.numeric(difftime(first_course_end_time, first_course_start_time, units = "days"))
  )

#-----------------------------
# 9) Join back to encounter-level data
#-----------------------------
pros_all_1st_infxn_abx <- pros_all_1st_infxn %>%
  left_join(encounter_abx_flag, by = "study_id") %>%
  left_join(abx_course_summary, by = "study_id") %>%
  mutate(
    abx_in_24h_after_2h_flag = if_else(is.na(abx_in_24h_after_2h_flag), 0L, abx_in_24h_after_2h_flag)
  )

#-----------------------------
# 10) Split back out if desired
#-----------------------------
pros_no_abx_1st_infxn_abx <- pros_all_1st_infxn_abx %>%
  filter(source_df == "pros_no_abx") %>%
  select(-source_df)

pros_yes_abx_1st_infxn_abx <- pros_all_1st_infxn_abx %>%
  filter(source_df == "pros_yes_abx") %>%
  select(-source_df)

#-----------------------------
# QC
#-----------------------------
pros_all_1st_infxn_abx %>%
  count(abx_in_24h_after_2h_flag)

all_abx_by_encounter %>%
  mutate(hours_from_admit = as.numeric(difftime(taken_time, picu_adm_date_time, units = "hours"))) %>%
  summarise(
    min_hours = min(hours_from_admit, na.rm = TRUE),
    median_hours = median(hours_from_admit, na.rm = TRUE),
    max_hours = max(hours_from_admit, na.rm = TRUE)
  )

window_hits %>%
  select(study_id, csn, picu_adm_date_time, taken_time, GENERIC_NAME, in_target_window) %>%
  arrange(study_id, taken_time) %>%
  slice_head(n = 20)

#### Ok now look at antibiotic duration for different subgroups #######


#-----------------------------
# Add cohort labels and cohort-specific prediction flags
#-----------------------------
pros_no_abx_for_summary <- pros_no_abx_1st_infxn_abx %>%
  dplyr::mutate(
    cohort = "Prospective abx-unexposed",
    pred_neg_by_model = model_score <= 0.05
  )

pros_yes_abx_for_summary <- pros_yes_abx_1st_infxn_abx %>%
  dplyr::mutate(
    cohort = "Prospective abx-exposed",
    pred_neg_by_model = model_score <= 0.074
  )

pros_combined_for_summary <- dplyr::bind_rows(
  pros_no_abx_for_summary,
  pros_yes_abx_for_summary
)

#-----------------------------
# Helper function for manuscript-style summary
#-----------------------------
make_manuscript_summary <- function(df, cohort_label) {

  sbi_neg <- df %>%
    dplyr::filter(sbi_present == 0)

  sbi_neg_abx <- sbi_neg %>%
    dplyr::filter(abx_in_24h_after_2h_flag == 1)

  sbi_neg_abx_pred_neg <- sbi_neg_abx %>%
    dplyr::filter(pred_neg_by_model)

  n_sbi_neg <- nrow(sbi_neg)
  n_sbi_neg_abx <- nrow(sbi_neg_abx)
  n_sbi_neg_abx_pred_neg <- nrow(sbi_neg_abx_pred_neg)

  tibble::tibble(
    cohort = cohort_label,

    `SBI-negative patients given antibiotics after PICU+2h` =
      paste0(
        n_sbi_neg_abx, "/", n_sbi_neg,
        " (", round(100 * n_sbi_neg_abx / n_sbi_neg, 1), "%)"
      ),

    `Antibiotic course duration, days` =
      paste0(
        round(stats::median(sbi_neg_abx$first_course_duration_days, na.rm = TRUE), 2),
        " [",
        round(stats::quantile(sbi_neg_abx$first_course_duration_days, probs = 0.25, na.rm = TRUE), 2),
        ", ",
        round(stats::quantile(sbi_neg_abx$first_course_duration_days, probs = 0.75, na.rm = TRUE), 2),
        "]"
      ),

    `Among SBI-negative patients given antibiotics, predicted SBI-negative by model` =
      paste0(
        n_sbi_neg_abx_pred_neg, "/", n_sbi_neg_abx,
        " (", round(100 * n_sbi_neg_abx_pred_neg / n_sbi_neg_abx, 1), "%)"
      ),

    `Antibiotic course duration in model-predicted SBI-negative group, days` =
      paste0(
        round(stats::median(sbi_neg_abx_pred_neg$first_course_duration_days, na.rm = TRUE), 2),
        " [",
        round(stats::quantile(sbi_neg_abx_pred_neg$first_course_duration_days, probs = 0.25, na.rm = TRUE), 2),
        ", ",
        round(stats::quantile(sbi_neg_abx_pred_neg$first_course_duration_days, probs = 0.75, na.rm = TRUE), 2),
        "]"
      )
  )
}

#-----------------------------
# Create manuscript-style summary table including overall row
#-----------------------------
pros_summary_manuscript <- dplyr::bind_rows(
  make_manuscript_summary(
    df = pros_no_abx_for_summary,
    cohort_label = "Prospective abx-unexposed"
  ),
  make_manuscript_summary(
    df = pros_yes_abx_for_summary,
    cohort_label = "Prospective abx-exposed"
  ),
  make_manuscript_summary(
    df = pros_combined_for_summary,
    cohort_label = "Overall"
  )
)

print(pros_summary_manuscript)


##### Now do abx analysis by subgroups #######
library(dplyr)
library(tibble)

# Assumes age is already in YEARS.
# If age is in months, replace age_years = age with age_years = age / 12

pros_no_abx_for_subgroups <- pros_no_abx_1st_infxn_abx %>%
  dplyr::mutate(
    cohort = "Prospective abx-unexposed",
    age_years = age
  )

pros_yes_abx_for_subgroups <- pros_yes_abx_1st_infxn_abx %>%
  dplyr::mutate(
    cohort = "Prospective abx-exposed",
    age_years = age
  )

pros_combined_for_subgroups <- dplyr::bind_rows(
  pros_no_abx_for_subgroups,
  pros_yes_abx_for_subgroups
)

#-----------------------------
# Helpers
#-----------------------------
fmt_n_pct <- function(num, den, digits = 1) {
  if (is.na(den) || den == 0) return(NA_character_)
  paste0(num, "/", den, " (", round(100 * num / den, digits), "%)")
}

fmt_median_iqr <- function(x, digits = 2) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)

  paste0(
    round(stats::median(x), digits),
    " [",
    round(as.numeric(stats::quantile(x, probs = 0.25, na.rm = TRUE)), digits),
    ", ",
    round(as.numeric(stats::quantile(x, probs = 0.75, na.rm = TRUE)), digits),
    "]"
  )
}

make_nonprediction_row <- function(df, subgroup_label) {
  sbi_neg <- df %>%
    dplyr::filter(sbi_present == 0)

  sbi_neg_abx <- sbi_neg %>%
    dplyr::filter(abx_in_24h_after_2h_flag == 1)

  tibble::tibble(
    subgroup = subgroup_label,
    n_sbi_negative = nrow(sbi_neg),
    `SBI-negative patients given antibiotics after PICU+2h` =
      fmt_n_pct(nrow(sbi_neg_abx), nrow(sbi_neg)),
    `Antibiotic course duration, days` =
      fmt_median_iqr(sbi_neg_abx$first_course_duration_days)
  )
}

#-----------------------------
# Create subgroup table
#-----------------------------
pros_subset_summary_manuscript <- dplyr::bind_rows(
  make_nonprediction_row(
    pros_combined_for_subgroups,
    "Overall"
  ),
  make_nonprediction_row(
    pros_combined_for_subgroups %>% dplyr::filter(pccc == 0),
    "PCCC = 0"
  ),
  make_nonprediction_row(
    pros_combined_for_subgroups %>% dplyr::filter(pccc == 1),
    "PCCC = 1"
  ),
  make_nonprediction_row(
    pros_combined_for_subgroups %>% dplyr::filter(malignancy_pccc == 0),
    "Malignancy PCCC = 0"
  ),
  make_nonprediction_row(
    pros_combined_for_subgroups %>% dplyr::filter(malignancy_pccc == 1),
    "Malignancy PCCC = 1"
  ),
  make_nonprediction_row(
    pros_combined_for_subgroups %>%
      dplyr::filter(age_years >= 3/12, age_years < 6/12),
    "Age >=3 months to <6 months"
  ),
  make_nonprediction_row(
    pros_combined_for_subgroups %>%
      dplyr::filter(age_years >= 6/12, age_years < 1),
    "Age >=6 months to <1 year"
  ),
  make_nonprediction_row(
    pros_combined_for_subgroups %>%
      dplyr::filter(age_years >= 1, age_years < 5),
    "Age >=1 year to <5 years"
  ),
  make_nonprediction_row(
    pros_combined_for_subgroups %>%
      dplyr::filter(age_years >= 5),
    "Age >=5 years"
  )
)

print(pros_subset_summary_manuscript)

### Now make stacked bar charts for SBI- patients and those given abx
library(dplyr)
library(ggplot2)
library(scales)
library(forcats)

#-----------------------------------------
# Create plotting dataframe for stacked bars
#-----------------------------------------
make_abx_prop_plot_df <- function(df) {

  sbi_neg_df <- df %>%
    dplyr::filter(sbi_present == 0)

  subgroup_definitions <- list(
    "Overall" = sbi_neg_df,
    "PCCC = 0" = sbi_neg_df %>% dplyr::filter(pccc == 0),
    "PCCC = 1" = sbi_neg_df %>% dplyr::filter(pccc == 1),
    "Malignancy PCCC = 0" = sbi_neg_df %>% dplyr::filter(malignancy_pccc == 0),
    "Malignancy PCCC = 1" = sbi_neg_df %>% dplyr::filter(malignancy_pccc == 1),
    "Age >=3 months to <6 months" = sbi_neg_df %>% dplyr::filter(age_years >= 3/12, age_years < 6/12),
    "Age >=6 months to <1 year" = sbi_neg_df %>% dplyr::filter(age_years >= 6/12, age_years < 1),
    "Age >=1 year to <5 years" = sbi_neg_df %>% dplyr::filter(age_years >= 1, age_years < 5),
    "Age >=5 years" = sbi_neg_df %>% dplyr::filter(age_years >= 5)
  )

  plot_df <- purrr::imap_dfr(
    subgroup_definitions,
    function(sub_df, subgroup_name) {

      den <- nrow(sub_df)

      n_no_abx <- sum(sub_df$abx_in_24h_after_2h_flag == 0, na.rm = TRUE)
      n_yes_abx <- sum(sub_df$abx_in_24h_after_2h_flag == 1, na.rm = TRUE)

      tibble::tibble(
        subgroup = subgroup_name,
        abx_group = c("No antibiotics after PICU+2h", "Received antibiotics after PICU+2h"),
        n = c(n_no_abx, n_yes_abx),
        total = den
      )
    }
  ) %>%
    dplyr::mutate(
      prop = dplyr::if_else(total > 0, n / total, NA_real_),
      label_n = paste0(n, "/", total),
      subgroup = factor(
        subgroup,
        levels = c(
          "Overall",
          "PCCC = 0",
          "PCCC = 1",
          "Malignancy PCCC = 0",
          "Malignancy PCCC = 1",
          "Age >=3 months to <6 months",
          "Age >=6 months to <1 year",
          "Age >=1 year to <5 years",
          "Age >=5 years"
        )
      ),
      abx_group = factor(
        abx_group,
        levels = c(
          "No antibiotics after PICU+2h",
          "Received antibiotics after PICU+2h"
        )
      )
    )

  return(plot_df)
}

#-----------------------------------------
# Build plotting dataframe
#-----------------------------------------
abx_prop_plot_df <- make_abx_prop_plot_df(pros_combined_for_subgroups)

print(abx_prop_plot_df)

#-----------------------------------------
# Stacked percentage bar chart
#-----------------------------------------
p_abx_subgroup_stacked <- ggplot(
  abx_prop_plot_df,
  aes(x = subgroup, y = prop, fill = abx_group)
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = label_n),
    position = position_stack(vjust = 0.5),
    size = 3.8,
    fontface = "bold",
    color = "white"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  scale_fill_manual(
    values = c(
      "No antibiotics after PICU+2h" = "blue",
      "Received antibiotics after PICU+2h" = "red"
    ),
    name = NULL
  ) +
  labs(
    title = "Antibiotic Use Among SBI-Negative Encounters \nBetween 2-26 Hours after PICU Admission by Subgroup",
    x = NULL,
    y = "Percent of SBI-negative encounters"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 11, face = "bold", angle = 35, hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_blank()
  )

p_abx_subgroup_stacked


abx_prop_plot_df <- abx_prop_plot_df %>%
  dplyr::mutate(
    label_both = paste0(label_n, "\n", scales::percent(prop, accuracy = 0.1))
  )

p_abx_subgroup_stacked_2 <- ggplot(
  abx_prop_plot_df,
  aes(x = subgroup, y = prop, fill = abx_group)
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = label_both),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    fontface = "bold",
    color = "white",
    lineheight = 0.9
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  scale_fill_manual(
    values = c(
      "No antibiotics after PICU+2h" = "blue",
      "Received antibiotics after PICU+2h" = "red"
    ),
    name = NULL
  ) +
  labs(
    title = "SBI-Negative Encounters: Antibiotic Use After PICU+2h by Subgroup",
    x = NULL,
    y = "Percent of SBI-negative encounters"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 11, face = "bold", angle = 35, hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_blank()
  )

p_abx_subgroup_stacked_2

# Create plotting dataframe for duration boxplots
# Uses SBI-negative encounters who received antibiotics
# after the +2h PICU timepoint
#-------------------------------------------------
make_abx_duration_plot_df <- function(df) {

  sbi_neg_df <- df %>%
    dplyr::filter(
      sbi_present == 0,
      abx_in_24h_after_2h_flag == 1,
      !is.na(first_course_duration_days),
      first_course_duration_days <= 7
    )

  subgroup_definitions <- list(
    "Overall" = sbi_neg_df,
    "PCCC = 0" = sbi_neg_df %>% dplyr::filter(as.character(pccc) == "0"),
    "PCCC = 1" = sbi_neg_df %>% dplyr::filter(as.character(pccc) == "1"),
    "Malignancy = 0" = sbi_neg_df %>% dplyr::filter(as.character(malignancy_pccc) == "0"),
    "Malignancy = 1" = sbi_neg_df %>% dplyr::filter(as.character(malignancy_pccc) == "1"),
    "Age >=3 months to <6 months" = sbi_neg_df %>% dplyr::filter(age_years >= 3/12, age_years < 6/12),
    "Age >=6 months to <1 year" = sbi_neg_df %>% dplyr::filter(age_years >= 6/12, age_years < 1),
    "Age >=1 year to <5 years" = sbi_neg_df %>% dplyr::filter(age_years >= 1, age_years < 5),
    "Age >=5 years" = sbi_neg_df %>% dplyr::filter(age_years >= 5)
  )

  plot_df <- purrr::imap_dfr(
    subgroup_definitions,
    function(sub_df, subgroup_name) {
      sub_df %>%
        dplyr::transmute(
          subgroup = subgroup_name,
          first_course_duration_days = first_course_duration_days
        )
    }
  ) %>%
    dplyr::mutate(
      subgroup = factor(
        subgroup,
        levels = c(
          "Overall",
          "PCCC = 0",
          "PCCC = 1",
          "Malignancy = 0",
          "Malignancy = 1",
          "Age >=3 months to <6 months",
          "Age >=6 months to <1 year",
          "Age >=1 year to <5 years",
          "Age >=5 years"
        )
      )
    )

  return(plot_df)
}

abx_duration_plot_df <- make_abx_duration_plot_df(pros_combined_for_subgroups)

#-------------------------------------------------
# Build N labels for each subgroup
#-------------------------------------------------
n_labels_df <- abx_duration_plot_df %>%
  dplyr::group_by(subgroup) %>%
  dplyr::summarise(
    n = dplyr::n(),
    ymax = max(first_course_duration_days, na.rm = TRUE),
    .groups = "drop"
  )

overall_ymax <- max(abx_duration_plot_df$first_course_duration_days, na.rm = TRUE)

n_labels_df <- n_labels_df %>%
  dplyr::mutate(
    label = paste0("N = ", n),
    label_y = ymax + 0.06 * overall_ymax
  )

#-------------------------------------------------
# Box and whisker plot
#-------------------------------------------------
p_abx_duration_subgroups <- ggplot(
  abx_duration_plot_df,
  aes(x = subgroup, y = first_course_duration_days)
) +
  geom_boxplot(
    width = 0.7,
    linewidth = 0.8,
    outlier.size = 2
  ) +
  geom_text(
    data = n_labels_df,
    aes(x = subgroup, y = label_y, label = label),
    inherit.aes = FALSE,
    size = 3.8,
    fontface = "bold"
  ) +
  scale_x_discrete(labels = \(x) stringr::str_wrap(x, width = 18)) +
  labs(
    title = "Duration of First Antibiotic Course Among SBI-Negative Encounters",
    subtitle = "Includes only SBI-negative encounters receiving antibiotics after PICU +2h",
    x = NULL,
    y = "Antibiotic duration (days)"
  ) +
  coord_cartesian(
    ylim = c(
      min(abx_duration_plot_df$first_course_duration_days, na.rm = TRUE),
      max(n_labels_df$label_y, na.rm = TRUE) * 1.05
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 11, face = "bold", angle = 35, hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold")
  )

p_abx_duration_subgroups


#### Determine p values using Mood's Median test ####

#--------------------------------------------
# 1. Pairwise Mood's median test helper
#--------------------------------------------
run_moods_median_pair <- function(data, group_a, group_b) {

  dat <- data %>%
    dplyr::filter(subgroup %in% c(group_a, group_b)) %>%
    dplyr::filter(!is.na(first_course_duration_days)) %>%
    dplyr::mutate(subgroup = droplevels(subgroup))

  if (nrow(dat) == 0) {
    stop("No rows found for the requested subgroup comparison.")
  }

  # Pooled median for the two groups being compared
  pooled_median <- median(dat$first_course_duration_days, na.rm = TRUE)

  # Mood's median test typically uses above median vs at/below median
  dat <- dat %>%
    dplyr::mutate(above_median = ifelse(first_course_duration_days > pooled_median,
                                        "Above pooled median",
                                        "At or below pooled median"))

  tab <- table(dat$subgroup, dat$above_median)

  # Expected counts to decide chi-square vs Fisher
  chisq_obj <- suppressWarnings(chisq.test(tab, correct = FALSE))
  expected_min <- min(chisq_obj$expected)

  if (expected_min < 5) {
    test_obj <- fisher.test(tab)
    test_used <- "Fisher's exact test on Mood's median table"
    p_value <- test_obj$p.value
    statistic <- NA_real_
  } else {
    test_obj <- chisq.test(tab, correct = FALSE)
    test_used <- "Chi-square test on Mood's median table"
    p_value <- test_obj$p.value
    statistic <- unname(test_obj$statistic)
  }

  group_summary <- dat %>%
    dplyr::group_by(subgroup) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median_duration_days = median(first_course_duration_days, na.rm = TRUE),
      .groups = "drop"
    )

  tibble::tibble(
    comparison = paste0(group_a, " vs ", group_b),
    group_1 = group_a,
    group_2 = group_b,
    pooled_median_days = pooled_median,
    n_group_1 = group_summary$n[group_summary$subgroup == group_a],
    n_group_2 = group_summary$n[group_summary$subgroup == group_b],
    median_group_1_days = group_summary$median_duration_days[group_summary$subgroup == group_a],
    median_group_2_days = group_summary$median_duration_days[group_summary$subgroup == group_b],
    test_used = test_used,
    statistic = statistic,
    p_value = p_value
  )
}

#--------------------------------------------
# 2. Define the comparisons you want
#--------------------------------------------
comparisons_to_run <- tribble(
  ~group_a,                       ~group_b,
  "PCCC = 0",                     "PCCC = 1",
  "Malignancy = 0",               "Malignancy = 1",
  "Age >=3 months to <6 months",  "Age >=5 years",
  "Age >=6 months to <1 year",    "Age >=5 years",
  "Age >=1 year to <5 years",     "Age >=5 years"
)

#--------------------------------------------
# 3. Run the pairwise tests
#--------------------------------------------
moods_pairwise_results <- purrr::pmap_dfr(
  comparisons_to_run,
  ~ run_moods_median_pair(
    data = abx_duration_plot_df,
    group_a = ..1,
    group_b = ..2
  )
)

#--------------------------------------------
# 4. Add Holm-adjusted p values
#--------------------------------------------
moods_pairwise_results <- moods_pairwise_results %>%
  dplyr::mutate(
    p_value_holm = p.adjust(p_value, method = "holm")
  )

moods_pairwise_results

# One more to determine variability in age groups abx duration #
run_moods_median_omnibus <- function(data, groups) {

  dat <- data %>%
    dplyr::filter(subgroup %in% groups) %>%
    dplyr::filter(!is.na(first_course_duration_days)) %>%
    dplyr::mutate(subgroup = droplevels(subgroup))

  pooled_median <- median(dat$first_course_duration_days, na.rm = TRUE)

  dat <- dat %>%
    dplyr::mutate(above_median = ifelse(first_course_duration_days > pooled_median,
                                        "Above pooled median",
                                        "At or below pooled median"))

  tab <- table(dat$subgroup, dat$above_median)

  chisq_obj <- suppressWarnings(chisq.test(tab, correct = FALSE))
  expected_min <- min(chisq_obj$expected)

  if (expected_min < 5) {
    test_obj <- fisher.test(tab)
    test_used <- "Fisher's exact test"
    statistic <- NA_real_
    p_value <- test_obj$p.value
  } else {
    test_obj <- chisq.test(tab, correct = FALSE)
    test_used <- "Chi-square test"
    statistic <- unname(test_obj$statistic)
    p_value <- test_obj$p.value
  }

  tibble::tibble(
    groups_tested = paste(groups, collapse = " | "),
    pooled_median_days = pooled_median,
    test_used = test_used,
    statistic = statistic,
    df = ifelse(test_used == "Chi-square test",
                unname(test_obj$parameter),
                NA_real_),
    p_value = p_value
  )
}

age_groups <- c(
  "Age >=3 months to <6 months",
  "Age >=6 months to <1 year",
  "Age >=1 year to <5 years",
  "Age >=5 years"
)

age_moods_omnibus <- run_moods_median_omnibus(
  data = abx_duration_plot_df,
  groups = age_groups
)

age_moods_omnibus

#### Create table for manuscript ####
moods_pairwise_pretty <- moods_pairwise_results %>%
  dplyr::transmute(
    comparison,
    `Median 1 (days)` = round(median_group_1_days, 2),
    `Median 2 (days)` = round(median_group_2_days, 2),
    `Pooled median (days)` = round(pooled_median_days, 2),
    `N group 1` = n_group_1,
    `N group 2` = n_group_2,
    `Raw p` = signif(p_value, 3),
    `Holm-adjusted p` = signif(p_value_holm, 3)
  )

moods_pairwise_pretty

## Perform comparisons of antibiotic use in SBI-negative children

#--------------------------------------------------
# 1. Start from patient-level SBI-negative data
#--------------------------------------------------
abx_subgroup_df <- pros_combined_for_subgroups %>%
  dplyr::filter(sbi_present == 0) %>%
  dplyr::mutate(
    abx_after_2h = as.integer(abx_in_24h_after_2h_flag),

    # Make subgroup variables explicit / stable
    pccc_grp = factor(as.character(pccc), levels = c("0", "1"),
                      labels = c("PCCC = 0", "PCCC = 1")),

    malignancy_grp = factor(as.character(malignancy_pccc), levels = c("0", "1"),
                            labels = c("Malignancy PCCC = 0", "Malignancy PCCC = 1")),

    age_group = dplyr::case_when(
      age >= (3/12) & age < (6/12) ~ "Age >=3 months to <6 months",
      age >= (6/12) & age < 1      ~ "Age >=6 months to <1 year",
      age >= 1 & age < 5           ~ "Age >=1 year to <5 years",
      age >= 5                     ~ "Age >=5 years",
      TRUE                         ~ NA_character_
    ),
    age_group = factor(
      age_group,
      levels = c(
        "Age >=3 months to <6 months",
        "Age >=6 months to <1 year",
        "Age >=1 year to <5 years",
        "Age >=5 years"
      )
    )
  )

# Optional checks
table(abx_subgroup_df$pccc_grp, useNA = "ifany")
table(abx_subgroup_df$malignancy_grp, useNA = "ifany")
table(abx_subgroup_df$age_group, useNA = "ifany")

#--------------------------------------------------
# 2. Helper: test difference in proportions
#    Uses chi-square unless small expected counts
#--------------------------------------------------
run_prop_test <- function(data, group_var, outcome_var = "abx_after_2h") {

  dat <- data %>%
    dplyr::filter(!is.na(.data[[group_var]]), !is.na(.data[[outcome_var]])) %>%
    dplyr::mutate(
      .group = factor(.data[[group_var]]),
      .outcome = factor(.data[[outcome_var]], levels = c(0, 1))
    )

  tab <- table(dat$.group, dat$.outcome)

  chisq_obj <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  expected_min <- min(chisq_obj$expected)

  if (expected_min < 5) {
    test_obj <- stats::fisher.test(tab)
    test_used <- "Fisher's exact test"
    statistic <- NA_real_
    df <- NA_real_
    p_value <- test_obj$p.value
  } else {
    test_obj <- stats::chisq.test(tab, correct = FALSE)
    test_used <- "Chi-square test"
    statistic <- unname(test_obj$statistic)
    df <- unname(test_obj$parameter)
    p_value <- test_obj$p.value
  }

  # Build per-group counts/proportions
  group_levels <- rownames(tab)

  group_summary <- tibble::tibble(
    group = group_levels,
    n_no_abx = as.integer(tab[, "0"]),
    n_yes_abx = as.integer(tab[, "1"]),
    total = n_no_abx + n_yes_abx,
    prop_yes_abx = n_yes_abx / total
  )

  list(
    table = tab,
    summary = group_summary,
    result = tibble::tibble(
      grouping_variable = group_var,
      test_used = test_used,
      statistic = statistic,
      df = df,
      p_value = p_value
    )
  )
}

#--------------------------------------------------
# 3. PCCC: 0 vs 1
#--------------------------------------------------
pccc_test <- run_prop_test(
  data = abx_subgroup_df,
  group_var = "pccc_grp"
)

pccc_test$table
pccc_test$summary
pccc_test$result

#--------------------------------------------------
# 4. Malignancy PCCC: 0 vs 1
#--------------------------------------------------
malignancy_test <- run_prop_test(
  data = abx_subgroup_df,
  group_var = "malignancy_grp"
)

malignancy_test$table
malignancy_test$summary
malignancy_test$result

#--------------------------------------------------
# 5. Age group: omnibus test across all 4 groups
#--------------------------------------------------
age_omnibus_test <- run_prop_test(
  data = abx_subgroup_df %>% dplyr::filter(!is.na(age_group)),
  group_var = "age_group"
)

age_omnibus_test$table
age_omnibus_test$summary
age_omnibus_test$result

#--------------------------------------------------
# 6. Age group: pairwise vs reference (Age >=5 years)
#    with Holm correction
#--------------------------------------------------
age_ref <- "Age >=5 years"

age_comparisons <- tibble::tribble(
  ~comparison_group,
  "Age >=3 months to <6 months",
  "Age >=6 months to <1 year",
  "Age >=1 year to <5 years"
)

run_age_vs_ref <- function(comp_group, ref_group = age_ref, data = abx_subgroup_df) {

  dat <- data %>%
    dplyr::filter(age_group %in% c(comp_group, ref_group)) %>%
    dplyr::mutate(
      age_group = factor(age_group, levels = c(comp_group, ref_group))
    )

  test_out <- run_prop_test(dat, "age_group")

  sum_df <- test_out$summary

  tibble::tibble(
    comparison = paste0(comp_group, " vs ", ref_group),

    n_yes_comp = sum_df$n_yes_abx[sum_df$group == comp_group],
    total_comp = sum_df$total[sum_df$group == comp_group],
    prop_yes_comp = sum_df$prop_yes_abx[sum_df$group == comp_group],

    n_yes_ref = sum_df$n_yes_abx[sum_df$group == ref_group],
    total_ref = sum_df$total[sum_df$group == ref_group],
    prop_yes_ref = sum_df$prop_yes_abx[sum_df$group == ref_group],

    test_used = test_out$result$test_used,
    statistic = test_out$result$statistic,
    df = test_out$result$df,
    p_value = test_out$result$p_value
  )
}

age_pairwise_results <- purrr::map_dfr(
  age_comparisons$comparison_group,
  run_age_vs_ref
) %>%
  dplyr::mutate(
    p_value_holm = p.adjust(p_value, method = "holm")
  )

age_pairwise_results

### Get pretty table for manuscript ###
pccc_result_pretty <- pccc_test$summary %>%
  dplyr::mutate(
    prop_label = paste0(n_yes_abx, "/", total, " (",
                        scales::percent(prop_yes_abx, accuracy = 0.1), ")")
  ) %>%
  dplyr::select(group, prop_label) %>%
  tidyr::pivot_wider(names_from = group, values_from = prop_label) %>%
  dplyr::mutate(
    comparison = "PCCC = 0 vs PCCC = 1",
    p_value = pccc_test$result$p_value,
    test_used = pccc_test$result$test_used
  ) %>%
  dplyr::select(comparison, everything())

malignancy_result_pretty <- malignancy_test$summary %>%
  dplyr::mutate(
    prop_label = paste0(n_yes_abx, "/", total, " (",
                        scales::percent(prop_yes_abx, accuracy = 0.1), ")")
  ) %>%
  dplyr::select(group, prop_label) %>%
  tidyr::pivot_wider(names_from = group, values_from = prop_label) %>%
  dplyr::mutate(
    comparison = "Malignancy PCCC = 0 vs 1",
    p_value = malignancy_test$result$p_value,
    test_used = malignancy_test$result$test_used
  ) %>%
  dplyr::select(comparison, everything())

age_pairwise_pretty <- age_pairwise_results %>%
  dplyr::transmute(
    comparison,
    `Comparison group` = paste0(n_yes_comp, "/", total_comp, " (",
                                scales::percent(prop_yes_comp, accuracy = 0.1), ")"),
    `Age >=5 years` = paste0(n_yes_ref, "/", total_ref, " (",
                             scales::percent(prop_yes_ref, accuracy = 0.1), ")"),
    test_used,
    p_value,
    p_value_holm
  )

pccc_result_pretty
malignancy_result_pretty
# age_omnibus_test$result
age_pairwise_pretty

### Now combine and make pretty in a single table ###
library(dplyr)
library(tibble)
library(flextable)
library(officer)
library(stringr)

## Helper to format p-values for manuscript tables
fmt_p <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "",
    x < 0.001 ~ "<0.001",
    TRUE ~ formatC(x, format = "f", digits = 3)
  )
}

## -----------------------------
## 1. Standardize each table
## -----------------------------

pccc_tbl <- pccc_result_pretty %>%
  dplyr::transmute(
    Characteristic = "Pediatric complex chronic condition",
    Comparison = comparison,
    `Group 1` = `PCCC = 0`,
    `Group 2` = `PCCC = 1`,
    `P value` = fmt_p(p_value),
    `Adjusted P value` = "",
    `Statistical test` = test_used
  )

malignancy_tbl <- malignancy_result_pretty %>%
  dplyr::transmute(
    Characteristic = "Malignancy-related pediatric complex chronic condition",
    Comparison = comparison,
    `Group 1` = `Malignancy PCCC = 0`,
    `Group 2` = `Malignancy PCCC = 1`,
    `P value` = fmt_p(p_value),
    `Adjusted P value` = "",
    `Statistical test` = test_used
  )

age_tbl <- age_pairwise_pretty %>%
  dplyr::transmute(
    Characteristic = "Age group",
    Comparison = comparison,
    `Group 1` = `Comparison group`,
    `Group 2` = `Age >=5 years`,
    `P value` = fmt_p(p_value),
    `Adjusted P value` = fmt_p(p_value_holm),
    `Statistical test` = test_used
  )

## -----------------------------
## 2. Combine into one table
## -----------------------------

combined_results_tbl <- dplyr::bind_rows(
  pccc_tbl,
  malignancy_tbl,
  age_tbl
) %>%
  dplyr::mutate(
    Characteristic = factor(
      Characteristic,
      levels = c(
        "Pediatric complex chronic condition",
        "Malignancy-related pediatric complex chronic condition",
        "Age group"
      )
    )
  ) %>%
  dplyr::arrange(Characteristic)

combined_results_tbl

ft_results <- combined_results_tbl %>%
  flextable::flextable() %>%
  flextable::set_header_labels(
    Characteristic = "Characteristic",
    Comparison = "Comparison",
    `Group 1` = "Group 1",
    `Group 2` = "Group 2",
    `P value` = "P value",
    `Adjusted P value` = "Holm-adjusted P value",
    `Statistical test` = "Test"
  ) %>%
  flextable::merge_v(j = "Characteristic") %>%
  flextable::valign(j = "Characteristic", valign = "top") %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "left", part = "all") %>%
  flextable::autofit() %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::theme_booktabs() %>%
  flextable::set_caption(
    caption = "Association of patient characteristics with antibiotic exposure after the model prediction timepoint"
  )

ft_results

source(file = "/phi/sbi/sbi_blake/aim_1_paper_materials/create_new_pros_model.R")
