# Script to evaluate changes over time in input availability, input values, and patient mix. To be used for prospective validation paper looking at reasons why the model had poor
# performance when evaluated prospectively.


# Call setup script to import needed filepaths
source("/phi/sbi/sbi_blake/aim_1_paper_materials/setup_aim_1.R")

# Load additional scripts
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
library(pROC)
library(caret)
library(iml)
library(rsample)
library(purrr)
library(glmnet)
library(caret)
library(tidyverse)
library(survey)
library(forcats)
library(gtsummary)
library(flextable)
library(officer)
library(tidyr)
library(ggplot2)
library(scales)
library(stringr)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(pROC)
library(yardstick)
library(ggplot2)
library(patchwork)
library(pROC)      # ROC + AUC
library(PRROC)     # PR curve + AUPRC
library(WeightIt)
library(cobalt)
library(emmeans)



# Run script to generate models and rf_df and rf_df_abx
source(file = "/phi/sbi/sbi_blake/aim_1_paper_materials/retro_models_to_2026.R")

# Load in similar dataframe for the prospective model data
pros_all <- read_csv(file = "/phi/sbi/sbi_blake/pros_all_just_b4_modeling_1_15_26_all_models.csv")

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
pt_sex <- df_demographics %>% dplyr::select(MRN, SEX, ETHNICITY)
pt_sex <- df_demographics %>% dplyr::mutate(sex = as.factor(df_demographics$SEX), mrn = as.character(pt_sex$MRN), ethnicity = as.factor(pt_sex$ETHNICITY))

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


vitals_full_file_path <- "/phi/sbi/sbi_data/ten_yr_data/vitals_20220415.csv"
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
uniq_vps_mrn_picu_adm_case_id <- arrange(.data = uniq_vps_mrn_picu_adm_case_id, picu_adm_date_time)
uniq_epic_mrn_hosp_admit <- uniq_epic_mrn_hosp_admit %>% mutate(case_id = as.integer(0))

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

# Not get number of vital sign values for each case_id
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
micro_data_retro <- read_csv(file = "/phi/sbi/sbi_data/ten_yr_data/micro_20220128.csv")

# Filter to only those cultures I want

proc_list <- c("MENINGITIS ENCEPHALITIS PANEL (MEP)", "MENINGITIS ENCEPHALITIS PANEL CONDITIONAL (MEP)")
comp_list <- c("CULTURE, ANAEROBIC BLOOD", "CULTURE, BLOOD", "CULTURE, BLOOD ANAEROBIC", "CULTURE, CSF", "CULTURE, URINE", "CULTURE, URINE", "CULTURE, URINE QUANT",
               "MEP CONDITIONAL", "MEP RESULTS", "RESULTS FOR MEP")
micro_data_retro_cultures <- micro_data_retro %>% filter(PROC_NAME %in% proc_list | COMPONENT_NAME %in% comp_list) # 266k -> ~19k rows

# Now need to join to case_id via the MRN and time the culture was taken
slim_mrn_retro <- full_retro_name_fix %>% dplyr::select(case_id, mrn, picu_adm_date_time) %>% distinct()

# Add culture category
# Assign tests to culture categories
b_cat <- c("CULTURE, ANAEROBIC BLOOD", "CULTURE, BLOOD", "CULTURE, BLOOD ANAEROBIC")
u_cat <- c("CULTURE, URINE", "CULTURE, URINE", "CULTURE, URINE QUANT")
c_cat <- c( "CULTURE, CSF", "MEP CONDITIONAL", "MEP RESULTS", "RESULTS FOR MEP", "MENINGITIS ENCEPHALITIS PANEL (MEP)")

micro_data_retro_cultures <- micro_data_retro_cultures %>% mutate(cult_cat = case_when(
  COMPONENT_NAME %in% b_cat ~ "blood_cx",
  COMPONENT_NAME %in% u_cat ~ "urine_cx",
  (COMPONENT_NAME %in% c_cat | PROC_NAME %in% c_cat) ~ "csf_cx",
  TRUE ~ "none"
))

micro_slim_retro <- micro_data_retro_cultures %>% dplyr::select(MRN, cult_cat, TAKEN_TIME)

# Compare culture times to PICU admission times and identify patients who had each culture obtained
micro_slim_retro <- micro_slim_retro %>% rename(mrn = MRN)
micro_slim_retro$mrn <- as.character(micro_slim_retro$mrn)

# Need to add in mrn to the main dfs
retro_full_w_vps_w_mrn <- retro_full_w_vps %>% left_join(slim_mrn_retro %>% dplyr::select(case_id, mrn) %>% distinct(), by = c("study_id" = "case_id"))
retro_full_w_vps_w_mrn <- retro_full_w_vps_w_mrn %>% relocate(study_id, mrn, picu_adm_date_time, sbi_present, model_score, adm_month, adm_quarter, adm_year)

# Ensure micro times correct
micro_slim_retro$TAKEN_TIME <- force_tz(micro_slim_retro$TAKEN_TIME, tz = "America/Denver")

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
          inherits(micro_dt$TAKEN_TIME,   "POSIXct"))

# Create a stable row id so we can aggregate back per admission row
retro_dt[, row_id := .I]

# Index micro for speed (doesn't reorder like setkey)
setindex(micro_dt, mrn, TAKEN_TIME)

# Non-equi join: cultures within window for each admission row
matched <- micro_dt[
  retro_dt,
  on = .(mrn,
         TAKEN_TIME >= window_start,
         TAKEN_TIME <= window_end),
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
          inherits(micro_dt$TAKEN_TIME,   "POSIXct"))

# Create a stable row id so we can aggregate back per admission row
retro_dt_abx[, row_id := .I]

# Index micro for speed (doesn't reorder like setkey)
setindex(micro_dt, mrn, TAKEN_TIME)

# Non-equi join: cultures within window for each admission row
matched_abx <- micro_dt[
  retro_dt_abx,
  on = .(mrn,
         TAKEN_TIME >= window_start,
         TAKEN_TIME <= window_end),
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
micro_data_retro$TAKEN_TIME <- force_tz(micro_data_retro$TAKEN_TIME, tz = "America/Denver")
micro_data_retro <- micro_data_retro %>% rename(mrn = MRN)
micro_data_retro$mrn <- as.character(micro_data_retro$mrn)


# 1) Create a skinny CRP event table (unique per MRN + time)
crp_events <- micro_data_retro %>%
  filter(PROC_NAME == "C-REACTIVE PROTEIN") %>%
  transmute(
    mrn = as.character(mrn),
    taken_time = TAKEN_TIME
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
  filter(PROC_NAME == "PROCALCITONIN") %>%
  transmute(
    mrn = as.character(mrn),
    taken_time = TAKEN_TIME
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
no_crp_micro_retro <- micro_data_retro %>% filter(!(PROC_NAME %in% c("C-REACTIVE PROTEIN", "CHOLESTEROL, PLEURAL FLUID", "GLUCOSE, CSF", "LDH, PLEURAL FLUID",
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
  micro_dt[, t   := as.POSIXct(TAKEN_TIME, tz = "America/Denver")]
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
glimpse(labs_pros)
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


# Calculate score of retro model oin prospective data
source(file = "/phi/sbi/sbi_blake/aim_1_paper_materials/retro_model_to_pros_data.R")
source(file = "/phi/sbi/sbi_blake/aim_1_paper_materials/retro_model_to_pros_susp_infxn.R")



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

# Prospective: unexposed (pros_no_abx_1st_infxn has sbi_present + model_score)
# NOTE: your pros sbi_present is 0/1
pros_no_abx <- pros_no_abx_1st_infxn %>%
  transmute(
    scenario  = "Pros • Abx-",
    truth_num = as.integer(sbi_present),
    truth     = factor(if_else(truth_num == 1L, "yes", "no"), levels = c("yes","no")),
    score     = to_prob01(model_score)
  )

# Prospective: exposed
pros_yes_abx <- pros_yes_abx_1st_infxn %>%
  transmute(
    scenario  = "Pros • Abx+",
    truth_num = as.integer(sbi_present),
    truth     = factor(if_else(truth_num == 1L, "yes", "no"), levels = c("yes","no")),
    score     = to_prob01(model_score)
  )

dat_all <- bind_rows(reto_yes_abx = retro_yes_abx,
                     reto_no_abx  = retro_no_abx,
                     pros_yes_abx = pros_yes_abx,
                     pros_no_abx  = pros_no_abx) %>%
  filter(!is.na(truth_num), !is.na(score), is.finite(score)) %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c("Retro • Abx+", "Retro • Abx-", "Pros • Abx+", "Pros • Abx-")
    )
  )

dat_list <- split(dat_all, dat_all$scenario)

# ----------------------------
# 2) Curves + metrics
# ----------------------------

roc_curve_df <- function(df_one) {
  if (dplyr::n_distinct(df_one$truth_num) < 2) {
    return(tibble(scenario = df_one$scenario[1], fpr = NA_real_, tpr = NA_real_))
  }
  roc_obj <- pROC::roc(df_one$truth_num, df_one$score, quiet = TRUE, direction = "<")
  tibble(
    scenario = df_one$scenario[1],
    fpr      = 1 - roc_obj$specificities,
    tpr      = roc_obj$sensitivities
  )
}

auroc_one <- function(df_one) {
  roc_obj <- pROC::roc(df_one$truth_num, df_one$score, quiet = TRUE, direction = "<")
  as.numeric(pROC::auc(roc_obj))
}

pr_curve_df <- function(df_one) {
  if (dplyr::n_distinct(df_one$truth_num) < 2) {
    return(tibble(scenario = df_one$scenario[1], recall = NA_real_, precision = NA_real_))
  }

  pc <- yardstick::pr_curve(df_one, truth, score, event_level = "first")

  # robust column extraction across yardstick versions
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
  yardstick::pr_auc(df_one, truth, score, event_level = "first") %>%
    pull(.estimate) %>%
    as.numeric()
}

roc_df <- map_dfr(dat_list, roc_curve_df)
pr_df  <- map_dfr(dat_list, pr_curve_df)

metrics_df <- map_dfr(dat_list, ~ tibble(
  scenario   = .x$scenario[1],
  n          = nrow(.x),
  prevalence = mean(.x$truth_num == 1L),
  AUROC      = auroc_one(.x),
  AUPRC      = auprc_one(.x)
)) %>%
  mutate(
    scenario = factor(scenario, levels = levels(dat_all$scenario)),
    label = paste0(
      "AUROC = ", sprintf("%.3f", AUROC), "\n",
      "AUPRC = ", sprintf("%.3f", AUPRC), "\n",
      "Prev = ", sprintf("%.3f", prevalence)
    )
  )

# ----------------------------
# 3) Plots
# ----------------------------

roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_line(linewidth = 1) +
  facet_wrap(~ scenario, ncol = 2) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "ROC Curves (Retrospective TEST + Prospective)",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  geom_text(
    data = metrics_df,
    aes(x = 0.62, y = 0.10, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.5
  ) +
  theme_bw()

pr_plot <- ggplot(pr_df, aes(x = recall, y = precision)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ scenario, ncol = 2) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Precision–Recall Curves (Retrospective TEST + Prospective)",
    x = "Recall",
    y = "Precision"
  ) +
  geom_text(
    data = metrics_df,
    aes(x = 0.62, y = 0.10, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.5
  ) +
  theme_bw()

print(roc_plot)
print(pr_plot)
metrics_df



# Ensure correct scale for model scores
rescale_to_unit <- function(x) {
  if (is.numeric(x) && max(x, na.rm = TRUE) > 1) {
    x / 100
  } else {
    x
  }
}

retro_no_abx_1st_infxn$model_score <-
  rescale_to_unit(retro_no_abx_1st_infxn$model_score)

retro_yes_abx_1st_infxn$model_score <-
  rescale_to_unit(retro_yes_abx_1st_infxn$model_score)

pros_no_abx_1st_infxn$model_score <-
  rescale_to_unit(pros_no_abx_1st_infxn$model_score)

pros_yes_abx_1st_infxn$model_score <-
  rescale_to_unit(pros_yes_abx_1st_infxn$model_score)


# ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ----------------------------
# Inputs: my 4 evaluation dataframes
# (these are the prospective/retrospective encounter-level eval frames)
# ----------------------------
dfs_eval <- list(
  "Retro • Abx-" = retro_no_abx_1st_infxn,
  "Retro • Abx+" = retro_yes_abx_1st_infxn,
  "Pros • Abx-"  = pros_no_abx_1st_infxn,
  "Pros • Abx+"  = pros_yes_abx_1st_infxn
)

score_col <- "model_score"
y_col     <- "sbi_present"

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

ruleout_stats_one <- function(df, scenario, score_col, y_col, threshold) {
  sc <- parse_scenario(scenario)

  dat <- df %>%
    transmute(
      score = as.numeric(.data[[score_col]]),
      truth = to_truth01(.data[[y_col]])
    ) %>%
    filter(!is.na(score), !is.na(truth), truth %in% c(0L, 1L))

  n_total <- nrow(dat)
  n_sbi0  <- sum(dat$truth == 0L)
  n_sbi1  <- sum(dat$truth == 1L)

  low <- dat$score < threshold

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

quad_df <- purrr::imap_dfr(dfs_eval, ~{
  abx_lbl <- ifelse(grepl("Abx\\+", .y), "Abx+", "Abx-")
  thr <- thr_tbl$threshold[thr_tbl$abx == abx_lbl][1]
  ruleout_stats_one(.x, scenario = .y, score_col = score_col, y_col = y_col, threshold = thr)
})

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
      TRUE ~ paste0("No patients < thr\nlow-risk n=0")
    )
  )

seg_df <- quad_df %>%
  select(abx, epoch, x_plot, y_plot) %>%
  tidyr::pivot_wider(names_from = epoch, values_from = c(x_plot, y_plot)) %>%
  filter(!is.na(x_plot_Retro), !is.na(x_plot_Pros)) %>%
  transmute(
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
    alpha = 0.6
  ) +

  # reference lines (edit/remove if you prefer)
  geom_hline(yintercept = 0.95, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 0.05, linetype = 2, alpha = 0.5) +

  # points: hollow + lower alpha when n_low==0
  geom_point(
    aes(shape = has_low, alpha = has_low),
    size = 4,
    stroke = 1.1
  ) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45), guide = "none") +

  geom_text_repel(
    aes(label = note),
    size = 3.2,
    show.legend = FALSE,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.alpha = 0.5,
    seed = 1
  ) +

  facet_wrap(~abx, ncol = 2) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +

  # IMPORTANT: zoom without dropping points
  coord_cartesian(ylim = c(y_lo, y_hi)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.5)) +

  labs(
    title = "Clinical rule-out performance at fixed SBI risk thresholds",
    subtitle = paste0(
      "Y-axis: NPV among low-risk; X-axis: fraction of SBI-negative encounters captured as low-risk\n",
      "Hollow points indicate zero low-risk predictions at the specified threshold."
    ),
    x = "Proportion of SBI-negative encounters classified as low risk",
    y = "Negative predictive value (NPV)",
    color = "Epoch"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

p_quadrant




# Now plot the distribution of scores
# ----------------------------
# Inputs
# ----------------------------
dfs_eval <- list(
  "Retro • Abx-" = retro_no_abx_1st_infxn,
  "Retro • Abx+" = retro_yes_abx_1st_infxn,
  "Pros • Abx-"  = pros_no_abx_1st_infxn,
  "Pros • Abx+"  = pros_yes_abx_1st_infxn
)

score_col <- "model_score"   # numeric 0-1
y_col     <- "sbi_present"   # not used here, but kept for later just in case

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
    pct_below = mean(score < threshold),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "N=", n, "\n",
      "% < thr=", percent(pct_below, accuracy = 0.1)
    )
  )

# ----------------------------
# Plot: histogram + density + vertical threshold line
# ----------------------------
p_dist <-
  ggplot(score_df, aes(x = score)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 25, alpha = 0.35, linewidth = 0.2) +
  geom_density(linewidth = 1) +
  geom_vline(aes(xintercept = threshold),
             linetype = "dashed", linewidth = 1) +
  geom_text(
    data = annot_df,
    aes(x = 0.98, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1, vjust = 1.2,
    size = 3.3
  ) +
  facet_grid(epoch ~ abx) +
  scale_x_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "Distribution of predicted SBI probability by epoch and antibiotic exposure",
    subtitle = "Dashed line = fixed rule-out threshold (Abx-: 0.05, Abx+: 0.074)",
    x = "Predicted SBI probability",
    y = "Density"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
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

p_dist_by_truth <-
  ggplot(score_truth_df, aes(x = score, fill = truth)) +
  geom_density(
    alpha = 0.35,                 # <- transparency so overlap looks purple
    linewidth = 0.2,
    adjust = 1                    # tweak if you want more/less smoothing
  ) +
  geom_vline(
    data = distinct(score_truth_df, epoch, abx, threshold),
    aes(xintercept = threshold),
    linetype = "dashed", linewidth = 1,
    inherit.aes = FALSE
  ) +
  facet_grid(epoch ~ abx) +
  scale_x_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_fill_manual(values = c("SBI-" = "blue", "SBI+" = "red")) +
  labs(
    title = "Distribution of predicted SBI probability by epoch, antibiotic exposure, and true SBI status",
    subtitle = "Dashed line = fixed rule-out threshold (Abx-: 0.05, Abx+: 0.074)",
    x = "Predicted SBI probability",
    y = "Density",
    fill = "True outcome"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
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


#Standardize inputs and create low_risk flag
prep_epoch_df <- function(df, epoch_label, cutoff_prob) {

  # scalar decision, computed once
  needs_pct_to_prob <- max(df$model_score, na.rm = TRUE) > 1

  df %>%
    dplyr::mutate(
      epoch = epoch_label,
      A = dplyr::if_else(epoch == "pros", 1L, 0L),

      model_prob = if (needs_pct_to_prob) model_score / 100 else model_score,

      low_risk = dplyr::if_else(!is.na(model_prob) & model_prob <= cutoff_prob, 1L, 0L),

      is_female = as.factor(is_female),
      race = as.factor(race),
      ethnicity = as.factor(ethnicity),
      pccc = as.factor(pccc),
      malignancy_pccc = as.factor(malignancy_pccc)
    )
}


# Build defensible covariate set
get_ps_covariates <- function(df, outcome_vars = character()) {
  # hard excludes (never in PS model)
  exclude <- c(
    "study_id","mrn","picu_adm_date_time",
    "epoch","A",
    "sbi_present",
    "model_score","model_prob","low_risk",
    outcome_vars
  )

  # candidate pool: everything else
  cand <- setdiff(names(df), exclude)

  # strongly recommended baseline set if present
  force_in <- intersect(
    c("age","is_female","race","ethnicity","pccc","malignancy_pccc",
      "los_before_icu_days","imv_at_picu_adm"),
    cand
  )

  # include preicu_* one-hot columns if present
  preicu_cols <- grep("^preicu_", cand, value = TRUE)

  # final covariate list
  unique(c(force_in, preicu_cols))
  print(unique(c(force_in, preicu_cols)))
}

get_ps_covariates_case_mix <- function(df) {

  exclude <- c(
    "study_id","mrn","picu_adm_date_time",
    "epoch","A",
    "sbi_present",
    "model_score","model_prob","low_risk"
  )

  cand <- setdiff(names(df), exclude)

  # explicitly remove "management/ascertainment" fields from case-mix PS
  mgmt_like <- intersect(
    c("crp_pres","pct_pres","cxr_pres","any_micro_pres",
      "blood_y_n","urine_y_n","csf_y_n",
      "suspected_infection"),
    cand
  )

  cand <- setdiff(cand, mgmt_like)

  # include demographic/comorbidity and pre-ICU location (if present)
  force_in <- intersect(
    c("age","is_female","race","ethnicity","pccc","malignancy_pccc",
      "los_before_icu_days","imv_at_picu_adm"),
    cand
  )

  # include preicu one-hot columns
  preicu_cols <- grep("^preicu_", cand, value = TRUE)

  # include a broad set of baseline phys summaries if present
  phys_like <- grep("^(hr_|sbp_|dbp_|rr_|temp_|o2sat_|fio2_|lactate_|hematocrit_|pco2_)", cand, value = TRUE)

  unique(c(force_in, preicu_cols, phys_like))
}

make_weights_pros_to_retro <- function(df, ps_trim = c(0.01, 0.99), weight_cap = 20) {

  covars <- get_ps_covariates_case_mix(df)
  ps_form <- as.formula(paste("A ~", paste(covars, collapse = " + ")))

  # ATT with treated = Pros (A=1) — NOTE: no stabilize arg
  wfit <- WeightIt::weightit(
    formula  = ps_form,
    data     = df,
    method   = "ps",
    estimand = "ATT"
  )

  out <- df %>%
    dplyr::mutate(
      w_raw = wfit$weights,
      ps    = wfit$ps
    ) %>%
    dplyr::filter(!is.na(ps), ps >= ps_trim[1], ps <= ps_trim[2]) %>%
    dplyr::mutate(w = pmin(w_raw, weight_cap))

  bal <- cobalt::bal.tab(wfit, un = TRUE, m.threshold = 0.1)
  print(bal)

  list(df_w = out, wfit = wfit, balance = bal, covariates = covars)
}


w_mean <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE)

ruleout_metrics <- function(df, w_col = NULL) {

  w <- if (is.null(w_col)) rep(1, nrow(df)) else df[[w_col]]

  # basic quantities
  y <- df$sbi_present
  lr <- df$low_risk

  # weights only where both defined
  ok <- !is.na(y) & !is.na(lr) & !is.na(w)
  y <- y[ok]; lr <- lr[ok]; w <- w[ok]

  # low-risk yield
  p_low_all <- w_mean(lr == 1, w)

  # yield among SBI-negative
  neg <- (y == 0)
  p_low_among_neg <- if (sum(w[neg]) > 0) w_mean(lr[neg] == 1, w[neg]) else NA_real_

  # NPV among low-risk (P(Y=0 | low_risk=1))
  lr1 <- (lr == 1)
  npv <- if (sum(w[lr1]) > 0) w_mean(y[lr1] == 0, w[lr1]) else NA_real_

  # false negatives among low-risk (Y=1 & low_risk=1)
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


# Helper function for running weightint plus diagnostics + one outcome model
run_iptw_epoch_analysis <- function(df, outcome, estimand = "ATE",
                                    ps_trim = c(0.01, 0.99),
                                    weight_cap = NULL) {

  stopifnot(outcome %in% names(df))

  covars <- get_ps_covariates(df, outcome_vars = outcome)

  ps_form <- as.formula(
    paste("A ~", paste(covars, collapse = " + "))
  )

  # propensity weighting (logistic PS by default)
  wfit <- if (estimand == "ATE") {
    weightit(
      formula  = ps_form,
      data     = df,
      method   = "ps",
      estimand = estimand,
      stabilize = TRUE
    )
  } else {
    weightit(
      formula  = ps_form,
      data     = df,
      method   = "ps",
      estimand = estimand
    )
  }


  out <- df %>%
    mutate(
      w = wfit$weights,
      ps = wfit$ps
    ) %>%
    # trim on propensity score for overlap (optional but usually helpful)
    filter(!is.na(ps), ps >= ps_trim[1], ps <= ps_trim[2])

  # optional truncation/capping of extreme weights
  if (!is.null(weight_cap)) {
    out <- out %>% mutate(w = pmin(w, weight_cap))
  }

  # balance diagnostics
  bal <- bal.tab(wfit, un = TRUE, m.threshold = 0.1)
  print(bal)

  # quick love plot (run interactively)
  # love.plot(bal, threshold = 0.1, abs = TRUE, var.order = "unadjusted")

  # survey design + weighted model
  des <- svydesign(ids = ~1, data = out, weights = ~w)

  # binary outcome model: effect of epoch (A)
  # (quasibinomial is stable with weights)
  fit <- svyglm(as.formula(paste(outcome, "~ A")), design = des, family = quasibinomial())

  # marginal weighted risks by epoch (on probability scale)
  emm <- emmeans(fit, ~ A, type = "response")

  # TRUE risk difference (pros - retro), on probability scale
  rd <- contrast(
    emm,
    method = list("RD_pros_minus_retro" = c(-1, 1)),
    type = "response"
  )

  # Risk ratio computed manually from the marginal risks
  emm_df <- as.data.frame(emm)
  prob_col <- intersect(c("prob", "response"), names(emm_df))[1]

  p0 <- emm_df[[prob_col]][emm_df$A == 0]
  p1 <- emm_df[[prob_col]][emm_df$A == 1]
  rr <- p1 / p0


  list(
    df = out,
    covariates = covars,
    weightit = wfit,
    balance = bal,
    model = fit,
    marginal_risks = emm,
    risk_difference = rd,
    risk_ratio = rr
  )
}

make_weights_pros_to_retro <- function(df,
                                       estimand = "ATT",
                                       ps_trim = c(0.01, 0.99),
                                       weight_cap = NULL,
                                       covariate_fun = get_ps_covariates) {

  # Create A_retro so that "treated" = Retro (keeps w=1),
  # and Pros gets weighted to look like Retro.
  df2 <- df %>%
    mutate(
      A_retro = if_else(epoch == "retro", 1L, 0L)
    )

  covars <- covariate_fun(df2)

  ps_form <- as.formula(paste("A_retro ~", paste(covars, collapse = " + ")))

  wfit <- weightit(
    formula  = ps_form,
    data     = df2,
    method   = "ps",
    estimand = estimand
  )

  df_w <- df2 %>%
    mutate(
      ps = wfit$ps,
      w  = wfit$weights
    ) %>%
    filter(!is.na(ps), ps >= ps_trim[1], ps <= ps_trim[2])

  if (!is.null(weight_cap)) {
    df_w <- df_w %>% mutate(w = pmin(w, weight_cap))
  }

  # diagnostics
  bal <- bal.tab(wfit, un = TRUE, m.threshold = 0.1)
  print(bal)

  list(
    df_w = df_w,
    weightit = wfit,
    balance = bal,
    covariates = covars
  )
}


# Now run on my two strata
df_noabx <- bind_rows(
  prep_epoch_df(retro_no_abx_1st_infxn, epoch_label = "retro", cutoff_prob = 0.05),
  prep_epoch_df(pros_no_abx_1st_infxn,     epoch_label = "pros",  cutoff_prob = 0.05)
)

# Fix sbi_present class
retro_yes_abx_1st_infxn <- retro_yes_abx_1st_infxn %>% mutate(sbi_present =
                            ifelse(sbi_present == "yes", yes = 1, no = 0))

# Now run on my two strata
df_yesabx <- bind_rows(
  prep_epoch_df(retro_yes_abx_1st_infxn, epoch_label = "retro", cutoff_prob = 0.074),
  prep_epoch_df(pros_yes_abx_1st_infxn,     epoch_label = "pros",  cutoff_prob = 0.074)
)


# 1) create weights
w_noabx <- make_weights_pros_to_retro(df_noabx, ps_trim = c(0.01, 0.99), weight_cap = 20)
df_noabx_w <- w_noabx$df_w

retro_unw <- df_noabx_w %>% filter(A == 0) %>% ruleout_metrics()

pros_unw  <- df_noabx_w %>% filter(A == 1) %>% ruleout_metrics()

# This is now Pros weighted to look like Retro
pros_wtd  <- df_noabx_w %>%
  filter(A == 1) %>%
  ruleout_metrics(w_col = "w")

bind_rows(
  retro_unw %>% mutate(epoch = "Retro", weighting = "Unweighted"),
  pros_unw  %>% mutate(epoch = "Pros",  weighting = "Unweighted"),
  pros_wtd  %>% mutate(epoch = "Pros",  weighting = "Weighted to Retro case-mix")
)

w_yesabx <- make_weights_pros_to_retro(df_yesabx, ps_trim = c(0.01, 0.99), weight_cap = 20)
df_yesabx_w <- w_yesabx$df_w

retro_unw <- df_yesabx_w %>% filter(A == 0) %>% ruleout_metrics()
pros_unw  <- df_yesabx_w %>% filter(A == 1) %>% ruleout_metrics()

pros_wtd  <- df_yesabx_w %>%
  filter(A == 1) %>%
  ruleout_metrics(w_col = "w")

bind_rows(
  retro_unw %>% mutate(epoch = "Retro", weighting = "Unweighted"),
  pros_unw  %>% mutate(epoch = "Pros",  weighting = "Unweighted"),
  pros_wtd  %>% mutate(epoch = "Pros",  weighting = "Weighted to Retro case-mix")
)

library(pROC)

weighted_auroc <- function(df, truth_col = "sbi_present", score_col = "model_prob", w_col = "w") {
  d <- df %>% filter(!is.na(.data[[truth_col]]), !is.na(.data[[score_col]]), !is.na(.data[[w_col]]))
  r <- pROC::roc(response = d[[truth_col]], predictor = d[[score_col]],
                 weights = d[[w_col]], quiet = TRUE, direction = "<")
  as.numeric(pROC::auc(r))
}

# Abx− Pros weighted to Retro
df_noabx_w %>% filter(A == 1) %>% weighted_auroc(truth_col="sbi_present", score_col="model_prob", w_col="w")

# Abx+ Pros weighted to Retro
df_yesabx_w %>% filter(A == 1) %>% weighted_auroc(truth_col="sbi_present", score_col="model_prob", w_col="w")

library(PRROC)

weighted_auprc <- function(df, truth_col="sbi_present", score_col="model_prob", w_col="w") {
  d <- df %>% filter(!is.na(.data[[truth_col]]), !is.na(.data[[score_col]]), !is.na(.data[[w_col]]))
  scores_pos <- d[[score_col]][d[[truth_col]] == 1]
  scores_neg <- d[[score_col]][d[[truth_col]] == 0]
  w_pos <- d[[w_col]][d[[truth_col]] == 1]
  w_neg <- d[[w_col]][d[[truth_col]] == 0]
  PRROC::pr.curve(scores.class0 = scores_pos, scores.class1 = scores_neg,
                  weights.class0 = w_pos, weights.class1 = w_neg,
                  curve = FALSE)$auc.integral
}

weighted_demo_pros_auroc_no_abx <- df_noabx_w %>% filter(A == 1) %>% weighted_auroc(truth_col="sbi_present", score_col="model_prob", w_col="w")
weighted_demo_pros_auroc_yes_abx <-  df_yesabx_w %>% filter(A == 1) %>% weighted_auroc(truth_col="sbi_present", score_col="model_prob", w_col="w")

weighted_demo_pros_auprc_no_abx <- df_noabx_w %>% filter(A == 1) %>% weighted_auprc(truth_col="sbi_present", score_col="model_prob", w_col="w")
weighted_demo_pros_auprc_yes_abx <- df_yesabx_w %>% filter(A == 1) %>% weighted_auprc(truth_col="sbi_present", score_col="model_prob", w_col="w")

weighted_demo_pros_auroc_no_abx
weighted_demo_pros_auroc_yes_abx
weighted_demo_pros_auprc_no_abx
weighted_demo_pros_auprc_yes_abx


#### Plot weighted and unweighted pros charts
# ============================================================
# Weighted vs unweighted Pros ROC + PR (2-panel)
# Assumes you already created:
#   df_noabx_w  (from w_noabx$df_w)  with columns: epoch, A, sbi_present (0/1), model_prob (0-1), w
#   df_yesabx_w (from w_yesabx$df_w) with same columns
# And you want Pros only (epoch=="pros" or A==1), comparing:
#   - Unweighted (w = 1)
#   - Weighted (w = IPTW weights)
# ============================================================

library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(pROC)
library(PRROC)
library(patchwork)

# ----------------------------
# Helper: build ROC points (weighted or unweighted)
# ----------------------------
# ------------------------------------------------------------
# Helper: ROC curve points
# ------------------------------------------------------------
roc_points_prroc <- function(df, truth_col="sbi_present", score_col="model_prob", w_col=NULL) {
  d <- df %>%
    filter(!is.na(.data[[truth_col]]), !is.na(.data[[score_col]])) %>%
    mutate(
      y = as.integer(.data[[truth_col]]),
      s = as.numeric(.data[[score_col]]),
      w = if (is.null(w_col)) rep(1, n()) else as.numeric(.data[[w_col]])
    ) %>%
    filter(is.finite(s), !is.na(w), w > 0)

  scores_pos <- d$s[d$y == 1]
  scores_neg <- d$s[d$y == 0]
  w_pos <- d$w[d$y == 1]
  w_neg <- d$w[d$y == 0]

  roc <- if (is.null(w_col)) {
    PRROC::roc.curve(scores.class0=scores_pos, scores.class1=scores_neg, curve=TRUE)
  } else {
    PRROC::roc.curve(scores.class0=scores_pos, scores.class1=scores_neg,
                     weights.class0=w_pos, weights.class1=w_neg, curve=TRUE)
  }

  tibble(
    fpr = roc$curve[,1],
    tpr = roc$curve[,2],
    auc = roc$auc
  )
}

pr_points_prroc <- function(df, truth_col="sbi_present", score_col="model_prob", w_col=NULL) {
  d <- df %>%
    filter(!is.na(.data[[truth_col]]), !is.na(.data[[score_col]])) %>%
    mutate(
      y = as.integer(.data[[truth_col]]),
      s = as.numeric(.data[[score_col]]),
      w = if (is.null(w_col)) rep(1, n()) else as.numeric(.data[[w_col]])
    ) %>%
    filter(is.finite(s), !is.na(w), w > 0)

  scores_pos <- d$s[d$y == 1]
  scores_neg <- d$s[d$y == 0]
  w_pos <- d$w[d$y == 1]
  w_neg <- d$w[d$y == 0]

  pr <- if (is.null(w_col)) {
    PRROC::pr.curve(scores.class0=scores_pos, scores.class1=scores_neg, curve=TRUE)
  } else {
    PRROC::pr.curve(scores.class0=scores_pos, scores.class1=scores_neg,
                    weights.class0=w_pos, weights.class1=w_neg, curve=TRUE)
  }

  tibble(
    recall = pr$curve[,1],
    precision = pr$curve[,2],
    auc = pr$auc.integral
  )
}

plot_weighted_vs_unweighted_pros <- function(df_w, title_suffix="", w_col="w") {
  d_pros <- df_w %>% filter(epoch=="pros", A==1)

  roc_unw <- roc_points_prroc(d_pros, w_col=NULL) %>% mutate(kind="Unweighted")
  roc_w   <- roc_points_prroc(d_pros, w_col=w_col) %>% mutate(kind="Weighted")
  roc_df  <- bind_rows(roc_unw, roc_w) %>%
    mutate(kind = factor(kind, levels=c("Unweighted","Weighted")))

  pr_unw <- pr_points_prroc(d_pros, w_col=NULL) %>% mutate(kind="Unweighted")
  pr_w   <- pr_points_prroc(d_pros, w_col=w_col) %>% mutate(kind="Weighted")
  pr_df  <- bind_rows(pr_unw, pr_w) %>%
    mutate(kind = factor(kind, levels=c("Unweighted","Weighted")))

  # annotations
  roc_lab <- tibble(
    kind = factor(c("Unweighted","Weighted"), levels=c("Unweighted","Weighted")),
    label = sprintf("AUROC=%.3f", c(unique(roc_unw$auc)[1], unique(roc_w$auc)[1]))
  )
  pr_lab <- tibble(
    kind = factor(c("Unweighted","Weighted"), levels=c("Unweighted","Weighted")),
    label = sprintf("AUPRC=%.3f", c(unique(pr_unw$auc)[1], unique(pr_w$auc)[1]))
  )

  lt_vals <- c("Unweighted"="solid", "Weighted"="dashed")

  p_roc <- ggplot(roc_df, aes(fpr, tpr, linetype=kind)) +
    geom_abline(intercept=0, slope=1, linetype=2) +
    geom_line(linewidth=1) +
    coord_equal(xlim=c(0,1), ylim=c(0,1)) +
    scale_linetype_manual(values=lt_vals, limits=names(lt_vals), drop=FALSE) +
    theme_bw() +
    labs(title=paste0("Prospective ROC (", title_suffix, ")"),
         x="False positive rate", y="True positive rate", linetype="Pros weighting") +
    geom_text(data=roc_lab, aes(x=0.65, y=0.10, label=label),
              inherit.aes=FALSE, hjust=0, size=3.5)

  p_pr <- ggplot(pr_df, aes(recall, precision, linetype=kind)) +
    geom_line(linewidth=1) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    scale_linetype_manual(values=lt_vals, limits=names(lt_vals), drop=FALSE) +
    theme_bw() +
    labs(title=paste0("Prospective PR (", title_suffix, ")"),
         x="Recall", y="Precision", linetype="Pros weighting") +
    geom_text(data=pr_lab, aes(x=0.65, y=0.10, label=label),
              inherit.aes=FALSE, hjust=0, size=3.5)

  p_roc | p_pr
}

p_yesabx <- plot_weighted_vs_unweighted_pros(df_yesabx_w, title_suffix="Abx+")
p_noabx  <- plot_weighted_vs_unweighted_pros(df_noabx_w,  title_suffix="Abx-")

p_yesabx
p_noabx




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

  # Ensure numeric 0/1 (some may be dbl already; just be safe)
  df <- df %>%
    mutate(across(all_of(pre_cols), ~ as.numeric(.)))

  # For each row, find which preicu_* column is "1"
  loc <- apply(df[, pre_cols, drop = FALSE], 1, function(x) {
    idx <- which(x == 1)
    if (length(idx) == 0) return(NA_character_)          # none flagged
    if (length(idx) > 1) return("Multiple/Conflicting")  # >1 flagged (rare but possible)
    pre_cols[idx]
  })

  df %>%
    mutate(
      preicu_location = loc,
      preicu_location = str_remove(preicu_location, paste0("^", prefix)),
      preicu_location = fct_explicit_na(as.factor(preicu_location), na_level = "Missing")
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
      # make is_female prettier if it's 0/1 stored as factor
      is_female = as.character(is_female),
      is_female = case_when(
        is_female %in% c("1", "TRUE", "T") ~ "Female",
        is_female %in% c("0", "FALSE", "F") ~ "Male",
        TRUE ~ NA_character_
      ),
      is_female = fct_explicit_na(as.factor(is_female), na_level = "Missing")
    )

  # Variables to include (PS-model variables + preicu location)
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

  # Keep only vars that exist in df (defensive)
  vars <- intersect(vars, names(df))

  tab <- df %>%
    dplyr::select(all_of(c("epoch", vars))) %>%
    tbl_summary(
      by = epoch,
      statistic = list(
        all_continuous() ~ "{median} [{p25}, {p75}]",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits = all_continuous() ~ 2,
      missing = "ifany"
    ) %>%
    modify_header(label ~ "**Characteristic**") %>%
    modify_spanning_header(all_stat_cols() ~ "**Epoch**") %>%
    bold_labels()

  if (!is.null(title)) {
    tab <- tab %>% modify_caption(paste0("**", title, "**"))
  }

  tab
}

# -----------------------------
# 3) Create the two Table 1s
# -----------------------------
table1_noabx  <- make_table1_epoch(df_noabx,  title = "Table 1. Baseline characteristics by epoch — No antibiotics stratum")
table1_yesabx <- make_table1_epoch(df_yesabx, title = "Table 1. Baseline characteristics by epoch — Antibiotics stratum")


#### Ok now we look at missingness and distribution of key inputs between retro and pros.
# =============================================================================
# Dataset shift analysis (retro vs prospective), stratified by abx exposure
# NOW:
#   - DOES NOT write CSVs or save PNGs (those lines are commented out)
#   - Returns objects into the R environment
#   - Prints the shift tables and prints the plots to the active device
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(readr)
})

# ----------------------------
# USER SETTINGS
# ----------------------------

EXCLUDE_COLS <- c(
  "study_id", "mrn", "picu_adm_date_time",
  "sbi_present", "model_score", "epoch",
  "adm_month", "adm_quarter", "adm_year", "adm_ym", "adm_yq",
  "pna_1_0", "blood", "cns", "urine", "ever_cx_neg_sepsis"
)

FORCE_CATEGORICAL <- c("race", "ethnicity", "is_female")

TOP_N_SMD_PLOT        <- 30
TOP_N_AVAIL_PLOT      <- 30
TOP_N_DENSITY_FEATURE <- 6

# OUT_DIR <- "shift_outputs"
# dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

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
# Plot helpers (RETURN ggplot objects + print them)
# ----------------------------
plot_density_overlays <- function(df_retro, df_pros, shift_tbl, stratum_name) {
  top_vars <- shift_tbl %>%
    filter(type == "continuous", is.finite(abs_smd)) %>%
    slice_head(n = TOP_N_DENSITY_FEATURE) %>%
    pull(variable)

  if (length(top_vars) == 0) return(NULL)

  long <- bind_rows(
    df_retro %>% select(all_of(top_vars)) %>% mutate(epoch2 = "retro"),
    df_pros  %>% select(all_of(top_vars)) %>% mutate(epoch2 = "prospective")
  ) %>%
    pivot_longer(cols = all_of(top_vars), names_to = "variable", values_to = "value") %>%
    filter(is.finite(value))

  ggplot(long, aes(x = value, group = epoch2)) +
    geom_density(aes(linetype = epoch2), linewidth = 0.7, adjust = 1) +
    facet_wrap(~ variable, scales = "free", ncol = 2) +
    labs(
      title = paste0("Density overlays (top shifted continuous predictors) — ", stratum_name),
      x = NULL, y = "Density", linetype = "Epoch"
    ) +
    theme_bw(base_size = 12)
}

plot_availability_diff <- function(shift_tbl, stratum_name, min_abs_diff = 0.01) {

  top <- shift_tbl %>%
    mutate(abs_avail_diff = abs(avail_diff)) %>%
    filter(is.finite(abs_avail_diff), abs_avail_diff >= min_abs_diff) %>%
    arrange(desc(abs_avail_diff)) %>%
    slice_head(n = TOP_N_AVAIL_PLOT) %>%
    mutate(variable = factor(variable, levels = rev(variable)))

  if (nrow(top) == 0) {
    message("No availability diffs >= ", min_abs_diff, " in ", stratum_name)
    return(NULL)
  }

  ggplot(top, aes(x = variable, y = avail_diff)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste0("Availability difference (prospective − retro) — ", stratum_name),
      subtitle = paste0("Showing variables with |diff| ≥ ", min_abs_diff),
      x = NULL,
      y = "Difference in availability"
    ) +
    theme_bw(base_size = 12)
}


plot_lollipop_abs_smd <- function(shift_tbl, stratum_name) {
  top <- shift_tbl %>%
    filter(is.finite(abs_smd)) %>%
    slice_head(n = TOP_N_SMD_PLOT) %>%
    mutate(variable = factor(variable, levels = rev(variable)))

  ggplot(top, aes(x = variable, y = abs_smd)) +
    geom_segment(aes(x = variable, xend = variable, y = 0, yend = abs_smd), linewidth = 0.6) +
    geom_point(size = 2) +
    coord_flip() +
    labs(
      title = paste0("Most shifted predictors by |SMD| — ", stratum_name),
      x = NULL,
      y = "|Standardized Mean Difference|"
    ) +
    theme_bw(base_size = 12)
}

# ----------------------------
# Run for each abx stratum (STORE objects + PRINT)
# ----------------------------
run_stratum <- function(df_retro, df_pros, stratum_name) {

  shift_tbl <- build_shift_table(df_retro, df_pros, stratum_name)

  # ---- previously saved to CSV (now commented out) ----
  # write_csv(
  #   shift_tbl %>%
  #     select(stratum, variable, type, retro_summary, pros_summary,
  #            retro_avail, pros_avail, avail_diff, smd, abs_smd, ks_stat),
  #   file.path(OUT_DIR, paste0(prefix, "_shift_table.csv"))
  # )

  # Create plots (objects returned)
  p_density <- plot_density_overlays(df_retro, df_pros, shift_tbl, stratum_name)
  p_avail   <- plot_availability_diff(shift_tbl, stratum_name)
  p_smd     <- plot_lollipop_abs_smd(shift_tbl, stratum_name)

  # ---- previously saved to PNG (now commented out) ----
  # if (!is.null(p_density)) ggsave(file.path(OUT_DIR, paste0(prefix, "_density_overlays.png")), p_density, width = 10, height = 8, dpi = 300)
  # ggsave(file.path(OUT_DIR, paste0(prefix, "_availability_diff.png")), p_avail, width = 10, height = 8, dpi = 300)
  # ggsave(file.path(OUT_DIR, paste0(prefix, "_abs_smd_lollipop.png")), p_smd, width = 10, height = 8, dpi = 300)

  # Print table + figures
  cat("\n============================================================\n")
  cat("SHIFT TABLE — ", stratum_name, "\n", sep = "")
  cat("============================================================\n")
  print(
    shift_tbl %>%
      select(variable, type, retro_summary, pros_summary,
             retro_avail, pros_avail, avail_diff, smd, abs_smd, ks_stat) %>%
      slice_head(n = 50)
  )
  cat("\n(Note: only top 50 rows printed; full table is in the returned object.)\n")

  cat("\n============================================================\n")
  cat("PLOTS — ", stratum_name, "\n", sep = "")
  cat("============================================================\n")

  if (!is.null(p_density)) print(p_density)
  print(p_avail)
  print(p_smd)

  # Return everything so it remains in the environment
  list(
    shift_table = shift_tbl,
    p_density   = p_density,
    p_avail     = p_avail,
    p_smd       = p_smd
  )
}

if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  grDevices::dev.new()  # opens a new graphics device window
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

### Now perform propensity weighting analysis to determine if measurement frequency
# and availability explains difference
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(WeightIt)
  library(cobalt)
  library(patchwork)
  library(scales)
  library(PRROC)
})

get_ps_covariates_case_mix_plus_measurement <- function(df) {

  # Never include these in the PS
  exclude <- c(
    "study_id","mrn","picu_adm_date_time",
    "epoch","A","A_retro",
    "sbi_present",
    "model_score","model_prob","low_risk"
  )

  cand <- setdiff(names(df), exclude)

  # baseline demographics/comorbidity and pre-ICU location (same as before)
  force_in <- intersect(
    c("age","is_female","race","ethnicity","pccc","malignancy_pccc",
      "los_before_icu_days","imv_at_picu_adm"),
    cand
  )

  preicu_cols <- grep("^preicu_", cand, value = TRUE)

  # measurement availability flags
  present_cols <- grep("(_pres|_present)$", cand, value = TRUE)

  # measurement frequency/count columns (adjust patterns to your naming conventions)
  count_cols <- unique(c(
    grep("^n_", cand, value = TRUE),          # e.g., n_hr_rows, n_temp_rows
    grep("_count$", cand, value = TRUE),      # e.g., fio2_count, o2sat_count
    grep("^number_", cand, value = TRUE)      # e.g., number_fio2, number_sat
  ))

  # optional: if you want baseline physiology summaries too (often helpful)
  phys_like <- grep("^(hr_|sbp_|dbp_|rr_|temp_|o2sat_|fio2_|lactate_|hematocrit_|pco2_)", cand, value = TRUE)

  covars <- unique(c(force_in, preicu_cols, phys_like, present_cols, count_cols))

  # drop any that became excluded by overlap of names
  covars <- setdiff(covars, exclude)

  covars
}

library(dplyr)
library(WeightIt)
library(cobalt)

make_weights_pros_to_retro_ebal <- function(df,
                                            ps_trim = NULL,          # usually not needed for ebal
                                            weight_cap = NULL,
                                            covariate_fun,
                                            moments = 1) {

  df2 <- df %>%
    mutate(A_retro = if_else(epoch == "retro", 1L, 0L))

  covars <- covariate_fun(df2)

  # entropy balancing weights Pros (A_retro=0) to match Retro (A_retro=1)
  ps_form <- as.formula(paste("A_retro ~", paste(covars, collapse = " + ")))

  wfit <- weightit(
    formula  = ps_form,
    data     = df2,
    method   = "ebal",
    estimand = "ATT",
    moments  = moments
  )

  df_w <- df2 %>%
    mutate(
      w = wfit$weights
    )

  # optional cap
  if (!is.null(weight_cap)) {
    df_w <- df_w %>% mutate(w = pmin(w, weight_cap))
  }

  # balance diagnostics
  bal <- bal.tab(wfit, un = TRUE, m.threshold = 0.1)
  print(bal)

  list(df_w = df_w, weightit = wfit, balance = bal, covariates = covars)
}


# Weighted AUC via weighted Mann–Whitney
weighted_auc <- function(truth01, score, w) {
  ok <- is.finite(score) & !is.na(truth01) & !is.na(w)
  truth01 <- truth01[ok]; score <- score[ok]; w <- w[ok]

  pos <- truth01 == 1
  neg <- truth01 == 0
  if (!any(pos) || !any(neg)) return(NA_real_)

  s_pos <- score[pos]; w_pos <- w[pos]
  s_neg <- score[neg]; w_neg <- w[neg]

  # pairwise compare: I(s_pos > s_neg) + 0.5 I(tie)
  # do it efficiently by sorting negatives
  ord_neg <- order(s_neg)
  s_neg <- s_neg[ord_neg]; w_neg <- w_neg[ord_neg]
  cum_w_neg <- cumsum(w_neg)
  total_w_neg <- sum(w_neg)

  # for each positive score, find weight of negatives below it and equal to it
  below_w <- function(x) {
    idx <- findInterval(x, s_neg, left.open = TRUE)  # strictly < x
    if (idx <= 0) return(0)
    cum_w_neg[idx]
  }
  equal_w <- function(x) {
    lo <- match(x, s_neg)
    if (is.na(lo)) return(0)
    # sum weights where s_neg == x
    sum(w_neg[s_neg == x])
  }

  contrib <- 0
  for (i in seq_along(s_pos)) {
    b <- below_w(s_pos[i])
    e <- equal_w(s_pos[i])
    contrib <- contrib + w_pos[i] * (b + 0.5 * e)
  }

  contrib / (sum(w_pos) * total_w_neg)
}

roc_points_weighted <- function(df, truth_col="sbi_present", score_col="model_prob", w_col=NULL) {
  d <- df %>%
    transmute(
      truth = as.integer(.data[[truth_col]]),
      score = as.numeric(.data[[score_col]]),
      w = if (is.null(w_col)) 1 else as.numeric(.data[[w_col]])
    ) %>%
    filter(!is.na(truth), is.finite(score), is.finite(w))

  thr <- sort(unique(d$score), decreasing = TRUE)
  if (length(thr) < 2) return(tibble(fpr=numeric(), tpr=numeric()))

  pos_w_total <- sum(d$w[d$truth == 1])
  neg_w_total <- sum(d$w[d$truth == 0])
  if (pos_w_total == 0 || neg_w_total == 0) return(tibble(fpr=numeric(), tpr=numeric()))

  out <- lapply(thr, function(t) {
    pred_pos <- d$score >= t
    tp_w <- sum(d$w[pred_pos & d$truth == 1])
    fp_w <- sum(d$w[pred_pos & d$truth == 0])
    tpr <- tp_w / pos_w_total
    fpr <- fp_w / neg_w_total
    c(fpr=fpr, tpr=tpr)
  })

  out <- do.call(rbind, out)
  tibble(fpr = out[, "fpr"], tpr = out[, "tpr"])
}

auprc_and_points_weighted <- function(df, truth_col="sbi_present", score_col="model_prob", w_col=NULL) {
  d <- df %>%
    transmute(
      truth = as.integer(.data[[truth_col]]),
      score = as.numeric(.data[[score_col]]),
      w = if (is.null(w_col)) 1 else as.numeric(.data[[w_col]])
    ) %>%
    filter(!is.na(truth), is.finite(score), is.finite(w))

  scores_pos <- d$score[d$truth == 1]
  scores_neg <- d$score[d$truth == 0]
  w_pos <- d$w[d$truth == 1]
  w_neg <- d$w[d$truth == 0]

  if (length(scores_pos) == 0 || length(scores_neg) == 0) {
    return(list(auprc = NA_real_, curve = tibble(recall=numeric(), precision=numeric())))
  }

  pr <- PRROC::pr.curve(
    scores.class0 = scores_pos,
    scores.class1 = scores_neg,
    weights.class0 = w_pos,
    weights.class1 = w_neg,
    curve = TRUE
  )

  curve <- as.data.frame(pr$curve)
  # PRROC curve columns: recall (x) and precision (y) are typically V1/V2; be robust:
  recall <- curve[[1]]
  precision <- curve[[2]]

  list(
    auprc = pr$auc.integral,
    curve = tibble(recall = recall, precision = precision)
  )
}

plot_weighted_vs_unweighted <- function(df_w, title_suffix = "") {

  d_pros <- df_w %>% filter(epoch == "pros")  # these are the units being reweighted
  stopifnot(nrow(d_pros) > 0)

  # ROC points + AUROC
  roc_unw <- roc_points_weighted(d_pros, w_col = NULL) %>% mutate(kind="Unweighted")
  roc_w   <- roc_points_weighted(d_pros, w_col = "w")  %>% mutate(kind="Weighted")

  auc_unw <- weighted_auc(d_pros$sbi_present, d_pros$model_prob, rep(1, nrow(d_pros)))
  auc_w   <- weighted_auc(d_pros$sbi_present, d_pros$model_prob, d_pros$w)

  roc_df <- bind_rows(roc_unw, roc_w)

  # PR points + AUPRC
  pr_unw <- auprc_and_points_weighted(d_pros, w_col = NULL)
  pr_w   <- auprc_and_points_weighted(d_pros, w_col = "w")

  pr_df <- bind_rows(
    pr_unw$curve %>% mutate(kind="Unweighted"),
    pr_w$curve   %>% mutate(kind="Weighted")
  )

  auprc_unw <- pr_unw$auprc
  auprc_w   <- pr_w$auprc

  # ---- plots ----
  p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr, linetype = kind)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    geom_line(linewidth = 1) +
    coord_equal(xlim = c(0,1), ylim = c(0,1)) +
    labs(
      title = paste0("Prospective ROC (", title_suffix, ")"),
      subtitle = paste0(
        "AUROC unweighted=", sprintf("%.3f", auc_unw),
        " | weighted=", sprintf("%.3f", auc_w)
      ),
      x = "False positive rate",
      y = "True positive rate",
      linetype = "Pros weighting"
    ) +
    theme_bw()

  p_pr <- ggplot(pr_df, aes(x = recall, y = precision, linetype = kind)) +
    geom_line(linewidth = 1) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    labs(
      title = paste0("Prospective PR (", title_suffix, ")"),
      subtitle = paste0(
        "AUPRC unweighted=", sprintf("%.3f", auprc_unw),
        " | weighted=", sprintf("%.3f", auprc_w)
      ),
      x = "Recall",
      y = "Precision",
      linetype = "Pros weighting"
    ) +
    theme_bw()

  p_roc + p_pr
}

check_weights <- function(df_w) {
  df_w %>%
    group_by(epoch) %>%
    summarise(
      n = n(),
      w_mean = mean(w),
      w_sd   = sd(w),
      w_min  = min(w),
      w_max  = max(w),
      prop_unique = n_distinct(round(w, 8))/n(),
      .groups="drop"
    )
}


# Weight Pros -> Retro using measurement+frequency+baseline covariates
w_noabx_meas_ebal <- make_weights_pros_to_retro_ebal(
  df_noabx,
  covariate_fun = get_ps_covariates_case_mix_plus_measurement,
  weight_cap = 20,
  moments = 1
)
df_noabx_w_meas_ebal <- w_noabx_meas_ebal$df_w
check_weights(df_noabx_w_meas_ebal)

w_yesabx_meas_ebal <- make_weights_pros_to_retro_ebal(
  df_yesabx,
  covariate_fun = get_ps_covariates_case_mix_plus_measurement,
  weight_cap = 20,
  moments = 1
)
df_yesabx_w_meas_ebal <- w_yesabx_meas_ebal$df_w
check_weights(df_yesabx_w_meas_ebal)




