# ---------------------------------------------------------------------------- #
# file:    setup_aim_1.R
#
# authors: Blake Martin   <blake.martin@childrenscolorado.org>
#          Tell Bennett   <tell.bennett@ucdenver.edu>
#
# This file is sourced by other scripts
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# Options and library loads
# ---------------------------------------------------------------------------- #
options(stringsAsFactors = FALSE)
options(dplyr.max_print = 10)
options(dplyr.min_print =  5)

######

options("digits" = 15)
options("scipen" = 10000)

######
set.seed(32312)

#######
library(tidyverse)
library(tidyr)
library(readr)
library(dplyr)
library(readxl)
library(purrr)
library(chron)
library(stringr)
library(stringi)
library(stats)
library(tibble)
library(robustbase) #For colMedians function
library(DescTools)
library(rriskDistributions)
library(tigerstats)
library(fst)
library(rpart)
library(ranger)
library(caTools)
library(caret)
library(pROC)
library(pccc)
library(stratifyR)
library(glmnet)
library(regclass)
library(devtools)
library(recipes)
library(ggpubr)
library(e1071)
library(mltools)
library(data.table)
library(kernlab)
library(MLmetrics)
library(OptimalCutpoints)
library(qwraps2)


# Set User path to repo and data
if (Sys.info()["user"] == "tdbennett") {
  options(user_path_to_repo = '~/sbi_prediction/',
          user_data_path = '/media/sf_C_DRIVE/sbi_data/')
} else {
  options(user_path_to_repo = "/home/martinbl/sbi_blake/",
          user_data_path = "/phi/sbi/sbi_data/ten_yr_data/",
          user_norms_path = "/phi/sbi/sbi_data/")
}

# Guard/fallback for optional user_norms_path configuration
# Some user setups may not define user_norms_path explicitly; in that case
# fall back to user_data_path so norm-based path construction still works.
user_norms_base_path <- getOption("user_norms_path")
if (is.null(user_norms_base_path) || !nzchar(user_norms_base_path)) {
  user_norms_base_path <- getOption("user_data_path")
  warning("user_norms_path not set; falling back to user_data_path for norm files.")
}

# File Paths for all SBI Data files
abx_file_path <- paste0(getOption("user_data_path"), "anti_infectives.csv")
adm_dx_file_path <- paste0(getOption("user_data_path"), "basic_encounter_information.csv")
bp_norms_lo_boys <- paste0(user_norms_base_path, "Lo_bp_norms_boys.csv")
bp_norms_lo_girls <- paste0(user_norms_base_path, "Lo_bp_norms_girls.csv")
bp_norms_cdc_boys <- paste0(user_norms_base_path, "CDC_bp_norms_boys.csv")
bp_norms_cdc_girls <- paste0(user_norms_base_path, "CDC_bp_norms_girls.csv")
demographics_file_path <- paste0(getOption("user_data_path"), "demographics_20220128.csv")
extra_weights_file_path <- paste0(getOption("user_data_path"), "heights_and_weights.csv")

gcs_file_path <- paste0(getOption("user_data_path"), "gcs_scores.csv")

Gemelli_bp_under_one_boys <- paste0(user_norms_base_path, "Gemelli_bp_under_one_boys.csv")
Gemelli_bp_under_one_girls <- paste0(user_norms_base_path, "Gemelli_bp_under_one_girls.csv")

hr_norms <- paste0(user_norms_base_path, "hr_percentiles.CSV")

info_admit_epic <- paste0(getOption("user_data_path"), "admit_info.xlsx")
options(prospective_model_data_file_path = getOption(
  "prospective_model_data_file_path",
  "/phi/sbi/sbi_blake/pros_all_just_b4_modeling_1_15_26_all_models.csv"
))

# Need updated lab info
labs_1_file_path <- paste0(getOption("user_data_path"), "labs.csv")
labs_2_file_path <- paste0(getOption("user_data_path"), "labs_2_fixed_dates.csv")

length_cdc_infant_boys <- paste0(user_norms_base_path, "CDC_length_infant_boys.csv")
length_cdc_infant_girls <- paste0(user_norms_base_path, "CDC_length_infant_girls.csv")
length_cdc_boys <- paste0(user_norms_base_path, "CDC_length_boys.csv")
length_cdc_girls <- paste0(user_norms_base_path, "CDC_length_girls.csv")

# micro_comments_file_path <- paste0(getOption("user_data_path"), "micro_with_comments_fixed_dates.csv")
micro_file_path <- paste0(getOption("user_data_path"), "micro_20220128.csv")

PRISM_PIM_file_path_1 <- paste0(user_norms_base_path, "PRISM_PIM.csv")
PRISM_PIM_file_path_2 <- paste0(getOption("user_data_path"), "2018_to_2020_vps.xlsx")
prob_list_file_path <- paste0(getOption("user_data_path"), "problem_list.csv")
pt_origin_file_path <- paste0(user_norms_base_path, "pt_origin.csv")
pt_origin_file_path_10_yr <- paste0(getOption("user_data_path"), "vps_ten_yr.xlsx")

pt_sex_file_path <- paste0(getOption("user_data_path"), "pt_sex.csv")
race_file_path <- paste0(getOption("user_data_path"), "race_death.csv")
rr_norms <- paste0(user_norms_base_path, "rr_percentiles.CSV")
vitals_full_file_path <- paste0(getOption("user_data_path"), "vitals.csv")
vps_file_path_10_yr <- paste0(getOption("user_data_path"), "vps_ten_yr.xlsx")
vps_file_path <- paste0(user_norms_base_path, "vps_full_admit_list_fixed_dates.xlsx")

weight_cdc_infant_boys <- paste0(user_norms_base_path, "CDC_weight_infant_boys.csv")
weight_cdc_infant_girls <- paste0(user_norms_base_path, "CDC_weight_infant_girls.csv")
weight_cdc_boys <- paste0(user_norms_base_path, "CDC_weight_boys.csv")
weight_cdc_girls <- paste0(user_norms_base_path, "CDC_weight_girls.csv")


# Script file paths
abx_script <- paste0(getOption("user_path_to_repo"), "abx_10yr.R")
adm_dx_script <- paste0(getOption("user_path_to_repo"), "adm_dx.R")

admit_info_from_epic <- paste0(getOption("user_path_to_repo"), "admit_epic.R")

biomarker_models_script <- paste0(getOption("user_path_to_repo"), "biomarker_models_10yr.R")
box_plots_script <- paste0(getOption("user_path_to_repo"), "box_plots_10yr.R")
cx_neg_sepsis_script <- paste0(getOption("user_path_to_repo"), "cx_neg_sepsis_10yr.R")
death_case_ids_script <-  paste0(getOption("user_path_to_repo"), "death_case_ids_10yr.R")
dummy_script <- paste0(getOption("user_path_to_repo"), "dummy_variable_creation.R")
load_gcs_script <- paste0(getOption("user_path_to_repo"), "get_gcs.R")
impute_labs_script <- paste0(getOption("user_path_to_repo"), "impute_labs_10yr.R")
lab_tidying_script <- paste0(getOption("user_path_to_repo"), "lab_tidying_10yr.R")
line_and_vent_entry <- paste0(getOption("user_path_to_repo"), "lines_and_vent_data_entry.R")
load_lines_demo_table <- paste0(getOption("user_path_to_repo"), "load_vps_lines_and_demo_table_10yr.R")
make_labs_factors_script <- paste0(getOption("user_path_to_repo"), "make_labs_factors_10yr.R")
micro_script <- paste0(getOption("user_path_to_repo"), "micro_tidying_10yr.R")
mod_p_val_script <- paste0(getOption("user_path_to_repo"), "model_p_value_10yr.R")
models_no_pna <- paste0(getOption("user_path_to_repo"), "models_no_pna.R")

models_script <- paste0(getOption("user_path_to_repo"), "models_10yr.R")
models_script_revision <- paste0(getOption("user_path_to_repo"), "models_revised.R")

pen_regr_script <- paste0(getOption("user_path_to_repo"), "pen_regr_script_10yr.R")
pen_regr_script_revision <- paste0(getOption("user_path_to_repo"), "pen_regr_script_revised_10yr.R")

PRISM_script <- paste0(getOption("user_path_to_repo"), "PRISM_PIM_10yr.R")
prob_r_script <- paste0(getOption("user_path_to_repo"), "prob_list_10yr.R")
pt_origin_script <- paste0(getOption("user_path_to_repo"), "pt_origin.R")
rf_from_crp_script <- paste0(getOption("user_path_to_repo"), "rf_from_crp.R")
roc_graph_script <- paste0(getOption("user_path_to_repo"), "roc_graph_10yr.R")
sbi_subgroups_script <- paste0(getOption("user_path_to_repo"), "sbi_subgroups_10yr.R")
sirs_script <- paste0(getOption("user_path_to_repo"), "sirs_10yr.R")

svm_script <- paste0(getOption("user_path_to_repo"), "svm_10yr.R")
svm_script_revision <- paste0(getOption("user_path_to_repo"), "svm_revised.R")

vitals_tidying_script <- paste0(getOption("user_path_to_repo"), "vitals_tidying_10yr.R")
vps_pna_eval_script <- paste0(getOption("user_path_to_repo"), "vps_pna_eval.R")
vps_pnas_script <- paste0(getOption("user_path_to_repo"), "vps_pnas_10yr.R")
vs_normalize_script <- paste0(getOption("user_path_to_repo"), "vs_normalize_10yr.R")