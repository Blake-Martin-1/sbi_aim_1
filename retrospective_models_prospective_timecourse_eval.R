# Evaluate retrospective RF models across all prospective hourly timepoints
#
# What this script does:
# 1) Loads retrospective RF models (rf_model, rf_model_abx) and predictor sets
# 2) Loads full prospective table (pros_all) with repeated timepoints
# 3) Cleans/projections data to model-ready inputs for each RF model type
# 4) Scores every hourly prospective timepoint (rounded hour; one row per encounter-hour)
# 5) Reports AUROC/AUPRC overall and within suspicion-of-infection subset
# 6) Generates an NPV-over-time plot (median + IQR from bootstrap) for each model

suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(DescTools)
  library(yardstick)
  library(pROC)
  library(data.table)
  library(lubridate)
})

# -----------------------------
# User-editable parameters
# -----------------------------
prospective_csv_path <- "/phi/sbi/sbi_blake/pros_all_just_b4_modeling_1_15_26_all_models.csv"
retro_models_script_path <- "/phi/sbi/sbi_blake/aim_1_paper_materials/retro_models_to_2026.R"
output_metrics_csv <- "prospective_timecourse_rf_metrics.csv"
output_predictions_csv <- "prospective_timecourse_rf_predictions.csv"
output_npv_plot <- "prospective_timecourse_npv_by_hour.png"
pros_micro_csv_path <- "/phi/sbi/prospective_data/Prospective/data/micro_export_pros_100125.csv"
pros_cxr_csv_path <- "/phi/sbi/prospective_data/Prospective/data/hosp_chest_xray_export_pros_100125.csv"

# Thresholds used for NPV calculations (can be adjusted)
threshold_no_abx <- 0.05
threshold_yes_abx <- 0.074

# Bootstrap controls for NPV-by-hour distribution summaries
n_boot <- 500
set.seed(2026)

# -----------------------------
# Utility functions
# -----------------------------

safe_factor01_to_yesno <- function(x) {
  if (is.factor(x) || is.character(x)) {
    x_chr <- as.character(x)
    out <- ifelse(x_chr %in% c("yes", "1", "TRUE", "true"), 1L, 0L)
    return(as.integer(out))
  }
  as.integer(ifelse(is.na(x), NA, x > 0))
}

truth_factor <- function(y01) {
  factor(ifelse(y01 == 1, "yes", "no"), levels = c("yes", "no"))
}

auroc_safe <- function(df, truth_col = "truth", prob_col = "p_yes") {
  dat <- df %>% dplyr::filter(!is.na(.data[[truth_col]]), !is.na(.data[[prob_col]]))
  if (nrow(dat) < 10 || dplyr::n_distinct(dat[[truth_col]]) < 2) return(NA_real_)

  # ensure levels match and event is "yes"
  dat[[truth_col]] <- factor(dat[[truth_col]], levels = c("yes", "no"))

  # IMPORTANT: pass estimate column positionally (no estimate =)
  yardstick::roc_auc(dat, !!rlang::sym(truth_col), !!rlang::sym(prob_col), event_level = "first") %>%
    dplyr::pull(.estimate)
}

auprc_safe <- function(df, truth_col = "truth", prob_col = "p_yes") {
  dat <- df %>% dplyr::filter(!is.na(.data[[truth_col]]), !is.na(.data[[prob_col]]))
  if (nrow(dat) < 10 || dplyr::n_distinct(dat[[truth_col]]) < 2) return(NA_real_)

  dat[[truth_col]] <- factor(dat[[truth_col]], levels = c("yes", "no"))

  yardstick::pr_auc(dat, !!rlang::sym(truth_col), !!rlang::sym(prob_col), event_level = "first") %>%
    dplyr::pull(.estimate)
}

npv_at_threshold_df <- function(df, threshold) {
  dat <- df %>% filter(!is.na(truth), !is.na(p_yes))
  if (nrow(dat) == 0) return(NA_real_)

  pred_neg <- dat$p_yes < threshold
  tn <- sum(pred_neg & dat$truth == "no")
  fn <- sum(pred_neg & dat$truth == "yes")

  if ((tn + fn) == 0) return(NA_real_)
  tn / (tn + fn)
}

bootstrap_hourly_npv <- function(df, threshold, n_boot = 200) {
  df <- df %>% filter(!is.na(hour_rounded), !is.na(truth), !is.na(p_yes))
  hrs <- sort(unique(df$hour_rounded))

  purrr::map_dfr(hrs, function(hh) {
    d_h <- df %>% filter(hour_rounded == hh)
    if (nrow(d_h) < 20) {
      return(tibble(
        hour_rounded = hh,
        n_rows = nrow(d_h),
        npv_median = NA_real_,
        npv_q25 = NA_real_,
        npv_q75 = NA_real_
      ))
    }

    boot_vals <- replicate(n_boot, {
      idx <- sample.int(nrow(d_h), size = nrow(d_h), replace = TRUE)
      npv_at_threshold_df(d_h[idx, ], threshold = threshold)
    })

    tibble(
      hour_rounded = hh,
      n_rows = nrow(d_h),
      npv_median = median(boot_vals, na.rm = TRUE),
      npv_q25 = quantile(boot_vals, probs = 0.25, na.rm = TRUE),
      npv_q75 = quantile(boot_vals, probs = 0.75, na.rm = TRUE)
    )
  })
}

derive_suspicion_from_micro_cxr <- function(pros_df, micro_csv_path, cxr_csv_path) {
  encounter_dt <- as.data.table(
    pros_df %>%
      dplyr::select(study_id, pat_mrn_id, picu_adm_date_time) %>% rename(mrn = pat_mrn_id) %>%
      distinct()
  )

  encounter_dt[, `:=`(
    mrn = as.character(mrn),
    picu_adm_date_time = force_tz(picu_adm_date_time, tz = "America/Denver"),
    suspected_infection = 0L,
    window_start = picu_adm_date_time - as.difftime(24, units = "hours"),
    window_end = picu_adm_date_time + as.difftime(24, units = "hours"),
    row_id = .I
  )]

  micro_dt <- as.data.table(
    readr::read_csv(micro_csv_path, show_col_types = FALSE) %>%
      mutate(
        pat_mrn_id = as.character(pat_mrn_id),
        specimn_taken_time = force_tz(as.POSIXct(specimn_taken_time), tz = "America/Denver")
      ) %>%
      rename(mrn = pat_mrn_id)
  )
  micro_dt <- micro_dt[!is.na(mrn) & !is.na(specimn_taken_time)]
  setindex(micro_dt, mrn, specimn_taken_time)

  micro_hits <- micro_dt[
    encounter_dt,
    on = .(mrn, specimn_taken_time >= window_start, specimn_taken_time <= window_end),
    nomatch = 0L,
    allow.cartesian = TRUE
  ][, .(suspected_infection = 1L), by = row_id]

  encounter_dt <- merge(encounter_dt, micro_hits, by = "row_id", all.x = TRUE, suffixes = c("", "_micro"))
  encounter_dt[!is.na(suspected_infection_micro), suspected_infection := 1L]
  encounter_dt[, suspected_infection_micro := NULL]

  cxr_dt <- as.data.table(
    readr::read_csv(cxr_csv_path, show_col_types = FALSE) %>%
      mutate(
        mrn = as.character(pat_mrn_id),
        chest_x_ray_order = force_tz(
          as.POSIXct(chest_x_ray_order, format = "%Y-%m-%d %H:%M:%OS"),
          tz = "America/Denver"
        )
      )
  )
  cxr_dt <- cxr_dt[chest_xray_yn == 1 & !is.na(chest_x_ray_order) & !is.na(mrn)]

  cxr_hit_ids <- unique(
    cxr_dt[
      encounter_dt,
      on = .(mrn, chest_x_ray_order >= window_start, chest_x_ray_order <= window_end),
      nomatch = 0L,
      .(study_id)
    ]$study_id
  )

  encounter_dt[study_id %in% cxr_hit_ids, suspected_infection := 1L]

  suspicion_map <- encounter_dt[, .(study_id, suspected_infection)]

  pros_df %>%
    left_join(as_tibble(suspicion_map), by = "study_id") %>%
    mutate(
      suspected_infection = replace_na(suspected_infection, 0L),
      suspicion_flag = suspected_infection
    )
}

# Attempt to identify a suspicion-of-infection column in prospective data
derive_suspicion_flag <- function(df) {
  candidates <- c(
    "any_micro_pres", "suspicion_of_infection", "suspicion_infection",
    "suspected_infection", "micro_sbi_1_0", "bcx_sent"
  )

  col_hit <- candidates[candidates %in% names(df)][1]
  if (is.na(col_hit)) {
    warning("No explicit suspicion column found; using all rows for both 'all' and 'suspicion' metrics.")
    return(df %>% mutate(suspicion_flag = 1L))
  }

  message("Using suspicion flag column: ", col_hit)
  df %>% mutate(suspicion_flag = safe_factor01_to_yesno(.data[[col_hit]]))
}

prepare_model_inputs <- function(df, predictors, is_abx_model = FALSE) {
  out <- df

  if (is_abx_model) {
    if ("hematocrit_mean" %in% names(out) && !"hematocrit_blood" %in% names(out)) {
      out <- out %>% rename(hematocrit_blood = hematocrit_mean)
    }
    if ("fio2_last" %in% names(out) && !"last_fio2" %in% names(out)) {
      out <- out %>% rename(last_fio2 = fio2_last)
    }
    if ("dbp_max" %in% names(out) && !"max_dbp" %in% names(out)) {
      out <- out %>% rename(max_dbp = dbp_max)
    }
  }

  missing_pred <- setdiff(predictors, names(out))
  if (length(missing_pred) > 0) {
    stop("Missing required predictor columns: ", paste(missing_pred, collapse = ", "))
  }

  out
}

score_model_dataset <- function(df_model, model_obj, predictors, threshold, model_label) {
  pred_prob <- predict(model_obj, newdata = df_model[, predictors, drop = FALSE], type = "prob")

  scored <- df_model %>%
    mutate(
      model_label = model_label,
      p_yes = pred_prob[, "yes"],
      truth = factor(ifelse(sbi_present == 1, "yes", "no"), levels = c("yes", "no")),
      threshold = threshold
    )

  all_metrics <- tibble(
    model_label = model_label,
    subset = "all_patients",
    n = nrow(scored),
    prevalence = mean(scored$truth == "yes", na.rm = TRUE),
    auroc = auroc_safe(scored),
    auprc = auprc_safe(scored),
    npv = npv_at_threshold_df(scored, threshold = threshold)
  )

  scored_susp <- scored %>% filter(suspicion_flag == 1)
  susp_metrics <- tibble(
    model_label = model_label,
    subset = "suspicion_of_infection",
    n = nrow(scored_susp),
    prevalence = mean(scored_susp$truth == "yes", na.rm = TRUE),
    auroc = auroc_safe(scored_susp),
    auprc = auprc_safe(scored_susp),
    npv = npv_at_threshold_df(scored_susp, threshold = threshold)
  )

  npv_hour_all <- bootstrap_hourly_npv(scored, threshold = threshold, n_boot = n_boot) %>%
    mutate(model_label = model_label, subset = "all_patients")

  npv_hour_susp <- bootstrap_hourly_npv(scored_susp, threshold = threshold, n_boot = n_boot) %>%
    mutate(model_label = model_label, subset = "suspicion_of_infection")

  list(
    scored = scored,
    metrics = bind_rows(all_metrics, susp_metrics),
    npv_hour = bind_rows(npv_hour_all, npv_hour_susp)
  )
}

# -----------------------------
# Main workflow
# -----------------------------

# Load retrospective model objects and predictor names:
# Expected objects from this source call:
#   rf_model, rf_model_abx, predictors, predictors_abx
source(retro_models_script_path)

# Load prospective full time-series data
pros_all <- readr::read_csv(prospective_csv_path, show_col_types = FALSE)

# Basic cleaning and hour binning
pros_all_clean <- pros_all %>%
  mutate(
    sbi_present = safe_factor01_to_yesno(sbi_present),
    hour_rounded = as.integer(round(hours_since_picu_adm, digits = 0))
  ) %>%
  filter(!is.na(sbi_present), !is.na(hour_rounded), hour_rounded >= 0, hour_rounded <= 24)

# Keep one row per encounter-hour-model_type, resolving close-time duplicates
# using smallest |hours - rounded_hour| and then latest score_time.
pros_all_clean <- pros_all_clean %>%
  mutate(hour_delta = abs(hours_since_picu_adm - hour_rounded)) %>%
  group_by(study_id, pat_enc_csn_id, model_type, hour_rounded) %>%
  arrange(hour_delta, desc(score_time), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(-hour_delta)

# Derive suspicion-of-infection for every encounter in pros_all_clean using
# micro orders and chest X-ray orders within +/- 24 hours of PICU admission.
pros_all_clean <- derive_suspicion_from_micro_cxr(
  pros_df = pros_all_clean,
  micro_csv_path = pros_micro_csv_path,
  cxr_csv_path = pros_cxr_csv_path
)

# Split by prospective model channel
pros_no_abx <- pros_all_clean %>% filter(model_type == "RF_no_abx")
pros_yes_abx <- pros_all_clean %>% filter(model_type == "RF_yes_abx")

# Align predictor names/shape
# For the abx model, this mirrors prior rename logic used in the codebase.
predictors_abx_use <- predictors_abx
predictors_abx_use[predictors_abx_use == "hematocrit_mean"] <- "hematocrit_blood"

pros_no_abx <- prepare_model_inputs(pros_no_abx, predictors = predictors, is_abx_model = FALSE)
pros_yes_abx <- prepare_model_inputs(pros_yes_abx, predictors = predictors_abx_use, is_abx_model = TRUE)

# Ensure factor fields expected by the trained models have proper type if present
for (nm in c("pccc", "malignancy_pccc", "race")) {
  if (nm %in% names(pros_no_abx)) pros_no_abx[[nm]] <- as.factor(pros_no_abx[[nm]])
  if (nm %in% names(pros_yes_abx)) pros_yes_abx[[nm]] <- as.factor(pros_yes_abx[[nm]])
}

# Score both models across all timepoints
res_no_abx <- score_model_dataset(
  df_model = pros_no_abx,
  model_obj = rf_model,
  predictors = predictors,
  threshold = threshold_no_abx,
  model_label = "RF_no_abx"
)

res_yes_abx <- score_model_dataset(
  df_model = pros_yes_abx,
  model_obj = rf_model_abx,
  predictors = predictors_abx_use,
  threshold = threshold_yes_abx,
  model_label = "RF_yes_abx"
)

metrics_out <- bind_rows(res_no_abx$metrics, res_yes_abx$metrics)
preds_out <- bind_rows(res_no_abx$scored, res_yes_abx$scored)
npv_hour_out <- bind_rows(res_no_abx$npv_hour, res_yes_abx$npv_hour)

# Save outputs
readr::write_csv(metrics_out, output_metrics_csv)
readr::write_csv(preds_out, output_predictions_csv)

# Plot NPV median/IQR over time
p_npv <- npv_hour_out %>%
  ggplot(aes(x = hour_rounded, y = npv_median, color = model_label, fill = model_label)) +
  geom_ribbon(aes(ymin = npv_q25, ymax = npv_q75), alpha = 0.20, color = NA) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~subset) +
  scale_x_continuous(breaks = 0:24) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Prospective NPV Over Time for Retrospective RF Models",
    subtitle = "Median and IQR from bootstrap resampling by rounded hour, by subset",
    x = "Hours since PICU admission (rounded)",
    y = "Negative Predictive Value (NPV)",
    color = "Model",
    fill = "Model"
  ) +
  theme_bw(base_size = 12)

print(p_npv)

ggsave(
  output_npv_plot,
  p_npv,
  width = 11,
  height = 6,
  dpi = 300
)

# Console summary
message("\n=== Overall performance summary ===")
print(metrics_out)
message("\nWrote: ", output_metrics_csv)
message("Wrote: ", output_predictions_csv)
message("Wrote: ", output_npv_plot)
