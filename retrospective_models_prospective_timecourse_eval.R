# retrospective_models_prospective_timecourse_eval.R
# Script to prepare all of the prospective dataset for input into the two random forest models. Specifically we will need to retain all
# timepoints, not just the point 2 hours after PICU admission. No need to recreate all of the 'suspicion of infection' code since we
# are not going to look at that subgroup at these other timepoints.

source(file = "retro_models_to_2026.R")

pros_all <- read_csv(file = prospective_model_data_file_path)

#------------------------------------------------------------
# Build hourly cohorts using the same timing logic as aim_1_main.R
#
# In aim_1_main.R, the "2-hour" prediction is:
#   hours_since_picu_adm >= 2 and < 4
# followed by the earliest available row per study_id.
#
# Here, for each target hour h, we use:
#   hours_since_picu_adm >= h and < h + 2
# followed by the earliest available row per study_id and hour.
#
# Therefore, picu_hour == 2 in this script now matches the
# aim_1_main.R 2-hour analytic timepoint.
#------------------------------------------------------------

make_aim1_timing_timecourse <- function(df, abx_value, target_hours = 1:24) {
  purrr::map_dfr(
    target_hours,
    function(hr) {
      df %>%
        dplyr::filter(
          abx_exp == abx_value,
          hours_since_picu_adm >= hr,
          hours_since_picu_adm < hr + 2
        ) %>%
        dplyr::arrange(study_id, hours_since_picu_adm, score_time) %>%
        dplyr::group_by(study_id) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(picu_hour = hr)
    }
  )
}

# Filter to RF model outputs
pros_rf_ue <- pros_all %>%
  dplyr::filter(model_type == "RF_no_abx")

pros_rf_ae <- pros_all %>%
  dplyr::filter(model_type == "RF_yes_abx")

# Build timecourse rows using aim_1_main.R-style timing
pros_24hr_rf_ue <- make_aim1_timing_timecourse(
  df = pros_rf_ue,
  abx_value = 0,
  target_hours = 1:24
)

pros_24hr_rf_ae <- make_aim1_timing_timecourse(
  df = pros_rf_ae,
  abx_value = 1,
  target_hours = 1:24
)


#------------------------------------------------------------
# Select required columns for no-abx model
#------------------------------------------------------------

thresh_rf_no_abx <- 0.05

required_no_abx <- c(
  "study_id", "picu_hour", "hours_since_picu_adm", "score_time",
  "picu_adm_date_time", "hsp_account_id",
  "sbi_present", "model_score", predictors
)

missing_from_source_no_abx <- setdiff(required_no_abx, colnames(pros_24hr_rf_ue))
missing_from_source_no_abx

pros_rf_same_24 <- pros_24hr_rf_ue %>%
  dplyr::select(dplyr::all_of(required_no_abx)) %>%
  dplyr::relocate(study_id, picu_adm_date_time, picu_hour, hours_since_picu_adm, score_time, sbi_present, model_score) %>%
  dplyr::mutate(epoch = "prospective")


#------------------------------------------------------------
# Select required columns for yes-abx model
#------------------------------------------------------------

thresh_rf_yes_abx <- 0.074

pros_24hr_rf_ae_renamed <- pros_24hr_rf_ae %>%
  dplyr::rename(
    last_fio2 = fio2_last,
    max_dbp = dbp_max
  )

predictors_abx[predictors_abx == "hematocrit_mean"] <- "hematocrit_blood"

pros_24hr_rf_ae_renamed <- pros_24hr_rf_ae_renamed %>%
  dplyr::rename(
    hematocrit_blood = hematocrit_mean
  )

required_yes_abx <- c(
  "study_id", "picu_hour", "hours_since_picu_adm", "score_time",
  "picu_adm_date_time", "hsp_account_id",
  "sbi_present", "model_score", predictors_abx
)

missing_from_source_yes_abx <- setdiff(required_yes_abx, colnames(pros_24hr_rf_ae_renamed))
missing_from_source_yes_abx

pros_rf_same_abx_24 <- pros_24hr_rf_ae_renamed %>%
  dplyr::select(dplyr::all_of(required_yes_abx)) %>%
  dplyr::relocate(study_id, picu_adm_date_time, picu_hour, hours_since_picu_adm, score_time, sbi_present, model_score) %>%
  dplyr::mutate(epoch = "prospective")


# Names used later
pros_24_no <- pros_rf_same_24
pros_24_yes <- pros_rf_same_abx_24

#------------------------------------------------------------
# Match hourly cohorts to the exact prospective encounter cohorts
# used in aim_1_main.R
#------------------------------------------------------------

required_aim1_objects <- c("pros_no_abx_1st_infxn", "pros_yes_temp")

missing_aim1_objects <- required_aim1_objects[
  !vapply(required_aim1_objects, exists, logical(1))
]

if (length(missing_aim1_objects) > 0) {
  stop(
    "The following aim_1_main.R cohort objects are missing: ",
    paste(missing_aim1_objects, collapse = ", "),
    ". Run or source the code that creates these objects before this script."
  )
}

quad_ids_no_abx <- pros_no_abx_1st_infxn %>%
  dplyr::distinct(study_id)

quad_ids_yes_abx <- pros_yes_temp %>%
  dplyr::distinct(study_id)

pros_24_no_matched <- pros_24_no %>%
  dplyr::semi_join(quad_ids_no_abx, by = "study_id")

pros_24_yes_matched <- pros_24_yes %>%
  dplyr::semi_join(quad_ids_yes_abx, by = "study_id")

# Optional QC: show how many encounters/rows are included
matched_qc <- dplyr::bind_rows(
  pros_24_no_matched %>%
    dplyr::summarise(
      group = "No antibiotic exposure before PICU",
      n_rows = dplyr::n(),
      n_encounters = dplyr::n_distinct(study_id),
      n_hour2_rows = sum(picu_hour == 2, na.rm = TRUE),
      n_hour2_encounters = dplyr::n_distinct(study_id[picu_hour == 2])
    ),
  pros_24_yes_matched %>%
    dplyr::summarise(
      group = "Antibiotic exposure before PICU",
      n_rows = dplyr::n(),
      n_encounters = dplyr::n_distinct(study_id),
      n_hour2_rows = sum(picu_hour == 2, na.rm = TRUE),
      n_hour2_encounters = dplyr::n_distinct(study_id[picu_hour == 2])
    )
)

matched_qc

#------------------------------------------------------------
# Recalculate predictions on matched hourly cohorts
#------------------------------------------------------------

validate_predictor_schema(pros_24_no_matched, predictors, "pros_24_no_matched")

rf_pred_prob_pros_24_matched <- predict(
  rf_model,
  pros_24_no_matched[, dplyr::all_of(predictors)],
  type = "prob"
)

validate_predictor_schema(pros_24_yes_matched, predictors_abx, "pros_24_yes_matched")

rf_pred_prob_pros_24_abx_matched <- predict(
  rf_model_abx,
  pros_24_yes_matched[, dplyr::all_of(predictors_abx)],
  type = "prob"
)

pred_df_pros_no_abx_24_matched <- pros_24_no_matched %>%
  dplyr::transmute(
    rowIndex = dplyr::row_number(),
    study_id = study_id,
    picu_hour = picu_hour,
    hours_since_picu_adm = hours_since_picu_adm,
    score_time = score_time,
    sbi_present = sbi_present,
    pred_prob_yes = rf_pred_prob_pros_24_matched[, "yes"]
  )

pred_df_pros_yes_abx_24_matched <- pros_24_yes_matched %>%
  dplyr::transmute(
    rowIndex = dplyr::row_number(),
    study_id = study_id,
    picu_hour = picu_hour,
    hours_since_picu_adm = hours_since_picu_adm,
    score_time = score_time,
    sbi_present = sbi_present,
    pred_prob_yes = rf_pred_prob_pros_24_abx_matched[, "yes"]
  )

#------------------------------------------------------------
# QC: confirm picu_hour == 2 matches aim_1_main.R timing/cohort
#------------------------------------------------------------

qc_hour2_no <- pred_df_pros_no_abx_24_matched %>%
  dplyr::filter(picu_hour == 2) %>%
  dplyr::summarise(
    group = "No antibiotic exposure before PICU",
    n_rows = dplyr::n(),
    n_encounters = dplyr::n_distinct(study_id),
    min_hours_since_picu_adm = min(hours_since_picu_adm, na.rm = TRUE),
    max_hours_since_picu_adm = max(hours_since_picu_adm, na.rm = TRUE)
  )

qc_hour2_yes <- pred_df_pros_yes_abx_24_matched %>%
  dplyr::filter(picu_hour == 2) %>%
  dplyr::summarise(
    group = "Antibiotic exposure before PICU",
    n_rows = dplyr::n(),
    n_encounters = dplyr::n_distinct(study_id),
    min_hours_since_picu_adm = min(hours_since_picu_adm, na.rm = TRUE),
    max_hours_since_picu_adm = max(hours_since_picu_adm, na.rm = TRUE)
  )

dplyr::bind_rows(qc_hour2_no, qc_hour2_yes)

#### Now make code to look at NPV over time ###
# Helper function

make_npv_by_hour <- function(df, threshold, group_label) {
  out <- df %>%
    mutate(
      pred_negative = pred_prob_yes <= threshold
    ) %>%
    group_by(picu_hour) %>%
    summarise(
      n_total = n(),
      n_pred_negative = sum(pred_negative, na.rm = TRUE),
      n_true_negative = sum(pred_negative & sbi_present == 0, na.rm = TRUE),
      npv = ifelse(n_pred_negative > 0, n_true_negative / n_pred_negative, NA_real_),
      .groups = "drop"
    )

  # Add Wilson CIs only where denominator > 0
  ci_df <- binom::binom.confint(
    x = out$n_true_negative,
    n = out$n_pred_negative,
    methods = "wilson"
  ) %>%
    dplyr::select(lower, upper)

  out %>%
    bind_cols(ci_df) %>%
    rename(
      ci_low = lower,
      ci_high = upper
    ) %>%
    mutate(
      group = group_label,
      threshold = threshold
    )
}

# ----------------------------
# Create summary data for both groups
# ----------------------------
npv_no_abx <- make_npv_by_hour(
  df = pred_df_pros_no_abx_24_matched,
  threshold = thresh_rf_no_abx,
  group_label = "No antibiotic exposure before PICU"
)

npv_yes_abx <- make_npv_by_hour(
  df = pred_df_pros_yes_abx_24_matched,
  threshold = thresh_rf_yes_abx,
  group_label = "Antibiotic exposure before PICU"
)

npv_plot_df <- bind_rows(npv_no_abx, npv_yes_abx) %>%
  mutate(
    facet_label = paste0(group, "\nThreshold for predicted SBI-negative: ≤ ", threshold)
  )

# ----------------------------
# Plot
# ----------------------------

library(dplyr)
library(ggplot2)
library(scales)

npv_plot_df <- npv_plot_df %>%
  mutate(
    tn_label = ifelse(
      n_pred_negative > 0,
      paste0(n_true_negative, "/", n_pred_negative),
      NA_character_
    )
  )

line_df <- npv_plot_df %>%
  filter(!is.na(npv)) %>%
  arrange(facet_label, picu_hour)

p_npv_facet <- ggplot(npv_plot_df, aes(x = picu_hour, y = npv)) +
  geom_ribbon(
    data = line_df,
    aes(x = picu_hour, ymin = ci_low, ymax = ci_high, group = facet_label),
    inherit.aes = FALSE,
    fill = "lightblue",
    alpha = 0.35
  ) +
  geom_line(
    data = line_df,
    aes(group = facet_label),
    linewidth = 1.2,
    color = "darkblue"
  ) +
  geom_point(
    data = line_df,
    size = 2.8,
    color = "darkblue"
  ) +
  geom_text(
    data = line_df,
    aes(label = tn_label),
    na.rm = TRUE,
    vjust = -0.9,
    size = 3.5,
    color = "darkblue"
  ) +
  facet_wrap(~ facet_label, ncol = 1) +
  scale_x_continuous(breaks = 1:24) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(ylim = c(0, 1.08)) +
  labs(
    title = "Negative predictive value over PICU time",
    subtitle = "Point labels show true negatives / total predicted negatives",
    x = "PICU hour",
    y = "Negative predictive value"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 13),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    strip.text = element_text(face = "bold", size = 14)
  )


p_npv_facet


### Repeat but for AUROC and AUPRC by hour ###
library(dplyr)
library(ggplot2)
library(scales)
library(pROC)
library(yardstick)
library(purrr)

#------------------------------------------------------------
# Metric helper functions matching aim_1_main.R
#------------------------------------------------------------

to_truth01 <- function(x) {
  if (is.logical(x)) {
    return(as.integer(x))
  }

  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.character(x)) {
    x_clean <- tolower(trimws(x))
    return(dplyr::case_when(
      x_clean %in% c("1", "yes", "y", "true", "t", "sbi+", "sbi_positive", "positive") ~ 1L,
      x_clean %in% c("0", "no", "n", "false", "f", "sbi-", "sbi_negative", "negative") ~ 0L,
      TRUE ~ NA_integer_
    ))
  }

  x_num <- suppressWarnings(as.numeric(x))
  dplyr::if_else(is.na(x_num), NA_integer_, as.integer(x_num))
}

to_prob01 <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  x_num[!is.finite(x_num)] <- NA_real_
  x_num
}

prep_metric_df <- function(df) {
  df %>%
    dplyr::transmute(
      sbi_present_num = to_truth01(sbi_present),
      score_yes = to_prob01(pred_prob_yes)
    ) %>%
    dplyr::filter(
      !is.na(sbi_present_num),
      !is.na(score_yes),
      is.finite(score_yes)
    ) %>%
    dplyr::mutate(
      # 1 = SBI-negative case/event
      # 0 = SBI-positive control/non-event
      case_num = dplyr::if_else(sbi_present_num == 0L, 1L, 0L),

      # yardstick event class is first level = "yes"
      case = factor(
        dplyr::if_else(case_num == 1L, "yes", "no"),
        levels = c("yes", "no")
      ),

      # Incoming model prediction is Pr(SBI-positive), so invert
      # to get Pr(SBI-negative), matching the event/case class.
      score_case = 1 - score_yes
    )
}

calc_auroc_aim1 <- function(df) {
  metric_df <- prep_metric_df(df)

  if (dplyr::n_distinct(metric_df$case_num) < 2) {
    return(NA_real_)
  }

  roc_obj <- pROC::roc(
    response = metric_df$case_num,
    predictor = metric_df$score_case,
    quiet = TRUE,
    direction = "<"
  )

  as.numeric(pROC::auc(roc_obj))
}

calc_auprc_aim1 <- function(df) {
  metric_df <- prep_metric_df(df)

  if (dplyr::n_distinct(metric_df$case_num) < 2) {
    return(NA_real_)
  }

  yardstick::pr_auc(
    metric_df,
    case,
    score_case,
    event_level = "first"
  ) %>%
    dplyr::pull(.estimate) %>%
    as.numeric()
}

bootstrap_metric_ci_aim1 <- function(df, metric_fun, n_boot = 1000) {
  metric_df <- prep_metric_df(df)

  if (dplyr::n_distinct(metric_df$case_num) < 2) {
    return(c(NA_real_, NA_real_))
  }

  boot_vals <- replicate(n_boot, {
    idx <- sample.int(nrow(df), size = nrow(df), replace = TRUE)
    boot_df <- df[idx, , drop = FALSE]

    boot_metric_df <- prep_metric_df(boot_df)

    if (dplyr::n_distinct(boot_metric_df$case_num) < 2) {
      return(NA_real_)
    }

    metric_fun(boot_df)
  })

  boot_vals <- boot_vals[is.finite(boot_vals)]

  if (length(boot_vals) < 20) {
    return(c(NA_real_, NA_real_))
  }

  as.numeric(stats::quantile(
    boot_vals,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  ))
}

make_metric_by_hour_aim1 <- function(df, group_label, metric = c("auroc", "auprc"), n_boot = 1000) {
  metric <- match.arg(metric)

  metric_fun <- switch(
    metric,
    "auroc" = calc_auroc_aim1,
    "auprc" = calc_auprc_aim1
  )

  df %>%
    dplyr::group_by(picu_hour) %>%
    dplyr::group_modify(~{
      dat <- .x

      metric_df <- prep_metric_df(dat)

      estimate <- metric_fun(dat)
      ci <- bootstrap_metric_ci_aim1(
        df = dat,
        metric_fun = metric_fun,
        n_boot = n_boot
      )

      tibble::tibble(
        n_total = nrow(metric_df),
        n_sbi_pos = sum(metric_df$case_num == 0L, na.rm = TRUE),
        n_sbi_neg = sum(metric_df$case_num == 1L, na.rm = TRUE),
        prevalence_sbi_negative = mean(metric_df$case_num == 1L, na.rm = TRUE),
        metric_value = estimate,
        ci_low = ci[1],
        ci_high = ci[2]
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      group = group_label,
      metric = toupper(metric)
    )
}

#------------------------------------------------------------
# Create hourly AUROC summaries using matched cohorts
#------------------------------------------------------------

auroc_no_abx <- make_metric_by_hour_aim1(
  df = pred_df_pros_no_abx_24_matched,
  group_label = "No antibiotic exposure before PICU",
  metric = "auroc",
  n_boot = 1000
)

auroc_yes_abx <- make_metric_by_hour_aim1(
  df = pred_df_pros_yes_abx_24_matched,
  group_label = "Antibiotic exposure before PICU",
  metric = "auroc",
  n_boot = 1000
)

auroc_plot_df <- dplyr::bind_rows(auroc_no_abx, auroc_yes_abx) %>%
  dplyr::mutate(
    facet_label = group
  )

#------------------------------------------------------------
# Create hourly AUPRC summaries using matched cohorts
#------------------------------------------------------------

auprc_no_abx <- make_metric_by_hour_aim1(
  df = pred_df_pros_no_abx_24_matched,
  group_label = "No antibiotic exposure before PICU",
  metric = "auprc",
  n_boot = 1000
)

auprc_yes_abx <- make_metric_by_hour_aim1(
  df = pred_df_pros_yes_abx_24_matched,
  group_label = "Antibiotic exposure before PICU",
  metric = "auprc",
  n_boot = 1000
)

auprc_plot_df <- dplyr::bind_rows(auprc_no_abx, auprc_yes_abx) %>%
  dplyr::mutate(
    facet_label = group
  )

#------------------------------------------------------------
# Add facet labels
#------------------------------------------------------------

auroc_plot_df <- auroc_plot_df %>%
  dplyr::mutate(
    facet_label = group
  )

auprc_plot_df <- auprc_plot_df %>%
  dplyr::mutate(
    facet_label = group
  )


#------------------------------------------------------------
# Keep only non-missing metric rows for lines/ribbons
#------------------------------------------------------------

auroc_line_df <- auroc_plot_df %>%
  dplyr::filter(!is.na(metric_value)) %>%
  dplyr::arrange(facet_label, picu_hour)

auprc_line_df <- auprc_plot_df %>%
  dplyr::filter(!is.na(metric_value)) %>%
  dplyr::arrange(facet_label, picu_hour)

# point labels with metrics

auroc_line_df <- auroc_line_df %>%
  dplyr::mutate(
    metric_label = sprintf("%.2f", metric_value)
  )

auprc_line_df <- auprc_line_df %>%
  dplyr::mutate(
    metric_label = sprintf("%.2f", metric_value)
  )

#------------------------------------------------------------
# AUROC plot
#------------------------------------------------------------

p_auroc_facet <- ggplot(auroc_plot_df, aes(x = picu_hour, y = metric_value)) +
  geom_ribbon(
    data = auroc_line_df,
    aes(
      x = picu_hour,
      ymin = ci_low,
      ymax = ci_high,
      group = facet_label
    ),
    inherit.aes = FALSE,
    fill = "lightblue",
    alpha = 0.35
  ) +
  geom_line(
    data = auroc_line_df,
    aes(group = facet_label),
    linewidth = 1.2,
    color = "darkblue"
  ) +
  geom_point(
    data = auroc_line_df,
    size = 2.8,
    color = "darkblue"
  ) +
  facet_wrap(~ facet_label, ncol = 1) +
  scale_x_continuous(breaks = 1:24) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(
    title = "AUROC over PICU time",
    subtitle = "Shaded ribbons show 95% bootstrap confidence intervals",
    x = "PICU hour",
    y = "AUROC"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 13),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    strip.text = element_text(face = "bold", size = 14)
  )

#------------------------------------------------------------
# AUPRC plot
#------------------------------------------------------------

#------------------------------------------------------------
# SBI-negative prevalence reference line for AUPRC plot
# Because AUPRC event/case = SBI-negative, the no-skill baseline
# is the prevalence of SBI-negative encounters.
#------------------------------------------------------------

#------------------------------------------------------------
# SBI-negative prevalence reference line for AUPRC plot
# Uses the exact same matched cohorts used for AUPRC.
#------------------------------------------------------------

auprc_prev_df <- dplyr::bind_rows(
  pred_df_pros_no_abx_24_matched %>%
    dplyr::distinct(study_id, sbi_present) %>%
    dplyr::mutate(facet_label = "No antibiotic exposure before PICU"),

  pred_df_pros_yes_abx_24_matched %>%
    dplyr::distinct(study_id, sbi_present) %>%
    dplyr::mutate(facet_label = "Antibiotic exposure before PICU")
) %>%
  dplyr::mutate(
    sbi_present_num = to_truth01(sbi_present),
    sbi_negative = sbi_present_num == 0L
  ) %>%
  dplyr::group_by(facet_label) %>%
  dplyr::summarise(
    n_encounters = dplyr::n(),
    n_sbi_neg = sum(sbi_negative, na.rm = TRUE),
    sbi_neg_prevalence = n_sbi_neg / n_encounters,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    prevalence_label = paste0(
      "SBI\u2212 prevalence = ",
      scales::percent(sbi_neg_prevalence, accuracy = 1)
    ),
    label_x = 1.3,
    label_y = pmax(sbi_neg_prevalence - 0.035, 0.02)
  )

p_auprc_facet <- ggplot(auprc_plot_df, aes(x = picu_hour, y = metric_value)) +
  geom_hline(
    data = auprc_prev_df,
    aes(yintercept = sbi_neg_prevalence),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.8,
    color = "gray35"
  ) +
  geom_text(
    data = auprc_prev_df,
    aes(
      x = label_x,
      y = label_y,
      label = prevalence_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.8,
    color = "gray25"
  ) +
  geom_ribbon(
    data = auprc_line_df,
    aes(
      x = picu_hour,
      ymin = ci_low,
      ymax = ci_high,
      group = facet_label
    ),
    inherit.aes = FALSE,
    fill = "lightblue",
    alpha = 0.35
  ) +
  geom_line(
    data = auprc_line_df,
    aes(group = facet_label),
    linewidth = 1.2,
    color = "darkblue"
  ) +
  geom_point(
    data = auprc_line_df,
    size = 2.8,
    color = "darkblue"
  ) +
  facet_wrap(~ facet_label, ncol = 1) +
  scale_x_continuous(breaks = 1:24) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(
    title = "AUPRC over PICU time",
    subtitle = "Shaded ribbons show 95% bootstrap confidence intervals; dashed line shows SBI\u2212 prevalence",
    x = "PICU hour",
    y = "AUPRC"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 13),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    strip.text = element_text(face = "bold", size = 14)
  )
# Add labels #
p_auroc_facet <- p_auroc_facet +
  geom_text(
    data = auroc_line_df,
    aes(label = metric_label),
    na.rm = TRUE,
    vjust = -0.8,
    size = 3.5,
    color = "darkblue"
  )

p_auprc_facet <- p_auprc_facet +
  geom_text(
    data = auprc_line_df,
    aes(label = metric_label),
    na.rm = TRUE,
    vjust = -0.8,
    size = 3.5,
    color = "darkblue"
  )

p_auroc_facet
p_auprc_facet
