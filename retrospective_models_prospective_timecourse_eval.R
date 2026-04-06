# Script to prepare all of the prospective dataset for input into the two random forest models. Specifically we will need to retain all
# timepoints, not just the point 2 hours after PICU admission. No need to recreate all of the 'suspicion of infection' code since we
# are not going to look at that subgroup at these other timepoints.

pros_all <- read_csv(file = prospective_model_data_file_path)

# Filter to just the antibiotic unexposed random forest model
pros_rf_ue <- pros_all %>% filter(model_type == "RF_no_abx")

# Filter to just the 24hr after picu admission
pros_24hr_rf_ue <- pros_rf_ue %>% filter(hours_since_picu_adm >= 0) %>% filter(hours_since_picu_adm <= 24)

pros_24hr_rf_ue <- pros_24hr_rf_ue %>%
  mutate(
    picu_hour = floor(hours_since_picu_adm)
  ) %>%
  arrange(study_id, picu_hour, score_time) %>%
  group_by(study_id, picu_hour) %>%
  slice(1) %>%
  ungroup()

#Now select only the same input and output columns as present in the retrospective dataset, make sure names and format match, scales match, etc
thresh_rf_no_abx <- 0.05

required_no_abx <- c(
  "study_id", "picu_hour", "picu_adm_date_time", "hsp_account_id",
  "sbi_present", "model_score", predictors
)

missing_from_source_no_abx <- setdiff(required_no_abx, colnames(pros_24hr_rf_ue))
missing_from_source_no_abx

pros_rf_same_24 <- pros_24hr_rf_ue %>%
  dplyr::select(dplyr::all_of(required_no_abx))

pros_rf_same_24 <- pros_rf_same_24 %>%
  dplyr::relocate(sbi_present, model_score)

pros_rf_same_24$epoch <- "prospective"


# Filter to just the antibiotic unexposed random forest model
pros_rf_ae <- pros_all %>% filter(model_type == "RF_yes_abx")

# Filter to just the 2hr mark
pros_24hr_rf_ae <- pros_rf_ae %>% filter(hours_since_picu_adm >= 0) %>% filter(hours_since_picu_adm <= 24)

pros_24hr_rf_ae <- pros_24hr_rf_ae %>%
  mutate(
    picu_hour = floor(hours_since_picu_adm)
  ) %>%
  arrange(study_id, picu_hour, score_time) %>%
  group_by(study_id, picu_hour) %>%
  slice(1) %>%
  ungroup()

#Now select only the same input and output columns as present in the retrospective dataset, make sure names and format match, scales match, etc
thresh_rf_yes_abx <- 0.074 #

pros_24hr_rf_ae_renamed <- pros_24hr_rf_ae %>%
  dplyr::rename(
    last_fio2 = fio2_last,
    max_dbp = dbp_max
  )

retro_rf_same_abx_renamed <- retro_rf_same_abx %>%
  dplyr::rename(
    hematocrit_mean = hematocrit_blood
  )

required_yes_abx <- c(
  "study_id", "picu_hour", "picu_adm_date_time", "hsp_account_id",
  "sbi_present", "model_score", predictors_abx
)

missing_from_source_yes_abx <- setdiff(required_yes_abx, colnames(pros_24hr_rf_ae_renamed))
missing_from_source_yes_abx

pros_24hr_rf_ae_renamed <- pros_24hr_rf_ae_renamed %>% rename(hematocrit_blood = hematocrit_mean)

pros_rf_same_abx_24 <- pros_24hr_rf_ae_renamed %>%
  dplyr::select(dplyr::all_of(required_yes_abx)) %>%
  dplyr::relocate(sbi_present, model_score)

pros_rf_same_abx_24$epoch <- "prospective"


## Catch up with df names used in aim_1_main (i.e. because we didn't do any of the micro stuff)
pros_no_abx_final_24 <- pros_rf_same_24
pros_yes_abx_final_24 <- pros_rf_same_abx_24

pros_no_abx_final_24 <- pros_no_abx_final_24 %>% relocate(study_id, picu_adm_date_time, picu_hour, sbi_present, model_score)
pros_yes_abx_final_24 <- pros_yes_abx_final_24 %>% relocate(study_id, picu_adm_date_time, picu_hour, sbi_present, model_score)

pros_24_no <- pros_no_abx_final_24
pros_24_yes <- pros_yes_abx_final_24

# Now port in some of the code from retro_model_to_pros_data.R to repeat with all hours of each study_id
# Now test the no_abx model on prospective cohort at all hours
validate_predictor_schema(pros_24_no, predictors, "pros_24_no")
rf_pred_prob_pros_24 <- predict(rf_model, pros_24_no[, all_of(c(predictors))], type = "prob")
pros_auroc_no_abx_24 <- colAUC(rf_pred_prob_pros_24, pros_24_no$sbi_present, plotROC = FALSE) # AUC of 0.736 in prospective test set
pros_auroc_no_abx_24 # overall AUROC 0.715

# Create dataframe to view probabilities and outcomes
pred_df_pros_no_abx_24 <- pros_24_no %>%
  transmute(
    rowIndex = row_number(),
    picu_hour = picu_hour,
    sbi_present = sbi_present,
    pred_prob_yes = rf_pred_prob_pros_24[, "yes"]
  )


# Now test the yes_abx model on prospective 24hr AE cohort
# First fix predictors_abx to match names in pros df
predictors_abx[predictors_abx == "hematocrit_mean"] <- "hematocrit_blood"
pros_yes_temp_24 <- pros_24_yes

validate_predictor_schema(pros_yes_temp_24, predictors_abx, "pros_24_yes")
rf_pred_prob_pros_24_abx <- predict(rf_model_abx, pros_yes_temp_24[, all_of(c(predictors_abx))], type = "prob")
pros_auroc_yes_abx_24 <- colAUC(rf_pred_prob_pros_24_abx, pros_24_yes$sbi_present, plotROC = FALSE)
pros_auroc_yes_abx_24 # overall AUROC 0.735

# Create dataframe to view probabilities and outcomes
pred_df_pros_yes_abx_24 <- pros_24_yes %>%
  transmute(
    rowIndex = row_number(),
    picu_hour = picu_hour,
    sbi_present = sbi_present,
    pred_prob_yes = rf_pred_prob_pros_24_abx[, "yes"]
  )


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
  df = pred_df_pros_no_abx_24,
  threshold = 0.05,
  group_label = "No antibiotic exposure before PICU"
)

npv_yes_abx <- make_npv_by_hour(
  df = pred_df_pros_yes_abx_24,
  threshold = 0.074,
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
  scale_x_continuous(breaks = 0:24) +
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


#####################################################################################################
# ----------------------------
# Match hourly cohorts to the CURRENT p_quadrant cohorts to ensure hourly map works
# ----------------------------

quad_ids_no_abx <- pros_no_abx_1st_infxn %>%
  dplyr::distinct(study_id)

# IMPORTANT:
# Current p_quadrant code appears to use pros_yes_temp for Pros • Abx+
# not pros_yes_abx_1st_infxn
quad_ids_yes_abx <- pros_yes_temp %>%
  dplyr::distinct(study_id)

pros_24_no_matched <- pros_24_no %>%
  dplyr::semi_join(quad_ids_no_abx, by = "study_id")

pros_24_yes_matched <- pros_24_yes %>%
  dplyr::semi_join(quad_ids_yes_abx, by = "study_id")

# Recalculate predictions on the matched hourly cohorts
rf_pred_prob_pros_24_matched <- predict(
  rf_model,
  pros_24_no_matched[, dplyr::all_of(predictors)],
  type = "prob"
)

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
    sbi_present = sbi_present,
    pred_prob_yes = rf_pred_prob_pros_24_matched[, "yes"]
  )

pred_df_pros_yes_abx_24_matched <- pros_24_yes_matched %>%
  dplyr::transmute(
    rowIndex = dplyr::row_number(),
    study_id = study_id,
    picu_hour = picu_hour,
    sbi_present = sbi_present,
    pred_prob_yes = rf_pred_prob_pros_24_abx_matched[, "yes"]
  )

# Rebuild hourly NPV summaries
npv_no_abx_matched <- make_npv_by_hour(
  df = pred_df_pros_no_abx_24_matched,
  threshold = 0.05,
  group_label = "No antibiotic exposure before PICU"
)

npv_yes_abx_matched <- make_npv_by_hour(
  df = pred_df_pros_yes_abx_24_matched,
  threshold = 0.074,
  group_label = "Antibiotic exposure before PICU"
)

npv_plot_df_matched <- dplyr::bind_rows(
  npv_no_abx_matched,
  npv_yes_abx_matched
) %>%
  dplyr::mutate(
    facet_label = paste0(group, "\nThreshold for predicted SBI-negative: ≤ ", threshold),
    tn_label = dplyr::if_else(
      n_pred_negative > 0,
      paste0(n_true_negative, "/", n_pred_negative),
      NA_character_
    )
  )

npv_plot_df_matched <- npv_plot_df_matched %>%
  dplyr::mutate(
    npv_plot = dplyr::if_else(is.na(npv), 0, npv),
    npv_missing = is.na(npv)
  )



line_df_matched <- npv_plot_df_matched %>%
  dplyr::filter(!is.na(npv)) %>%
  dplyr::arrange(facet_label, picu_hour)

line_df_matched <- npv_plot_df_matched %>%
  dplyr::filter(!is.na(npv)) %>%
  dplyr::arrange(facet_label, picu_hour)

p_npv_facet_matched <- ggplot2::ggplot(npv_plot_df_matched, ggplot2::aes(x = picu_hour, y = npv)) +
  ggplot2::geom_ribbon(
    data = line_df_matched,
    ggplot2::aes(x = picu_hour, ymin = ci_low, ymax = ci_high, group = facet_label),
    inherit.aes = FALSE,
    fill = "lightblue",
    alpha = 0.35
  ) +
  ggplot2::geom_line(
    data = line_df_matched,
    ggplot2::aes(group = facet_label),
    linewidth = 1.2,
    color = "darkblue"
  ) +
  geom_point(
    data = npv_plot_df_matched,
    aes(y = npv_plot, shape = npv_missing),
    size = 2.8,
    color = "darkblue",
    fill = "white"
  ) +
  ggplot2::geom_text(
    data = npv_plot_df_matched %>% dplyr::filter(npv_missing),
    ggplot2::aes(y = npv_plot, label = "0/0"),
    vjust = -0.9,
    size = 3.2,
    color = "darkblue"
  ) +
  scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 1),
    guide = "none"
  ) +
  ggplot2::geom_text(
    data = line_df_matched,
    ggplot2::aes(label = tn_label),
    na.rm = TRUE,
    vjust = -0.9,
    size = 3.5,
    color = "darkblue"
  ) +
  ggplot2::facet_wrap(~ facet_label, ncol = 1) +
  ggplot2::scale_x_continuous(breaks = 0:24) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::coord_cartesian(ylim = c(0, 1.08)) +
  ggplot2::labs(
    title = "Negative predictive value over PICU time",
    x = "PICU hour",
    y = "Negative predictive value"
  ) +
  ggplot2::theme_bw(base_size = 16)

p_npv_facet_matched

##########################################################################################





### Repeat but for AUROC and AUPRC by hour ###
library(dplyr)
library(ggplot2)
library(scales)
library(pROC)
library(PRROC)
library(purrr)

#------------------------------------------------------------
# Metric functions
#------------------------------------------------------------

calc_auroc <- function(truth, score) {
  truth_num <- as.integer(truth)

  if (length(unique(truth_num[!is.na(truth_num)])) < 2) {
    return(NA_real_)
  }

  roc_obj <- pROC::roc(
    response = truth_num,
    predictor = score,
    quiet = TRUE,
    direction = "<"
  )

  as.numeric(pROC::auc(roc_obj))
}

calc_auprc <- function(truth, score) {
  truth_num <- as.integer(truth)

  if (length(unique(truth_num[!is.na(truth_num)])) < 2) {
    return(NA_real_)
  }

  pos_scores <- score[truth_num == 1]
  neg_scores <- score[truth_num == 0]

  if (length(pos_scores) == 0 || length(neg_scores) == 0) {
    return(NA_real_)
  }

  pr_obj <- PRROC::pr.curve(
    scores.class0 = pos_scores,
    scores.class1 = neg_scores,
    curve = FALSE
  )

  as.numeric(pr_obj$auc.integral)
}

#------------------------------------------------------------
# Bootstrap CI helper
#------------------------------------------------------------

bootstrap_metric_ci <- function(df, metric_fun, n_boot = 1000, conf_level = 0.95) {

  # Need both classes present
  if (length(unique(df$sbi_present[!is.na(df$sbi_present)])) < 2) {
    return(tibble::tibble(
      estimate = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_
    ))
  }

  estimate <- metric_fun(df$sbi_present, df$pred_prob_yes)

  boot_vals <- replicate(n_boot, {
    idx <- sample.int(n = nrow(df), size = nrow(df), replace = TRUE)
    boot_df <- df[idx, , drop = FALSE]
    metric_fun(boot_df$sbi_present, boot_df$pred_prob_yes)
  })

  boot_vals <- boot_vals[!is.na(boot_vals)]

  alpha <- 1 - conf_level

  if (length(boot_vals) < 20) {
    return(tibble::tibble(
      estimate = estimate,
      ci_low = NA_real_,
      ci_high = NA_real_
    ))
  }

  tibble::tibble(
    estimate = estimate,
    ci_low = as.numeric(stats::quantile(boot_vals, probs = alpha / 2, na.rm = TRUE)),
    ci_high = as.numeric(stats::quantile(boot_vals, probs = 1 - alpha / 2, na.rm = TRUE))
  )
}

#------------------------------------------------------------
# Generic summary-by-hour function
#------------------------------------------------------------

make_metric_by_hour <- function(df, group_label, metric = c("auroc", "auprc"), n_boot = 1000) {

  metric <- match.arg(metric)

  metric_fun <- switch(
    metric,
    "auroc" = calc_auroc,
    "auprc" = calc_auprc
  )

  out <- df %>%
    dplyr::group_by(picu_hour) %>%
    dplyr::group_modify(~{
      dat <- .x

      ci_res <- bootstrap_metric_ci(
        df = dat,
        metric_fun = metric_fun,
        n_boot = n_boot,
        conf_level = 0.95
      )

      tibble::tibble(
        n_total = nrow(dat),
        n_sbi_pos = sum(dat$sbi_present == 1, na.rm = TRUE),
        n_sbi_neg = sum(dat$sbi_present == 0, na.rm = TRUE),
        metric_value = ci_res$estimate,
        ci_low = ci_res$ci_low,
        ci_high = ci_res$ci_high
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      group = group_label,
      metric = toupper(metric)
    )

  out
}

#------------------------------------------------------------
# Create hourly AUROC summaries
#------------------------------------------------------------

auroc_no_abx <- make_metric_by_hour(
  df = pred_df_pros_no_abx_24,
  group_label = "No antibiotic exposure before PICU",
  metric = "auroc",
  n_boot = 1000
)

auroc_yes_abx <- make_metric_by_hour(
  df = pred_df_pros_yes_abx_24,
  group_label = "Antibiotic exposure before PICU",
  metric = "auroc",
  n_boot = 1000
)

auroc_plot_df <- dplyr::bind_rows(auroc_no_abx, auroc_yes_abx)

#------------------------------------------------------------
# Create hourly AUPRC summaries
#------------------------------------------------------------

auprc_no_abx <- make_metric_by_hour(
  df = pred_df_pros_no_abx_24,
  group_label = "No antibiotic exposure before PICU",
  metric = "auprc",
  n_boot = 1000
)

auprc_yes_abx <- make_metric_by_hour(
  df = pred_df_pros_yes_abx_24,
  group_label = "Antibiotic exposure before PICU",
  metric = "auprc",
  n_boot = 1000
)

auprc_plot_df <- dplyr::bind_rows(auprc_no_abx, auprc_yes_abx)

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

# Omit hour 0
auroc_plot_df_no0 <- auroc_plot_df %>%
  dplyr::filter(picu_hour != 0)

auprc_plot_df_no0 <- auprc_plot_df %>%
  dplyr::filter(picu_hour != 0)

auroc_line_df_no0 <- auroc_line_df %>%
  dplyr::filter(picu_hour != 0)

auprc_line_df_no0 <- auprc_line_df %>%
  dplyr::filter(picu_hour != 0)

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
  scale_x_continuous(breaks = 0:24) +
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

p_auprc_facet <- ggplot(auprc_plot_df, aes(x = picu_hour, y = metric_value)) +
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
  scale_x_continuous(breaks = 0:24) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(
    title = "AUPRC over PICU time",
    subtitle = "Shaded ribbons show 95% bootstrap confidence intervals",
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

# If we want hour 0 plotted
p_auroc_facet
p_auprc_facet

# Plot the non 0 hour plots
p_auroc_facet_no0 <- ggplot(auroc_plot_df_no0, aes(x = picu_hour, y = metric_value)) +
  geom_ribbon(
    data = auroc_line_df_no0,
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
    data = auroc_line_df_no0,
    aes(group = facet_label),
    linewidth = 1.2,
    color = "darkblue"
  ) +
  geom_point(
    data = auroc_line_df_no0,
    size = 2.8,
    color = "darkblue"
  ) +
  geom_text(
    data = auroc_line_df_no0,
    aes(label = metric_label),
    na.rm = TRUE,
    vjust = -0.8,
    size = 3.5,
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

p_auprc_facet_no0 <- ggplot(auprc_plot_df_no0, aes(x = picu_hour, y = metric_value)) +
  geom_ribbon(
    data = auprc_line_df_no0,
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
    data = auprc_line_df_no0,
    aes(group = facet_label),
    linewidth = 1.2,
    color = "darkblue"
  ) +
  geom_point(
    data = auprc_line_df_no0,
    size = 2.8,
    color = "darkblue"
  ) +
  geom_text(
    data = auprc_line_df_no0,
    aes(label = metric_label),
    na.rm = TRUE,
    vjust = -0.8,
    size = 3.5,
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
    subtitle = "Shaded ribbons show 95% bootstrap confidence intervals",
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

p_auroc_facet_no0
p_auprc_facet_no0