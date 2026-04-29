##### Probability based decision policy exploration ######

### =========================================================
### Sequential SBI policy grid search
### - tune on rf_valid_df
### - lock best policy
### - evaluate once on rf_test_df
### =========================================================

### Packages
library(dplyr)
library(tidyr)
library(tibble)
library(future)
library(furrr)

n_workers <- max(1L, as.integer(parallelly::availableCores()) - 1L)

future::plan(
  future::multisession,
  workers = n_workers
)

message("Using ", n_workers, " parallel workers")


### ---------------------------------------------------------
### 1) Prep data
###    Conservative choice:
###    - expand each study_id to every clock hour 0:24
###    - if an hour is missing, it remains NA
###    - missing hours break multi-hour rules
### ---------------------------------------------------------

prep_seq_data <- function(df) {

  base_df <- df %>%
    dplyr::select(study_id, sbi_present, hours_since_picu_adm, rf_prob) %>%
    dplyr::mutate(
      study_id = as.character(study_id),
      sbi_present = as.integer(sbi_present),
      hours_since_picu_adm = as.integer(hours_since_picu_adm),
      rf_prob = as.numeric(rf_prob)
    ) %>%
    dplyr::filter(
      !is.na(study_id),
      !is.na(sbi_present),
      !is.na(hours_since_picu_adm),
      hours_since_picu_adm >= 0,
      hours_since_picu_adm <= 24
    ) %>%
    dplyr::group_by(study_id, hours_since_picu_adm) %>%
    dplyr::summarise(
      sbi_present = dplyr::first(sbi_present),
      rf_prob = mean(rf_prob, na.rm = TRUE),
      .groups = "drop"
    )

  truth_df <- base_df %>%
    dplyr::group_by(study_id) %>%
    dplyr::summarise(
      sbi_present = dplyr::first(sbi_present),
      .groups = "drop"
    )

  seq_df <- tidyr::expand_grid(
    study_id = unique(base_df$study_id),
    hours_since_picu_adm = 0:24
  ) %>%
    dplyr::left_join(
      base_df %>% dplyr::select(study_id, hours_since_picu_adm, rf_prob),
      by = c("study_id", "hours_since_picu_adm")
    ) %>%
    dplyr::left_join(
      truth_df,
      by = "study_id"
    ) %>%
    dplyr::arrange(study_id, hours_since_picu_adm)

  seq_df
}

rf_valid_seq <- prep_seq_data(rf_valid_df)
rf_test_seq  <- prep_seq_data(rf_test_df)

build_patient_list <- function(seq_df) {
  split_df <- split(seq_df, seq_df$study_id)

  lapply(
    split_df,
    function(pt_df) {
      pt_df <- pt_df %>% dplyr::arrange(hours_since_picu_adm)

      list(
        study_id = pt_df$study_id[[1]],
        sbi_present = as.integer(pt_df$sbi_present[[1]]),
        probs = as.numeric(pt_df$rf_prob),
        hours = as.integer(pt_df$hours_since_picu_adm)
      )
    }
  )
}

valid_patients <- build_patient_list(rf_valid_seq)
test_patients  <- build_patient_list(rf_test_seq)

### ---------------------------------------------------------
### 2) Rule helper functions
### ---------------------------------------------------------

low_rule_met <- function(probs, i, params) {

  low_rule <- params$low_rule[[1]]
  low_threshold <- params$low_threshold[[1]]

  if (low_rule == "single_low") {
    current_val <- probs[i]

    if (is.na(current_val)) {
      return(FALSE)
    }

    return(current_val <= low_threshold)
  }

  if (low_rule == "consecutive") {
    low_k <- params$low_k[[1]]

    if (is.na(low_k) || i < low_k) {
      return(FALSE)
    }

    window_vals <- probs[(i - low_k + 1):i]

    if (any(is.na(window_vals))) {
      return(FALSE)
    }

    return(all(window_vals <= low_threshold))
  }

  if (low_rule == "m_of_n") {
    low_m <- params$low_m[[1]]
    low_n <- params$low_n[[1]]

    if (any(is.na(c(low_m, low_n))) || i < low_n) {
      return(FALSE)
    }

    window_vals <- probs[(i - low_n + 1):i]

    if (any(is.na(window_vals))) {
      return(FALSE)
    }

    return(sum(window_vals <= low_threshold) >= low_m)
  }

  stop("Unknown low_rule.")
}

high_rule_met <- function(probs, i, params) {

  high_rule <- params$high_rule[[1]]
  high_threshold <- params$high_threshold[[1]]

  if (high_rule == "single_high") {
    current_val <- probs[i]

    if (is.na(current_val)) {
      return(FALSE)
    }

    return(current_val >= high_threshold)
  }

  if (high_rule == "m_of_n_high") {
    high_m <- params$high_m[[1]]
    high_n <- params$high_n[[1]]

    if (any(is.na(c(high_m, high_n))) || i < high_n) {
      return(FALSE)
    }

    window_vals <- probs[(i - high_n + 1):i]

    if (any(is.na(window_vals))) {
      return(FALSE)
    }

    return(sum(window_vals >= high_threshold) >= high_m)
  }

  if (high_rule == "consecutive_high") {
    high_k <- params$high_k[[1]]

    if (is.na(high_k) || i < high_k) {
      return(FALSE)
    }

    window_vals <- probs[(i - high_k + 1):i]

    if (any(is.na(window_vals))) {
      return(FALSE)
    }

    return(all(window_vals >= high_threshold))
  }

  stop("Unknown high_rule.")
}

### ---------------------------------------------------------
### 3) Apply one policy to one patient
###    High-risk stop rule takes precedence over rule-out
### ---------------------------------------------------------

apply_policy_one_patient <- function(pt, params) {

  probs <- pt$probs
  hours <- pt$hours

  final_state <- "indeterminate"
  decision_hour <- NA_integer_

  for (i in seq_along(probs)) {

    ## First: stop rule for high-risk / not eligible
    if (high_rule_met(probs = probs, i = i, params = params)) {
      final_state <- "not_eligible"
      decision_hour <- hours[i]
      break
    }

    ## Second: rule-out if at/after minimum hour
    if (hours[i] >= params$min_ruleout_hour[[1]]) {
      if (low_rule_met(probs = probs, i = i, params = params)) {
        final_state <- "ruled_out"
        decision_hour <- hours[i]
        break
      }
    }
  }

  list(
    study_id = pt$study_id,
    sbi_present = pt$sbi_present,
    final_state = final_state,
    decision_hour = decision_hour
  )
}

apply_policy_dataset <- function(patient_list, params) {
  n <- length(patient_list)

  study_id <- character(n)
  sbi_present <- integer(n)
  final_state <- character(n)
  decision_hour <- rep(NA_integer_, n)

  for (j in seq_len(n)) {
    res_j <- apply_policy_one_patient(
      pt = patient_list[[j]],
      params = params
    )

    study_id[[j]] <- res_j$study_id
    sbi_present[[j]] <- res_j$sbi_present
    final_state[[j]] <- res_j$final_state
    decision_hour[[j]] <- res_j$decision_hour
  }

  data.frame(
    study_id = study_id,
    sbi_present = sbi_present,
    final_state = final_state,
    decision_hour = decision_hour,
    stringsAsFactors = FALSE
  )
}

### ---------------------------------------------------------
### 4) Summarize one policy at the encounter level
### ---------------------------------------------------------

summarise_policy <- function(decisions_df, params) {

  n_patients <- nrow(decisions_df)
  sbi_present <- decisions_df$sbi_present
  final_state <- decisions_df$final_state
  decision_hour <- decisions_df$decision_hour

  is_true_neg <- sbi_present == 0
  is_true_pos <- sbi_present == 1
  is_ruled_out <- final_state == "ruled_out"
  is_not_eligible <- final_state == "not_eligible"
  is_indeterminate <- final_state == "indeterminate"

  n_true_neg <- sum(is_true_neg, na.rm = TRUE)
  n_true_pos <- sum(is_true_pos, na.rm = TRUE)

  n_ruled_out <- sum(is_ruled_out, na.rm = TRUE)
  n_ruled_out_true_neg <- sum(is_ruled_out & is_true_neg, na.rm = TRUE)
  n_ruled_out_false_neg <- sum(is_ruled_out & is_true_pos, na.rm = TRUE)

  if (n_ruled_out > 0) {
    npv <- n_ruled_out_true_neg / n_ruled_out

    ## Exact 95% CI for NPV
    npv_ci <- stats::binom.test(
      x = n_ruled_out_true_neg,
      n = n_ruled_out,
      conf.level = 0.95
    )$conf.int

    npv_lower <- npv_ci[1]
    npv_upper <- npv_ci[2]
  } else {
    npv <- NA_real_
    npv_lower <- NA_real_
    npv_upper <- NA_real_
  }

  n_not_eligible <- sum(is_not_eligible, na.rm = TRUE)

  if (n_not_eligible > 0) {
    ppv_not_eligible <- mean(sbi_present[is_not_eligible] == 1, na.rm = TRUE)
  } else {
    ppv_not_eligible <- NA_real_
  }

  tibble::tibble(
    policy_id = params$policy_id[[1]],
    min_ruleout_hour = params$min_ruleout_hour[[1]],
    low_threshold = params$low_threshold[[1]],
    low_rule = params$low_rule[[1]],
    low_k = params$low_k[[1]],
    low_m = params$low_m[[1]],
    low_n = params$low_n[[1]],
    high_threshold = params$high_threshold[[1]],
    high_rule = params$high_rule[[1]],
    high_k = params$high_k[[1]],
    high_m = params$high_m[[1]],
    high_n = params$high_n[[1]],

    n_patients = n_patients,
    n_true_neg = n_true_neg,
    n_true_pos = n_true_pos,

    n_ruled_out = n_ruled_out,
    n_ruled_out_true_neg = n_ruled_out_true_neg,
    n_ruled_out_false_neg = n_ruled_out_false_neg,

    npv = npv,
    npv_lower = npv_lower,
    npv_upper = npv_upper,

    prop_all_patients_ruled_out =
      ifelse(n_patients > 0, n_ruled_out / n_patients, NA_real_),

    prop_true_neg_ruled_out =
      ifelse(n_true_neg > 0, n_ruled_out_true_neg / n_true_neg, NA_real_),

    prop_true_pos_ruleout_error =
      ifelse(n_true_pos > 0, n_ruled_out_false_neg / n_true_pos, NA_real_),

    median_hour_ruleout_all =
      ifelse(n_ruled_out > 0,
             stats::median(decision_hour[is_ruled_out], na.rm = TRUE),
             NA_real_),

    median_hour_ruleout_true_neg =
      ifelse(sum(is_ruled_out & is_true_neg, na.rm = TRUE) > 0,
             stats::median(
               decision_hour[is_ruled_out & is_true_neg],
               na.rm = TRUE
             ),
             NA_real_),

    n_not_eligible = n_not_eligible,

    prop_true_pos_not_eligible =
      ifelse(n_true_pos > 0,
             sum(is_not_eligible & is_true_pos, na.rm = TRUE) / n_true_pos,
             NA_real_),

    prop_true_neg_not_eligible =
      ifelse(n_true_neg > 0,
             sum(is_not_eligible & is_true_neg, na.rm = TRUE) / n_true_neg,
             NA_real_),

    ppv_not_eligible = ppv_not_eligible,

    prop_indeterminate =
      ifelse(n_patients > 0, sum(is_indeterminate, na.rm = TRUE) / n_patients, NA_real_),

    prop_true_neg_indeterminate =
      ifelse(n_true_neg > 0,
             sum(is_indeterminate & is_true_neg, na.rm = TRUE) / n_true_neg,
             NA_real_),

    prop_true_pos_indeterminate =
      ifelse(n_true_pos > 0,
             sum(is_indeterminate & is_true_pos, na.rm = TRUE) / n_true_pos,
             NA_real_)
  )
}

### ---------------------------------------------------------
### 5) Build the policy grid
###    Lower thresholds centered around 0.12
### ---------------------------------------------------------

low_rule_grid <- dplyr::bind_rows(
  tibble::tibble(
    low_rule = "single_low",
    low_k = NA_integer_,
    low_m = NA_integer_,
    low_n = NA_integer_
  ),
  tibble::tibble(
    low_rule = "consecutive",
    low_k = c(2L, 3L),
    low_m = NA_integer_,
    low_n = NA_integer_
  ),
  tibble::tibble(
    low_rule = "m_of_n",
    low_k = NA_integer_,
    low_m = c(2L, 2L, 3L),
    low_n = c(3L, 4L, 4L)
  )
)

high_rule_grid <- dplyr::bind_rows(
  tibble::tibble(
    high_rule = "single_high",
    high_k = NA_integer_,
    high_m = NA_integer_,
    high_n = NA_integer_
  ),
  tibble::tibble(
    high_rule = "m_of_n_high",
    high_k = NA_integer_,
    high_m = c(2L, 2L, 3L),
    high_n = c(3L, 4L, 4L)
  ),
  tibble::tibble(
    high_rule = "consecutive_high",
    high_k = c(2L, 3L),
    high_m = NA_integer_,
    high_n = NA_integer_
  )
)

# Match the legacy search space from old_decision_policy_code.R exactly
min_ruleout_hour_grid <- 1
low_threshold_grid <- seq(0.05, 0.20, by = 0.01)
high_threshold_grid <- seq(0.25, 0.90, by = 0.05)

policy_grid <- tidyr::crossing(
  min_ruleout_hour = min_ruleout_hour_grid,
  low_threshold = low_threshold_grid,
  high_threshold = high_threshold_grid,
  low_rule_grid,
  high_rule_grid
) %>%
  dplyr::filter(high_threshold > low_threshold) %>%
  dplyr::mutate(policy_id = dplyr::row_number())

### ---------------------------------------------------------
### 6) Run grid search on VALIDATION set
### ---------------------------------------------------------

run_grid_search <- function(patient_list, policy_grid) {

  furrr::future_map_dfr(
    seq_len(nrow(policy_grid)),
    function(i) {

      if (i %% 100 == 0) {
        message("Completed ", i, " / ", nrow(policy_grid), " policies")
      }

      params_i <- policy_grid[i, , drop = FALSE]

      decisions_i <- apply_policy_dataset(
        patient_list = patient_list,
        params = params_i
      )

      summarise_policy(
        decisions_df = decisions_i,
        params = params_i
      )
    },
    .options = furrr::furrr_options(seed = TRUE)
  )
}

valid_grid_results <- run_grid_search(
  patient_list = valid_patients,
  policy_grid = policy_grid
)

### ---------------------------------------------------------
### ---------------------------------------------------------
### 7) Rank policies
###    Primary constraint: point estimate NPV >= 0.95
###    Then maximize rule-out among true SBI-negative patients
###    Then favor earlier rule-out
### ---------------------------------------------------------

valid_grid_ranked <- valid_grid_results %>%
  dplyr::mutate(
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
  ) %>%
  dplyr::arrange(
    dplyr::desc(meets_npv_target),
    dplyr::desc(coverage_score),
    speed_score,
    dplyr::desc(npv),
    dplyr::desc(n_ruled_out),
    dplyr::desc(prop_true_pos_not_eligible)
  )

top_valid_policies <- valid_grid_ranked %>%
  dplyr::select(
    policy_id,
    meets_npv_target,
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
    npv,
    npv_lower,
    npv_upper,
    n_ruled_out,
    n_ruled_out_true_neg,
    n_ruled_out_false_neg,
    prop_true_neg_ruled_out,
    median_hour_ruleout_true_neg,
    prop_true_pos_not_eligible,
    prop_true_neg_not_eligible,
    prop_indeterminate
  ) %>%
  dplyr::slice_head(n = 25)

best_policy <- valid_grid_ranked %>%
  dplyr::slice(1)

top_valid_policies
best_policy

### ---------------------------------------------------------
### 8) Apply the single winning policy to VALIDATION and TEST
### ---------------------------------------------------------

validation_decisions_best <- apply_policy_dataset(
  patient_list = valid_patients,
  params = best_policy
)

test_decisions_best <- apply_policy_dataset(
  patient_list = test_patients,
  params = best_policy
)

validation_best_summary <- summarise_policy(
  decisions_df = validation_decisions_best,
  params = best_policy
) %>%
  dplyr::mutate(dataset = "validation")

test_best_summary <- summarise_policy(
  decisions_df = test_decisions_best,
  params = best_policy
) %>%
  dplyr::mutate(dataset = "test")

best_policy_performance <- dplyr::bind_rows(
  validation_best_summary,
  test_best_summary
)

best_policy_performance

### ---------------------------------------------------------
### 9) Optional: inspect encounter-level decisions from best policy
### ---------------------------------------------------------

validation_decisions_best %>%
  dplyr::count(final_state, sbi_present)

test_decisions_best %>%
  dplyr::count(final_state, sbi_present)

### ---------------------------------------------------------
### 10) Optional: merge best-policy decision back to patient-hour data
###     so you can inspect example trajectories
### ---------------------------------------------------------

valid_trajectories_with_decision <- rf_valid_seq %>%
  dplyr::left_join(
    validation_decisions_best,
    by = c("study_id", "sbi_present")
  )

test_trajectories_with_decision <- rf_test_seq %>%
  dplyr::left_join(
    test_decisions_best,
    by = c("study_id", "sbi_present")
  )

#### Now create figure for manuscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

## -----------------------------
## Prep labels for final states
## -----------------------------
plot_decisions <- test_decisions_best %>%
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
panel_a_df <- plot_decisions %>%
  dplyr::count(truth_label, final_state_label, name = "n") %>%
  dplyr::group_by(truth_label) %>%
  dplyr::mutate(
    prop = n / sum(n),
    label = paste0(scales::percent(prop, accuracy = 1), "\n(n=", n, ")")
  ) %>%
  dplyr::ungroup()

p_a <- ggplot2::ggplot(
  panel_a_df,
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
    y = "Patients within SBI-negative group",
    fill = NULL,
    title = "Policy disposition by true SBI status"
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
## For ruled_out patients, decision_hour is the hour of rule-out.
## For others, they contribute 0 to cumulative rule-out.

cum_df <- plot_decisions %>%
  dplyr::filter(!is.na(truth_label)) %>%
  dplyr::select(study_id, truth_label, final_state, decision_hour)

hours_df <- tidyr::expand_grid(
  study_id = unique(cum_df$study_id),
  hour = 0:24
) %>%
  dplyr::left_join(
    cum_df,
    by = "study_id"
  ) %>%
  dplyr::mutate(
    cum_ruled_out = dplyr::case_when(
      final_state == "ruled_out" & !is.na(decision_hour) & decision_hour <= hour ~ 1,
      TRUE ~ 0
    )
  )

panel_b_df <- hours_df %>%
  dplyr::group_by(truth_label, hour) %>%
  dplyr::summarise(
    prop_cum_ruled_out = mean(cum_ruled_out, na.rm = TRUE),
    .groups = "drop"
  )

p_b <- ggplot2::ggplot(
  panel_b_df,
  ggplot2::aes(x = hour, y = prop_cum_ruled_out)
) +
  ggplot2::geom_area(
    data = panel_b_df %>% dplyr::filter(truth_label == "SBI-negative"),
    fill = "blue",
    alpha = 0.18
  ) +
  ggplot2::geom_area(
    data = panel_b_df %>% dplyr::filter(truth_label == "SBI-positive"),
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
    title = "Cumulative rule-out over time"
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

## ----------------------------------------
## Combine panels
## ----------------------------------------
policy_figure <- p_a + p_b + patchwork::plot_layout(widths = c(1, 1.15))

policy_figure

future::plan(future::sequential)
