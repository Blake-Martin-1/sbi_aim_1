source("plot_save_helpers.R")
#### Apply new nmodel to prospective results from Apr 2026 ####

updated_25 <- readr::read_csv(
  file = file.path(sbi_blake_phi_path, "model_events_parsed_rf.csv"),
  na = c("", "NULL"),
  col_types = readr::cols(.default = readr::col_character()),
  progress = FALSE,
  lazy = FALSE
)

## Drop columns ending in _RAW
updated_25_noraw <- updated_25[, !grepl("_RAW$", names(updated_25)), drop = FALSE]
updated_25_noraw <- updated_25_noraw[, !grepl("_TS_LEN$", names(updated_25_noraw)), drop = FALSE]
updated_25_noraw <- updated_25_noraw[, !grepl("_IS_COMPLETE$", names(updated_25_noraw)), drop = FALSE]

# Remove remainder of NA columns
all_null_cols <- names(updated_25_noraw)[
  sapply(updated_25_noraw, function(x) all(is.na(x)))
]

all_null_cols
length(all_null_cols)

## Create a mapping table so you can inspect old -> new names
name_map <- data.frame(
  old_name = names(updated_25_noraw),
  new_name = tolower(names(updated_25_noraw)),
  in_predictors = tolower(names(updated_25_noraw)) %in% predictors,
  stringsAsFactors = FALSE
)

## Apply the rename
updated_25_renamed <- updated_25_noraw
names(updated_25_renamed) <- name_map$new_name

## Optional: move predictor columns to the front, keep all others after them
updated_25_renamed <- updated_25_renamed %>%
  dplyr::select(dplyr::any_of(predictors), dplyr::everything())

## Optional checks
missing_predictors <- setdiff(predictors, names(updated_25_renamed))
extra_columns <- setdiff(names(updated_25_renamed), predictors)

## View results
name_map %>% dplyr::arrange(dplyr::desc(in_predictors), new_name)
missing_predictors
extra_columns

# Fix score_time column
updated_25_renamed <- updated_25_renamed %>%
  dplyr::mutate(
    picu_adm_date_time = lubridate::ymd_hms(
      icu_start_instant_str,
      tz = "America/Denver"
    )
  )

# Create score time via picu_adm_date_time and time since PICU admission
updated_25_renamed <- updated_25_renamed %>%
  dplyr::mutate(
    hours_since_icu_num = as.numeric(hours_since_icu),
    score_time = picu_adm_date_time + lubridate::dhours(hours_since_icu_num)
  )

# REmove duplicate columns and rows without PICU admission times
updated_25_renamed <- updated_25_renamed %>% dplyr::select(-age_y, -pccc_y) %>% rename(age = age_x, pccc = pccc_x)
updated_25_renamed <- updated_25_renamed %>% filter(!is.na(picu_adm_date_time))
updated_25_renamed <- updated_25_renamed %>% filter(!is.na(hours_since_icu_num))
updated_25_renamed <- updated_25_renamed %>% filter(!is.na(pccc))

## Ok now process data to make sure format matches the dataframe used by the models ##
pros_25 <- updated_25_renamed

# Lower case of column names
colnames(pros_25) <- str_to_lower(colnames(pros_25))

# Rearrange important cols to front
pros_25 <- pros_25 %>% relocate(csn, picu_adm_date_time, score_time, hours_since_icu_num, score_value)
pros_25 <- pros_25 %>% dplyr::arrange(csn, picu_adm_date_time, hours_since_icu_num)

# Change name to identify which model generated which scores
pros_25 <- pros_25 %>% rename(auto_score = score_value) #auto score is that generated in Epic

# Now ensure all the classes of the predictor columns match up
for (col in predictors) {

  target_class <- class(test_df[[col]])

  if (inherits(test_df[[col]], "factor")) {

    pros_25[[col]] <- factor(
      pros_25[[col]],
      levels = levels(test_df[[col]]),
      ordered = is.ordered(test_df[[col]])
    )

  } else if (inherits(test_df[[col]], "POSIXct")) {

    pros_25[[col]] <- as.POSIXct(
      pros_25[[col]],
      tz = attr(test_df[[col]], "tzone")
    )

  } else if (inherits(test_df[[col]], "Date")) {

    pros_25[[col]] <- as.Date(pros_25[[col]])

  } else if (inherits(test_df[[col]], "integer")) {

    pros_25[[col]] <- as.integer(pros_25[[col]])

  } else if (inherits(test_df[[col]], "numeric")) {

    pros_25[[col]] <- as.numeric(pros_25[[col]])

  } else if (inherits(test_df[[col]], "logical")) {

    pros_25[[col]] <- as.logical(pros_25[[col]])

  } else if (inherits(test_df[[col]], "character")) {

    pros_25[[col]] <- as.character(pros_25[[col]])

  }
}

# Now generate the scores using local model version
newdata <- pros_25 %>%
  dplyr::select(dplyr::all_of(predictors))

offline_prob <- predict(
  rf_tune,
  newdata = newdata,
  type = "prob"
)

offline_scores <- offline_prob[, "pos"]

comp_scores <- pros_25 %>% mutate(local_score = offline_scores)
comp_scores$auto_score <- as.numeric(comp_scores$auto_score) / 100

# make histogram plots of the two score types
comp_scores_long <- comp_scores %>%
  dplyr::select(auto_score, local_score) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "score_type",
    values_to = "score"
  ) %>%
  dplyr::mutate(
    score_type = dplyr::recode(
      score_type,
      auto_score = "Epic auto score",
      local_score = "Local model score"
    )
  )

p_score_density_comparison <- ggplot(comp_scores_long, aes(x = score, color = score_type, fill = score_type)) +
  geom_density(alpha = 0.25, linewidth = 1.1, na.rm = TRUE) +
  scale_color_manual(
    values = c(
      "Local model score" = "red",
      "Epic auto score" = "blue"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Local model score" = "red",
      "Epic auto score" = "blue"
    )
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "Distribution of Epic Auto Scores vs Local Model Scores",
    x = "Predicted probability of SBI",
    y = "Density",
    color = "Score source",
    fill = "Score source"
  ) +
  theme_minimal(base_size = 14)

p_score_density_comparison
save_aim1_plot(p_score_density_comparison, "epic_auto_vs_local_model_score_density.tiff")

# Now perform analysis of why and when scores differ
hour_col <- "hours_since_icu_num"

score_diff_df <- comp_scores %>%
  dplyr::select(
    csn,
    hour = dplyr::all_of(hour_col),
    auto_score,
    local_score
  ) %>%
  dplyr::filter(
    !is.na(csn),
    !is.na(hour),
    !is.na(auto_score),
    !is.na(local_score)
  ) %>%
  dplyr::mutate(
    hour = as.numeric(hour),
    auto_score = as.numeric(auto_score),
    local_score = as.numeric(local_score),

    # Main difference metrics
    diff = local_score - auto_score,
    abs_diff = abs(diff),
    sq_diff = diff^2,

    # Bland-Altman x-axis
    mean_score = (local_score + auto_score) / 2,

    # Direction
    local_higher = local_score > auto_score,
    auto_higher = auto_score > local_score
  )

# Calculate overall summary
overall_diff_summary <- score_diff_df %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_csn = dplyr::n_distinct(csn),

    mean_auto_score = mean(auto_score),
    mean_local_score = mean(local_score),

    mean_diff = mean(diff),
    sd_diff = sd(diff),
    median_diff = median(diff),
    q1_diff = stats::quantile(diff, 0.25),
    q3_diff = stats::quantile(diff, 0.75),
    min_diff = min(diff),
    max_diff = max(diff),

    mean_abs_diff = mean(abs_diff),
    median_abs_diff = median(abs_diff),
    q90_abs_diff = stats::quantile(abs_diff, 0.90),
    q95_abs_diff = stats::quantile(abs_diff, 0.95),
    q99_abs_diff = stats::quantile(abs_diff, 0.99),
    max_abs_diff = max(abs_diff),

    rmse = sqrt(mean(sq_diff)),

    pearson_cor = stats::cor(auto_score, local_score, method = "pearson"),
    spearman_cor = stats::cor(auto_score, local_score, method = "spearman"),

    pct_abs_diff_gt_0_001 = mean(abs_diff > 0.001) * 100,
    pct_abs_diff_gt_0_005 = mean(abs_diff > 0.005) * 100,
    pct_abs_diff_gt_0_01 = mean(abs_diff > 0.01) * 100,
    pct_abs_diff_gt_0_05 = mean(abs_diff > 0.05) * 100,
    pct_abs_diff_gt_0_10 = mean(abs_diff > 0.10) * 100,

    pct_local_higher = mean(local_higher) * 100,
    pct_auto_higher = mean(auto_higher) * 100
  )

overall_diff_summary

# Now analyze hour specific differences
diff_by_hour <- score_diff_df %>%
  dplyr::group_by(hour) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_csn = dplyr::n_distinct(csn),

    mean_auto_score = mean(auto_score),
    mean_local_score = mean(local_score),

    mean_diff = mean(diff),
    sd_diff = sd(diff),
    se_diff = sd(diff) / sqrt(dplyr::n()),
    lower_95_mean_diff = mean_diff - 1.96 * se_diff,
    upper_95_mean_diff = mean_diff + 1.96 * se_diff,

    median_diff = median(diff),
    q1_diff = stats::quantile(diff, 0.25),
    q3_diff = stats::quantile(diff, 0.75),

    mean_abs_diff = mean(abs_diff),
    median_abs_diff = median(abs_diff),
    q90_abs_diff = stats::quantile(abs_diff, 0.90),
    q95_abs_diff = stats::quantile(abs_diff, 0.95),

    rmse = sqrt(mean(sq_diff)),

    pct_abs_diff_gt_0_01 = mean(abs_diff > 0.01) * 100,
    pct_abs_diff_gt_0_05 = mean(abs_diff > 0.05) * 100,
    pct_abs_diff_gt_0_10 = mean(abs_diff > 0.10) * 100,

    pct_local_higher = mean(local_higher) * 100,
    pct_auto_higher = mean(auto_higher) * 100,

    .groups = "drop"
  ) %>%
  dplyr::arrange(hour)

diff_by_hour

# Now look at encounter level differences
diff_by_csn <- score_diff_df %>%
  dplyr::group_by(csn) %>%
  dplyr::summarise(
    n_scores = dplyr::n(),
    first_hour = min(hour),
    last_hour = max(hour),

    mean_auto_score = mean(auto_score),
    mean_local_score = mean(local_score),

    mean_diff = mean(diff),
    median_diff = median(diff),

    mean_abs_diff = mean(abs_diff),
    median_abs_diff = median(abs_diff),
    max_abs_diff = max(abs_diff),

    rmse = sqrt(mean(sq_diff)),

    pct_local_higher = mean(local_higher) * 100,
    pct_auto_higher = mean(auto_higher) * 100,

    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(mean_abs_diff))

diff_by_csn

# Now identify hour of disagreement by encounter
max_diff_per_csn <- score_diff_df %>%
  dplyr::group_by(csn) %>%
  dplyr::slice_max(order_by = abs_diff, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    csn,
    hour,
    auto_score,
    local_score,
    diff,
    abs_diff
  ) %>%
  dplyr::arrange(dplyr::desc(abs_diff))

max_diff_per_csn

# Summarixze when these disagreements occur
max_diff_hour_summary <- max_diff_per_csn %>%
  dplyr::count(hour, name = "n_csn_with_max_diff") %>%
  dplyr::mutate(
    pct_csn_with_max_diff = 100 * n_csn_with_max_diff / sum(n_csn_with_max_diff)
  ) %>%
  dplyr::arrange(hour)

max_diff_hour_summary

# Plot differences
p_score_diff_density <- ggplot(score_diff_df, aes(x = diff)) +
  geom_density(fill = "gray80", color = "black", alpha = 0.6, na.rm = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Distribution of Local Score Minus Epic Auto Score",
    x = "Local score - Epic auto score",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

p_score_diff_density
save_aim1_plot(p_score_diff_density, "local_minus_epic_auto_score_difference_density.tiff")

# Absolute mean difference plot
p_abs_score_diff_density <- ggplot(score_diff_df, aes(x = abs_diff)) +
  geom_density(fill = "gray80", color = "black", alpha = 0.6, na.rm = TRUE) +
  labs(
    title = "Distribution of Absolute Differences Between Scores",
    x = "|Local score - Epic auto score|",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

p_abs_score_diff_density
save_aim1_plot(p_abs_score_diff_density, "absolute_score_difference_density.tiff")

# Pot of mean differences by hour
p_mean_diff_by_hour <- ggplot(diff_by_hour, aes(x = hour, y = mean_diff)) +
  geom_ribbon(
    aes(ymin = lower_95_mean_diff, ymax = upper_95_mean_diff),
    alpha = 0.2
  ) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = sort(unique(diff_by_hour$hour))) +
  labs(
    title = "Mean Difference Between Local and Epic Auto Scores by PICU Hour",
    x = "Hours since PICU admission",
    y = "Mean difference: local score - Epic auto score"
  ) +
  theme_minimal(base_size = 14)

p_mean_diff_by_hour
save_aim1_plot(p_mean_diff_by_hour, "mean_local_minus_epic_auto_score_difference_by_picu_hour.tiff")

# Plot boxplots of differences by hour
p_diff_boxplot_by_hour <- ggplot(score_diff_df, aes(x = factor(hour), y = diff)) +
  geom_boxplot(outlier.alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Distribution of Score Differences by PICU Hour",
    x = "Hours since PICU admission",
    y = "Local score - Epic auto score"
  ) +
  theme_minimal(base_size = 14)

p_diff_boxplot_by_hour
save_aim1_plot(p_diff_boxplot_by_hour, "local_minus_epic_auto_score_difference_boxplot_by_picu_hour.tiff")

# Bland altman plot
ba_limits <- score_diff_df %>%
  dplyr::summarise(
    mean_diff = mean(diff),
    sd_diff = sd(diff),
    lower_limit = mean_diff - 1.96 * sd_diff,
    upper_limit = mean_diff + 1.96 * sd_diff
  )

p_bland_altman_scores <- ggplot(score_diff_df, aes(x = mean_score, y = diff)) +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_hline(
    yintercept = ba_limits$mean_diff,
    linewidth = 1
  ) +
  geom_hline(
    yintercept = ba_limits$lower_limit,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = ba_limits$upper_limit,
    linetype = "dashed"
  ) +
  labs(
    title = "Bland-Altman Plot: Local Score vs Epic Auto Score",
    x = "Mean of local and Epic auto scores",
    y = "Local score - Epic auto score"
  ) +
  theme_minimal(base_size = 14)

p_bland_altman_scores
save_aim1_plot(p_bland_altman_scores, "bland_altman_local_vs_epic_auto_scores.tiff")

# Scatterplot
p_auto_vs_local_scatter <- ggplot(score_diff_df, aes(x = auto_score, y = local_score)) +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Epic Auto Score vs Local Model Score",
    x = "Epic auto score",
    y = "Local model score"
  ) +
  theme_minimal(base_size = 14)

p_auto_vs_local_scatter
save_aim1_plot(p_auto_vs_local_scatter, "epic_auto_vs_local_model_score_scatter.tiff")

# Single summary list of difference measures
score_comparison_results <- list(
  overall_diff_summary = overall_diff_summary,
  diff_by_hour = diff_by_hour,
  diff_by_csn = diff_by_csn,
  max_diff_per_csn = max_diff_per_csn,
  max_diff_hour_summary = max_diff_hour_summary,
  threshold_summary = threshold_summary
)

score_comparison_results


# Make hourly boxplots of differences
score_diff_hourly <- comp_scores %>%
  dplyr::mutate(
    hour_bin = cut(
      hours_since_icu_num,
      breaks = seq(0, 24, by = 1),
      include.lowest = TRUE,
      right = TRUE,
      labels = paste0(1:24)
    ),
    score_diff = local_score - auto_score,
    abs_score_diff = abs(score_diff)
  ) %>%
  dplyr::filter(
    !is.na(hour_bin),
    !is.na(auto_score),
    !is.na(local_score)
  )

p_hourly_score_diff_boxplot <- ggplot(score_diff_hourly, aes(x = hour_bin, y = score_diff)) +
  geom_boxplot(outlier.alpha = 0.15, outlier.size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Difference Between Local and Epic Auto Scores by PICU Hour",
    x = "PICU hour bin",
    y = "Local score - Epic auto score"
  ) +
  theme_minimal(base_size = 14)

p_hourly_score_diff_boxplot
save_aim1_plot(p_hourly_score_diff_boxplot, "hourly_local_minus_epic_auto_score_difference_boxplot.tiff")

# Make a table behind this plot for hourly differences
hourly_diff_summary <- score_diff_hourly %>%
  dplyr::group_by(hour_bin) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_csn = dplyr::n_distinct(csn),
    median_diff = median(score_diff, na.rm = TRUE),
    q1_diff = stats::quantile(score_diff, 0.25, na.rm = TRUE),
    q3_diff = stats::quantile(score_diff, 0.75, na.rm = TRUE),
    mean_abs_diff = mean(abs_score_diff, na.rm = TRUE),
    median_abs_diff = median(abs_score_diff, na.rm = TRUE),
    q95_abs_diff = stats::quantile(abs_score_diff, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

hourly_diff_summary


