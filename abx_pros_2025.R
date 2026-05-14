# Script to load in and tidy all antibiotics to determine antibiotic exposure and administration in the PICU as well as duration

abx_raw <- read.csv(file = "/phi/sbi/prospective_data/Prospective/data/antinfective_export_pros_100125_updated.csv")

abx_df <- abx_raw %>% mutate(order_inst = as.POSIXct(order_inst, tz = "UTC"), taken_time = as.POSIXct(taken_time, tz = "UTC"))
abx_df <- abx_df %>% mutate(order_inst = force_tz(order_inst, tzone = "America/Denver"),
                              taken_time = force_tz(taken_time, tzone = "America/Denver"))

# Fix classes and variable names
abx_df$PAT_MRN_ID <- as.character(abx_df$PAT_MRN_ID)
abx_df <- abx_df %>% rename(mrn = PAT_MRN_ID, hsp_account_id = HSP_ACCOUNT_ID)
colnames(abx_df) <- str_to_lower(colnames(abx_df))

# Convert pros_all and abx_df to data.tables
setDT(pros_all)
setDT(abx_df)

# Ensure that 'score_time' and 'taken_time' are in POSIXct datetime format
pros_all[, score_time := as.POSIXct(score_time)]
abx_df[, taken_time := as.POSIXct(taken_time)]

# Create 'start' and 'end' times for the 48-hour window in pros_all
pros_all[, start := score_time - 48*3600]  # 48 hours before score_time
pros_all[, end := score_time]

# Create 'start' and 'end' times in abx_df (both are 'taken_time' since it's a point event)
abx_df[, start := taken_time]
abx_df[, end := taken_time]

# Set keys on 'pat_enc_csn_id', 'start', and 'end' for both data.tables
setkey(pros_all, pat_enc_csn_id, start, end)
setkey(abx_df, pat_enc_csn_id, start, end)

# Add a row identifier to pros_all for easy indexing
pros_all[, rowid := .I]

# Perform a non-equi join to find overlaps where antibiotics were administered within the 48-hour window
abx_df$pat_enc_csn_id <- as.character(abx_df$pat_enc_csn_id) # change to character to allow comparison
overlaps <- foverlaps(abx_df, pros_all, type = "any", nomatch = 0L)

# Initialize 'abx_exp' column to 0
pros_all[, abx_exp := 0L]

# Identify rows in pros_all where overlaps occurred and set 'abx_exp' to 1
rows_with_abx <- unique(overlaps$rowid)
pros_all[rows_with_abx, abx_exp := 1L]

# Now abx_exposure information is captured.
# Add in info on which medication is an antibiotics
drug_class_table <- read.csv(file = "/phi/sbi/sbi_blake/aim_1_paper_materials/unique_antimicrobials_classified.csv")
drug_class_table <- drug_class_table %>% dplyr::select(-source_index)

# Join in info using medication_name
abx_df <- abx_df %>% left_join(drug_class_table, by = "medication_name")

########################### Now will move on to identifying duration of antibiotics ######################################################
# Define a function to calculate antibiotic duration after the score
calculate_abx_duration <- function(pat_id, score_time, picu_adm_time, abx_df) {
  # Step 1: Filter antibiotics for the given patient (pat_enc_csn_id) and only after PICU admission
  abx_filtered <- abx_df %>%
    filter(pat_enc_csn_id == pat_id & taken_time >= picu_adm_time) %>%
    arrange(taken_time) # Sort by time

  # Step 2: Check if any antibiotics are found for this patient after PICU admission
  if (nrow(abx_filtered) == 0) {
    return(0) # No antibiotics for this patient after PICU admission
  }

  # Step 3: Further filter antibiotics to only those given after the score time
  abx_after_score <- abx_filtered %>%
    filter(taken_time >= score_time)

  # Step 4: Check if any antibiotics are found after the score time
  if (nrow(abx_after_score) == 0) {
    return(0) # No antibiotics after the score time
  }

  # Step 5: Handle the case where there's only one dose after score time
  if (nrow(abx_after_score) == 1) {
    # If there's only one dose, consider the duration as 1 day (24 hours)
    return(1) # Single dose => duration is 1 day
  }

  # Step 6: Identify continuous antibiotics administration (within 26-hour gaps)
  start_time <- abx_after_score$taken_time[1] # Start with the first dose time
  last_time <- start_time

  # Loop through remaining antibiotic doses to find continuous administrations
  for (i in 2:nrow(abx_after_score)) {
    current_time <- abx_after_score$taken_time[i]

    # Check if the time difference between doses is less than 26 hours
    if (difftime(current_time, last_time, units = "hours") <= 26) {
      last_time <- current_time
    } else {
      break # If gap is greater than 26 hours, break the loop (end of continuous course)
    }
  }

  # Step 7: Calculate duration (in days) from the first dose to 24 hours after the last dose
  end_time <- last_time + hours(24) # Add 24 hours to the final dose time
  abx_duration <- as.numeric(difftime(end_time, start_time, units = "days"))

  return(abx_duration)
}

# Step 8: Apply this function to each row of the pros_all dataframe
pros_all <- pros_all %>%
  rowwise() %>%
  mutate(abx_duration_after_score = calculate_abx_duration(
    pat_enc_csn_id, score_time, picu_adm_date_time, abx_df
  )) %>%
  ungroup()


# # # Store resulting dataframes for future use
# write.fst(x = pros_all, path = "/phi/sbi/sbi_blake/pros_all_2025_validation_after_abx_5_8_25.fst")
# write.fst(x = abx_df, path = "/phi/sbi/sbi_blake/abx_df_2025_validation_after_abx_5_8_25.fst")


# Skip all prior code and read in the fst files for pros_all at this step and the abd_df dataframe
# pros_all <- read.fst(path = "/phi/sbi/sbi_blake/pros_all_2025_validation_after_abx_5_8_25.fst")
# abx_df <- read.fst(path = "/phi/sbi/sbi_blake/abx_df_2025_validation_after_abx_5_8_25.fst")




