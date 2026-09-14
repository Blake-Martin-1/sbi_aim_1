# Utilities for summarizing prospective SBI types at the PICU-encounter level.

prospective_sbi_source_lookup <- tibble::tribble(
  ~specimen_source, ~broad_category,
  "Abdomen", "abdomen",
  "Abdominal Wall", "abdomen",
  "Abscess", "other",
  "Ankle", "msk",
  "Anus", "gastrointestinal",
  "Arm", "msk",
  "Axilla", "msk",
  "Back", "msk",
  "Biliary Tract", "gastrointestinal",
  "Blood, Arterial", "blood",
  "Blood, Capillary", "blood",
  "Blood, Line", "blood",
  "Blood, Unknown", "blood",
  "Blood, Venous", "blood",
  "Brain", "cns",
  "Bronchus", "respiratory",
  "Buccal", "ent",
  "Buttock", "msk",
  "Catheter", "unknown",
  "Catheter, Arterial Line", "blood",
  "Catheter, Direct Intracardiac", "blood",
  "Catheter, HD/Apheresis", "blood",
  "Catheter, Non-tunneled", "blood",
  "Catheter, PICC", "blood",
  "Catheter, Port", "blood",
  "Catheter, Tunneled", "blood",
  "Cellulitis", "msk",
  "Cervix", "gu",
  "Chest", "msk",
  "Colon", "gastrointestinal",
  "Colon, Rectum", "gastrointestinal",
  "Colon, Sigmoid", "gastrointestinal",
  "Cornea", "other",
  "Craniotomy", "cns",
  "Dialysis Effluent", "blood",
  "Duodenum", "gastrointestinal",
  "Ear", "other",
  "Eye", "other",
  "Face", "msk",
  "Femur", "msk",
  "Finger", "msk",
  "Foot", "msk",
  "Forearm", "msk",
  "G-tube", "gastrointestinal",
  "Gallbladder", "gastrointestinal",
  "Gastric", "gastrointestinal",
  "Groin", "msk",
  "Hand", "msk",
  "Head", "cns",
  "Heart", "mediastinum",
  "Hip", "msk",
  "JP Drain", "unknown",
  "Knee", "msk",
  "Leg", "msk",
  "Lip", "ent",
  "Liver", "gastrointestinal",
  "Lumbar Puncture", "cns",
  "Lung", "respiratory",
  "Mandible", "msk",
  "Maxilla", "msk",
  "Mediastinum", "mediastinum",
  "Mouth", "ent",
  "Multiple Sources - Specify in Comments", "unknown",
  "Naris", "respiratory",
  "Nasopharyngeal Swab", "respiratory",
  "Nasopharynx", "respiratory",
  "Neck", "msk",
  "Nose", "respiratory",
  "NULL", "unknown",
  "Oral Lesion", "ent",
  "Ostomy", "gastrointestinal",
  "Other - Specify in Comments", "unknown",
  "Pelvis", "gu",
  "Penis", "gu",
  "Pericardium", "mediastinum",
  "Perineum", "msk",
  "Peripheral, PIV Start", "blood",
  "Peripheral, Venipuncture", "blood",
  "Peritoneum", "gastrointestinal",
  "Pleura", "respiratory",
  "Pleural Space", "respiratory",
  "Scalp", "msk",
  "Sella", "cns",
  "Shoulder", "msk",
  "Shunt", "cns",
  "Sinus", "cns",
  "Skull", "cns",
  "Small Intestine", "gastrointestinal",
  "Stomach, Body and Antrum", "gastrointestinal",
  "Surface Lesion", "unknown",
  "Testis", "gu",
  "Thigh", "msk",
  "Throat", "ent",
  "Tibia", "msk",
  "Toe", "msk",
  "Tongue", "ent",
  "Trachea", "respiratory",
  "Urine, 1st Void", "gu",
  "Urine, Bagged", "gu",
  "Urine, Clean Catch", "gu",
  "Urine, Cotton Ball", "gu",
  "Urine, Indwelling Catheter", "gu",
  "Urine, Nephrostomy", "gu",
  "Urine, Open Bladder", "gu",
  "Urine, Straight Catheter", "gu",
  "Urine, Suprapubic", "gu",
  "Vagina", "gu",
  "Vulva", "gu",
  "Wrist", "msk"
)

normalize_specimen_source <- function(x) {
  x |>
    stringr::str_replace_all("\u00a0", " ") |>
    stringr::str_to_lower() |>
    # The prospective microbiology extract uses snake_case values (for
    # example, peripheral_venipuncture), while the source dictionary was
    # supplied as display labels (Peripheral, Venipuncture). Normalize both
    # forms to the same punctuation-independent key before joining.
    stringr::str_replace_all("[^[:alnum:]]+", "_") |>
    stringr::str_replace_all("^_+|_+$", "")
}

prospective_sbi_source_lookup <- prospective_sbi_source_lookup |>
  dplyr::mutate(specimen_source_key = normalize_specimen_source(specimen_source)) |>
  dplyr::select(specimen_source_key, broad_category)

is_positive_sbi_flag <- function(x) {
  tolower(trimws(as.character(x))) %in% c("1", "true", "yes", "y")
}

# Each study_id contributes at most once to a category, even when multiple
# organisms/specimens from that category were recorded or when it appears in
# both antibiotic-exposure cohorts. A patient may contribute to more than one
# category. Antibiotic exposure is intentionally not part of this summary.
summarize_prospective_sbi_types <- function(
    sbi_micro,
    abx_unexposed,
    abx_exposed,
    window_hours = 24) {
  required_micro <- c("mrn", "specimen_source", "time_obtained")
  required_cohort <- c(
    "study_id", "mrn", "picu_adm_date_time", "ever_cx_neg_sepsis",
    "pna_1_0"
  )

  if (!all(required_micro %in% names(sbi_micro))) {
    stop("sbi_micro is missing: ", paste(setdiff(required_micro, names(sbi_micro)), collapse = ", "))
  }
  if (!all(required_cohort %in% names(abx_unexposed)) ||
      !all(required_cohort %in% names(abx_exposed))) {
    stop("Both cohort data frames must contain study_id, mrn, and picu_adm_date_time")
  }

  encounters <- dplyr::bind_rows(abx_unexposed, abx_exposed) |>
    dplyr::select(dplyr::all_of(required_cohort)) |>
    dplyr::mutate(mrn = as.character(mrn)) |>
    dplyr::distinct()

  total_patients <- dplyr::n_distinct(encounters$study_id, na.rm = TRUE)

  microbiology_sbi_types <- sbi_micro |>
    dplyr::transmute(
      mrn = as.character(mrn),
      specimen_source = as.character(specimen_source),
      specimen_source_key = normalize_specimen_source(specimen_source),
      time_obtained
    ) |>
    dplyr::filter(!is.na(mrn), !is.na(time_obtained)) |>
    dplyr::inner_join(encounters, by = "mrn", relationship = "many-to-many") |>
    dplyr::filter(
      time_obtained >= picu_adm_date_time - lubridate::hours(window_hours),
      time_obtained <= picu_adm_date_time + lubridate::hours(window_hours)
    ) |>
    dplyr::left_join(prospective_sbi_source_lookup, by = "specimen_source_key") |>
    dplyr::mutate(
      broad_category = dplyr::coalesce(broad_category, "unknown"),
      specimen_source = dplyr::coalesce(specimen_source, "NULL")
    ) |>
    dplyr::distinct(study_id, broad_category, specimen_source)

  non_microbiology_sbi_types <- dplyr::bind_rows(
    encounters |>
      dplyr::filter(is_positive_sbi_flag(ever_cx_neg_sepsis)) |>
      dplyr::transmute(
        study_id,
        broad_category = "culture_negative_sepsis",
        specimen_source = NA_character_
      ),
    encounters |>
      dplyr::filter(is_positive_sbi_flag(pna_1_0)) |>
      dplyr::transmute(
        study_id,
        broad_category = "bacterial_pneumonia",
        specimen_source = "VPS pna_1_0"
      )
  )

  encounter_sbi_types <- dplyr::bind_rows(
    microbiology_sbi_types,
    non_microbiology_sbi_types
  ) |>
    dplyr::distinct(study_id, broad_category, specimen_source)

  all_categories <- sort(unique(c(
    prospective_sbi_source_lookup$broad_category,
    "culture_negative_sepsis", "bacterial_pneumonia"
  )))
  sbi_type_summary <- encounter_sbi_types |>
    dplyr::distinct(study_id, broad_category) |>
    dplyr::count(broad_category, name = "patients_with_sbi") |>
    tidyr::complete(
      broad_category = all_categories,
      fill = list(patients_with_sbi = 0L)
    ) |>
    dplyr::mutate(
      total_patients = total_patients,
      proportion_of_all_patients = patients_with_sbi / total_patients
    ) |>
    dplyr::arrange(broad_category)

  list(
    encounter_sbi_types = encounter_sbi_types,
    sbi_type_summary = sbi_type_summary
  )
}
