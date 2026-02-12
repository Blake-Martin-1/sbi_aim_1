# Proposed maintenance tasks

## 1) Typo fix task
**Issue found:** The top-level description in `aim_1_main.R` says "put availability" instead of "input availability".

**Proposed task:**
- Update the opening comment to read "input availability" so the script description is accurate and professional.

**Why it matters:**
- This is a visible typo in the primary analysis script header and can create confusion for readers.

## 2) Bug fix task
**Issue found:** In `setup_aim_1.R`, `user_norms_path` is only set in the `else` branch, but later file paths always use `getOption("user_norms_path")`.

**Proposed task:**
- Set `user_norms_path` in both user branches (including `tdbennett`) or add a guarded fallback.
- Add an early assertion that required options (`user_data_path`, `user_norms_path`) are non-NULL before constructing file paths.

**Why it matters:**
- For users in the `tdbennett` branch, norm-file paths may be built from `NULL`, leading to file read failures.

## 3) Comment/documentation discrepancy task
**Issue found:** The header in `setup_aim_1.R` says `file: setup.R`, but the actual filename is `setup_aim_1.R`.

**Proposed task:**
- Update the header metadata to match the real filename.
- Optionally add a brief note that this setup script is specific to AIM 1.

**Why it matters:**
- Mismatched header metadata makes navigation and maintenance harder, especially when debugging source chains.

## 4) Test improvement task
**Issue found:** There are no automated checks around predictor schema alignment before scoring prospective data in `retro_model_to_pros_data.R`.

**Proposed task:**
- Add a small `testthat` test (or validation helper) that verifies:
  - all required predictor columns exist before `predict()`,
  - renamed columns map correctly for both no-abx and yes-abx paths,
  - no unexpected extra/missing columns after alignment.

**Why it matters:**
- This script currently relies on manual renaming and assumptions; schema drift will fail late and be hard to diagnose.
