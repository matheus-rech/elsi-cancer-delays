# 03_elsi_regional.R
# Compute ELSI-Brazil Wave 2 regional summary statistics for ecological linkage
# Output: analysis/data/elsi_regional.rds

library(dplyr)
library(readr)
library(srvyr)

# ── 1. Load data ──────────────────────────────────────────────────────────────
csv_path <- "/Users/matheusrech/Downloads/fingertips-br/files/elsi_native_addons/app/elsi_data.csv"

cat("Reading ELSI data...\n")
elsi <- read_csv(csv_path, show_col_types = FALSE)
cat(sprintf("Loaded: %d rows × %d cols\n", nrow(elsi), ncol(elsi)))

# ── 2. Create binary indicators ───────────────────────────────────────────────
elsi <- elsi |>
  mutate(
    frail_binary = as.integer(frailty_status == "Frail"),
    dep_binary   = as.integer(!is.na(depression_score) & depression_score >= 12),
    adl_any      = as.integer(!is.na(adl_disability)   & adl_disability >= 1),
    soc_low      = as.integer(!is.na(social_score)     & social_score <= 1),
    multi3       = as.integer(!is.na(n_chronic)        & n_chronic >= 3),
    # ensure numeric survey vars
    cancer_dx_n  = as.integer(cancer_dx == 1),
    hosp_n       = as.integer(hospitalized == 1),
    plan_n       = as.integer(has_health_plan == 1)
  )

cat("Binary indicators created.\n")

# ── 3. Survey design ──────────────────────────────────────────────────────────
# Drop rows with missing design variables
design_vars <- c("upa", "estrato", "weight", "region_label")
elsi_clean <- elsi |>
  filter(
    !is.na(upa), !is.na(estrato), !is.na(weight),
    !is.na(region_label), weight > 0
  )
cat(sprintf("Rows after dropping missing design vars: %d\n", nrow(elsi_clean)))

options(survey.lonely.psu = "adjust")

svy <- elsi_clean |>
  as_survey_design(
    ids    = upa,
    strata = estrato,
    weights = weight,
    nest   = TRUE
  )

cat("Survey design set up.\n")

# ── 4. Regional summary ───────────────────────────────────────────────────────
cat("Computing regional survey-weighted estimates...\n")

regional <- svy |>
  group_by(region_label) |>
  summarise(
    # Prevalences (proportions)
    cancer_prev   = survey_mean(cancer_dx_n,  na.rm = TRUE, vartype = "ci"),
    frail_prev    = survey_mean(frail_binary, na.rm = TRUE, vartype = "ci"),
    dep_prev      = survey_mean(dep_binary,   na.rm = TRUE, vartype = "ci"),
    adl_prev      = survey_mean(adl_any,      na.rm = TRUE, vartype = "ci"),
    soc_low_prev  = survey_mean(soc_low,      na.rm = TRUE, vartype = "ci"),
    hosp_prev     = survey_mean(hosp_n,       na.rm = TRUE, vartype = "ci"),
    plan_prev     = survey_mean(plan_n,       na.rm = TRUE, vartype = "ci"),
    # Continuous means
    mean_fi       = survey_mean(frailty_index,  na.rm = TRUE, vartype = "ci"),
    mean_social   = survey_mean(social_score,   na.rm = TRUE, vartype = "ci"),
    mean_income   = survey_mean(income_pc,      na.rm = TRUE, vartype = "ci"),
    mean_age      = survey_mean(age,            na.rm = TRUE, vartype = "ci"),
    mean_educ     = survey_mean(education_yrs,  na.rm = TRUE, vartype = "ci"),
    # Sample n (unweighted)
    n = unweighted(n())
  )

cat("Done.\n")

# ── 5. Save ───────────────────────────────────────────────────────────────────
out_dir <- "/Users/matheusrech/Downloads/fingertips-br/files/elsi_native_addons/analysis/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, "elsi_regional.rds")
saveRDS(regional, out_path)
cat(sprintf("Saved: %s\n", out_path))

# ── 6. Print table ────────────────────────────────────────────────────────────
cat("\n── ELSI Regional Summary (5 regions) ──\n")
print(regional, n = Inf, width = 120)

# Compact view of key prevalences
cat("\n── Key prevalences by region ──\n")
regional |>
  select(region_label, n,
         cancer_prev, frail_prev, dep_prev, adl_prev,
         hosp_prev, plan_prev, mean_fi, mean_income) |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  print(n = Inf, width = 120)
