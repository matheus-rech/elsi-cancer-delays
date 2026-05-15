#!/usr/bin/env Rscript
# 00_prepare_data.R
# Download ELSI DTA from GitHub Release → process → elsi_data.csv + variable_labels.csv
# Preserves Stata variable labels as a data dictionary
#
# Data source: ELSI-Brazil (English Stata 13 files)
# Release: https://github.com/matheus-rech/elsi-cancer-delays/releases/tag/v0.1.0-data

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
})

REPO    <- "matheus-rech/elsi-cancer-delays"
TAG     <- "v0.1.0-data"
DTA_URL <- paste0("https://github.com/", REPO, "/releases/download/", TAG, "/")
DATA_DIR <- "analysis/data"
APP_DIR  <- "app"

dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(APP_DIR,  showWarnings = FALSE, recursive = TRUE)

# ── 1. Download DTA files ─────────────────────────────────────────────────────
download_if_missing <- function(filename, destdir) {
  dest <- file.path(destdir, filename)
  if (file.exists(dest)) {
    cat(sprintf("[%s] Cached (%.1f MB)\n", filename, file.size(dest) / 1024^2))
    return(dest)
  }
  url <- paste0(DTA_URL, filename)
  cat(sprintf("[%s] Downloading from GitHub Release...\n", filename))
  download.file(url, dest, mode = "wb", quiet = FALSE)
  cat(sprintf("[%s] Done (%.1f MB)\n", filename, file.size(dest) / 1024^2))
  dest
}

wave2_path    <- download_if_missing("ELSI.English.2nd.wave.stata13.dta", DATA_DIR)
baseline_path <- download_if_missing("ELSI.English.baseline.stata13.dta", DATA_DIR)

# ── 2. Read DTA (preserving labels) ──────────────────────────────────────────
cat("\nReading Wave 2 DTA...\n")
d <- read_dta(wave2_path)
cat(sprintf("  %d rows x %d cols\n", nrow(d), ncol(d)))

# ── 3. Extract variable labels as data dictionary ─────────────────────────────
labels_df <- tibble(
  variable = names(d),
  label = sapply(names(d), function(v) {
    lbl <- attr(d[[v]], "label")
    if (is.null(lbl)) NA_character_ else as.character(lbl)
  }, USE.NAMES = FALSE),
  class = sapply(d, function(x) class(x)[1])
)
labels_path <- file.path(DATA_DIR, "variable_labels.csv")
write.csv(labels_df, labels_path, row.names = FALSE)
cat(sprintf("Variable dictionary: %d variables → %s\n", nrow(labels_df), labels_path))

# ── 4. Map English DTA names to analysis names ───────────────────────────────
# English DTA uses sex2/age2/calibrated_weight; questionnaire codes (n60, p70, etc.) are identical
# Also has sex/age as unnamed vars — use sex2/age2 (respondent-specific)

out <- d %>%
  transmute(
    # --- Demographics & SES ---
    sex           = as.integer(sex2),
    age           = as.numeric(age2),
    region        = as.integer(region),
    zona          = as.integer(area),        # 1=urban, 2=rural (English DTA: 'area')
    marital       = as.integer(e7),
    race          = as.integer(e9),
    education_yrs = as.numeric(e22),
    n_children    = as.numeric(e11),
    income_pc     = as.numeric(if ("rendadompc" %in% names(d)) rendadompc
                               else if ("percapita_income" %in% names(d)) percapita_income
                               else calibrated_weight * NA),  # fallback

    # Survey design
    weight        = as.numeric(calibrated_weight),
    estrato       = as.numeric(stratum),
    upa           = as.numeric(psu),

    # --- Self-rated health ---
    self_health   = as.integer(n1),

    # --- Cancer ---
    cancer_dx       = as.integer(n60),
    cancer_age_dx   = as.numeric(n60_1),
    cancer_tx_2y    = as.integer(n60_3),
    cancer_recur    = as.integer(n60_4),
    tx_chemo        = as.integer(n60_51),
    tx_surgery      = as.integer(n60_52),
    tx_radiation    = as.integer(n60_53),
    tx_symptom_med  = as.integer(n60_54),
    tx_other        = as.integer(n60_57),
    cancer_status   = as.integer(n60_7),

    # --- Chronic conditions ---
    hypertension    = as.integer(n28),
    diabetes        = as.integer(n35),
    cholesterol     = as.integer(n44),
    heart_attack    = as.integer(n46),
    angina          = as.integer(n48),
    heart_failure   = as.integer(n50),
    asthma          = as.integer(n54),
    copd            = as.integer(n55),
    arthritis       = as.integer(n56),
    osteoporosis    = as.integer(n57),
    back_problem    = as.integer(n58),
    depression_dx   = as.integer(n59),
    kidney_disease  = as.integer(n61),
    parkinson       = as.integer(n62),
    alzheimer       = as.integer(n63),
    cabg_stent      = as.integer(n66),

    # --- Falls ---
    fall_12m        = as.integer(n18),
    n_falls         = as.numeric(n19),
    hip_fracture    = as.integer(n21),

    # --- Fried frailty components ---
    weight_loss_3m  = as.integer(n69),
    exhaustion_1    = as.integer(n72),
    exhaustion_2    = as.integer(n73),
    grip_1          = as.numeric(mf27),
    grip_2          = as.numeric(mf28),
    grip_3          = as.numeric(mf29),
    gait_time_1     = as.numeric(mf35s),
    gait_time_2     = as.numeric(mf38s),
    balance_side    = as.numeric(mf30),
    balance_semi    = as.numeric(mf31),
    balance_tandem  = as.numeric(mf32),

    # --- Anthropometry & vitals ---
    height_cm       = as.numeric(mf13),
    weight_kg       = as.numeric(mf22),
    sbp             = as.numeric(mf9),
    dbp             = as.numeric(mf10),
    pulse           = as.numeric(mf10_1),
    waist_cm        = as.numeric(mf17m),
    hip_cm          = as.numeric(mf19m),

    # --- Social participation ---
    soc_letters     = as.integer(p70),
    soc_visit       = as.integer(p71),
    soc_invite      = as.integer(p72),
    soc_go_out      = as.integer(p73),
    soc_organized   = as.integer(p74),
    soc_civic       = as.integer(p75),
    soc_internet    = as.integer(p76),
    soc_drive       = as.integer(p77),
    soc_games       = as.integer(p78),
    soc_hobby       = as.integer(p79),
    soc_short_trip  = as.integer(p80),
    soc_long_trip   = as.integer(p81),
    soc_volunteer   = as.integer(p82),

    # --- Depression (CES-D 8) ---
    dep_depressed        = as.integer(r2),
    dep_everything_effort = as.integer(r3),
    dep_sleep_bad        = as.integer(r4),
    dep_happy            = as.integer(r5),
    dep_lonely           = as.integer(r6),
    dep_enjoyed_life     = as.integer(r7),
    dep_sad              = as.integer(r8),
    dep_could_not_go     = as.integer(r9),

    # --- Physical activity ---
    vigorous_activity = as.integer(l1),
    moderate_activity = as.integer(l2),
    light_activity    = as.integer(l3),

    # --- Cognitive ---
    memory_self  = as.integer(q3),
    memory_trend = as.integer(q4),

    # --- Nagi items ---
    nagi_jog         = as.integer(p5),
    nagi_walk_1km    = as.integer(p6),
    nagi_walk_100m   = as.integer(p7),
    nagi_stairs_many = as.integer(p8),
    nagi_stairs_one  = as.integer(p9),
    nagi_sit_2h      = as.integer(p10),
    nagi_stand_up    = as.integer(p10_1),
    nagi_bend        = as.integer(p12),
    nagi_arms_up     = as.integer(p13),
    nagi_push        = as.integer(p14),
    nagi_lift_5kg    = as.integer(p15),
    nagi_pick_coin   = as.integer(p16),

    # --- ADLs ---
    adl_cross_room  = as.integer(p37),
    adl_dress       = as.integer(p40),
    adl_bath        = as.integer(p43),
    adl_eat         = as.integer(p46),
    adl_bed         = as.integer(p49),
    adl_toilet      = as.integer(p55),

    # --- IADLs ---
    iadl_hygiene    = as.integer(p17),
    iadl_meals      = as.integer(p20),
    iadl_money      = as.integer(p22),
    iadl_transport  = as.integer(p24),
    iadl_shopping   = as.integer(p26),
    iadl_phone      = as.integer(p28),
    iadl_meds       = as.integer(p30),
    iadl_light_hw   = as.integer(p33),
    iadl_heavy_hw   = as.integer(p35),

    # --- Incontinence ---
    incontinence    = as.integer(p58),

    # --- Healthcare utilization ---
    has_health_plan    = as.integer(u1),
    doctor_visits      = as.numeric(u6),
    specialist_visits  = as.numeric(u8),
    has_regular_service = as.integer(u10),
    easy_appointment   = as.integer(u12),
    appt_24h           = as.integer(u13),
    hospitalized       = as.integer(u42),
    n_hospitalizations = as.numeric(u43),
    hospital_days      = as.numeric(u44),
    emergency_care     = as.integer(u49),

    # --- Lifestyle ---
    ever_smoked     = as.integer(l30_0),
    current_smoker  = as.integer(l30),
    cigs_per_day    = as.numeric(l30_1),
    alcohol_freq    = as.integer(l24),
    alcohol_days_wk = as.numeric(l25),
    alcohol_doses   = as.numeric(l26),
    binge_male      = as.integer(l28),
    binge_female    = as.integer(l29),
    veg_days_wk     = as.numeric(l15),

    # --- Sleep & other ---
    sleep_quality   = as.integer(n74),
    sleep_meds      = as.integer(n75),
    flu_vaccine     = as.integer(n67),
    phys_health_bad_days  = as.numeric(n2),
    mental_health_bad_days = as.numeric(n3),
    bedridden       = as.integer(p1),
    wheelchair      = as.integer(p2),
    uses_walking_aid = as.integer(p3)
  )

cat("Selected:", nrow(out), "rows x", ncol(out), "cols\n")

# ── 5. Clean missing value codes ─────────────────────────────────────────────
clean_89 <- function(x) ifelse(x %in% c(8, 9, 88, 99, 888, 999), NA_real_, as.numeric(x))

vars_to_clean <- c(
  "cancer_age_dx","cancer_tx_2y","cancer_recur","cancer_status",
  "tx_chemo","tx_surgery","tx_radiation","tx_symptom_med","tx_other",
  "hypertension","diabetes","cholesterol","heart_attack","angina",
  "heart_failure","asthma","copd","arthritis","osteoporosis","back_problem",
  "depression_dx","kidney_disease","parkinson","alzheimer","cabg_stent",
  "fall_12m","n_falls","hip_fracture",
  "weight_loss_3m","exhaustion_1","exhaustion_2",
  "self_health","memory_self","memory_trend",
  "has_health_plan","hospitalized","n_hospitalizations","emergency_care",
  "easy_appointment","appt_24h","has_regular_service",
  "ever_smoked","current_smoker","alcohol_freq",
  "binge_male","binge_female","sleep_quality","sleep_meds","flu_vaccine",
  "bedridden","wheelchair","uses_walking_aid","incontinence",
  "education_yrs","race","marital"
)
for (v in vars_to_clean) {
  if (v %in% names(out)) out[[v]] <- clean_89(out[[v]])
}

soc_binary <- c("soc_letters","soc_visit","soc_invite","soc_go_out",
                "soc_organized","soc_civic","soc_internet","soc_drive",
                "soc_games","soc_hobby","soc_short_trip","soc_long_trip",
                "soc_volunteer")
for (v in soc_binary) out[[v]] <- ifelse(out[[v]] %in% c(8, 9), NA_integer_, out[[v]])

func_vars <- c(grep("^nagi_", names(out), value = TRUE),
               grep("^adl_",  names(out), value = TRUE),
               grep("^iadl_", names(out), value = TRUE))
for (v in func_vars) out[[v]] <- ifelse(out[[v]] > 3 | out[[v]] < 0, NA_integer_, as.integer(out[[v]]))

dep_neg <- c("dep_depressed","dep_everything_effort","dep_sleep_bad",
             "dep_lonely","dep_sad","dep_could_not_go")
dep_pos <- c("dep_happy","dep_enjoyed_life")
for (v in c(dep_neg, dep_pos)) out[[v]] <- ifelse(out[[v]] %in% c(8, 9), NA_integer_, out[[v]])

for (v in c("vigorous_activity","moderate_activity","light_activity")) {
  out[[v]] <- ifelse(out[[v]] %in% c(8, 9), NA_integer_, out[[v]])
}

# ── 6. Derive labels ─────────────────────────────────────────────────────────
recode_yn <- function(x) case_when(x == 0 ~ "No", x == 1 ~ "Yes", TRUE ~ NA_character_)

out <- out %>% mutate(
  sex_label = case_when(sex == 0 ~ "Female", sex == 1 ~ "Male"),
  region_label = case_when(
    region == 1 ~ "North", region == 2 ~ "Northeast",
    region == 3 ~ "Southeast", region == 4 ~ "South",
    region == 5 ~ "Central-West"),
  zona_label = case_when(zona == 1 ~ "Urban", zona == 2 ~ "Rural"),
  marital_label = case_when(
    marital == 1 ~ "Single", marital == 2 ~ "Married/Partner",
    marital == 3 ~ "Divorced/Separated", marital == 4 ~ "Widowed"),
  race_label = case_when(
    race == 1 ~ "White", race == 2 ~ "Black", race == 3 ~ "Mixed/Pardo",
    race == 4 ~ "Asian", race == 5 ~ "Indigenous"),
  age_group = case_when(
    age < 60 ~ "50-59", age < 70 ~ "60-69", age < 80 ~ "70-79", age >= 80 ~ "80+"),
  cancer_label = recode_yn(cancer_dx),
  self_health_label = case_when(
    self_health == 1 ~ "Very good", self_health == 2 ~ "Good",
    self_health == 3 ~ "Fair", self_health == 4 ~ "Poor",
    self_health == 5 ~ "Very poor"),
  cancer_status_label = case_when(
    cancer_status == 1 ~ "In treatment", cancer_status == 2 ~ "Controlled/remission",
    cancer_status == 3 ~ "Getting worse"),
  education_cat = case_when(
    is.na(education_yrs) ~ NA_character_,
    education_yrs == 0 ~ "None",
    education_yrs <= 4 ~ "1-4 years",
    education_yrs <= 8 ~ "5-8 years",
    education_yrs <= 11 ~ "9-11 years",
    TRUE ~ "12+ years"),
  income_quartile = case_when(
    is.na(income_pc) ~ NA_character_,
    income_pc <= quantile(income_pc, 0.25, na.rm = TRUE) ~ "Q1 (lowest)",
    income_pc <= quantile(income_pc, 0.50, na.rm = TRUE) ~ "Q2",
    income_pc <= quantile(income_pc, 0.75, na.rm = TRUE) ~ "Q3",
    TRUE ~ "Q4 (highest)"),
  hypertension_label  = recode_yn(hypertension),
  diabetes_label      = recode_yn(diabetes),
  cholesterol_label   = recode_yn(cholesterol),
  heart_failure_label = recode_yn(heart_failure),
  arthritis_label     = recode_yn(arthritis),
  smoking_label = case_when(
    current_smoker == 1 ~ "Current", ever_smoked == 1 ~ "Former", TRUE ~ "Never")
)

# ── 7. Composite scores ──────────────────────────────────────────────────────
out$social_score <- rowSums(out[, soc_binary] == 1, na.rm = TRUE)

out$depression_score <- rowSums(
  cbind(out[, dep_neg], 4 - out[, dep_pos]), na.rm = FALSE)

chronic_cols <- c("hypertension","diabetes","cholesterol","heart_attack","angina",
                  "heart_failure","asthma","copd","arthritis","osteoporosis",
                  "back_problem","depression_dx","kidney_disease","cancer_dx")
out$n_chronic <- rowSums(out[, chronic_cols] == 1, na.rm = TRUE)

adl_cols  <- grep("^adl_",  names(out), value = TRUE)
iadl_cols <- grep("^iadl_", names(out), value = TRUE)
nagi_cols <- grep("^nagi_", names(out), value = TRUE)
out$adl_disability    <- rowSums(out[, adl_cols]  >= 1, na.rm = TRUE)
out$iadl_disability   <- rowSums(out[, iadl_cols] >= 1, na.rm = TRUE)
out$nagi_limitations  <- rowSums(out[, nagi_cols] >= 1, na.rm = TRUE)

out$grip_strength <- pmax(out$grip_1, out$grip_2, out$grip_3, na.rm = TRUE)
out$grip_strength <- ifelse(out$grip_strength > 90 | out$grip_strength <= 0, NA_real_, out$grip_strength)

out$gait_time <- pmin(out$gait_time_1, out$gait_time_2, na.rm = TRUE)
out$gait_time <- ifelse(out$gait_time <= 0 | out$gait_time > 60, NA_real_, out$gait_time)
out$gait_speed <- ifelse(!is.na(out$gait_time) & out$gait_time > 0, 3 / out$gait_time, NA_real_)

out$bmi <- ifelse(
  !is.na(out$height_cm) & out$height_cm > 100 & out$height_cm < 220 &
  !is.na(out$weight_kg) & out$weight_kg > 25 & out$weight_kg < 200,
  out$weight_kg / (out$height_cm / 100)^2, NA_real_)
out$bmi_cat <- case_when(
  is.na(out$bmi) ~ NA_character_,
  out$bmi < 18.5 ~ "Underweight", out$bmi < 25 ~ "Normal",
  out$bmi < 30 ~ "Overweight", TRUE ~ "Obese")

out$whr <- ifelse(
  !is.na(out$waist_cm) & out$waist_cm > 40 & !is.na(out$hip_cm) & out$hip_cm > 40,
  out$waist_cm / out$hip_cm, NA_real_)

# ── 8. Fried Frailty Phenotype ────────────────────────────────────────────────
out$fried_wt <- as.integer(out$weight_loss_3m == 1)
out$fried_exhaust <- as.integer(
  (!is.na(out$exhaustion_1) & out$exhaustion_1 >= 3) |
  (!is.na(out$exhaustion_2) & out$exhaustion_2 >= 3))
out$fried_activity <- as.integer(
  (!is.na(out$vigorous_activity) & out$vigorous_activity >= 4) &
  (!is.na(out$moderate_activity) & out$moderate_activity >= 4) &
  (!is.na(out$light_activity)    & out$light_activity >= 4))

out <- out %>%
  group_by(sex) %>%
  mutate(
    gait_p80 = quantile(gait_time, 0.80, na.rm = TRUE),
    fried_gait = as.integer(!is.na(gait_time) & gait_time >= gait_p80)) %>%
  ungroup() %>% select(-gait_p80)

out <- out %>%
  group_by(sex) %>%
  mutate(
    grip_p20 = quantile(grip_strength, 0.20, na.rm = TRUE),
    fried_grip = as.integer(!is.na(grip_strength) & grip_strength <= grip_p20)) %>%
  ungroup() %>% select(-grip_p20)

fried_cols <- c("fried_wt","fried_exhaust","fried_activity","fried_gait","fried_grip")
out$fried_score <- rowSums(out[, fried_cols], na.rm = TRUE)
fried_available <- rowSums(!is.na(out[, fried_cols]))
out$fried_score <- ifelse(fried_available >= 3, out$fried_score, NA_integer_)

out$frailty_status <- case_when(
  is.na(out$fried_score) ~ NA_character_,
  out$fried_score == 0 ~ "Robust",
  out$fried_score <= 2 ~ "Pre-frail",
  TRUE ~ "Frail")

# ── 9. Deficit Accumulation Frailty Index ─────────────────────────────────────
deficit_items <- out %>% transmute(
  d01 = as.integer(hypertension == 1), d02 = as.integer(diabetes == 1),
  d03 = as.integer(cholesterol == 1),  d04 = as.integer(heart_attack == 1),
  d05 = as.integer(angina == 1),       d06 = as.integer(heart_failure == 1),
  d07 = as.integer(asthma == 1),       d08 = as.integer(copd == 1),
  d09 = as.integer(arthritis == 1),    d10 = as.integer(osteoporosis == 1),
  d11 = as.integer(back_problem == 1), d12 = as.integer(depression_dx == 1),
  d13 = as.integer(kidney_disease == 1), d14 = as.integer(cancer_dx == 1),
  d15 = as.integer(self_health >= 4),
  d16 = as.integer(!is.na(depression_score) & depression_score >= 12),
  d17 = as.integer(adl_disability >= 1), d18 = as.integer(iadl_disability >= 1),
  d19 = as.integer(fall_12m == 1),     d20 = as.integer(hospitalized == 1),
  d21 = as.integer(weight_loss_3m == 1), d22 = as.integer(incontinence == 1),
  d23 = as.integer(!is.na(memory_self) & memory_self >= 4),
  d24 = as.integer(fried_activity == 1),
  d25 = as.integer(!is.na(bmi) & (bmi < 18.5 | bmi >= 30)),
  d26 = as.integer(!is.na(grip_strength) & fried_grip == 1),
  d27 = as.integer(!is.na(gait_speed) & fried_gait == 1),
  d28 = as.integer(nagi_walk_100m >= 1), d29 = as.integer(nagi_stairs_one >= 1),
  d30 = as.integer(nagi_bend >= 1),    d31 = as.integer(nagi_lift_5kg >= 1),
  d32 = as.integer(!is.na(sleep_quality) & sleep_quality >= 3),
  d33 = as.integer(current_smoker == 1), d34 = as.integer(bedridden == 1),
  d35 = as.integer(uses_walking_aid == 1))

n_deficits  <- rowSums(deficit_items, na.rm = TRUE)
n_available <- rowSums(!is.na(deficit_items))
out$frailty_index <- ifelse(n_available >= 20, n_deficits / n_available, NA_real_)

# ── 10. QoL proxy ────────────────────────────────────────────────────────────
out$qol_proxy <- with(out, {
  sh   <- ifelse(!is.na(self_health) & self_health <= 5, (5 - self_health) / 4 * 100, NA_real_)
  dep  <- ifelse(!is.na(depression_score), (1 - depression_score / 24) * 100, NA_real_)
  soc  <- social_score / 13 * 100
  nagi <- ifelse(!is.na(nagi_limitations), (1 - nagi_limitations / 12) * 100, NA_real_)
  rowMeans(cbind(sh, dep, soc, nagi), na.rm = TRUE)
})

# ── 11. Final derived variables ───────────────────────────────────────────────
out$late_diagnosis <- case_when(
  is.na(out$cancer_age_dx) | out$cancer_age_dx > 200 ~ NA_character_,
  out$cancer_age_dx >= 70 ~ "Late (\u226570)",
  out$cancer_age_dx >= 50 ~ "Middle (50-69)",
  TRUE ~ "Early (<50)")

out$binge_drinking <- case_when(
  out$sex == 1 & out$binge_male == 1 ~ 1L, out$sex == 0 & out$binge_female == 1 ~ 1L,
  out$sex == 1 & out$binge_male == 0 ~ 0L, out$sex == 0 & out$binge_female == 0 ~ 0L,
  TRUE ~ NA_integer_)

out$multimorbidity <- case_when(
  out$n_chronic >= 3 ~ "3+", out$n_chronic == 2 ~ "2",
  out$n_chronic == 1 ~ "1", out$n_chronic == 0 ~ "0")

# ── 12. Write output ─────────────────────────────────────────────────────────
csv_path <- file.path(APP_DIR, "elsi_data.csv")
write.csv(out, csv_path, row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("\nWritten: %s (%.1f MB)\n", csv_path, file.size(csv_path) / 1024^2))
cat(sprintf("Dimensions: %d rows x %d cols\n", nrow(out), ncol(out)))

# Summary
cat("\n=== ELSI DATA SUMMARY ===\n")
cat(sprintf("Cancer prevalence: %.1f%%\n", 100 * mean(out$cancer_dx == 1, na.rm = TRUE)))
cat(sprintf("Mean age: %.1f\n", mean(out$age, na.rm = TRUE)))
cat(sprintf("Female: %.1f%%\n", 100 * mean(out$sex == 0, na.rm = TRUE)))
cat(sprintf("Frail: %.1f%%\n", 100 * mean(out$frailty_status == "Frail", na.rm = TRUE)))
cat(sprintf("Mean FI: %.3f\n", mean(out$frailty_index, na.rm = TRUE)))
cat(sprintf("Hospitalized 12m: %.1f%%\n", 100 * mean(out$hospitalized == 1, na.rm = TRUE)))
