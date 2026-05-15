## 06_table1.R
## Table 1 — Cohort characteristics by Lei 12.732 compliance
## Input:  analysis/data/apac_cohort.rds
## Output: analysis/tables/table1_characteristics.html

suppressPackageStartupMessages({
  library(dplyr)
  library(gtsummary)
  library(gt)
})

IN_FILE  <- "analysis/data/apac_cohort.rds"
OUT_HTML <- "analysis/tables/table1_characteristics.html"

# ── 1. Read cohort ────────────────────────────────────────────────────────────
message("Reading cohort …")
cohort <- readRDS(IN_FILE)
message(sprintf("  %s rows × %d cols", format(nrow(cohort), big.mark = ","), ncol(cohort)))

# ── 2. Derive display variables ───────────────────────────────────────────────
message("Deriving display variables …")

df <- cohort %>%
  mutate(
    # 1. Compliance variable
    compliance = factor(
      ifelse(compliant_60d, "\u226460 days", ">60 days"),
      levels = c("\u226460 days", ">60 days")
    ),

    # 2. Race label (SUS APAC codes: "01"=Branca, "02"=Preta, "03"=Parda)
    race_label = case_when(
      race == "01" ~ "White",
      race == "02" ~ "Black",
      race == "03" ~ "Mixed",
      TRUE         ~ "Other/Unknown"
    ),

    # Ensure factors with sensible ordering
    sex_label    = factor(sex_label,    levels = c("Female", "Male")),
    race_label   = factor(race_label,   levels = c("White", "Black", "Mixed", "Other/Unknown")),
    age_grp      = factor(age_grp,      levels = c("50-59", "60-69", "70-79", "80+")),
    cancer_group = factor(cancer_group, levels = c("Breast", "Prostate", "Colorectal",
                                                    "Lung", "Stomach", "Cervical",
                                                    "Leukemia", "Lymphoma", "Other")),
    stage_label  = factor(stage_label,  levels = c("In situ", "I", "II", "III", "IV")),
    region       = factor(region,       levels = c("North", "Northeast", "Southeast",
                                                    "South", "Central-West")),
    # death is already 0/1 integer in cohort; keep as logical for yes/no display
    death        = as.logical(death)
  ) %>%
  select(
    compliance,
    age,
    sex_label,
    race_label,
    age_grp,
    cancer_group,
    stage_label,
    region,
    ttt,
    death
  )

# ── 3. Rename columns to display names ────────────────────────────────────────
df <- df %>%
  rename(
    `Age (years)`              = age,
    `Sex`                      = sex_label,
    `Race/ethnicity`           = race_label,
    `Age group`                = age_grp,
    `Cancer type`              = cancer_group,
    `Stage`                    = stage_label,
    `Region`                   = region,
    `Time to treatment (days)` = ttt,
    `In-hospital death`        = death
  )

message(sprintf("  Analysis dataset: %s rows × %d cols", format(nrow(df), big.mark = ","), ncol(df)))

# ── 4. Build gtsummary table ───────────────────────────────────────────────────
message("Building gtsummary table …")

tbl <- df %>%
  tbl_summary(
    by      = compliance,
    missing = "no",
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous()  ~ 1,
      all_categorical() ~ c(0, 1)
    )
  ) %>%
  add_overall() %>%
  add_p() %>%
  bold_p() %>%
  bold_labels() %>%
  modify_spanning_header(
    all_stat_cols(stat_0 = FALSE) ~ "**Lei 12.732 Compliance**"
  ) %>%
  modify_caption(
    "**Table 1.** Cohort characteristics by Lei 12.732 compliance (time to treatment \u226460 days)"
  )

# ── 5. Save as HTML ───────────────────────────────────────────────────────────
message(sprintf("Saving HTML → %s", OUT_HTML))
dir.create(dirname(OUT_HTML), showWarnings = FALSE, recursive = TRUE)

tbl %>%
  as_gt() %>%
  gt::gtsave(OUT_HTML)

message("Done. Table saved to: ", OUT_HTML)

# Quick console preview
print(tbl)
