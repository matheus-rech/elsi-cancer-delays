# 05_cox_regression.R
# Cox proportional hazards regression — predictors of time to first chemotherapy
# Output: analysis/figures/fig5_cox_forest.pdf
#         analysis/tables/table2_cox.html

suppressPackageStartupMessages({
  for (pkg in c("gtsummary", "gt")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  library(survival)
  library(broom)
  library(ggplot2)
  library(dplyr)
  library(gtsummary)
  library(gt)
})

# ── paths ──────────────────────────────────────────────────────────────────────
BASE  <- "analysis"
DATA  <- file.path(BASE, "data", "apac_cohort.rds")
FIGS  <- file.path(BASE, "figures")
TABS  <- file.path(BASE, "tables")
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)
dir.create(TABS, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load & prepare modeling dataset ────────────────────────────────────────
d_raw <- readRDS(DATA)

# race mapping: 01=White, 02=Black, 03=Mixed, 04=Asian, 05=Indigenous; else NA
race_map <- c("01" = "White", "02" = "Black", "03" = "Mixed",
              "04" = "Asian",  "05" = "Indigenous")

d_cox <- d_raw |>
  mutate(
    time365   = pmin(ttt, 365),
    event365  = as.integer(ttt <= 365),
    race_label = race_map[race]          # NA for any race not in map (e.g. "99")
  ) |>
  filter(
    !is.na(sex_label),
    !is.na(age_grp),
    !is.na(region),
    !is.na(race_label),                  # drops unmapped race codes
    cancer_group != "Other",
    !is.na(cancer_group),
    stage_label %in% c("I", "II", "III", "IV"),
    !is.na(stage_label)
  )

cat(sprintf("Modeling dataset: %d patients (excluded %d)\n",
            nrow(d_cox), nrow(d_raw) - nrow(d_cox)))
cat("Event rate:", round(mean(d_cox$event365) * 100, 1), "%\n")

# ── 2. Set reference levels ────────────────────────────────────────────────────
d_cox <- d_cox |>
  mutate(
    sex_label    = relevel(factor(sex_label),    ref = "Female"),
    age_grp      = relevel(factor(age_grp),      ref = "50-59"),
    race_label   = relevel(factor(race_label),   ref = "White"),
    cancer_group = relevel(factor(cancer_group), ref = "Breast"),
    stage_label  = relevel(factor(stage_label),  ref = "I"),
    region       = relevel(factor(region),        ref = "Southeast")
  )

cat("\nFactor levels:\n")
cat("  sex_label   :", levels(d_cox$sex_label),   "\n")
cat("  age_grp     :", levels(d_cox$age_grp),     "\n")
cat("  race_label  :", levels(d_cox$race_label),  "\n")
cat("  cancer_group:", levels(d_cox$cancer_group),"\n")
cat("  stage_label :", levels(d_cox$stage_label), "\n")
cat("  region      :", levels(d_cox$region),      "\n")

# ── 3. Fit Cox model ───────────────────────────────────────────────────────────
cox_fit <- coxph(
  Surv(time365, event365) ~ age_grp + sex_label + race_label +
    cancer_group + stage_label + region,
  data = d_cox
)

cat("\n\n══════════════════════════════════════════════════════\n")
cat("Cox model summary\n")
cat("══════════════════════════════════════════════════════\n")
print(summary(cox_fit))

cat("\n\n══════════════════════════════════════════════════════\n")
cat("Schoenfeld test (proportional hazards assumption)\n")
cat("══════════════════════════════════════════════════════\n")
ph_test <- cox.zph(cox_fit)
print(ph_test)

# ── 4. Forest plot (ggplot2) ───────────────────────────────────────────────────
hr_df <- tidy(cox_fit, conf.int = TRUE, exponentiate = TRUE) |>
  filter(term != "(Intercept)") |>
  mutate(
    # Pretty labels
    label = case_when(
      term == "age_grp60-69"              ~ "Age 60-69 vs 50-59",
      term == "age_grp70-79"              ~ "Age 70-79 vs 50-59",
      term == "age_grp80+"               ~ "Age 80+ vs 50-59",
      term == "sex_labelMale"             ~ "Male vs Female",
      term == "race_labelBlack"           ~ "Race: Black vs White",
      term == "race_labelMixed"           ~ "Race: Mixed vs White",
      term == "race_labelAsian"           ~ "Race: Asian vs White",
      term == "race_labelIndigenous"      ~ "Race: Indigenous vs White",
      term == "cancer_groupCervical"      ~ "Cancer: Cervical vs Breast",
      term == "cancer_groupColorectal"    ~ "Cancer: Colorectal vs Breast",
      term == "cancer_groupLeukemia"      ~ "Cancer: Leukemia vs Breast",
      term == "cancer_groupLung"          ~ "Cancer: Lung vs Breast",
      term == "cancer_groupLymphoma"      ~ "Cancer: Lymphoma vs Breast",
      term == "cancer_groupProstate"      ~ "Cancer: Prostate vs Breast",
      term == "cancer_groupStomach"       ~ "Cancer: Stomach vs Breast",
      term == "stage_labelII"             ~ "Stage II vs I",
      term == "stage_labelIII"            ~ "Stage III vs I",
      term == "stage_labelIV"             ~ "Stage IV vs I",
      term == "regionNorth"               ~ "Region: North vs Southeast",
      term == "regionNortheast"           ~ "Region: Northeast vs Southeast",
      term == "regionSouth"               ~ "Region: South vs Southeast",
      term == "regionCentral-West"        ~ "Region: Central-West vs Southeast",
      TRUE ~ term
    ),
    # Group ordering for plot
    group = case_when(
      grepl("^Age", label)    ~ "Age",
      grepl("^Male", label)   ~ "Sex",
      grepl("^Race", label)   ~ "Race",
      grepl("^Cancer", label) ~ "Cancer type",
      grepl("^Stage", label)  ~ "Stage",
      grepl("^Region", label) ~ "Region",
      TRUE ~ "Other"
    ),
    group = factor(group, levels = c("Age", "Sex", "Race", "Cancer type", "Stage", "Region")),
    sig   = p.value < 0.05
  ) |>
  arrange(group, label)

# Ordered factor for y-axis
hr_df$label <- factor(hr_df$label, levels = rev(hr_df$label))

forest_plot <- ggplot(hr_df, aes(x = estimate, y = label, colour = group)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.25, linewidth = 0.6, alpha = 0.8) +
  geom_point(aes(shape = sig), size = 2.8) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2, 3),
    labels = c("0.5", "0.75", "1", "1.25", "1.5", "2", "3")
  ) +
  scale_colour_brewer(palette = "Dark2", name = "Variable group") +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16),
                     labels = c("FALSE" = "p ≥ 0.05", "TRUE" = "p < 0.05"),
                     name = "Significance") +
  labs(
    title    = "Predictors of Time to First Chemotherapy",
    subtitle = "Cox proportional hazards model — hazard ratios (95% CI)",
    x        = "Hazard Ratio (log scale)",
    y        = NULL,
    caption  = sprintf("n = %s patients | Reference: Female, 50-59, White, Breast, Stage I, Southeast",
                       format(nrow(d_cox), big.mark = ","))
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(colour = "grey40", size = 10),
    legend.position = "right",
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.4),
    axis.text.y   = element_text(size = 9)
  )

fig5_path <- file.path(FIGS, "fig5_cox_forest.pdf")
ggsave(fig5_path, plot = forest_plot, width = 10, height = 8,
       units = "in", dpi = 300, device = "pdf")
cat("\nSaved forest plot:", fig5_path, "\n")

# ── 5. Formatted table via gtsummary ──────────────────────────────────────────
cox_tbl <- tbl_regression(cox_fit, exponentiate = TRUE) |>
  bold_p() |>
  bold_labels()

tab2_path <- file.path(TABS, "table2_cox.html")
gt::gtsave(as_gt(cox_tbl), filename = tab2_path)
cat("Saved gtsummary table:", tab2_path, "\n")

# ── 6. Print key HR findings ───────────────────────────────────────────────────
cat("\n\n══════════════════════════════════════════════════════\n")
cat("Key hazard ratio findings (HR [95% CI], p-value)\n")
cat("══════════════════════════════════════════════════════\n")
hr_df |>
  mutate(across(c(estimate, conf.low, conf.high), \(x) round(x, 2)),
         p.value = round(p.value, 4)) |>
  select(label, HR = estimate, CI_low = conf.low, CI_high = conf.high, p = p.value) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\nDONE — Task 5 complete.\n")
