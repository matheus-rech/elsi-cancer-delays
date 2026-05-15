# 09_competing_risks.R
# Competing risks analysis — treatment vs death-before-treatment
# Output: analysis/figures/fig8_competing_risks.pdf
#         analysis/tables/table3_finegray.html

suppressPackageStartupMessages({
  for (pkg in c("tidycmprsk", "ggsurvfit", "gtsummary", "gt", "patchwork", "aod")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  library(survival)
  library(tidycmprsk)
  library(ggsurvfit)
  library(gtsummary)
  library(gt)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# ── paths ──────────────────────────────────────────────────────────────────────
BASE  <- "analysis"
DATA  <- file.path(BASE, "data", "apac_cohort.rds")
FIGS  <- file.path(BASE, "figures")
TABS  <- file.path(BASE, "tables")
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)
dir.create(TABS, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load & build competing risk status ─────────────────────────────────────
# Event codes:
#   0 = censored at 365 days (ttt > 365 and no in-hospital death)
#   1 = treated within 365 days (received chemotherapy — primary event)
#   2 = died before treatment (ttt > 365 AND death == 1 — competing event)
#
# Conservative rule: if ttt <= 365 → treated regardless of death flag (they
# received an APAC-registered chemo cycle). If ttt > 365 AND death == 1 →
# died before receiving outpatient chemo. Otherwise censored.

d_raw <- readRDS(DATA)

d_fg <- d_raw |>
  mutate(
    time365  = pmin(ttt, 365),
    cr_status = case_when(
      ttt <= 365              ~ 1L,   # treated (primary event)
      ttt >  365 & death == 1 ~ 2L,   # died before treatment (competing)
      TRUE                    ~ 0L    # censored
    ),
    cr_status = factor(cr_status, levels = c(0, 1, 2),
                       labels = c("Censored", "Treated", "Died before treatment"))
  ) |>
  filter(
    !is.na(sex_label),
    !is.na(age_grp),
    !is.na(region),
    !is.na(cancer_group),
    cancer_group != "Other",
    stage_label %in% c("I", "II", "III", "IV"),
    !is.na(stage_label)
  )

cat(sprintf("Competing risks dataset: %d patients (excluded %d)\n",
            nrow(d_fg), nrow(d_raw) - nrow(d_fg)))
cat("Event distribution:\n")
print(table(d_fg$cr_status))
cat("Event proportions:\n")
print(round(prop.table(table(d_fg$cr_status)) * 100, 1))

# ── 2. Set reference levels ────────────────────────────────────────────────────
d_fg <- d_fg |>
  mutate(
    sex_label    = relevel(factor(sex_label),    ref = "Female"),
    age_grp      = relevel(factor(age_grp),      ref = "50-59"),
    cancer_group = relevel(factor(cancer_group), ref = "Breast"),
    stage_label  = relevel(factor(stage_label),  ref = "I"),
    region       = relevel(factor(region),        ref = "Southeast")
  )

cat("\nFactor reference levels confirmed:\n")
cat("  sex_label   :", levels(d_fg$sex_label)[1],    "\n")
cat("  age_grp     :", levels(d_fg$age_grp)[1],      "\n")
cat("  cancer_group:", levels(d_fg$cancer_group)[1], "\n")
cat("  stage_label :", levels(d_fg$stage_label)[1],  "\n")
cat("  region      :", levels(d_fg$region)[1],       "\n")

# ── 3. Cumulative incidence — overall ─────────────────────────────────────────
cat("\nFitting overall cumulative incidence...\n")

cif_overall <- cuminc(Surv(time365, cr_status) ~ 1, data = d_fg)

n_total    <- nrow(d_fg)
n_treated  <- sum(d_fg$cr_status == "Treated")
n_died_bt  <- sum(d_fg$cr_status == "Died before treatment")
n_censored <- sum(d_fg$cr_status == "Censored")

p_overall <- cif_overall |>
  ggcuminc(outcome = c("Treated", "Died before treatment")) +
  geom_vline(xintercept = 60, linetype = "dashed",
             colour = "grey30", linewidth = 0.6) +
  annotate("text", x = 63, y = 0.03,
           label = "Lei 12.732\n(60 days)",
           hjust = 0, size = 3, colour = "grey30") +
  scale_x_continuous(
    limits = c(0, 365),
    breaks = seq(0, 365, by = 60),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_colour_manual(
    values = c("Treated" = "#2E86AB", "Died before treatment" = "#E74C3C"),
    name   = "Event"
  ) +
  scale_fill_manual(
    values = c("Treated" = "#2E86AB", "Died before treatment" = "#E74C3C"),
    name   = "Event"
  ) +
  add_confidence_interval() +
  add_risktable(
    risktable_stats = "{n.risk}",
    theme = theme_risktable_default(axis.text.y.size = 9)
  ) +
  labs(
    title    = "Cumulative Incidence of Treatment vs Death Before Treatment",
    subtitle = sprintf(
      "n = %s | Treated: %s | Died before treatment: %s | Censored: %s",
      format(n_total,    big.mark = ","),
      format(n_treated,  big.mark = ","),
      format(n_died_bt,  big.mark = ","),
      format(n_censored, big.mark = ",")
    ),
    x = "Days from diagnosis",
    y = "Cumulative incidence"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(colour = "grey40", size = 9),
    legend.position = "right"
  )

# ── 4. Cumulative incidence of treatment by region ─────────────────────────────
cat("Fitting cumulative incidence by region...\n")

region_levels  <- c("North", "Northeast", "Southeast", "South", "Central-West")
region_palette <- c(
  "North"        = "#E74C3C",
  "Northeast"    = "#F39C12",
  "Southeast"    = "#3498DB",
  "South"        = "#27AE60",
  "Central-West" = "#8E44AD"
)

d_region <- d_fg |>
  filter(region %in% region_levels) |>
  mutate(region = factor(region, levels = region_levels))

cif_region <- cuminc(Surv(time365, cr_status) ~ region, data = d_region)

p_region <- cif_region |>
  ggcuminc(outcome = "Treated") +
  geom_vline(xintercept = 60, linetype = "dashed",
             colour = "grey30", linewidth = 0.6) +
  annotate("text", x = 63, y = 0.03,
           label = "Lei 12.732\n(60 days)",
           hjust = 0, size = 3, colour = "grey30") +
  scale_x_continuous(
    limits = c(0, 365),
    breaks = seq(0, 365, by = 60),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_colour_manual(values = region_palette, name = "Region") +
  scale_fill_manual(values  = region_palette, name = "Region") +
  add_confidence_interval() +
  add_risktable(
    risktable_stats = "{n.risk}",
    theme = theme_risktable_default(axis.text.y.size = 8)
  ) +
  labs(
    title    = "Cumulative Incidence of Treatment by Region",
    subtitle = "Accounting for death before treatment as competing event",
    x        = "Days from diagnosis",
    y        = "Cumulative incidence of treatment"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(colour = "grey40", size = 9),
    legend.position = "right"
  )

# ── 5. Fine-Gray subdistribution hazard model ──────────────────────────────────
FG_CACHE <- file.path(BASE, "data", "fg_fit_cache.rds")

if (file.exists(FG_CACHE)) {
  cat("\nLoading cached Fine-Gray model...\n")
  fg_fit <- readRDS(FG_CACHE)
} else {
  cat("\nFitting Fine-Gray model (this may take 10-20 min on large data)...\n")
  fg_fit <- crr(
    Surv(time365, cr_status) ~ age_grp + sex_label + cancer_group + stage_label + region,
    data = d_fg
  )
  saveRDS(fg_fit, FG_CACHE)
  cat("Fine-Gray model cached to:", FG_CACHE, "\n")
}

cat("\n\n══════════════════════════════════════════════════════\n")
cat("Fine-Gray model summary (outcome: treated within 365 days)\n")
cat("══════════════════════════════════════════════════════\n")
print(summary(fg_fit))

# ── 6. Fine-Gray table (gtsummary) ────────────────────────────────────────────
cat("\nBuilding Fine-Gray table...\n")

fg_tbl <- tbl_regression(
  fg_fit,
  exponentiate = TRUE,
  label = list(
    age_grp      ~ "Age group",
    sex_label    ~ "Sex",
    cancer_group ~ "Cancer type",
    stage_label  ~ "Stage",
    region       ~ "Region"
  )
) |>
  bold_p() |>
  bold_labels() |>
  modify_header(estimate = "**SHR (95% CI)**") |>
  modify_caption(
    sprintf(
      "**Table 3. Fine-Gray Subdistribution Hazard Ratios for Time to Treatment** (n = %s; competing event: death before treatment; reference: Female, 50-59, Breast, Stage I, Southeast)",
      format(nrow(d_fg), big.mark = ",")
    )
  ) |>
  add_global_p()

tab3_path <- file.path(TABS, "table3_finegray.html")
gt::gtsave(as_gt(fg_tbl), filename = tab3_path)
cat("Saved Fine-Gray table:", tab3_path, "\n")

# ── 7. Combine plots and save figure ──────────────────────────────────────────
cat("\nCombining plots and saving figure...\n")

# patchwork cannot directly combine ggsurvfit objects with risk tables;
# save each panel separately within one PDF using layout on multiple pages.

fig8_path <- file.path(FIGS, "fig8_competing_risks.pdf")

pdf(fig8_path, width = 12, height = 14)

# Page 1 — overall cumulative incidence
print(p_overall)

# Page 2 — by region
print(p_region)

dev.off()

cat("Saved competing risks figure:", fig8_path, "\n")

# ── 8. Print key SHR findings ─────────────────────────────────────────────────
cat("\n\n══════════════════════════════════════════════════════\n")
cat("Key subdistribution hazard ratio findings (SHR [95% CI], p-value)\n")
cat("══════════════════════════════════════════════════════\n")

shr_df <- broom::tidy(fg_fit, conf.int = TRUE, exponentiate = TRUE) |>
  mutate(
    label = case_when(
      term == "age_grp60-69"              ~ "Age 60-69 vs 50-59",
      term == "age_grp70-79"              ~ "Age 70-79 vs 50-59",
      term == "age_grp80+"               ~ "Age 80+ vs 50-59",
      term == "sex_labelMale"             ~ "Male vs Female",
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
    )
  )

shr_df |>
  mutate(across(c(estimate, conf.low, conf.high), \(x) round(x, 2)),
         p.value = round(p.value, 4)) |>
  select(label, SHR = estimate, CI_low = conf.low, CI_high = conf.high, p = p.value) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\nDONE — Task 9 complete.\n")
cat(sprintf("  Figure : %s\n", fig8_path))
cat(sprintf("  Table  : %s\n", tab3_path))
