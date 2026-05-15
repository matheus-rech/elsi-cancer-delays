# Cancer Treatment Delays in Older Brazilians — Publication Pipeline

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a reproducible analysis pipeline and publication-quality figures/tables for a paper on cancer treatment delays among older adults in Brazil, combining ELSI-Brazil survey data with DATASUS APAC-ONCO administrative records.

**Architecture:** Two complementary data sources analyzed independently then linked ecologically. APAC-QT (chemotherapy) microdatasus downloads → individual-level KM + Cox. ELSI Wave 2 → survey-weighted population health context. Ecological join at the UF (state) level for correlation analysis. All outputs in a Quarto document with reproducible R code.

**Tech Stack:** R (microdatasus, survival, survminer, srvyr, ggplot2, gtsummary, geobr, sf, ranger, quarto)

**Target journal:** BMJ Global Health / Lancet Regional Health Americas / Int J Cancer

---

### Task 1: Download APAC-QT data for 5 representative UFs (one per macro-region)

Representative UFs: SP (Southeast), BA (Northeast), RS (South), PA (North), GO (Central-West). Full year 2022. This gives geographic diversity without downloading all 27 UFs (which can be expanded later).

**Files:**
- Create: `analysis/01_download_apac.R`
- Create: `analysis/data/apac_qt_5uf_2022.rds`

- [ ] **Step 1: Create the download script**

```r
# analysis/01_download_apac.R
library(microdatasus)
library(dplyr)

ufs <- c("SP", "BA", "RS", "PA", "GO")
results <- list()

for (uf in ufs) {
  cat("Downloading SIA-AQ for", uf, "2022...\n")
  d <- fetch_datasus(
    year_start = 2022, month_start = 1,
    year_end = 2022, month_end = 12,
    uf = uf, information_system = "SIA-AQ",
    timeout = 600
  )
  d$uf_download <- uf
  results[[uf]] <- d
  cat("  ", nrow(d), "rows\n")
}

apac <- bind_rows(results)
cat("Total:", nrow(apac), "rows\n")
saveRDS(apac, "analysis/data/apac_qt_5uf_2022_raw.rds")
```

- [ ] **Step 2: Run the download**

Run: `Rscript analysis/01_download_apac.R`
Expected: ~1-2M rows across 5 UFs, saved as RDS (~50-100 MB)

- [ ] **Step 3: Commit**

```bash
git add analysis/01_download_apac.R
git commit -m "feat: download APAC-QT data for 5 representative UFs"
```

---

### Task 2: Build the APAC analytic cohort (deduplication, variable derivation)

Transform raw APAC records into a one-row-per-patient analytic dataset with time-to-treatment, demographics, cancer type, staging, and geographic variables.

**Files:**
- Create: `analysis/02_build_apac_cohort.R`
- Create: `analysis/data/apac_cohort.rds`

- [ ] **Step 1: Write the cohort-building script**

```r
# analysis/02_build_apac_cohort.R
library(dplyr)

raw <- readRDS("analysis/data/apac_qt_5uf_2022_raw.rds")

cohort <- raw %>%
  mutate(
    dt_dx = as.Date(AQ_DTIDEN, format = "%Y%m%d"),
    dt_tx = as.Date(AQ_DTINTR, format = "%Y%m%d"),
    age = as.numeric(AP_NUIDADE),
    sex = AP_SEXO,
    race = AP_RACACOR,
    mun_res = AP_MUNPCN,
    uf_res = substr(AP_MUNPCN, 1, 2),
    cid = AQ_CID10,
    stage = AQ_ESTADI,
    obito = as.integer(AP_OBITO),
    cns = AP_CNSPCN
  ) %>%
  filter(!is.na(dt_dx), !is.na(dt_tx), age >= 50) %>%
  # Deduplicate: first APAC per patient-cancer pair
  group_by(cns, cid) %>%
  arrange(dt_tx) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    ttt = as.numeric(dt_tx - dt_dx),
    sex_label = case_when(sex == "F" ~ "Female", sex == "M" ~ "Male"),
    age_grp = cut(age, c(50, 60, 70, 80, Inf), right = FALSE,
                  labels = c("50-59", "60-69", "70-79", "80+")),
    stage_label = case_when(
      stage == "0" ~ "In situ", stage == "1" ~ "I",
      stage == "2" ~ "II", stage == "3" ~ "III", stage == "4" ~ "IV"),
    cancer_group = case_when(
      grepl("^C50", cid) ~ "Breast",
      cid == "C61" ~ "Prostate",
      cid %in% c("C18","C180","C181","C182","C183","C184","C185",
                  "C186","C187","C188","C189","C19","C20") ~ "Colorectal",
      grepl("^C34", cid) ~ "Lung",
      grepl("^C16", cid) ~ "Stomach",
      grepl("^C53", cid) ~ "Cervical",
      grepl("^C92|^C91", cid) ~ "Leukemia",
      grepl("^C82|^C83|^C85", cid) ~ "Lymphoma",
      !is.na(cid) ~ "Other"),
    uf_label = case_when(
      uf_res == "35" ~ "SP", uf_res == "29" ~ "BA",
      uf_res == "43" ~ "RS", uf_res == "15" ~ "PA",
      uf_res == "52" ~ "GO", TRUE ~ uf_res),
    region = case_when(
      uf_label == "SP" ~ "Southeast", uf_label == "BA" ~ "Northeast",
      uf_label == "RS" ~ "South", uf_label == "PA" ~ "North",
      uf_label == "GO" ~ "Central-West"),
    compliant_60d = as.integer(ttt <= 60),
    compliant_30d = as.integer(ttt <= 30)
  ) %>%
  filter(ttt >= 0, ttt <= 1825)  # 0 to 5 years

cat("Cohort:", nrow(cohort), "unique patient-cancers\n")
saveRDS(cohort, "analysis/data/apac_cohort.rds")
```

- [ ] **Step 2: Run and verify**

Run: `Rscript analysis/02_build_apac_cohort.R`
Expected: ~200-400K unique elderly cancer cases across 5 UFs

- [ ] **Step 3: Commit**

```bash
git add analysis/02_build_apac_cohort.R
git commit -m "feat: build deduplicated APAC analytic cohort"
```

---

### Task 3: ELSI regional summary for ecological linkage

Compute survey-weighted regional health indicators from ELSI to correlate with APAC treatment delays at the UF/region level.

**Files:**
- Create: `analysis/03_elsi_regional.R`
- Create: `analysis/data/elsi_regional.rds`
- Read: `app/elsi_data.csv` (enhanced dataset from prep_data_enhanced.R)

- [ ] **Step 1: Write the ELSI regional summary script**

```r
# analysis/03_elsi_regional.R
library(dplyr)
library(srvyr)

elsi <- read.csv("app/elsi_data.csv")

svy <- elsi %>%
  mutate(
    frail_binary = as.integer(!is.na(frailty_status) & frailty_status == "Frail"),
    dep_binary = as.integer(!is.na(depression_score) & depression_score >= 12),
    adl_any = as.integer(adl_disability >= 1),
    soc_low = as.integer(social_score <= 1),
    multi3 = as.integer(n_chronic >= 3)
  ) %>%
  as_survey_design(ids = upa, strata = estrato, weights = weight, nest = TRUE)

regional <- svy %>%
  group_by(region_label) %>%
  summarise(
    n_elsi = n(),
    cancer_prev = survey_mean(cancer_dx == 1, na.rm = TRUE),
    frail_prev = survey_mean(frail_binary, na.rm = TRUE),
    mean_fi = survey_mean(frailty_index, na.rm = TRUE),
    dep_prev = survey_mean(dep_binary, na.rm = TRUE),
    adl_prev = survey_mean(adl_any, na.rm = TRUE),
    soc_low_prev = survey_mean(soc_low, na.rm = TRUE),
    mean_social = survey_mean(social_score, na.rm = TRUE),
    hosp_prev = survey_mean(hospitalized == 1, na.rm = TRUE),
    plan_prev = survey_mean(has_health_plan == 1, na.rm = TRUE),
    mean_income = survey_mean(income_pc, na.rm = TRUE),
    mean_age = survey_mean(age, na.rm = TRUE),
    mean_educ = survey_mean(education_yrs, na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(regional, "analysis/data/elsi_regional.rds")
cat("ELSI regional summary: 5 regions x", ncol(regional), "indicators\n")
print(regional)
```

- [ ] **Step 2: Run and verify**

Run: `Rscript analysis/03_elsi_regional.R`
Expected: 5-row data frame with survey-weighted indicators per macro-region

- [ ] **Step 3: Commit**

```bash
git add analysis/03_elsi_regional.R
git commit -m "feat: compute ELSI regional health indicators for ecological linkage"
```

---

### Task 4: Publication-quality KM figures

Generate 4 KM figures (overall, by cancer type, by stage, by region) with risk tables, using survminer and ggplot2, saved as 300dpi TIFF/PDF.

**Files:**
- Create: `analysis/04_km_figures.R`
- Create: `analysis/figures/fig1_km_overall.pdf`
- Create: `analysis/figures/fig2_km_cancer_type.pdf`
- Create: `analysis/figures/fig3_km_stage.pdf`
- Create: `analysis/figures/fig4_km_region.pdf`

- [ ] **Step 1: Install survminer if needed**

Run: `Rscript -e 'if (!requireNamespace("survminer")) install.packages("survminer")'`

- [ ] **Step 2: Write the KM figure script**

```r
# analysis/04_km_figures.R
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)

cohort <- readRDS("analysis/data/apac_cohort.rds")

# Cap at 365 days for visualization
cohort$time365 <- pmin(cohort$ttt, 365)
cohort$event365 <- as.integer(cohort$ttt <= 365)

dir.create("analysis/figures", showWarnings = FALSE, recursive = TRUE)

# Color palette
region_cols <- c("North"="#E74C3C", "Northeast"="#F39C12",
                 "Southeast"="#3498DB", "South"="#27AE60", "Central-West"="#8E44AD")

# Fig 1: Overall KM with 60-day reference line
fit1 <- survfit(Surv(time365, event365) ~ 1, data = cohort)
p1 <- ggsurvplot(fit1, data = cohort, fun = "event",
                 risk.table = TRUE, risk.table.col = "strata",
                 conf.int = TRUE, pval = FALSE,
                 xlab = "Days from diagnosis", ylab = "Cumulative proportion treated",
                 title = "Time to First Chemotherapy — Elderly (≥50 years), Brazil 2022",
                 subtitle = paste0("n = ", format(nrow(cohort), big.mark = ","),
                                   " | Median: ", median(cohort$ttt), " days"),
                 ggtheme = theme_minimal(base_size = 14),
                 break.x.by = 30, xlim = c(0, 365),
                 surv.scale = "percent")
p1$plot <- p1$plot +
  geom_vline(xintercept = 60, linetype = "dashed", color = "#1F3864", linewidth = 0.8) +
  annotate("text", x = 65, y = 0.05, label = "Lei 12.732\n(60 days)",
           hjust = 0, color = "#1F3864", size = 3.5, fontface = "italic")
ggsave("analysis/figures/fig1_km_overall.pdf", print(p1), width = 10, height = 8, dpi = 300)

# Fig 2: By cancer type (top 6)
top_cancers <- c("Breast","Prostate","Colorectal","Lung","Stomach","Cervical")
d2 <- cohort %>% filter(cancer_group %in% top_cancers)
fit2 <- survfit(Surv(time365, event365) ~ cancer_group, data = d2)
p2 <- ggsurvplot(fit2, data = d2, fun = "event",
                 risk.table = TRUE, pval = TRUE,
                 xlab = "Days from diagnosis", ylab = "Cumulative proportion treated",
                 title = "Time to Treatment by Cancer Type",
                 ggtheme = theme_minimal(base_size = 14),
                 break.x.by = 60, xlim = c(0, 365),
                 legend.title = "Cancer", surv.scale = "percent")
p2$plot <- p2$plot +
  geom_vline(xintercept = 60, linetype = "dashed", color = "#1F3864")
ggsave("analysis/figures/fig2_km_cancer_type.pdf", print(p2), width = 12, height = 10, dpi = 300)

# Fig 3: By stage
d3 <- cohort %>% filter(stage_label %in% c("I","II","III","IV"))
fit3 <- survfit(Surv(time365, event365) ~ stage_label, data = d3)
p3 <- ggsurvplot(fit3, data = d3, fun = "event",
                 risk.table = TRUE, pval = TRUE,
                 xlab = "Days from diagnosis", ylab = "Cumulative proportion treated",
                 title = "Time to Treatment by Staging at Diagnosis",
                 ggtheme = theme_minimal(base_size = 14),
                 palette = c("#27AE60","#F1C40F","#E67E22","#E74C3C"),
                 break.x.by = 60, xlim = c(0, 365),
                 legend.title = "Stage", surv.scale = "percent")
p3$plot <- p3$plot +
  geom_vline(xintercept = 60, linetype = "dashed", color = "#1F3864")
ggsave("analysis/figures/fig3_km_stage.pdf", print(p3), width = 12, height = 10, dpi = 300)

# Fig 4: By region
d4 <- cohort %>% filter(!is.na(region))
fit4 <- survfit(Surv(time365, event365) ~ region, data = d4)
p4 <- ggsurvplot(fit4, data = d4, fun = "event",
                 risk.table = TRUE, pval = TRUE,
                 xlab = "Days from diagnosis", ylab = "Cumulative proportion treated",
                 title = "Time to Treatment by Brazilian Macro-Region",
                 ggtheme = theme_minimal(base_size = 14),
                 palette = region_cols,
                 break.x.by = 60, xlim = c(0, 365),
                 legend.title = "Region", surv.scale = "percent")
p4$plot <- p4$plot +
  geom_vline(xintercept = 60, linetype = "dashed", color = "#1F3864")
ggsave("analysis/figures/fig4_km_region.pdf", print(p4), width = 12, height = 10, dpi = 300)

cat("All figures saved to analysis/figures/\n")
```

- [ ] **Step 3: Run and verify all 4 figures**

Run: `Rscript analysis/04_km_figures.R`
Expected: 4 PDF files in analysis/figures/

- [ ] **Step 4: Commit**

```bash
git add analysis/04_km_figures.R analysis/figures/
git commit -m "feat: publication-quality KM figures with risk tables"
```

---

### Task 5: Cox proportional hazards regression

Multivariable Cox model for predictors of treatment delay. Outcome = time to treatment. Covariates: age, sex, race, cancer type, stage, region. Output: forest plot + formatted table.

**Files:**
- Create: `analysis/05_cox_regression.R`
- Create: `analysis/figures/fig5_cox_forest.pdf`
- Create: `analysis/tables/table2_cox.html`

- [ ] **Step 1: Install gtsummary if needed**

Run: `Rscript -e 'if (!requireNamespace("gtsummary")) install.packages("gtsummary")'`

- [ ] **Step 2: Write the Cox regression script**

```r
# analysis/05_cox_regression.R
library(dplyr)
library(survival)
library(ggplot2)
library(broom)
library(gtsummary)

cohort <- readRDS("analysis/data/apac_cohort.rds")

# Prepare modeling dataset
d_cox <- cohort %>%
  filter(!is.na(sex_label), !is.na(age_grp), !is.na(region),
         cancer_group != "Other", !is.na(stage_label),
         stage_label != "In situ") %>%
  mutate(
    sex_label = factor(sex_label, levels = c("Female", "Male")),
    age_grp = factor(age_grp, levels = c("50-59","60-69","70-79","80+")),
    region = factor(region, levels = c("Southeast","South","Central-West","Northeast","North")),
    stage_label = factor(stage_label, levels = c("I","II","III","IV")),
    cancer_group = factor(cancer_group),
    race_label = case_when(
      race == "01" ~ "White", race == "02" ~ "Black",
      race == "03" ~ "Mixed", race == "04" ~ "Asian",
      race == "05" ~ "Indigenous", TRUE ~ NA_character_),
    race_label = factor(race_label, levels = c("White","Mixed","Black","Asian","Indigenous"))
  ) %>%
  filter(!is.na(race_label))

cat("Cox modeling dataset:", nrow(d_cox), "rows\n")

# Cap follow-up at 365 days
d_cox$time365 <- pmin(d_cox$ttt, 365)
d_cox$event365 <- as.integer(d_cox$ttt <= 365)

# Fit Cox model
# Note: in "time to treatment" analysis, event = receiving treatment
# Higher HR = faster treatment (shorter time to event)
cox_fit <- coxph(
  Surv(time365, event365) ~ age_grp + sex_label + race_label +
    cancer_group + stage_label + region,
  data = d_cox
)

cat("\n=== COX MODEL SUMMARY ===\n")
print(summary(cox_fit))

# Test proportional hazards assumption
cat("\n=== SCHOENFELD TEST ===\n")
ph_test <- cox.zph(cox_fit)
print(ph_test)

# Forest plot
tidy_cox <- tidy(cox_fit, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(term = gsub("sex_label|age_grp|race_label|cancer_group|stage_label|region", "", term))

dir.create("analysis/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("analysis/tables", showWarnings = FALSE, recursive = TRUE)

p_forest <- ggplot(tidy_cox, aes(estimate, reorder(term, estimate))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5, color = "#1F3864") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "#1F3864") +
  scale_x_log10() +
  labs(x = "Hazard Ratio (95% CI, log scale)",
       y = NULL,
       title = "Predictors of Time to First Chemotherapy",
       subtitle = "HR > 1 = faster treatment | HR < 1 = longer delay") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", color = "#1F3864"))

ggsave("analysis/figures/fig5_cox_forest.pdf", p_forest, width = 10, height = 8, dpi = 300)

# Formatted table via gtsummary
tbl <- tbl_regression(cox_fit, exponentiate = TRUE) %>%
  bold_p() %>%
  bold_labels()
gt::gtsave(as_gt(tbl), "analysis/tables/table2_cox.html")

cat("Forest plot and table saved.\n")
```

- [ ] **Step 3: Run and verify**

Run: `Rscript analysis/05_cox_regression.R`
Expected: Cox model with HR, 95% CI, p-values. Forest plot PDF. HTML table.

- [ ] **Step 4: Commit**

```bash
git add analysis/05_cox_regression.R analysis/figures/fig5_cox_forest.pdf analysis/tables/
git commit -m "feat: Cox regression for treatment delay predictors"
```

---

### Task 6: Table 1 — Cohort characteristics (gtsummary)

Standard "Table 1" for the paper: demographics, cancer characteristics, and treatment outcomes, stratified by compliance with Lei 12.732 (≤60 days vs >60 days).

**Files:**
- Create: `analysis/06_table1.R`
- Create: `analysis/tables/table1_characteristics.html`

- [ ] **Step 1: Write the Table 1 script**

```r
# analysis/06_table1.R
library(dplyr)
library(gtsummary)

cohort <- readRDS("analysis/data/apac_cohort.rds")

d_tbl <- cohort %>%
  mutate(
    compliance = ifelse(compliant_60d == 1, "≤60 days", ">60 days"),
    compliance = factor(compliance, levels = c("≤60 days", ">60 days")),
    race_label = case_when(
      race == "01" ~ "White", race == "02" ~ "Black",
      race == "03" ~ "Mixed", TRUE ~ "Other/Unknown")
  ) %>%
  select(
    compliance, age, sex_label, race_label, age_grp,
    cancer_group, stage_label, region, ttt, obito
  ) %>%
  rename(
    `Age (years)` = age,
    Sex = sex_label,
    `Race/ethnicity` = race_label,
    `Age group` = age_grp,
    `Cancer type` = cancer_group,
    Stage = stage_label,
    Region = region,
    `Time to treatment (days)` = ttt,
    `In-hospital death` = obito
  )

tbl1 <- tbl_summary(
  d_tbl,
  by = compliance,
  statistic = list(
    all_continuous() ~ "{median} ({p25}, {p75})",
    all_categorical() ~ "{n} ({p}%)"
  ),
  missing = "no"
) %>%
  add_overall() %>%
  add_p() %>%
  bold_p() %>%
  bold_labels() %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Lei 12.732 Compliance**")

dir.create("analysis/tables", showWarnings = FALSE, recursive = TRUE)
gt::gtsave(as_gt(tbl1), "analysis/tables/table1_characteristics.html")
cat("Table 1 saved.\n")
```

- [ ] **Step 2: Run and verify**

Run: `Rscript analysis/06_table1.R`
Expected: HTML table with demographics stratified by compliance

- [ ] **Step 3: Commit**

```bash
git add analysis/06_table1.R analysis/tables/table1_characteristics.html
git commit -m "feat: Table 1 cohort characteristics by Lei 12.732 compliance"
```

---

### Task 7: Ecological analysis — ELSI health indicators × APAC treatment delays

Correlate regional-level ELSI indicators (frailty, social isolation, insurance, income) with APAC treatment delays. Produce scatter plots and correlation matrix.

**Files:**
- Create: `analysis/07_ecological.R`
- Create: `analysis/figures/fig6_ecological.pdf`

- [ ] **Step 1: Write the ecological analysis script**

```r
# analysis/07_ecological.R
library(dplyr)
library(ggplot2)
library(patchwork)

elsi_reg <- readRDS("analysis/data/elsi_regional.rds")
cohort <- readRDS("analysis/data/apac_cohort.rds")

# APAC regional summary
apac_reg <- cohort %>%
  group_by(region) %>%
  summarise(
    n_apac = n(),
    median_ttt = median(ttt),
    pct_compliant = mean(compliant_60d) * 100,
    pct_over_180d = mean(ttt > 180) * 100,
    .groups = "drop"
  )

# Merge
merged <- elsi_reg %>%
  rename(region = region_label) %>%
  inner_join(apac_reg, by = "region")

cat("Merged ecological data:\n")
print(merged %>% select(region, frail_prev, plan_prev, mean_income,
                         soc_low_prev, median_ttt, pct_compliant))

# Scatter plots
make_scatter <- function(data, x, y, xlab, ylab, title) {
  ggplot(data, aes(.data[[x]], .data[[y]], label = region)) +
    geom_point(size = 4, color = "#E74C3C") +
    geom_text(vjust = -1, size = 3.5) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "#1F3864") +
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 12, color = "#1F3864"))
}

p1 <- make_scatter(merged, "plan_prev", "pct_compliant",
  "Health plan coverage (%)", "Lei 12.732 compliance (%)",
  "Insurance coverage vs treatment compliance")

p2 <- make_scatter(merged, "frail_prev", "median_ttt",
  "Frailty prevalence (%)", "Median time to treatment (days)",
  "Frailty burden vs treatment delay")

p3 <- make_scatter(merged, "mean_income", "pct_compliant",
  "Mean per capita income (R$)", "Lei 12.732 compliance (%)",
  "Income vs treatment compliance")

p4 <- make_scatter(merged, "soc_low_prev", "pct_over_180d",
  "Social isolation prevalence (%)", "Patients waiting >180 days (%)",
  "Social isolation vs extreme delay")

combined <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Ecological Analysis: ELSI Population Health Indicators vs APAC Treatment Delays",
    subtitle = "Each point = one Brazilian macro-region (n = 5)",
    theme = theme(plot.title = element_text(face = "bold", size = 16, color = "#1F3864")))

dir.create("analysis/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("analysis/figures/fig6_ecological.pdf", combined, width = 14, height = 10, dpi = 300)
cat("Ecological figure saved.\n")
```

- [ ] **Step 2: Run and verify**

Run: `Rscript analysis/07_ecological.R`
Expected: 4-panel scatter plot PDF

- [ ] **Step 3: Commit**

```bash
git add analysis/07_ecological.R analysis/figures/fig6_ecological.pdf
git commit -m "feat: ecological analysis linking ELSI indicators to APAC treatment delays"
```

---

### Task 8: Choropleth map — Regional treatment compliance

Map showing Lei 12.732 compliance (%) by macro-region alongside key ELSI indicators, using geobr shapefiles.

**Files:**
- Create: `analysis/08_choropleth.R`
- Create: `analysis/figures/fig7_choropleth.pdf`

- [ ] **Step 1: Write the choropleth script**

```r
# analysis/08_choropleth.R
library(dplyr)
library(ggplot2)
library(sf)
library(geobr)
library(patchwork)

elsi_reg <- readRDS("analysis/data/elsi_regional.rds")
cohort <- readRDS("analysis/data/apac_cohort.rds")

apac_reg <- cohort %>%
  group_by(region) %>%
  summarise(
    pct_compliant = mean(compliant_60d) * 100,
    median_ttt = median(ttt),
    .groups = "drop"
  )

merged <- elsi_reg %>%
  rename(region = region_label) %>%
  inner_join(apac_reg, by = "region")

# Load Brazil regions shapefile
br <- read_region(year = 2020, showProgress = FALSE)
br$region <- case_when(
  br$name_region == "Norte" ~ "North",
  br$name_region == "Nordeste" ~ "Northeast",
  br$name_region == "Sudeste" ~ "Southeast",
  br$name_region == "Sul" ~ "South",
  br$name_region == "Centro-Oeste" ~ "Central-West"
)
map_data <- br %>% left_join(merged, by = "region")

make_map <- function(data, fill_var, title, palette = "YlOrRd", direction = 1) {
  ggplot(data) +
    geom_sf(aes(fill = .data[[fill_var]]), color = "white", linewidth = 0.5) +
    geom_sf_text(aes(label = paste0(region, "\n", round(.data[[fill_var]], 1))),
                 size = 3, color = "black") +
    scale_fill_distiller(palette = palette, direction = direction, name = NULL) +
    labs(title = title) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5, color = "#1F3864"),
          legend.position = "bottom")
}

p1 <- make_map(map_data, "pct_compliant", "Lei 12.732 Compliance (%)", "RdYlGn", 1)
p2 <- make_map(map_data, "median_ttt", "Median Days to Treatment", "YlOrRd", 1)
p3 <- make_map(map_data, "frail_prev", "ELSI: Frailty Prevalence", "YlOrRd", 1)
p4 <- make_map(map_data, "plan_prev", "ELSI: Health Plan Coverage (%)", "RdYlGn", 1)

combined <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Cancer Treatment Delays and Population Health Vulnerability — Brazil, 2022",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)))

dir.create("analysis/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("analysis/figures/fig7_choropleth.pdf", combined, width = 14, height = 12, dpi = 300)
cat("Choropleth saved.\n")
```

- [ ] **Step 2: Run and verify**

Run: `Rscript analysis/08_choropleth.R`
Expected: 4-panel choropleth PDF

- [ ] **Step 3: Commit**

```bash
git add analysis/08_choropleth.R analysis/figures/fig7_choropleth.pdf
git commit -m "feat: regional choropleth maps combining APAC and ELSI data"
```

---

### Task 9: Integrate into Shiny app — "Time to Treatment" tab

Add a new tab to the existing Shiny app with interactive KM curves, compliance dashboard, and ecological analysis powered by the pre-computed APAC cohort data.

**Files:**
- Modify: `app/app.R` (add new nav_panel after "Regional Map" tab)
- Read: `app/apac_qt_sp2022.rds` (already exists from earlier work)

- [ ] **Step 1: Add the Time to Treatment tab to app.R**

Add a new `nav_panel("Time to Treatment", ...)` with:
- KM curve (plotly interactive) with 60-day reference line
- Compliance metrics (value boxes: % ≤30d, ≤60d, >60d, >180d)
- Stratification controls (cancer type, stage, age group, region)
- Data source note explaining APAC-QT/DATASUS origin

- [ ] **Step 2: Test the app**

Run: `Rscript -e 'shiny::runApp("app/", port=7878, launch.browser=TRUE)'`
Expected: New tab renders with KM curves and compliance metrics

- [ ] **Step 3: Commit**

```bash
git add app/app.R
git commit -m "feat: add Time to Treatment tab with APAC-derived KM analysis"
```
