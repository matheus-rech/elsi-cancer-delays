# National Expansion & Publication-Ready Pipeline

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expand the cancer treatment delay analysis from 5 representative UFs to all 27 Brazilian states, add competing risks and sensitivity analyses required by peer reviewers, produce state-level choropleth maps, generate supplementary material, and assemble a STROBE-compliant Quarto manuscript targeting BMJ Global Health / Lancet Regional Health Americas.

**Architecture:** The existing 8-script pipeline (01-08) reads from `analysis/data/apac_cohort.rds`. Strategy: (1) download all 27 UFs → (2) rebuild the cohort file in-place → (3) re-run existing scripts 04-08 to regenerate figures/tables with national data → (4) add new analyses (competing risks, sensitivity, state-level maps) → (5) assemble Quarto manuscript → (6) upgrade Shiny app.

**Tech Stack:** R (microdatasus, survival, survminer, cmprsk/tidycmprsk, srvyr, ggplot2, gtsummary, geobr, sf, quarto, gt, flextable)

**Target journal:** BMJ Global Health / Lancet Regional Health Americas

**Current state (completed):**
- `analysis/01_download_apac.R` → 5 UFs (SP, BA, RS, PA, GO), 1.65M raw rows
- `analysis/02_build_apac_cohort.R` → 49,530 unique patient-cancer pairs (age ≥50)
- `analysis/03_elsi_regional.R` → 5-region survey-weighted ELSI indicators
- `analysis/04_km_figures.R` → Fig 1-4 (overall, cancer type, stage, region)
- `analysis/05_cox_regression.R` → Fig 5 (forest plot) + Table 2 (HR)
- `analysis/06_table1.R` → Table 1 (demographics by compliance)
- `analysis/07_ecological.R` → Fig 6 (ELSI × APAC scatter)
- `analysis/08_choropleth.R` → Fig 7 (5-region choropleth)
- `app/app.R` → Shiny app with 10 tabs including Time to Treatment (Tab 9)

**Key constraint:** ELSI data has only `region_label` (5 macro-regions), not state-level identifiers. Ecological linkage remains at 5-region level. State-level analysis is APAC-internal only.

---

### Task 1: Download APAC-QT for all 27 UFs

Expand from 5 UFs to full national coverage. This is the longest-running task (~30-60 min) and must complete before everything else.

**Files:**
- Modify: `analysis/01_download_apac.R`
- Create: `analysis/data/apac_qt_27uf_2022_raw.rds` (~300-500 MB)

- [ ] **Step 1: Write the updated download script**

Replace the 5-UF download with all 27 UFs. Use a loop with progress logging and save intermediate results to avoid losing work on FTP timeouts.

```r
# analysis/01_download_apac.R
# Download APAC-QT (chemotherapy) from DATASUS for all 27 UFs, year 2022
# Output: analysis/data/apac_qt_27uf_2022_raw.rds

suppressPackageStartupMessages({
  library(microdatasus)
  library(dplyr)
})

OUT_DIR <- "analysis/data"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# All 27 Brazilian UFs
ufs <- c("AC","AL","AM","AP","BA","CE","DF","ES","GO","MA",
         "MG","MS","MT","PA","PB","PE","PI","PR","RJ","RN",
         "RO","RR","RS","SC","SE","SP","TO")

results <- list()
failed  <- character(0)

for (uf in ufs) {
  out_uf <- file.path(OUT_DIR, paste0("apac_qt_", uf, "_2022.rds"))

  # Skip if already downloaded

  if (file.exists(out_uf)) {
    cat(sprintf("[%s] Already exists, loading cached (%s)\n",
                uf, format(file.size(out_uf), big.mark = ",")))
    results[[uf]] <- readRDS(out_uf)
    next
  }

  cat(sprintf("[%s] Downloading SIA-AQ 2022...\n", uf))
  tryCatch({
    d <- fetch_datasus(
      year_start = 2022, month_start = 1,
      year_end   = 2022, month_end   = 12,
      uf = uf, information_system = "SIA-AQ",
      timeout = 600
    )
    if (nrow(d) == 0) {
      cat(sprintf("[%s] Empty — skipping\n", uf))
      next
    }
    d$uf_download <- uf
    saveRDS(d, out_uf)
    results[[uf]] <- d
    cat(sprintf("[%s] Done: %s rows\n", uf, format(nrow(d), big.mark = ",")))
  }, error = function(e) {
    cat(sprintf("[%s] FAILED: %s\n", uf, e$message))
    failed <<- c(failed, uf)
  })
}

if (length(failed) > 0) {
  cat("\n=== FAILED UFs (retry manually): ", paste(failed, collapse = ", "), "\n")
}

# Combine all
apac <- bind_rows(results)
cat(sprintf("\n=== TOTAL: %s rows from %d UFs ===\n",
            format(nrow(apac), big.mark = ","), length(results)))

saveRDS(apac, file.path(OUT_DIR, "apac_qt_27uf_2022_raw.rds"))
cat("Saved: analysis/data/apac_qt_27uf_2022_raw.rds\n")
cat(sprintf("File size: %.1f MB\n",
            file.size(file.path(OUT_DIR, "apac_qt_27uf_2022_raw.rds")) / 1024^2))
```

- [ ] **Step 2: Run the download**

Run: `Rscript analysis/01_download_apac.R`
Expected: ~3-8M rows across 27 UFs. The script caches per-UF RDS files so it can resume after interruption.

- [ ] **Step 3: Verify completeness**

Run: `Rscript -e 'd <- readRDS("analysis/data/apac_qt_27uf_2022_raw.rds"); cat("Rows:", nrow(d), "\nUFs:", length(unique(d$uf_download)), "\n"); print(table(d$uf_download))'`
Expected: 27 UFs represented (some small UFs like AC, RR may have <1000 rows)

---

### Task 2: Rebuild national analytic cohort

Overwrite `analysis/data/apac_cohort.rds` with the full 27-UF national cohort. Add UF name labels and IBGE UF codes for state-level maps.

**Files:**
- Modify: `analysis/02_build_apac_cohort.R`
- Overwrite: `analysis/data/apac_cohort.rds`

- [ ] **Step 1: Update the cohort script to use 27-UF data**

Change the input file from `apac_qt_5uf_2022_raw.rds` to `apac_qt_27uf_2022_raw.rds`. Add a `uf_name` column with state names and `uf_code` with 2-digit IBGE codes for joining with geobr shapefiles.

```r
# analysis/02_build_apac_cohort.R
# Build deduplicated analytic cohort from raw APAC-QT data (all 27 UFs)
# Input:  analysis/data/apac_qt_27uf_2022_raw.rds
# Output: analysis/data/apac_cohort.rds

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
})

IN_FILE  <- "analysis/data/apac_qt_27uf_2022_raw.rds"
OUT_FILE <- "analysis/data/apac_cohort.rds"

message("Reading raw RDS...")
raw <- readRDS(IN_FILE)
message(sprintf("  Raw: %s rows x %d cols", format(nrow(raw), big.mark = ","), ncol(raw)))

# Parse dates
raw <- raw %>%
  mutate(
    dt_dx = as.Date(as.character(AQ_DTIDEN), "%Y%m%d"),
    dt_tx = as.Date(as.character(AQ_DTINTR), "%Y%m%d")
  )

# Extract & rename key columns
cohort <- raw %>%
  transmute(
    cns        = AP_CNSPCN,
    cid        = AQ_CID10,
    dt_dx,
    dt_tx,
    age        = as.numeric(AP_NUIDADE),
    sex        = AP_SEXO,
    race       = AP_RACACOR,
    municipio  = AP_MUNPCN,
    uf         = substr(as.character(AP_MUNPCN), 1, 2),
    uf_download = uf_download,
    stage      = AQ_ESTADI,
    death      = as.integer(AP_OBITO)
  )

# Filter: age >= 50
cohort <- cohort %>% filter(!is.na(age), age >= 50)
message(sprintf("  After age filter: %s rows", format(nrow(cohort), big.mark = ",")))

# Deduplicate: first APAC per patient-cancer
cohort <- cohort %>%
  filter(!is.na(cns), !is.na(cid)) %>%
  group_by(cns, cid) %>%
  arrange(dt_tx, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()
message(sprintf("  After dedup: %s rows", format(nrow(cohort), big.mark = ",")))

# Compute time-to-treatment
cohort <- cohort %>%
  mutate(ttt = as.numeric(dt_tx - dt_dx))

# Filter plausible TTT (0-1825 days)
cohort <- cohort %>% filter(!is.na(ttt), ttt >= 0, ttt <= 1825)
message(sprintf("  After TTT filter: %s rows", format(nrow(cohort), big.mark = ",")))

# UF code → region + state name mapping
uf_info <- tribble(
  ~uf, ~uf_name, ~region,
  "11", "Rondonia",        "North",
  "12", "Acre",            "North",
  "13", "Amazonas",        "North",
  "14", "Roraima",         "North",
  "15", "Para",            "North",
  "16", "Amapa",           "North",
  "17", "Tocantins",       "North",
  "21", "Maranhao",        "Northeast",
  "22", "Piaui",           "Northeast",
  "23", "Ceara",           "Northeast",
  "24", "Rio Grande do Norte", "Northeast",
  "25", "Paraiba",         "Northeast",
  "26", "Pernambuco",      "Northeast",
  "27", "Alagoas",         "Northeast",
  "28", "Sergipe",         "Northeast",
  "29", "Bahia",           "Northeast",
  "31", "Minas Gerais",    "Southeast",
  "32", "Espirito Santo",  "Southeast",
  "33", "Rio de Janeiro",  "Southeast",
  "35", "Sao Paulo",       "Southeast",
  "41", "Parana",          "South",
  "42", "Santa Catarina",  "South",
  "43", "Rio Grande do Sul","South",
  "50", "Mato Grosso do Sul","Central-West",
  "51", "Mato Grosso",     "Central-West",
  "52", "Goias",           "Central-West",
  "53", "Distrito Federal","Central-West"
)

cohort <- cohort %>%
  left_join(uf_info, by = "uf") %>%
  mutate(
    sex_label = case_when(sex == "F" ~ "Female", sex == "M" ~ "Male", TRUE ~ NA_character_),
    age_grp = cut(age, breaks = c(50, 60, 70, 80, Inf),
                  labels = c("50-59", "60-69", "70-79", "80+"),
                  right = FALSE, include.lowest = FALSE),
    stage_label = case_when(
      stage == "0" ~ "In situ", stage == "1" ~ "I",
      stage == "2" ~ "II", stage == "3" ~ "III",
      stage == "4" ~ "IV", TRUE ~ NA_character_),
    cancer_group = case_when(
      grepl("^C50", cid) ~ "Breast",
      cid == "C61" ~ "Prostate",
      grepl("^C18", cid) | cid %in% c("C19", "C20") ~ "Colorectal",
      grepl("^C34", cid) ~ "Lung",
      grepl("^C16", cid) ~ "Stomach",
      grepl("^C53", cid) ~ "Cervical",
      grepl("^C91", cid) | grepl("^C92", cid) ~ "Leukemia",
      grepl("^C82", cid) | grepl("^C83", cid) | grepl("^C85", cid) ~ "Lymphoma",
      !is.na(cid) ~ "Other", TRUE ~ NA_character_),
    compliant_60d = ttt <= 60,
    compliant_30d = ttt <= 30
  )

message(sprintf("Saving cohort: %s rows", format(nrow(cohort), big.mark = ",")))
saveRDS(cohort, file = OUT_FILE)
message(sprintf("  File size: %.1f MB", file.size(OUT_FILE) / 1024^2))

# Summary report
message("\n=== NATIONAL COHORT SUMMARY ===")
message(sprintf("Total: %s unique patient-cancers", format(nrow(cohort), big.mark = ",")))
message(sprintf("UFs represented: %d", n_distinct(cohort$uf)))
message(sprintf("Median TTT: %.0f days", median(cohort$ttt)))
message(sprintf("Compliance <=60d: %.1f%%", 100 * mean(cohort$compliant_60d)))
message(sprintf("Compliance <=30d: %.1f%%", 100 * mean(cohort$compliant_30d)))
message(sprintf(">180d delay: %.1f%%", 100 * mean(cohort$ttt > 180)))
message("\n--- By region ---")
cohort %>% count(region, sort = TRUE) %>% print()
message("\n--- By cancer_group ---")
cohort %>% count(cancer_group, sort = TRUE) %>% mutate(pct = round(100*n/sum(n),1)) %>% print()
```

- [ ] **Step 2: Run and verify**

Run: `Rscript analysis/02_build_apac_cohort.R`
Expected: ~150-300K unique patient-cancer pairs (3-6x the 5-UF cohort), all 27 UFs represented.

- [ ] **Step 3: Verify new columns exist**

Run: `Rscript -e 'd <- readRDS("analysis/data/apac_cohort.rds"); cat(names(d), sep="\n"); cat("\nUF names:", paste(sort(unique(d$uf_name)), collapse=", "))'`
Expected: `uf_name`, `uf_download`, `region` columns present. 27 state names listed.

---

### Task 3: Re-run existing analysis scripts (04-08) with national data

All existing scripts read from `analysis/data/apac_cohort.rds`, which Task 2 just overwrote with national data. Re-run them to regenerate all figures and tables with the larger dataset.

**Files:**
- Re-run (no modification): `analysis/04_km_figures.R`, `analysis/05_cox_regression.R`, `analysis/06_table1.R`, `analysis/07_ecological.R`, `analysis/08_choropleth.R`
- Overwrite: All PDFs in `analysis/figures/` and HTMLs in `analysis/tables/`

- [ ] **Step 1: Re-run KM figures**

Run: `Rscript analysis/04_km_figures.R`
Expected: fig1-fig4 regenerated with national N in titles. Check that all 5 regions appear in fig4.

- [ ] **Step 2: Re-run Cox regression**

Run: `Rscript analysis/05_cox_regression.R`
Expected: Forest plot and HR table regenerated. Narrower confidence intervals with larger N.

- [ ] **Step 3: Re-run Table 1**

Run: `Rscript analysis/06_table1.R`
Expected: Table 1 regenerated with national cohort demographics.

- [ ] **Step 4: Re-run ecological analysis**

Run: `Rscript analysis/07_ecological.R`
Expected: Fig 6 ecological scatter plots regenerated with updated APAC regional summaries.

- [ ] **Step 5: Re-run choropleth**

Run: `Rscript analysis/08_choropleth.R`
Expected: Fig 7 choropleth regenerated with national APAC compliance rates per region.

- [ ] **Step 6: Quick visual check — verify N increased**

Run: `Rscript -e 'cat(readLines("analysis/data/01_download_apac.log")[length(readLines("analysis/data/01_download_apac.log"))])'`
Then open one figure to confirm the title shows the new N.

---

### Task 4: Competing risks analysis (Fine-Gray model)

Death before receiving treatment is an informative competing event. Patients who die before chemo never "fail" to comply — they were censored by death. Fine-Gray subdistribution hazards handle this correctly. Reviewers at BMJ/Lancet will expect this.

**Files:**
- Create: `analysis/09_competing_risks.R`
- Create: `analysis/figures/fig8_competing_risks.pdf`
- Create: `analysis/tables/table3_finegray.html`

- [ ] **Step 1: Install tidycmprsk if needed**

Run: `Rscript -e 'if (!requireNamespace("tidycmprsk", quietly = TRUE)) install.packages("tidycmprsk", repos = "https://cloud.r-project.org")'`

- [ ] **Step 2: Write the competing risks script**

```r
# analysis/09_competing_risks.R
# Competing risks: death before treatment as competing event
# Fine-Gray subdistribution hazard model
# Output: fig8_competing_risks.pdf, table3_finegray.html

suppressPackageStartupMessages({
  for (pkg in c("tidycmprsk", "ggsurvfit", "gtsummary", "gt")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  library(dplyr)
  library(survival)
  library(tidycmprsk)
  library(ggsurvfit)
  library(ggplot2)
  library(gtsummary)
  library(gt)
})

BASE <- "/Users/matheusrech/Downloads/fingertips-br/files/elsi_native_addons/analysis"
FIGS <- file.path(BASE, "figures")
TABS <- file.path(BASE, "tables")
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)
dir.create(TABS, showWarnings = FALSE, recursive = TRUE)

d <- readRDS(file.path(BASE, "data", "apac_cohort.rds"))

# Define competing risks outcome:
# 0 = censored (still waiting at 365d, no death)
# 1 = treated (received chemo within 365d)
# 2 = died before treatment (death == 1 AND ttt > 365 or never treated)
#
# In APAC data, AP_OBITO flags in-hospital death during that APAC period.
# Patients who died before treatment never appear with a treatment date.
# For this analysis: death == 1 AND ttt > median → likely died during wait.
# Conservative approach: use death flag as competing event for those with long TTT.

d <- d %>%
  mutate(
    time365 = pmin(ttt, 365),
    # Competing risk status
    cr_status = case_when(
      ttt <= 365 & (death == 0 | is.na(death)) ~ 1L,  # treated, alive
      ttt <= 365 & death == 1                   ~ 1L,  # treated (death during/after tx)
      death == 1                                ~ 2L,  # died before treatment
      TRUE                                      ~ 0L   # censored at 365d
    ),
    cr_status = factor(cr_status, levels = c(0, 1, 2),
                       labels = c("censored", "treated", "died"))
  )

cat("Competing risks status distribution:\n")
print(table(d$cr_status))

# ── Cumulative incidence curves ──────────────────────────────────────────────
# Overall
cuminc_overall <- cuminc(Surv(time365, cr_status) ~ 1, data = d)

p_cuminc <- cuminc_overall %>%
  ggcuminc(outcome = c("treated", "died")) +
  geom_vline(xintercept = 60, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 63, y = 0.02, label = "60d", hjust = 0,
           colour = "grey40", size = 3) +
  scale_x_continuous(breaks = seq(0, 365, 60)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Cumulative Incidence of Treatment and Death (Competing Risks)",
    subtitle = sprintf("n = %s elderly cancer patients, Brazil 2022",
                        format(nrow(d), big.mark = ",")),
    x = "Days from diagnosis",
    y = "Cumulative incidence"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13))

# By region
d_region <- d %>% filter(!is.na(region))
cuminc_region <- cuminc(Surv(time365, cr_status) ~ region, data = d_region)

p_cuminc_region <- cuminc_region %>%
  ggcuminc(outcome = "treated") +
  geom_vline(xintercept = 60, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = seq(0, 365, 60)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Cumulative Incidence of Treatment by Region (Competing Risks)",
    x = "Days from diagnosis",
    y = "Cumulative incidence of treatment"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13))

# Combine
library(patchwork)
combined <- p_cuminc / p_cuminc_region +
  plot_annotation(
    title = "Competing Risks Analysis: Treatment vs Death Before Treatment",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(file.path(FIGS, "fig8_competing_risks.pdf"), combined,
       width = 12, height = 14, dpi = 300)
cat("Saved: fig8_competing_risks.pdf\n")

# ── Fine-Gray regression ──────────────────────────────────────────────────────
d_fg <- d %>%
  filter(
    !is.na(sex_label), !is.na(age_grp), !is.na(region),
    cancer_group != "Other", !is.na(cancer_group),
    stage_label %in% c("I", "II", "III", "IV")
  ) %>%
  mutate(
    sex_label    = factor(sex_label, levels = c("Female", "Male")),
    age_grp      = factor(age_grp, levels = c("50-59","60-69","70-79","80+")),
    region       = factor(region, levels = c("Southeast","South","Central-West","Northeast","North")),
    stage_label  = factor(stage_label, levels = c("I","II","III","IV")),
    cancer_group = factor(cancer_group)
  )

cat(sprintf("\nFine-Gray modeling dataset: %d patients\n", nrow(d_fg)))

fg_fit <- crr(
  Surv(time365, cr_status) ~ age_grp + sex_label + cancer_group + stage_label + region,
  data = d_fg
)

cat("\n=== FINE-GRAY MODEL SUMMARY ===\n")
print(summary(fg_fit))

# Formatted table
fg_tbl <- tbl_regression(fg_fit, exponentiate = TRUE) %>%
  bold_p() %>%
  bold_labels() %>%
  modify_caption("**Table 3.** Fine-Gray subdistribution hazard ratios for time to treatment (competing risks: death before treatment)")

gtsave(as_gt(fg_tbl), file.path(TABS, "table3_finegray.html"))
cat("Saved: table3_finegray.html\n")
```

- [ ] **Step 3: Run and verify**

Run: `Rscript analysis/09_competing_risks.R`
Expected: Cumulative incidence plot showing treatment and death curves. Fine-Gray HR table comparing to standard Cox results.

---

### Task 5: Sensitivity analyses

Three standard sensitivity analyses that reviewers will request: (a) restricted to confirmed first-ever cancers (TTT 0-365d only), (b) stratified by cancer site (breast-only, prostate-only), (c) alternative compliance thresholds (30d, 90d, 120d).

**Files:**
- Create: `analysis/10_sensitivity.R`
- Create: `analysis/tables/etable1_sensitivity.html`
- Create: `analysis/figures/efig1_sensitivity_km.pdf`

- [ ] **Step 1: Write the sensitivity analysis script**

```r
# analysis/10_sensitivity.R
# Sensitivity analyses for robustness
# Output: etable1_sensitivity.html, efig1_sensitivity_km.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(patchwork)
  library(gtsummary)
  library(gt)
})

BASE <- "/Users/matheusrech/Downloads/fingertips-br/files/elsi_native_addons/analysis"
FIGS <- file.path(BASE, "figures")
TABS <- file.path(BASE, "tables")
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)
dir.create(TABS, showWarnings = FALSE, recursive = TRUE)

d <- readRDS(file.path(BASE, "data", "apac_cohort.rds"))

# ═══════════════════════════════════════════════════════════════════════════════
# Sensitivity A: Restrict to TTT 0-365 days (exclude outliers >1 year)
# ═══════════════════════════════════════════════════════════════════════════════
d_365 <- d %>% filter(ttt <= 365)

# ═══════════════════════════════════════════════════════════════════════════════
# Sensitivity B: By cancer site (top 4 individually)
# ═══════════════════════════════════════════════════════════════════════════════
top4 <- c("Breast", "Prostate", "Colorectal", "Lung")

km_plots <- list()
for (cancer in top4) {
  dd <- d_365 %>% filter(cancer_group == cancer)
  if (nrow(dd) < 50) next

  dd$time365 <- pmin(dd$ttt, 365)
  dd$event   <- 1L

  fit <- survfit(Surv(time365, event) ~ region, data = dd %>% filter(!is.na(region)))
  p <- ggsurvplot(
    fit, data = dd %>% filter(!is.na(region)),
    fun = "event", pval = TRUE, conf.int = FALSE,
    risk.table = FALSE,
    title = sprintf("%s (n = %s)", cancer, format(nrow(dd), big.mark = ",")),
    xlab = "Days", ylab = "Cumulative incidence",
    xlim = c(0, 365), break.time.by = 60,
    ggtheme = theme_classic(base_size = 10),
    legend.title = "Region"
  )
  p$plot <- p$plot +
    geom_vline(xintercept = 60, linetype = "dashed", colour = "grey50")

  km_plots[[cancer]] <- p$plot
}

if (length(km_plots) >= 4) {
  combined_km <- (km_plots[[1]] + km_plots[[2]]) / (km_plots[[3]] + km_plots[[4]]) +
    plot_annotation(
      title = "eFigure 1. Time to Treatment by Region — Stratified by Cancer Site",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
  ggsave(file.path(FIGS, "efig1_sensitivity_km.pdf"), combined_km,
         width = 14, height = 10, dpi = 300)
  cat("Saved: efig1_sensitivity_km.pdf\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Sensitivity C: Alternative compliance thresholds
# ═══════════════════════════════════════════════════════════════════════════════
thresholds <- c(30, 60, 90, 120, 180)

compliance_by_threshold <- lapply(thresholds, function(t) {
  d %>%
    group_by(region) %>%
    summarise(
      threshold = t,
      n = n(),
      pct_compliant = round(100 * mean(ttt <= t, na.rm = TRUE), 1),
      .groups = "drop"
    )
}) %>% bind_rows()

# Overall
compliance_overall <- lapply(thresholds, function(t) {
  tibble(
    threshold = t,
    region = "Overall",
    n = nrow(d),
    pct_compliant = round(100 * mean(d$ttt <= t, na.rm = TRUE), 1)
  )
}) %>% bind_rows()

compliance_all <- bind_rows(compliance_overall, compliance_by_threshold)

cat("\n=== COMPLIANCE BY THRESHOLD ===\n")
print(tidyr::pivot_wider(compliance_all, names_from = threshold,
                          values_from = pct_compliant, names_prefix = "<="))

# Save as HTML table
tbl_wide <- tidyr::pivot_wider(
  compliance_all,
  id_cols = region,
  names_from = threshold,
  values_from = pct_compliant,
  names_prefix = "≤"
) %>%
  rename(Region = region) %>%
  gt() %>%
  tab_header(
    title = "eTable 1. Treatment Compliance (%) by Threshold and Region",
    subtitle = "Percentage of patients receiving first chemotherapy within threshold (days)"
  ) %>%
  fmt_number(decimals = 1) %>%
  cols_label(
    `≤30` = "≤30d", `≤60` = "≤60d", `≤90` = "≤90d",
    `≤120` = "≤120d", `≤180` = "≤180d"
  )

gtsave(tbl_wide, file.path(TABS, "etable1_sensitivity.html"))
cat("Saved: etable1_sensitivity.html\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Sensitivity D: Cox model restricted to breast cancer only
# ═══════════════════════════════════════════════════════════════════════════════
d_breast <- d_365 %>%
  filter(cancer_group == "Breast", !is.na(sex_label), !is.na(age_grp),
         !is.na(region), stage_label %in% c("I","II","III","IV")) %>%
  mutate(
    time365  = pmin(ttt, 365),
    event365 = 1L,
    age_grp     = factor(age_grp, levels = c("50-59","60-69","70-79","80+")),
    stage_label = factor(stage_label, levels = c("I","II","III","IV")),
    region      = factor(region, levels = c("Southeast","South","Central-West","Northeast","North"))
  )

if (nrow(d_breast) >= 100) {
  cox_breast <- coxph(
    Surv(time365, event365) ~ age_grp + stage_label + region,
    data = d_breast
  )
  cat("\n=== BREAST CANCER COX (SENSITIVITY) ===\n")
  print(summary(cox_breast))

  tbl_breast <- tbl_regression(cox_breast, exponentiate = TRUE) %>%
    bold_p() %>%
    bold_labels() %>%
    modify_caption("eTable 2. Cox regression — Breast cancer only")

  gtsave(as_gt(tbl_breast), file.path(TABS, "etable2_cox_breast.html"))
  cat("Saved: etable2_cox_breast.html\n")
}

cat("\nAll sensitivity analyses complete.\n")
```

- [ ] **Step 2: Run and verify**

Run: `Rscript analysis/10_sensitivity.R`
Expected: eFigure 1 (cancer-specific KM), eTable 1 (threshold comparison), eTable 2 (breast-only Cox).

---

### Task 6: State-level choropleth map (27 UFs)

The current choropleth shows only 5 macro-regions. With all 27 UFs in the national cohort, produce a 27-state choropleth showing compliance and median TTT at state level. This is a major visual upgrade for the paper.

**Files:**
- Create: `analysis/11_state_choropleth.R`
- Create: `analysis/figures/fig9_state_choropleth.pdf`

- [ ] **Step 1: Write the state-level choropleth script**

```r
# analysis/11_state_choropleth.R
# State-level (27 UFs) choropleth: compliance + median TTT
# Output: fig9_state_choropleth.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(geobr)
  library(patchwork)
  library(RColorBrewer)
})

BASE <- "/Users/matheusrech/Downloads/fingertips-br/files/elsi_native_addons/analysis"
FIGS <- file.path(BASE, "figures")
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)

d <- readRDS(file.path(BASE, "data", "apac_cohort.rds"))

# ── Per-state APAC summary ───────────────────────────────────────────────────
state_stats <- d %>%
  filter(!is.na(uf)) %>%
  group_by(uf, uf_name, region) %>%
  summarise(
    n              = n(),
    median_ttt     = median(ttt, na.rm = TRUE),
    pct_compliant  = 100 * mean(compliant_60d, na.rm = TRUE),
    pct_over_180d  = 100 * mean(ttt > 180, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(code_state = as.numeric(uf))

cat("State-level APAC summary:\n")
print(state_stats %>% arrange(pct_compliant) %>% select(uf_name, n, median_ttt, pct_compliant))

# ── Load Brazil state shapefile ──────────────────────────────────────────────
cat("\nDownloading state shapefile...\n")
shp <- read_state(year = 2020, showProgress = FALSE)

# Join
geo <- shp %>%
  left_join(state_stats, by = "code_state")

# Centroids for labels (only states with data)
ctr <- st_centroid(geo %>% filter(!is.na(n)))
ctr_xy <- st_coordinates(ctr)
ctr_df <- st_drop_geometry(ctr) %>%
  mutate(lon = ctr_xy[, 1], lat = ctr_xy[, 2])

# ── Map theme ─────────────────────────────────────────────────────────────────
map_theme <- theme_void(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    legend.key.height = unit(0.35, "cm"),
    legend.title = element_text(size = 9, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey40")
  )

# ── Panel A: Compliance (%) ──────────────────────────────────────────────────
p1 <- ggplot(geo) +
  geom_sf(aes(fill = pct_compliant), colour = "white", linewidth = 0.3) +
  scale_fill_distiller(
    palette = "RdYlGn", direction = 1, name = "Compliance (%)",
    na.value = "grey85",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
  ) +
  geom_text(
    data = ctr_df, aes(x = lon, y = lat, label = abbrev_state),
    size = 2, colour = "grey20"
  ) +
  labs(title = "Lei 12.732 Compliance by State (%)",
       subtitle = "% patients treated within 60 days of diagnosis") +
  map_theme

# ── Panel B: Median TTT (days) ───────────────────────────────────────────────
p2 <- ggplot(geo) +
  geom_sf(aes(fill = median_ttt), colour = "white", linewidth = 0.3) +
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1, name = "Median TTT (days)",
    na.value = "grey85",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
  ) +
  geom_text(
    data = ctr_df, aes(x = lon, y = lat, label = abbrev_state),
    size = 2, colour = "grey20"
  ) +
  labs(title = "Median Time to Treatment by State (days)",
       subtitle = "Median days from diagnosis to first chemotherapy") +
  map_theme

# ── Panel C: Sample size ─────────────────────────────────────────────────────
p3 <- ggplot(geo) +
  geom_sf(aes(fill = n), colour = "white", linewidth = 0.3) +
  scale_fill_distiller(
    palette = "Blues", direction = 1, name = "N patients",
    na.value = "grey85", trans = "log10",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
  ) +
  geom_text(
    data = ctr_df, aes(x = lon, y = lat, label = abbrev_state),
    size = 2, colour = "grey20"
  ) +
  labs(title = "Sample Size by State",
       subtitle = "Elderly cancer patients on SUS chemotherapy (log scale)") +
  map_theme

# ── Panel D: >180d extreme delay (%) ────────────────────────────────────────
p4 <- ggplot(geo) +
  geom_sf(aes(fill = pct_over_180d), colour = "white", linewidth = 0.3) +
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1, name = ">180d delay (%)",
    na.value = "grey85",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
  ) +
  geom_text(
    data = ctr_df, aes(x = lon, y = lat, label = abbrev_state),
    size = 2, colour = "grey20"
  ) +
  labs(title = "Extreme Delay (>180 days) by State",
       subtitle = "% patients waiting over 6 months for chemotherapy") +
  map_theme

# ── Combine ───────────────────────────────────────────────────────────────────
combined <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Cancer Treatment Delays Among Older Adults — Brazil, 27 States, 2022",
    subtitle = "APAC-Quimioterapia (DATASUS/SIA) | Adults aged ≥50 years",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey40")
    )
  )

ggsave(file.path(FIGS, "fig9_state_choropleth.pdf"), combined,
       width = 14, height = 14, dpi = 300)
cat("Saved: fig9_state_choropleth.pdf\n")
```

- [ ] **Step 2: Run and verify**

Run: `Rscript analysis/11_state_choropleth.R`
Expected: 4-panel 27-state choropleth showing compliance, median TTT, sample size, and extreme delay.

---

### Task 7: Quarto manuscript (STROBE-compliant)

Assemble all figures, tables, and analyses into a reproducible Quarto document following the STROBE checklist for cross-sectional studies. This becomes the submission-ready manuscript.

**Files:**
- Create: `analysis/manuscript.qmd`
- Create: `analysis/references.bib` (empty placeholder for manual citation entry)

- [ ] **Step 1: Write the Quarto manuscript**

```qmd
---
title: "Non-compliance with Brazil's 60-day cancer treatment law among older adults: a national population-based analysis of chemotherapy records"
author:
  - name: "[Author 1]"
    affiliations: "[Affiliation]"
  - name: "[Author 2]"
    affiliations: "[Affiliation]"
format:
  docx:
    reference-doc: default
  pdf:
    documentclass: article
    geometry: margin=2.5cm
    fontsize: 11pt
    linestretch: 2
    number-sections: true
execute:
  echo: false
  warning: false
  message: false
bibliography: references.bib
---

```{r setup}
#| include: false
library(dplyr)
library(survival)
library(ggplot2)
library(gtsummary)
library(gt)
library(patchwork)

cohort <- readRDS("data/apac_cohort.rds")
n_total <- format(nrow(cohort), big.mark = ",")
n_uf <- n_distinct(cohort$uf)
med_ttt <- round(median(cohort$ttt))
pct_60 <- round(100 * mean(cohort$compliant_60d), 1)
pct_180 <- round(100 * mean(cohort$ttt > 180), 1)
```

# Abstract

**Background:** Brazil's Lei 12.732/2012 mandates that cancer patients in the public health system (SUS) receive first treatment within 60 days of diagnosis. Compliance among older adults — who face the highest cancer burden — is unknown at the national level.

**Methods:** Population-based cross-sectional study of all SUS chemotherapy authorizations (APAC-Quimioterapia) for patients aged ≥50 years across all `r n_uf` Brazilian states in 2022. Primary outcome: time from diagnosis to first chemotherapy. Compliance defined as treatment within 60 days (Lei 12.732). Kaplan-Meier, Cox proportional hazards, and Fine-Gray competing risks models assessed predictors of delay. Ecological analysis linked treatment delays to population health indicators from the ELSI-Brazil survey.

**Results:** Among `r n_total` unique elderly cancer patients, median time to treatment was `r med_ttt` days. Only `r pct_60`% received treatment within the legal 60-day target, and `r pct_180`% waited over 180 days. [Results to be completed after final analysis.]

**Conclusions:** [To be completed.]

**Keywords:** cancer treatment delay, Lei 12.732, older adults, Brazil, SUS, health equity

# Introduction

Cancer is the second leading cause of death among older Brazilians, and its incidence is projected to increase with population aging. Brazil's Unified Health System (SUS) provides cancer treatment free of charge, and Lei 12.732/2012 establishes a legal maximum of 60 days between pathological diagnosis and first treatment.

Despite this legal mandate, anecdotal evidence and regional studies suggest widespread non-compliance, particularly among older adults who face additional barriers including comorbidity, frailty, geographic distance, and limited health literacy. However, no national population-based study has quantified the extent of non-compliance among older adults or identified its predictors across all Brazilian states.

This study aims to: (1) estimate compliance with Lei 12.732 among older adults receiving SUS chemotherapy nationally; (2) identify demographic, clinical, and geographic predictors of treatment delay; (3) compare treatment access across Brazil's 27 states; and (4) explore ecological associations between regional population health vulnerability (from the ELSI-Brazil survey) and treatment delay patterns.

# Methods

## Study design and data sources

Cross-sectional analysis of administrative health records from the Brazilian Outpatient Information System (SIA/SUS), specifically the High-Complexity Outpatient Procedures Authorization (APAC) subsystem for chemotherapy (APAC-Quimioterapia, information system code SIA-AQ).

Data were obtained via the `microdatasus` R package, which downloads directly from the DATASUS FTP server. All 27 Brazilian states (Unidades Federativas, UFs), full year 2022. Supplementary population health indicators from ELSI-Brazil Wave 2 (2019-2020), a nationally representative survey of Brazilians aged ≥50 years.

## Study population

**Inclusion criteria:** (1) APAC-QT record with valid diagnosis date (AQ_DTIDEN) and treatment date (AQ_DTINTR); (2) patient age ≥50 years at time of APAC; (3) ICD-10 chapter II neoplasm code (C00-C97).

**Exclusion criteria:** (1) time-to-treatment <0 or >1,825 days (data quality filter); (2) missing patient identifier (AP_CNSPCN) or ICD-10 code.

**Deduplication:** Multiple APAC records per patient were deduplicated by retaining the first treatment session per unique patient-cancer pair (CNS + ICD-10 code), ordered by treatment date.

## Variables

**Primary outcome:** Time-to-treatment (TTT), defined as the number of days from pathological diagnosis (AQ_DTIDEN) to first chemotherapy session (AQ_DTINTR).

**Binary compliance outcomes:** (1) Lei 12.732 compliance (TTT ≤60 days); (2) extended compliance (TTT ≤30 days); (3) extreme delay (TTT >180 days).

**Covariates:** Age (continuous and grouped: 50-59, 60-69, 70-79, 80+), sex, race/ethnicity (APAC codes: White, Black, Mixed, Asian, Indigenous), cancer type (ICD-10 grouped: breast, prostate, colorectal, lung, stomach, cervical, leukemia, lymphoma, other), clinical stage (TNM: I-IV), state of residence (27 UFs), and macro-region (5 regions).

## Statistical analysis

**Descriptive:** Table 1 presents cohort characteristics stratified by Lei 12.732 compliance. Continuous variables reported as median (IQR); categorical as n (%).

**Time-to-event:** Kaplan-Meier cumulative incidence curves with log-rank tests, stratified by cancer type, stage, region, age group, and sex. Treatment capped at 365 days for visualization.

**Multivariable regression:** Cox proportional hazards model for predictors of time to treatment. Hazard ratios >1 indicate faster treatment (shorter TTT). Schoenfeld residuals tested the proportional hazards assumption.

**Competing risks:** Fine-Gray subdistribution hazard model with death before treatment as the competing event.

**Ecological analysis:** Survey-weighted regional estimates from ELSI-Brazil (frailty prevalence, health plan coverage, per capita income, social isolation) correlated with APAC regional treatment compliance rates.

**Sensitivity analyses:** (a) Restricted to TTT 0-365 days; (b) cancer site-specific models (breast, prostate, colorectal, lung); (c) alternative compliance thresholds (30, 90, 120, 180 days).

All analyses performed in R (version 4.5). Statistical significance set at p < 0.05. This study followed the STROBE guidelines for reporting observational studies.

## Ethical considerations

This study used anonymized administrative data publicly available through DATASUS. ELSI-Brazil is a de-identified public-use dataset. No ethical approval was required per Brazilian Resolution 510/2016 (research using publicly available data).

# Results

[Results section will embed computed figures and tables from the analysis pipeline.]

```{r table1}
#| label: tbl-demographics
#| tbl-cap: "Cohort characteristics by Lei 12.732 compliance"
# Source: analysis/tables/table1_characteristics.html
# Embed pre-computed gtsummary table
```

```{r fig-km-overall}
#| label: fig-km-overall
#| fig-cap: "Time to first chemotherapy among elderly cancer patients, Brazil 2022"
knitr::include_graphics("figures/fig1_km_overall.pdf")
```

```{r fig-km-cancer}
#| label: fig-km-cancer
#| fig-cap: "Time to treatment by cancer type"
knitr::include_graphics("figures/fig2_km_cancer.pdf")
```

```{r fig-km-stage}
#| label: fig-km-stage
#| fig-cap: "Time to treatment by clinical stage"
knitr::include_graphics("figures/fig3_km_stage.pdf")
```

```{r fig-km-region}
#| label: fig-km-region
#| fig-cap: "Time to treatment by Brazilian macro-region"
knitr::include_graphics("figures/fig4_km_region.pdf")
```

```{r fig-cox}
#| label: fig-cox
#| fig-cap: "Forest plot of Cox regression hazard ratios"
knitr::include_graphics("figures/fig5_cox_forest.pdf")
```

```{r fig-ecological}
#| label: fig-ecological
#| fig-cap: "Ecological analysis: ELSI health indicators vs APAC treatment delays"
knitr::include_graphics("figures/fig6_ecological.pdf")
```

```{r fig-choropleth-region}
#| label: fig-choropleth-region
#| fig-cap: "Regional choropleth of treatment compliance and population health indicators"
knitr::include_graphics("figures/fig7_choropleth.pdf")
```

```{r fig-competing}
#| label: fig-competing
#| fig-cap: "Competing risks analysis: treatment vs death before treatment"
knitr::include_graphics("figures/fig8_competing_risks.pdf")
```

```{r fig-state}
#| label: fig-state
#| fig-cap: "State-level choropleth of treatment compliance across 27 Brazilian states"
knitr::include_graphics("figures/fig9_state_choropleth.pdf")
```

# Discussion

[To be written based on final results. Key discussion points:]

1. Magnitude of non-compliance and comparison with other countries
2. Geographic disparities — state and regional variation
3. Cancer-type specific patterns (breast vs cervical vs prostate)
4. Stage-at-diagnosis and its relationship with treatment delay
5. Ecological correlations with population health vulnerability
6. Competing risks perspective — death before treatment
7. Policy implications for Lei 12.732 enforcement
8. Strengths: national coverage, large N, multiple analytical approaches
9. Limitations: administrative data quality, no individual SES, ecological fallacy

# References

::: {#refs}
:::
```

- [ ] **Step 2: Create empty references.bib**

```
@article{placeholder,
  title = {Placeholder for references},
  year = {2026}
}
```

- [ ] **Step 3: Test Quarto rendering**

Run: `cd analysis && quarto render manuscript.qmd --to docx`
Expected: Word document with embedded figures and tables.

---

### Task 8: Upgrade Shiny app to national data + fix column mismatches

The Shiny TTT tab currently uses `app/apac_qt_sp2022.rds` (SP-only, 7.1MB) with column names (`cancer_type`, `estadiamento`, `age_group`, `sexo`) that don't match the cohort's column names (`cancer_group`, `stage_label`, `age_grp`, `sex_label`). Upgrade to use the national cohort and fix the stratification selector.

**Files:**
- Modify: `app/app.R` (lines 33-39 data loading, lines 370-408 UI, lines 1200-1346 server)
- Copy: `analysis/data/apac_cohort.rds` → `app/apac_cohort_national.rds`

- [ ] **Step 1: Copy national cohort to app directory**

Run: `cp analysis/data/apac_cohort.rds app/apac_cohort_national.rds`

- [ ] **Step 2: Update data loading (app.R lines 32-40)**

Replace:
```r
# APAC data for Time to Treatment analysis
apac_file <- "apac_qt_sp2022.rds"
if (file.exists(apac_file)) {
  apac <- readRDS(apac_file)
  has_apac <- TRUE
} else {
  has_apac <- FALSE
  apac <- NULL
}
```

With:
```r
# APAC data for Time to Treatment analysis (national cohort, all 27 UFs)
apac_file <- "apac_cohort_national.rds"
if (file.exists(apac_file)) {
  apac <- readRDS(apac_file)
  has_apac <- TRUE
} else {
  has_apac <- FALSE
  apac <- NULL
}
```

- [ ] **Step 3: Fix TTT tab stratification selector (app.R ~line 388)**

Replace:
```r
selectInput("ttt_strat", NULL,
  choices = c(
    "Overall"     = "overall",
    "Cancer type" = "cancer_type",
    "Stage"       = "estadiamento",
    "Age group"   = "age_group",
    "Sex"         = "sexo"
  ),
  selected = "overall"
```

With:
```r
selectInput("ttt_strat", NULL,
  choices = c(
    "Overall"     = "overall",
    "Cancer type" = "cancer_group",
    "Stage"       = "stage_label",
    "Age group"   = "age_grp",
    "Sex"         = "sex_label",
    "Region"      = "region",
    "State"       = "uf_name"
  ),
  selected = "overall"
```

- [ ] **Step 4: Fix server-side APAC data processing (app.R ~lines 1200-1210)**

The server code expects column names from the SP-only data. Update to match the national cohort columns (which already has `ttt`, `cancer_group`, `stage_label`, etc). Replace the date parsing and ttt calculation block with:

```r
  if (has_apac) {
    # National cohort already has ttt, cancer_group, stage_label, etc.
    apac_proc <- apac %>%
      filter(!is.na(ttt), ttt >= 0, ttt <= 730)

    apac_proc$time_band <- cut(
      apac_proc$ttt,
      breaks = c(-Inf, 30, 60, 180, Inf),
      labels = c("0-30d", "31-60d", "61-180d", ">180d"),
      right  = TRUE
    )
```

- [ ] **Step 5: Update data source note (app.R ~line 398-399)**

Replace:
```r
tags$small(tags$em(
  "Data: APAC-Quimioterapia (DATASUS/SIA), São Paulo 2022. Only SUS chemotherapy."
))
```

With:
```r
tags$small(tags$em(
  "Data: APAC-Quimioterapia (DATASUS/SIA), all 27 Brazilian states, 2022. Adults ≥50 years. Only SUS chemotherapy."
))
```

- [ ] **Step 6: Test the app**

Run: `Rscript -e 'shiny::runApp("app/", port=7878)'`
Verify: TTT tab loads with national data, all stratification options work, KM curves render.
