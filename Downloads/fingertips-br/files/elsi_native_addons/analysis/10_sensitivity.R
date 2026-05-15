# 10_sensitivity.R
# Sensitivity analyses for the APAC cohort study
# Output:
#   analysis/figures/efig1_sensitivity_km.pdf   — 4-panel cancer-specific KM by region
#   analysis/tables/etable1_sensitivity.html    — compliance by threshold table
#   analysis/tables/etable2_cox_breast.html     — breast-only Cox regression table

suppressPackageStartupMessages({
  for (pkg in c("survminer", "patchwork", "gtsummary", "gt", "tidyr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  library(survival)
  library(survminer)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(gtsummary)
  library(gt)
})

# ── paths ──────────────────────────────────────────────────────────────────────
BASE <- "/Users/matheusrech/Downloads/fingertips-br/files/elsi_native_addons/analysis"
DATA <- file.path(BASE, "data", "apac_cohort.rds")
FIGS <- file.path(BASE, "figures")
TABS <- file.path(BASE, "tables")
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)
dir.create(TABS, showWarnings = FALSE, recursive = TRUE)

# ── shared constants ───────────────────────────────────────────────────────────
REGION_LEVELS  <- c("North", "Northeast", "Southeast", "South", "Central-West")
REGION_PALETTE <- c(
  "North"        = "#E74C3C",
  "Northeast"    = "#F39C12",
  "Southeast"    = "#3498DB",
  "South"        = "#27AE60",
  "Central-West" = "#8E44AD"
)

# ── load data ──────────────────────────────────────────────────────────────────
d_raw <- readRDS(DATA)

d <- d_raw |>
  mutate(
    time365  = pmin(ttt, 365),
    event365 = as.integer(ttt <= 365),
    region   = factor(region, levels = REGION_LEVELS)
  )

cat(sprintf("Full cohort: %d patients | median TTT = %.0f days\n",
            nrow(d), median(d$ttt, na.rm = TRUE)))

# ══════════════════════════════════════════════════════════════════════════════
# Sensitivity A — Restrict to TTT 0–365 days (exclude outliers)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Sensitivity A: TTT 0–365 days ──\n")

d_A <- d |> filter(ttt >= 0, ttt <= 365)

cat(sprintf(
  "After restricting to 0-365 days: %d patients (excluded %d; %.1f%%)\n",
  nrow(d_A),
  nrow(d) - nrow(d_A),
  100 * (nrow(d) - nrow(d_A)) / nrow(d)
))

cat(sprintf(
  "Restricted cohort: median TTT = %.0f days; mean = %.0f days\n",
  median(d_A$ttt, na.rm = TRUE),
  mean(d_A$ttt, na.rm = TRUE)
))

# ══════════════════════════════════════════════════════════════════════════════
# Sensitivity B — Cancer-specific KM by region (top 4)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Sensitivity B: Cancer-specific KM by region ──\n")

TOP4 <- c("Breast", "Prostate", "Colorectal", "Lung")

# Helper: build one KM panel (ggplot2 object) for a single cancer site
make_cancer_km_panel <- function(cancer, data) {

  dd <- data |>
    filter(cancer_group == cancer, region %in% REGION_LEVELS) |>
    mutate(region = droplevels(factor(region, levels = REGION_LEVELS)))

  n_cancer <- nrow(dd)
  if (n_cancer < 10) {
    # Fallback empty panel if too few observations
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = sprintf("%s\n(n < 10)", cancer), size = 5) +
        theme_void()
    )
  }

  fit <- survfit(Surv(time365, event) ~ region,
                 data = dd |> mutate(event = 1L))

  # Extract cumulative incidence from survfit (fun = "event")
  km_data <- surv_summary(fit, data = dd)
  km_data$cumhaz_pct <- (1 - km_data$surv) * 100   # cumulative incidence %

  # Map strata to region labels
  km_data$region_label <- sub("region=", "", km_data$strata)
  km_data$region_label <- factor(km_data$region_label, levels = REGION_LEVELS)

  palette_sub <- REGION_PALETTE[levels(km_data$region_label)]

  p <- ggplot(km_data, aes(x = time, y = cumhaz_pct,
                            colour = region_label,
                            group  = region_label)) +
    geom_step(linewidth = 0.7) +
    geom_vline(xintercept = 60, linetype = "dashed",
               colour = "grey30", linewidth = 0.5) +
    annotate("text", x = 63, y = 3,
             label = "60d", hjust = 0, size = 2.8, colour = "grey30") +
    scale_colour_manual(values = palette_sub, name = "Region") +
    scale_x_continuous(limits = c(0, 365),
                       breaks = c(0, 60, 120, 180, 240, 300, 365)) +
    scale_y_continuous(limits = c(0, 100),
                       labels = function(x) paste0(x, "%")) +
    labs(
      title    = cancer,
      subtitle = sprintf("n = %s", format(n_cancer, big.mark = ",")),
      x        = "Days to treatment",
      y        = "Cumulative incidence"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title    = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(colour = "grey50", size = 9),
      legend.position = "right",
      legend.key.size  = unit(0.45, "cm"),
      legend.text      = element_text(size = 8),
      axis.text.x      = element_text(size = 8, angle = 30, hjust = 1),
      axis.text.y      = element_text(size = 8)
    )

  p
}

panels <- lapply(TOP4, make_cancer_km_panel, data = d)
names(panels) <- TOP4

# Assemble 4-panel figure using patchwork
# Collect common legend from first panel with a legend
combined_km <- (panels[[1]] | panels[[2]]) /
               (panels[[3]] | panels[[4]]) +
  plot_annotation(
    title   = "eFigure 1 — Cumulative Incidence of Treatment by Region",
    subtitle = "Top 4 cancer sites, restricted to 0–365 days; 60-day reference (dashed)",
    caption  = "KM curves show cumulative incidence of treatment (1 – survival function)\nAll patients treated by definition (event = 1)",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "grey40", size = 10),
      plot.caption  = element_text(colour = "grey50", size = 8)
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

efig1_path <- file.path(FIGS, "efig1_sensitivity_km.pdf")
ggsave(efig1_path, plot = combined_km, width = 14, height = 10,
       units = "in", dpi = 300, device = "pdf")
cat("Saved:", efig1_path, "\n")

# ══════════════════════════════════════════════════════════════════════════════
# Sensitivity C — Compliance by alternative thresholds (30/60/90/120/180 days)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Sensitivity C: Compliance by threshold table ──\n")

THRESHOLDS <- c(30, 60, 90, 120, 180)

# Compute compliance for every region + overall at each threshold
compute_compliance <- function(data, thresholds) {

  # Per region
  region_rows <- data |>
    filter(!is.na(region)) |>
    group_by(region) |>
    summarise(
      n = n(),
      !!!setNames(
        lapply(thresholds, function(th) {
          expr(round(mean(ttt <= !!th, na.rm = TRUE) * 100, 1))
        }),
        paste0("d", thresholds)
      ),
      .groups = "drop"
    ) |>
    mutate(Group = as.character(region))

  # Overall row
  overall_row <- data |>
    summarise(
      n = n(),
      !!!setNames(
        lapply(thresholds, function(th) {
          expr(round(mean(ttt <= !!th, na.rm = TRUE) * 100, 1))
        }),
        paste0("d", thresholds)
      )
    ) |>
    mutate(Group = "Overall", region = NA_character_)

  bind_rows(region_rows, overall_row) |>
    select(Group, n, everything(), -any_of("region"))
}

comp_tbl <- compute_compliance(d, THRESHOLDS)

cat("\nCompliance table preview:\n")
print(comp_tbl)

# Format as gt table
col_labels <- setNames(
  paste0("≤ ", THRESHOLDS, "d (%)"),
  paste0("d", THRESHOLDS)
)

gt_comp <- comp_tbl |>
  gt(rowname_col = "Group") |>
  tab_header(
    title    = md("**eTable 1 — Compliance Rate by Time Threshold and Region**"),
    subtitle = md("Percentage of patients starting treatment within each threshold")
  ) |>
  cols_label(
    n   = md("**N**"),
    !!!lapply(col_labels, md)
  ) |>
  fmt_number(columns = "n", decimals = 0, use_seps = TRUE) |>
  fmt_number(columns = paste0("d", THRESHOLDS), decimals = 1) |>
  cols_align(align = "center", columns = -Group) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_stub(rows = "Overall")
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#F0F4F8"),
      cell_text(weight = "bold")
    ),
    locations = cells_stub(rows = "Overall")
  ) |>
  tab_style(
    style = cell_fill(color = "#F0F4F8"),
    locations = cells_body(rows = Group == "Overall")
  ) |>
  # Highlight the 60-day column (main study threshold)
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(columns = "d60")
  ) |>
  tab_style(
    style = cell_borders(sides = "left", color = "#2E86AB", weight = px(2)),
    locations = cells_body(columns = "d60")
  ) |>
  tab_style(
    style = cell_borders(sides = "left", color = "#2E86AB", weight = px(2)),
    locations = cells_column_labels(columns = "d60")
  ) |>
  tab_footnote(
    footnote = "Bold column = primary study threshold (Lei 12.732/2013, 60 days)",
    locations = cells_column_labels(columns = "d60")
  ) |>
  tab_footnote(
    footnote = "Compliance = proportion of patients whose TTT ≤ threshold",
    locations = cells_title("subtitle")
  ) |>
  opt_table_font(font = google_font("Source Sans Pro")) |>
  opt_row_striping()

etable1_path <- file.path(TABS, "etable1_sensitivity.html")
gtsave(gt_comp, filename = etable1_path)
cat("Saved:", etable1_path, "\n")

# ══════════════════════════════════════════════════════════════════════════════
# Sensitivity D — Breast cancer-only Cox regression
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Sensitivity D: Breast cancer-only Cox regression ──\n")

d_breast <- d |>
  filter(
    cancer_group == "Breast",
    !is.na(age_grp),
    !is.na(stage_label),
    !is.na(region),
    stage_label %in% c("In situ", "I", "II", "III", "IV")
  ) |>
  mutate(
    age_grp     = relevel(factor(age_grp),     ref = "50-59"),
    stage_label = relevel(factor(stage_label), ref = "I"),
    region      = relevel(factor(region),       ref = "Southeast")
  )

cat(sprintf("Breast cancer subset: %d patients\n", nrow(d_breast)))
cat(sprintf("Event rate (treated ≤365d): %.1f%%\n",
            mean(d_breast$event365, na.rm = TRUE) * 100))

cat("\nFactor reference levels:\n")
cat("  age_grp    :", levels(d_breast$age_grp)[1],     "\n")
cat("  stage_label:", levels(d_breast$stage_label)[1], "\n")
cat("  region     :", levels(d_breast$region)[1],      "\n")

cox_breast <- coxph(
  Surv(time365, event365) ~ age_grp + stage_label + region,
  data = d_breast
)

cat("\nBreast-only Cox model summary:\n")
print(summary(cox_breast))

# Schoenfeld test
ph_breast <- cox.zph(cox_breast)
cat("\nSchoenfeld proportional hazards test:\n")
print(ph_breast)

# gtsummary table
cox_breast_tbl <- tbl_regression(
  cox_breast,
  exponentiate = TRUE,
  label = list(
    age_grp     ~ "Age group",
    stage_label ~ "Clinical stage",
    region      ~ "Region"
  )
) |>
  bold_p(t = 0.05) |>
  bold_labels() |>
  modify_header(label = "**Characteristic**") |>
  modify_caption(
    caption = paste0(
      "**eTable 2 — Breast Cancer-Only Cox Proportional Hazards Model**  \n",
      sprintf("n = %s patients | Reference: Age 50-59, Stage I, Southeast",
              format(nrow(d_breast), big.mark = ","))
    )
  ) |>
  modify_footnote(
    estimate = "HR = hazard ratio; higher HR = faster time to treatment",
    ci       = "95% confidence interval"
  )

etable2_path <- file.path(TABS, "etable2_cox_breast.html")
gt::gtsave(as_gt(cox_breast_tbl), filename = etable2_path)
cat("Saved:", etable2_path, "\n")

# ── Session summary ────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════\n")
cat("Sensitivity analyses complete.\n")
cat("  Sensitivity A: descriptive restriction to TTT 0-365d\n")
cat(sprintf("  Sensitivity B: 4-panel KM -> %s\n",
            file.path("figures", "efig1_sensitivity_km.pdf")))
cat(sprintf("  Sensitivity C: threshold table -> %s\n",
            file.path("tables", "etable1_sensitivity.html")))
cat(sprintf("  Sensitivity D: breast Cox -> %s\n",
            file.path("tables", "etable2_cox_breast.html")))
cat("══════════════════════════════════════════════════════════\n")
