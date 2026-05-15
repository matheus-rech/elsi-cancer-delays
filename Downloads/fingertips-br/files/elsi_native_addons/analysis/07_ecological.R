# =============================================================================
# Task 7 — Ecological analysis: ELSI health indicators × APAC treatment delays
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# 1. Load ELSI regional summary
# -----------------------------------------------------------------------------
elsi <- readRDS("analysis/data/elsi_regional.rds")

# -----------------------------------------------------------------------------
# 2. Load APAC cohort and compute per-region summary
# -----------------------------------------------------------------------------
apac_raw <- readRDS("analysis/data/apac_cohort.rds")

apac <- apac_raw |>
  group_by(region) |>
  summarise(
    n_apac        = n(),
    median_ttt    = median(ttt, na.rm = TRUE),
    pct_compliant = mean(compliant_60d, na.rm = TRUE) * 100,
    pct_over_180d = mean(ttt > 180,    na.rm = TRUE) * 100,
    .groups = "drop"
  ) |>
  rename(region_label = region)   # align join key with elsi

# -----------------------------------------------------------------------------
# 3. Join
# -----------------------------------------------------------------------------
merged <- left_join(elsi, apac, by = "region_label")

cat("\n=== Merged ecological table ===\n")
print(
  merged |>
    select(region_label, n_apac, median_ttt, pct_compliant, pct_over_180d,
           plan_prev, frail_prev, mean_income, soc_low_prev)
)

# -----------------------------------------------------------------------------
# 4–6. Four scatter plots
# -----------------------------------------------------------------------------

plot_scatter <- function(data, xvar, yvar, xlabel, ylabel, title) {
  ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(size = 4, color = "#E74C3C") +
    geom_text(aes(label = region_label), vjust = -1, size = 3.5) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed",
                color = "#1F3864") +
    labs(x = xlabel, y = ylabel, title = title) +
    theme_minimal(base_size = 13)
}

p1 <- plot_scatter(
  merged,
  xvar   = "plan_prev",
  yvar   = "pct_compliant",
  xlabel = "Insurance coverage (plan_prev)",
  ylabel = "Treatment compliance (%)",
  title  = "Insurance coverage vs treatment compliance"
)

p2 <- plot_scatter(
  merged,
  xvar   = "frail_prev",
  yvar   = "median_ttt",
  xlabel = "Frailty prevalence (frail_prev)",
  ylabel = "Median treatment delay (days)",
  title  = "Frailty burden vs treatment delay"
)

p3 <- plot_scatter(
  merged,
  xvar   = "mean_income",
  yvar   = "pct_compliant",
  xlabel = "Mean household income (BRL)",
  ylabel = "Treatment compliance (%)",
  title  = "Income vs treatment compliance"
)

p4 <- plot_scatter(
  merged,
  xvar   = "soc_low_prev",
  yvar   = "pct_over_180d",
  xlabel = "Social isolation prevalence (soc_low_prev)",
  ylabel = "% patients with delay > 180 d",
  title  = "Social isolation vs extreme delay"
)

# -----------------------------------------------------------------------------
# 7. Combine and save
# -----------------------------------------------------------------------------
combined <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Ecological analysis: ELSI health indicators × APAC treatment delays",
    theme = theme(plot.title = element_text(size = 15, face = "bold"))
  )

dir.create("analysis/figures", showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = "analysis/figures/fig6_ecological.pdf",
  plot     = combined,
  width    = 14,
  height   = 10,
  dpi      = 300
)

cat("\nSaved: analysis/figures/fig6_ecological.pdf\n")
