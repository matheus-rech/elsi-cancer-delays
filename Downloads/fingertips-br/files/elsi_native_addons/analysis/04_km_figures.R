# 04_km_figures.R
# Publication-quality Kaplan-Meier figures (cumulative incidence of treatment)
# Output: analysis/figures/fig1_km_overall.pdf through fig4_km_region.pdf

suppressPackageStartupMessages({
  if (!requireNamespace("survminer", quietly = TRUE))
    install.packages("survminer", repos = "https://cloud.r-project.org", quiet = TRUE)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
})

# ── paths ──────────────────────────────────────────────────────────────────────
BASE  <- "analysis"
DATA  <- file.path(BASE, "data", "apac_cohort.rds")
FIGS  <- file.path(BASE, "figures")
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)

# ── load & prep ────────────────────────────────────────────────────────────────
d <- readRDS(DATA)

d <- d |>
  mutate(
    time365  = pmin(ttt, 365),
    event365 = as.integer(ttt <= 365)
  )

cat(sprintf("Cohort: %d patients | median TTT = %.0f days\n",
            nrow(d), median(d$ttt, na.rm = TRUE)))

# ── shared helpers ─────────────────────────────────────────────────────────────
vline_60 <- geom_vline(xintercept = 60, linetype = "dashed",
                       colour = "grey30", linewidth = 0.6)

annotate_60 <- function(p) {
  p$plot <- p$plot +
    vline_60 +
    annotate("text", x = 62, y = 0.02, label = "Lei 12.732\n(60 days)",
             hjust = 0, size = 3, colour = "grey30")
  p
}

save_km <- function(p, fname, w = 10, h = 9) {
  path <- file.path(FIGS, fname)
  pdf(path, width = w, height = h)
  print(p, newpage = FALSE)
  dev.off()
  cat("Saved:", path, "\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Fig 1 — Overall KM
# ═══════════════════════════════════════════════════════════════════════════════
n_total  <- nrow(d)
med_ttt  <- median(d$ttt, na.rm = TRUE)

fit1 <- survfit(Surv(time365, event365) ~ 1, data = d)

p1 <- ggsurvplot(
  fit1,
  data          = d,
  fun           = "event",
  conf.int      = TRUE,
  risk.table    = TRUE,
  risk.table.height = 0.22,
  xlab          = "Days from diagnosis to first treatment",
  ylab          = "Cumulative incidence of treatment",
  title         = sprintf(
    "Time to Treatment - Overall (n = %s; median = %d days)",
    format(n_total, big.mark = ","), round(med_ttt)
  ),
  xlim          = c(0, 365),
  break.time.by = 60,
  palette       = "#2E86AB",
  ggtheme       = theme_classic(base_size = 12),
  surv.median.line = "hv",
  legend        = "none"
)

p1 <- annotate_60(p1)
save_km(p1, "fig1_km_overall.pdf", w = 10, h = 9)

# ═══════════════════════════════════════════════════════════════════════════════
# Fig 2 — By cancer type (top 6)
# ═══════════════════════════════════════════════════════════════════════════════
top6 <- c("Breast", "Prostate", "Colorectal", "Lung", "Stomach", "Cervical")

d2 <- d |>
  filter(cancer_group %in% top6) |>
  mutate(cancer_group = factor(cancer_group, levels = top6))

fit2 <- survfit(Surv(time365, event365) ~ cancer_group, data = d2)

p2 <- ggsurvplot(
  fit2,
  data              = d2,
  fun               = "event",
  conf.int          = FALSE,
  risk.table        = TRUE,
  risk.table.height = 0.28,
  pval              = TRUE,
  pval.size         = 4,
  xlab              = "Days from diagnosis to first treatment",
  ylab              = "Cumulative incidence of treatment",
  title             = "Time to Treatment by Cancer Type",
  xlim              = c(0, 365),
  break.time.by     = 60,
  legend.title      = "Cancer type",
  legend.labs       = top6,
  palette           = c("#E74C3C", "#3498DB", "#27AE60",
                        "#F39C12", "#8E44AD", "#1ABC9C"),
  ggtheme           = theme_classic(base_size = 12)
)

p2 <- annotate_60(p2)
save_km(p2, "fig2_km_cancer.pdf", w = 12, h = 10)

# ═══════════════════════════════════════════════════════════════════════════════
# Fig 3 — By stage (I–IV only)
# ═══════════════════════════════════════════════════════════════════════════════
stage_levels <- c("I", "II", "III", "IV")

d3 <- d |>
  filter(stage_label %in% stage_levels) |>
  mutate(stage_label = factor(stage_label, levels = stage_levels))

fit3 <- survfit(Surv(time365, event365) ~ stage_label, data = d3)

p3 <- ggsurvplot(
  fit3,
  data              = d3,
  fun               = "event",
  conf.int          = FALSE,
  risk.table        = TRUE,
  risk.table.height = 0.26,
  pval              = TRUE,
  pval.size         = 4,
  xlab              = "Days from diagnosis to first treatment",
  ylab              = "Cumulative incidence of treatment",
  title             = "Time to Treatment by Clinical Stage",
  xlim              = c(0, 365),
  break.time.by     = 60,
  legend.title      = "Stage",
  legend.labs       = stage_levels,
  palette           = c("#F1C40F", "#E67E22", "#E74C3C", "#922B21"),
  ggtheme           = theme_classic(base_size = 12)
)

p3 <- annotate_60(p3)
save_km(p3, "fig3_km_stage.pdf", w = 12, h = 10)

# ═══════════════════════════════════════════════════════════════════════════════
# Fig 4 — By region
# ═══════════════════════════════════════════════════════════════════════════════
region_levels  <- c("North", "Northeast", "Southeast", "South", "Central-West")
region_palette <- c(
  "North"        = "#E74C3C",
  "Northeast"    = "#F39C12",
  "Southeast"    = "#3498DB",
  "South"        = "#27AE60",
  "Central-West" = "#8E44AD"
)

d4 <- d |>
  filter(region %in% region_levels) |>
  mutate(region = factor(region, levels = region_levels))

fit4 <- survfit(Surv(time365, event365) ~ region, data = d4)

p4 <- ggsurvplot(
  fit4,
  data              = d4,
  fun               = "event",
  conf.int          = FALSE,
  risk.table        = TRUE,
  risk.table.height = 0.28,
  pval              = TRUE,
  pval.size         = 4,
  xlab              = "Days from diagnosis to first treatment",
  ylab              = "Cumulative incidence of treatment",
  title             = "Time to Treatment by Region",
  xlim              = c(0, 365),
  break.time.by     = 60,
  legend.title      = "Region",
  legend.labs       = region_levels,
  palette           = unname(region_palette[region_levels]),
  ggtheme           = theme_classic(base_size = 12)
)

p4 <- annotate_60(p4)
save_km(p4, "fig4_km_region.pdf", w = 12, h = 10)

cat("\nAll 4 KM figures saved to:", FIGS, "\n")
