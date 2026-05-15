# 08_choropleth.R
# Choropleth maps — regional treatment compliance + ELSI indicators
# Output: analysis/figures/fig7_choropleth.pdf (14x12 in, 300 dpi)

suppressPackageStartupMessages({
  for (pkg in c("dplyr", "ggplot2", "sf", "geobr", "patchwork",
                "RColorBrewer", "ggrepel")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
})

BASE  <- "analysis"
FIGS  <- file.path(BASE, "figures")
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load ELSI regional data ────────────────────────────────────────────────
cat("Loading ELSI regional data...\n")
elsi_reg <- readRDS(file.path(BASE, "data", "elsi_regional.rds")) |>
  rename(region = region_label) |>
  mutate(
    frail_prev_pct = frail_prev * 100,
    plan_prev_pct  = plan_prev  * 100
  )

cat(sprintf("  ELSI regions: %s\n", paste(elsi_reg$region, collapse = ", ")))

# ── 2. Compute per-region compliance stats from APAC cohort ──────────────────
cat("Loading APAC cohort and computing regional stats...\n")
apac <- readRDS(file.path(BASE, "data", "apac_cohort.rds"))

apac_region <- apac |>
  filter(!is.na(region)) |>
  group_by(region) |>
  summarise(
    pct_compliant = 100 * mean(compliant_60d, na.rm = TRUE),
    median_ttt    = median(ttt, na.rm = TRUE),
    .groups = "drop"
  )

cat("  APAC per-region summary:\n")
print(apac_region)

# ── 3. Join ELSI + APAC ───────────────────────────────────────────────────────
merged <- left_join(elsi_reg, apac_region, by = "region")
cat("\nMerged data:\n")
print(merged |> select(region, pct_compliant, median_ttt, frail_prev_pct, plan_prev_pct))

# ── 4. Load Brazil region shapefile and translate names ──────────────────────
cat("\nDownloading Brazil region shapefile (geobr)...\n")
shp <- geobr::read_region(year = 2020, showProgress = FALSE)

# Translate Portuguese region names to English
name_map <- c(
  "Norte"        = "North",
  "Nordeste"     = "Northeast",
  "Sudeste"      = "Southeast",
  "Sul"          = "South",
  "Centro Oeste" = "Central-West",
  "Centro-Oeste" = "Central-West"
)

shp <- shp |>
  mutate(region = dplyr::recode(name_region, !!!name_map))

cat(sprintf("  Shapefile regions: %s\n", paste(shp$region, collapse = ", ")))

# ── 5. Spatial join ───────────────────────────────────────────────────────────
geo <- left_join(shp, merged, by = "region")

# Compute centroids for label placement
ctr_sf  <- sf::st_centroid(geo)
ctr_xy  <- sf::st_coordinates(ctr_sf)
centroids <- sf::st_drop_geometry(ctr_sf) |>
  mutate(lon = ctr_xy[, 1], lat = ctr_xy[, 2])

# ── 6. Shared theme ──────────────────────────────────────────────────────────
map_theme <- function() {
  theme_void(base_size = 11) +
    theme(
      legend.position  = "bottom",
      legend.key.width = unit(1.8, "cm"),
      legend.key.height = unit(0.35, "cm"),
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8),
      plot.title       = element_text(size = 11, face = "bold", hjust = 0.5,
                                      margin = margin(b = 4)),
      plot.subtitle    = element_text(size = 8, hjust = 0.5, colour = "grey40",
                                      margin = margin(b = 2)),
      plot.margin      = margin(4, 4, 4, 4)
    )
}

make_map <- function(geo_data, centroids_data, fill_var, fill_label,
                     palette, direction = 1, label_fmt = "%.1f%%",
                     subtitle = NULL) {
  ggplot(geo_data) +
    geom_sf(aes(fill = .data[[fill_var]]), colour = "white", linewidth = 0.5) +
    scale_fill_distiller(
      palette   = palette,
      direction = direction,
      name      = fill_label,
      guide     = guide_colorbar(title.position = "top", title.hjust = 0.5)
    ) +
    geom_text(
      data    = centroids_data,
      aes(x = lon, y = lat,
          label = sprintf(paste0(region, "\n", label_fmt),
                          .data[[fill_var]])),
      size    = 2.8,
      fontface = "bold",
      colour  = "grey20",
      lineheight = 1.1
    ) +
    labs(title = fill_label, subtitle = subtitle) +
    map_theme()
}

# ── 7. Build four choropleth panels ──────────────────────────────────────────
cat("\nBuilding choropleth maps...\n")

p1 <- make_map(
  geo_data      = geo,
  centroids_data = centroids,
  fill_var      = "pct_compliant",
  fill_label    = "Lei 12.732 Compliance (%)",
  palette       = "RdYlGn",
  direction     = 1,
  label_fmt     = "%.1f%%",
  subtitle      = "% treated within 60 days"
)

p2 <- make_map(
  geo_data      = geo,
  centroids_data = centroids,
  fill_var      = "median_ttt",
  fill_label    = "Median Days to Treatment",
  palette       = "YlOrRd",
  direction     = 1,
  label_fmt     = "%.0f d",
  subtitle      = "Median time-to-treatment (days)"
)

p3 <- make_map(
  geo_data      = geo,
  centroids_data = centroids,
  fill_var      = "frail_prev_pct",
  fill_label    = "ELSI: Frailty Prevalence",
  palette       = "YlOrRd",
  direction     = 1,
  label_fmt     = "%.1f%%",
  subtitle      = "Survey-weighted prevalence (%)"
)

p4 <- make_map(
  geo_data      = geo,
  centroids_data = centroids,
  fill_var      = "plan_prev_pct",
  fill_label    = "ELSI: Health Plan Coverage (%)",
  palette       = "RdYlGn",
  direction     = 1,
  label_fmt     = "%.1f%%",
  subtitle      = "Survey-weighted coverage (%)"
)

# ── 8. Combine with patchwork ─────────────────────────────────────────────────
cat("Combining panels with patchwork...\n")
combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title   = "Cancer Treatment Delays and Population Health Vulnerability \u2014 Brazil, 2022",
    theme   = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5,
                                margin = margin(b = 8))
    )
  )

# ── 9. Save ───────────────────────────────────────────────────────────────────
out_file <- file.path(FIGS, "fig7_choropleth.pdf")
cat(sprintf("Saving to: %s\n", out_file))

ggsave(
  filename = out_file,
  plot     = combined,
  width    = 14,
  height   = 12,
  units    = "in",
  dpi      = 300,
  device   = "pdf"
)

cat(sprintf("Done. File size: %.1f KB\n", file.size(out_file) / 1024))
