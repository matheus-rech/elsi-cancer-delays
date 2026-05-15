# 11_state_choropleth.R
# State-level choropleth maps (27 UFs) — treatment compliance + delays
# Output: analysis/figures/fig9_state_choropleth.pdf (14x14 in, 300 dpi)

suppressPackageStartupMessages({
  for (pkg in c("dplyr", "ggplot2", "sf", "geobr", "patchwork", "RColorBrewer")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
})

BASE <- "analysis"
FIGS <- file.path(BASE, "figures")
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load APAC cohort and compute per-state stats ──────────────────────────
cat("Loading APAC cohort...\n")
apac <- readRDS(file.path(BASE, "data", "apac_cohort.rds"))

cat(sprintf("  Rows: %d | Columns: %s\n", nrow(apac), paste(names(apac), collapse = ", ")))

apac_state <- apac |>
  filter(!is.na(uf), !is.na(ttt)) |>
  mutate(
    code_state   = as.numeric(uf),
    over_180d    = ttt > 180
  ) |>
  group_by(code_state) |>
  summarise(
    n             = n(),
    median_ttt    = median(ttt, na.rm = TRUE),
    pct_compliant = 100 * mean(compliant_60d, na.rm = TRUE),
    pct_over_180d = 100 * mean(over_180d, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(log10_n = log10(n))

cat("  Per-state summary:\n")
print(apac_state)

# ── 2. Load state shapefile ───────────────────────────────────────────────────
cat("\nDownloading Brazil state shapefile (geobr)...\n")
shp <- geobr::read_state(year = 2020, showProgress = FALSE)
cat(sprintf("  Shapefile rows: %d | code_state range: %d–%d\n",
            nrow(shp), min(shp$code_state), max(shp$code_state)))

# ── 3. Left join shapefile with per-state stats ───────────────────────────────
geo <- left_join(shp, apac_state, by = "code_state")

matched <- sum(!is.na(geo$n))
cat(sprintf("  States with data: %d / %d\n", matched, nrow(geo)))

# ── 4. Centroids for state abbreviation labels ────────────────────────────────
ctr_sf  <- sf::st_centroid(geo)
ctr_xy  <- sf::st_coordinates(ctr_sf)
centroids <- sf::st_drop_geometry(ctr_sf) |>
  mutate(lon = ctr_xy[, 1], lat = ctr_xy[, 2])

# ── 5. Shared map theme ───────────────────────────────────────────────────────
map_theme <- function() {
  theme_void(base_size = 10) +
    theme(
      legend.position   = "bottom",
      legend.key.width  = unit(1.6, "cm"),
      legend.key.height = unit(0.32, "cm"),
      legend.title      = element_text(size = 8, face = "bold"),
      legend.text       = element_text(size = 7),
      plot.title        = element_text(size = 10, face = "bold", hjust = 0.5,
                                       margin = margin(b = 3)),
      plot.subtitle     = element_text(size = 7.5, hjust = 0.5, colour = "grey40",
                                       margin = margin(b = 2)),
      plot.margin       = margin(4, 4, 4, 4)
    )
}

# ── 6. Panel builder ──────────────────────────────────────────────────────────
make_state_map <- function(geo_data, centroids_data,
                           fill_var, fill_label, subtitle,
                           palette, direction = 1,
                           log_scale = FALSE) {

  # Choose scale
  if (log_scale) {
    fill_scale <- scale_fill_distiller(
      palette   = palette,
      direction = direction,
      name      = fill_label,
      na.value  = "grey85",
      labels    = function(x) formatC(10^x, format = "fg", big.mark = ","),
      guide     = guide_colorbar(title.position = "top", title.hjust = 0.5)
    )
    label_col <- fill_var  # log10_n
  } else {
    fill_scale <- scale_fill_distiller(
      palette   = palette,
      direction = direction,
      name      = fill_label,
      na.value  = "grey85",
      guide     = guide_colorbar(title.position = "top", title.hjust = 0.5)
    )
    label_col <- fill_var
  }

  # Centroid labels: only draw where data exist
  ctr_with_data <- centroids_data |> filter(!is.na(.data[[label_col]]))

  ggplot(geo_data) +
    geom_sf(aes(fill = .data[[fill_var]]),
            colour = "white", linewidth = 0.3) +
    fill_scale +
    geom_text(
      data     = ctr_with_data,
      aes(x = lon, y = lat, label = abbrev_state),
      size     = 1.9,
      fontface = "bold",
      colour   = "grey15"
    ) +
    labs(title = fill_label, subtitle = subtitle) +
    map_theme()
}

# ── 7. Build four panels ──────────────────────────────────────────────────────
cat("\nBuilding choropleth panels...\n")

p1 <- make_state_map(
  geo_data       = geo,
  centroids_data = centroids,
  fill_var       = "pct_compliant",
  fill_label     = "A  Lei 12.732 Compliance (%)",
  subtitle       = "% treated within 60 days",
  palette        = "RdYlGn",
  direction      = 1,
  log_scale      = FALSE
)

p2 <- make_state_map(
  geo_data       = geo,
  centroids_data = centroids,
  fill_var       = "median_ttt",
  fill_label     = "B  Median Time to Treatment",
  subtitle       = "Median days from diagnosis to treatment",
  palette        = "YlOrRd",
  direction      = 1,
  log_scale      = FALSE
)

p3 <- make_state_map(
  geo_data       = geo,
  centroids_data = centroids,
  fill_var       = "log10_n",
  fill_label     = "C  Sample Size (N patients)",
  subtitle       = "Log\u2081\u2080 scale",
  palette        = "Blues",
  direction      = 1,
  log_scale      = TRUE
)

p4 <- make_state_map(
  geo_data       = geo,
  centroids_data = centroids,
  fill_var       = "pct_over_180d",
  fill_label     = "D  Extreme Delay (% >180 days)",
  subtitle       = "% with time-to-treatment exceeding 180 days",
  palette        = "YlOrRd",
  direction      = 1,
  log_scale      = FALSE
)

# ── 8. Combine with patchwork ─────────────────────────────────────────────────
cat("Combining panels...\n")
combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title   = "Cancer Treatment Timeliness by Brazilian State (27 UFs) \u2014 APAC 2022",
    caption = "Source: DATASUS/SIA-APAC. Lei 12.732/2012 mandates treatment within 60 days of diagnosis.",
    theme   = theme(
      plot.title   = element_text(size = 12, face = "bold", hjust = 0.5,
                                  margin = margin(b = 6)),
      plot.caption = element_text(size = 7, colour = "grey50", hjust = 0.5,
                                  margin = margin(t = 6))
    )
  )

# ── 9. Save ───────────────────────────────────────────────────────────────────
out_file <- file.path(FIGS, "fig9_state_choropleth.pdf")
cat(sprintf("Saving to: %s\n", out_file))

ggsave(
  filename = out_file,
  plot     = combined,
  width    = 14,
  height   = 14,
  units    = "in",
  dpi      = 300,
  device   = "pdf"
)

cat(sprintf("Done. File size: %.1f KB\n", file.size(out_file) / 1024))
