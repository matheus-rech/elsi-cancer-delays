## 01_download_apac.R
## Download APAC-QT (SIA-AQ) data for ALL 27 UFs — full year 2022
## Caches per-UF RDS files to resume after interruption

library(microdatasus)
library(dplyr)

# ── Configuration ────────────────────────────────────────────────────────────
UFS <- c("AC","AL","AM","AP","BA","CE","DF","ES","GO","MA",
         "MG","MS","MT","PA","PB","PE","PI","PR","RJ","RN",
         "RO","RR","RS","SC","SE","SP","TO")
YEAR    <- 2022
OUT_DIR <- "analysis/data"
OUT     <- file.path(OUT_DIR, "apac_qt_27uf_2022_raw.rds")
LOG     <- file.path(OUT_DIR, "01_download_apac_27uf.log")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Start log
sink(LOG, split = TRUE)
cat(sprintf("=== APAC-QT 27-UF Download — %s ===\n\n", Sys.time()))

# ── Download loop (with per-UF cache) ────────────────────────────────────────
results <- list()
failed  <- character(0)

for (uf in UFS) {
  cache_file <- file.path(OUT_DIR, paste0("apac_qt_", uf, "_2022.rds"))

  # Resume from cache

  if (file.exists(cache_file)) {
    cat(sprintf("[%s] Cached — loading %s\n", uf,
                format(file.size(cache_file), big.mark = ",")))
    results[[uf]] <- readRDS(cache_file)
    next
  }

  cat(sprintf("[%s] Downloading SIA-AQ %d...\n", uf, YEAR))
  tryCatch({
    df <- fetch_datasus(
      year_start         = YEAR,
      month_start        = 1,
      year_end           = YEAR,
      month_end          = 12,
      uf                 = uf,
      information_system = "SIA-AQ",
      timeout            = 600
    )

    if (nrow(df) == 0) {
      cat(sprintf("[%s] Empty result — skipping\n", uf))
      next
    }

    df$uf_download <- uf
    saveRDS(df, cache_file)
    results[[uf]] <- df
    cat(sprintf("[%s] Done: %s rows\n", uf, format(nrow(df), big.mark = ",")))

  }, error = function(e) {
    cat(sprintf("[%s] FAILED: %s\n", uf, e$message))
    failed <<- c(failed, uf)
  })
}

# ── Report failures ──────────────────────────────────────────────────────────
if (length(failed) > 0) {
  cat(sprintf("\n=== FAILED UFs (retry manually): %s\n", paste(failed, collapse = ", ")))
}

# ── Combine all ──────────────────────────────────────────────────────────────
cat("\nBinding all UFs...\n")
combined <- bind_rows(results)

cat(sprintf("\n=== TOTAL: %s rows from %d UFs ===\n",
            format(nrow(combined), big.mark = ","), length(results)))

# Per-UF breakdown
per_uf <- combined %>% count(uf_download, name = "n_rows") %>% arrange(desc(n_rows))
print(as.data.frame(per_uf))

# ── Save combined ────────────────────────────────────────────────────────────
saveRDS(combined, file = OUT)
cat(sprintf("\nSaved: %s\n", OUT))
cat(sprintf("File size: %.1f MB\n", file.size(OUT) / 1024^2))

# ── Verify ───────────────────────────────────────────────────────────────────
check <- readRDS(OUT)
stopifnot(nrow(check) == nrow(combined))
cat(sprintf("Readback OK — %s rows confirmed.\n", format(nrow(check), big.mark = ",")))

cat(sprintf("\n=== Completed at %s ===\n", Sys.time()))
sink()
