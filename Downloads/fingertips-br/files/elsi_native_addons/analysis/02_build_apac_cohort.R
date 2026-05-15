## 02_build_apac_cohort.R
## Build deduplicated analytic cohort from raw APAC-QT data
## Input:  analysis/data/apac_qt_5uf_2022_raw.rds  (1.65M rows, 75 cols)
## Output: analysis/data/apac_cohort.rds

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
})

IN_FILE  <- "analysis/data/apac_qt_27uf_2022_raw.rds"
OUT_FILE <- "analysis/data/apac_cohort.rds"

# ── 1. Read raw data ──────────────────────────────────────────────────────────
message("Reading raw RDS …")
raw <- readRDS(IN_FILE)
message(sprintf("  Raw: %s rows × %d cols", format(nrow(raw), big.mark = ","), ncol(raw)))

# ── 2. Parse dates ────────────────────────────────────────────────────────────
message("Parsing dates …")
raw <- raw %>%
  mutate(
    dt_dx = as.Date(as.character(AQ_DTIDEN), "%Y%m%d"),
    dt_tx = as.Date(as.character(AQ_DTINTR), "%Y%m%d")
  )

# ── 3. Extract & rename key columns ──────────────────────────────────────────
message("Extracting columns …")
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

# ── 4. Filter: age >= 50 ──────────────────────────────────────────────────────
message("Filtering age >= 50 …")
cohort <- cohort %>% filter(!is.na(age), age >= 50)
message(sprintf("  After age filter: %s rows", format(nrow(cohort), big.mark = ",")))

# ── 5. Deduplicate: first APAC per patient-cancer ─────────────────────────────
message("Deduplicating (first APAC per patient-cancer) …")
cohort <- cohort %>%
  filter(!is.na(cns), !is.na(cid)) %>%
  group_by(cns, cid) %>%
  arrange(dt_tx, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()
message(sprintf("  After dedup: %s rows", format(nrow(cohort), big.mark = ",")))

# ── 6. Compute time-to-treatment (days) ───────────────────────────────────────
message("Computing TTT …")
cohort <- cohort %>%
  mutate(ttt = as.numeric(dt_tx - dt_dx))

# ── 7. Filter: plausible TTT (0–1825 days = 0–5 years) ───────────────────────
message("Filtering plausible TTT (0–1825 days) …")
cohort <- cohort %>% filter(!is.na(ttt), ttt >= 0, ttt <= 1825)
message(sprintf("  After TTT filter: %s rows", format(nrow(cohort), big.mark = ",")))

# ── 8. Derive variables ───────────────────────────────────────────────────────
message("Deriving analysis variables …")

## UF code → state name mapping
uf_to_name <- c(
  "11" = "Rondonia",        "12" = "Acre",            "13" = "Amazonas",
  "14" = "Roraima",         "15" = "Para",            "16" = "Amapa",
  "17" = "Tocantins",       "21" = "Maranhao",        "22" = "Piaui",
  "23" = "Ceara",           "24" = "Rio Grande do Norte", "25" = "Paraiba",
  "26" = "Pernambuco",      "27" = "Alagoas",         "28" = "Sergipe",
  "29" = "Bahia",           "31" = "Minas Gerais",    "32" = "Espirito Santo",
  "33" = "Rio de Janeiro",  "35" = "Sao Paulo",       "41" = "Parana",
  "42" = "Santa Catarina",  "43" = "Rio Grande do Sul","50" = "Mato Grosso do Sul",
  "51" = "Mato Grosso",     "52" = "Goias",           "53" = "Distrito Federal"
)

## UF numeric code (first 2 chars of 6-char municipality code)
uf_to_region <- c(
  "11" = "North",       # RO
  "12" = "North",       # AC
  "13" = "North",       # AM
  "14" = "North",       # RR
  "15" = "North",       # PA
  "16" = "North",       # AP
  "17" = "North",       # TO
  "21" = "Northeast",   # MA
  "22" = "Northeast",   # PI
  "23" = "Northeast",   # CE
  "24" = "Northeast",   # RN
  "25" = "Northeast",   # PB
  "26" = "Northeast",   # PE
  "27" = "Northeast",   # AL
  "28" = "Northeast",   # SE
  "29" = "Northeast",   # BA
  "31" = "Southeast",   # MG
  "32" = "Southeast",   # ES
  "33" = "Southeast",   # RJ
  "35" = "Southeast",   # SP
  "41" = "South",       # PR
  "42" = "South",       # SC
  "43" = "South",       # RS
  "50" = "Central-West",# MS
  "51" = "Central-West",# MT
  "52" = "Central-West",# GO
  "53" = "Central-West" # DF
)

cohort <- cohort %>%
  mutate(
    # Sex label
    sex_label = case_when(
      sex == "F" ~ "Female",
      sex == "M" ~ "Male",
      TRUE       ~ NA_character_
    ),

    # Age group
    age_grp = cut(
      age,
      breaks = c(50, 60, 70, 80, Inf),
      labels = c("50-59", "60-69", "70-79", "80+"),
      right  = FALSE,
      include.lowest = FALSE
    ),

    # Stage label (AQ_ESTADI codes "0"–"4")
    stage_label = case_when(
      stage == "0" ~ "In situ",
      stage == "1" ~ "I",
      stage == "2" ~ "II",
      stage == "3" ~ "III",
      stage == "4" ~ "IV",
      TRUE         ~ NA_character_
    ),

    # Cancer group (based on CID-10 prefix)
    cancer_group = case_when(
      grepl("^C50", cid)            ~ "Breast",
      cid == "C61"                  ~ "Prostate",
      grepl("^C18", cid) | cid %in% c("C19", "C20") ~ "Colorectal",
      grepl("^C34", cid)            ~ "Lung",
      grepl("^C16", cid)            ~ "Stomach",
      grepl("^C53", cid)            ~ "Cervical",
      grepl("^C91", cid) | grepl("^C92", cid) ~ "Leukemia",
      grepl("^C82", cid) | grepl("^C83", cid) | grepl("^C85", cid) ~ "Lymphoma",
      !is.na(cid)                   ~ "Other",
      TRUE                          ~ NA_character_
    ),

    # State name and region (mapped from 2-digit UF IBGE code)
    uf_name = dplyr::recode(uf, !!!uf_to_name, .default = NA_character_),
    region = dplyr::recode(uf, !!!uf_to_region, .default = NA_character_),

    # Treatment compliance flags
    compliant_60d = ttt <= 60,
    compliant_30d = ttt <= 30
  )

# ── 9. Save ───────────────────────────────────────────────────────────────────
message(sprintf("Saving cohort to: %s", OUT_FILE))
saveRDS(cohort, file = OUT_FILE)
message(sprintf("  File size: %.1f MB", file.size(OUT_FILE) / 1024^2))

# ── 10. Summary report ────────────────────────────────────────────────────────
message("\n=== COHORT SUMMARY ===")
message(sprintf("Total rows: %s", format(nrow(cohort), big.mark = ",")))

message("\n--- Rows per UF (download state) ---")
uf_counts <- cohort %>%
  count(uf, region, sort = TRUE) %>%
  as.data.frame()
print(uf_counts)

message("\n--- Rows per cancer_group ---")
cg_counts <- cohort %>%
  count(cancer_group, sort = TRUE) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  as.data.frame()
print(cg_counts)

message("\n--- Time-to-treatment (days) summary ---")
print(summary(cohort$ttt))
message(sprintf("  SD: %.1f days", sd(cohort$ttt, na.rm = TRUE)))

message("\n--- Compliance rates ---")
message(sprintf("  <= 30 days: %s / %s  (%.1f%%)",
  format(sum(cohort$compliant_30d, na.rm = TRUE), big.mark = ","),
  format(nrow(cohort), big.mark = ","),
  100 * mean(cohort$compliant_30d, na.rm = TRUE)))
message(sprintf("  <= 60 days: %s / %s  (%.1f%%)",
  format(sum(cohort$compliant_60d, na.rm = TRUE), big.mark = ","),
  format(nrow(cohort), big.mark = ","),
  100 * mean(cohort$compliant_60d, na.rm = TRUE)))

message("\n--- Sex distribution ---")
print(table(cohort$sex_label, useNA = "always"))

message("\n--- Age group distribution ---")
print(table(cohort$age_grp, useNA = "always"))

message("\n--- Stage distribution ---")
print(table(cohort$stage_label, useNA = "always"))

message("\n--- Region distribution ---")
print(table(cohort$region, useNA = "always"))

message("\nDone. Cohort saved to: ", OUT_FILE)
