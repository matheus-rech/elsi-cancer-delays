# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

**ELSI-Brazil Cancer & Aging Research Platform** — a comprehensive analysis pipeline and interactive Shiny app combining:
1. **ELSI-Brazil Wave 2** survey data (nationally representative, adults aged 50+)
2. **DATASUS APAC-Quimioterapia** administrative records (SUS chemotherapy, all 27 Brazilian states, 2022)

The project produces publication-quality figures, tables, and a STROBE-compliant manuscript for a paper on cancer treatment delays among older Brazilians (Lei 12.732 compliance).

## Project structure

```
elsi_native_addons/
  app/
    app.R                     ← 10-tab Shiny app (ELSI + APAC Time to Treatment)
    elsi_data.csv             ← ELSI Wave 2 enhanced dataset
    apac_cohort_national.rds  ← National APAC cohort for Shiny TTT tab
  analysis/
    01_download_apac.R        ← Download APAC-QT from DATASUS (all 27 UFs)
    02_build_apac_cohort.R    ← Build deduplicated analytic cohort (age >=50)
    03_elsi_regional.R        ← Survey-weighted ELSI regional indicators
    04_km_figures.R           ← Kaplan-Meier figures (fig1-fig4)
    05_cox_regression.R       ← Cox proportional hazards (fig5, table2)
    06_table1.R               ← Table 1 demographics by compliance
    07_ecological.R           ← Ecological analysis ELSI x APAC (fig6)
    08_choropleth.R           ← 5-region choropleth (fig7)
    09_competing_risks.R      ← Fine-Gray competing risks (fig8, table3)
    10_sensitivity.R          ← Sensitivity analyses (efig1, etable1-2)
    11_state_choropleth.R     ← 27-state choropleth (fig9)
    manuscript.qmd            ← Quarto manuscript (STROBE structure)
    references.bib            ← Bibliography placeholder
    data/                     ← RDS data files (not in git — too large)
    figures/                  ← Publication PDFs (fig1-fig9, efig1)
    tables/                   ← HTML tables (table1-3, etable1-2)
  docs/superpowers/plans/     ← Implementation plans
  index.html                  ← PWA shell (legacy — needs Shinylive bundle)
  manifest.webmanifest        ← PWA manifest
  CLAUDE.md                   ← This file
```

## Key data files (not in git)

| File | Size | Contents |
|------|------|----------|
| `analysis/data/apac_qt_27uf_2022_raw.rds` | ~202 MB | Raw APAC download, 3.98M rows |
| `analysis/data/apac_cohort.rds` | ~2.4 MB | Analytic cohort, 110,050 patients |
| `analysis/data/elsi_regional.rds` | ~5 KB | 5-region survey-weighted indicators |
| `app/elsi_data.csv` | ~20 MB | ELSI Wave 2 enhanced dataset |
| `app/apac_cohort_national.rds` | ~2.4 MB | Copy of cohort for Shiny app |

## How to reproduce the analysis

```bash
# 1. Download APAC data (requires internet, ~30-60 min)
Rscript analysis/01_download_apac.R

# 2. Build cohort
Rscript analysis/02_build_apac_cohort.R

# 3. ELSI regional summary
Rscript analysis/03_elsi_regional.R

# 4-11. Generate all figures, tables, and analyses
for script in analysis/0{4,5,6,7,8}_*.R analysis/{09,10,11}_*.R; do
  Rscript "$script"
done

# Copy cohort to app directory
cp analysis/data/apac_cohort.rds app/apac_cohort_national.rds

# Render manuscript
cd analysis && quarto render manuscript.qmd --to docx
```

## R dependencies

`microdatasus`, `dplyr`, `tidyr`, `survival`, `survminer`, `tidycmprsk`, `ggsurvfit`,
`srvyr`, `survey`, `ggplot2`, `patchwork`, `gtsummary`, `gt`, `geobr`, `sf`,
`RColorBrewer`, `broom`, `ranger`, `shapviz`, `kernelshap`, `MatchIt`, `cluster`,
`factoextra`, `DT`, `plotly`, `bslib`, `shiny`, `leaflet`, `httr2`, `waiter`, `scales`

## Key findings (national cohort, 110,050 patients)

- **Median time to treatment:** 159 days (2.65x the 60-day legal target)
- **Lei 12.732 compliance:** Only 28.6% receive treatment within 60 days
- **>180 day delay:** ~45% of patients
- **Breast cancer:** 45.5% of cohort, among the slowest to treat
- **Regional disparity:** South/Central-West faster than Southeast (counterintuitively)
- **Stage IV paradox:** Treated slower than Stage I (systemic bottleneck, not clinical triage)

## Target journal

BMJ Global Health / Lancet Regional Health Americas

## Parent project

Lives under `fingertips-br/files/`. The parent is a Fingertips-style public health API for Brazil (Python/FastAPI + Supabase).
