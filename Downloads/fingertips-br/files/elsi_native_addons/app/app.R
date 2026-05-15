library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(tidyr)
library(srvyr)
library(survey)
library(plotly)
library(DT)
library(patchwork)
library(scales)
library(ranger)
library(shapviz)
library(kernelshap)
library(MatchIt)
library(cluster)
library(factoextra)
library(broom)
library(httr2)
library(sf)
library(geobr)
library(leaflet)
library(RColorBrewer)
library(waiter)
library(survival)

# ============================================================
# DATA
# ============================================================
df <- read.csv("elsi_data.csv", stringsAsFactors = FALSE)

# APAC data for Time to Treatment analysis (national cohort, all 27 UFs)
apac_file <- "apac_cohort_national.rds"
if (file.exists(apac_file)) {
  apac <- readRDS(apac_file)
  has_apac <- TRUE
} else {
  has_apac <- FALSE
  apac <- NULL
}

# Factor levels
df$age_group <- factor(df$age_group, levels = c("50-59","60-69","70-79","80+"))
df$self_health_label <- factor(df$self_health_label,
  levels = c("Very good","Good","Fair","Poor","Very poor"))
df$region_label <- factor(df$region_label,
  levels = c("North","Northeast","Southeast","South","Central-West"))
df$income_quartile <- factor(df$income_quartile,
  levels = c("Q1 (lowest)","Q2","Q3","Q4 (highest)"))
df$education_cat <- factor(df$education_cat,
  levels = c("None","1-4 years","5-8 years","9-11 years","12+ years"))
df$frailty_status <- factor(df$frailty_status,
  levels = c("Robust","Pre-frail","Frail"))
df$bmi_cat <- factor(df$bmi_cat,
  levels = c("Underweight","Normal","Overweight","Obese"))
df$multimorbidity <- factor(df$multimorbidity, levels = c("0","1","2","3+"))
df$race_label <- factor(df$race_label,
  levels = c("White","Black","Mixed/Pardo","Asian","Indigenous"))
df$smoking_label <- factor(df$smoking_label, levels = c("Never","Former","Current"))
df$late_diagnosis <- factor(df$late_diagnosis,
  levels = c("Early (<50)","Middle (50-69)","Late (≥70)"))

# Social activity labels
soc_vars <- c(
  soc_letters = "Letters/phone", soc_visit = "Visit friends",
  soc_invite = "Invite home", soc_go_out = "Public places",
  soc_organized = "Organized activities", soc_civic = "Civic associations",
  soc_internet = "Internet use", soc_drive = "Drove a car",
  soc_games = "Games with peers", soc_hobby = "Hobbies/crafts",
  soc_short_trip = "Short trips", soc_long_trip = "Long trips",
  soc_volunteer = "Volunteer work"
)

# ESMO theme
esmo_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "#1F3864"),
    plot.subtitle = element_text(color = "#555555", size = 12),
    axis.title = element_text(color = "#333333"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )
esmo_colors <- c("No" = "#4A90D9", "Yes" = "#E74C3C")
esmo_fill <- scale_fill_manual(values = esmo_colors, name = "Cancer")

# Pre-load Brazil region shapefile
tryCatch({
  br_regions <- read_region(year = 2020, showProgress = FALSE)
  br_regions$region_label <- case_when(
    br_regions$name_region == "Norte" ~ "North",
    br_regions$name_region == "Nordeste" ~ "Northeast",
    br_regions$name_region == "Sudeste" ~ "Southeast",
    br_regions$name_region == "Sul" ~ "South",
    br_regions$name_region == "Centro-Oeste" ~ "Central-West"
  )
  has_map <- TRUE
}, error = function(e) {
  has_map <<- FALSE
  br_regions <<- NULL
})

# ============================================================
# UI
# ============================================================
ui <- page_navbar(
  title = "ELSI-Brazil: Cancer & Aging Research Platform",
  id = "main_nav",
  theme = bs_theme(
    version = 5, bootswatch = "flatly",
    primary = "#1F3864", "navbar-bg" = "#1F3864"
  ),
  sidebar = sidebar(
    width = 280, title = "Global Filters",
    selectInput("filter_sex", "Sex",
      choices = c("All", "Female", "Male"), selected = "All"),
    selectInput("filter_age", "Age group",
      choices = c("All", "50-59", "60-69", "70-79", "80+"), selected = "All"),
    selectInput("filter_region", "Region",
      choices = c("All", "North", "Northeast", "Southeast", "South", "Central-West"),
      selected = "All"),
    selectInput("filter_race", "Race/ethnicity",
      choices = c("All", "White", "Black", "Mixed/Pardo"), selected = "All"),
    selectInput("filter_education", "Education",
      choices = c("All", "None", "1-4 years", "5-8 years", "9-11 years", "12+ years"),
      selected = "All"),
    hr(),
    tags$small(tags$em("ELSI-Brazil Wave 2 (n=9,949). Survey-weighted estimates."))
  ),

  # ---- Tab 1: Overview ----
  nav_panel("Overview", icon = icon("chart-bar"),
    layout_columns(col_widths = c(2, 2, 2, 2, 2, 2),
      value_box("N (weighted)", textOutput("vb_n"), theme = "primary"),
      value_box("Cancer prev.", textOutput("vb_cancer"), theme = "danger"),
      value_box("Mean age", textOutput("vb_age"), theme = "info"),
      value_box("Frail %", textOutput("vb_frail"), theme = "warning"),
      value_box("ADL disability", textOutput("vb_adl"), theme = "secondary"),
      value_box("Hospitalized", textOutput("vb_hosp"), theme = "dark")
    ),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Cancer Prevalence by Age & Sex"),
        plotOutput("plot_cancer_age_sex", height = "350px")),
      card(card_header("Cancer Prevalence by Region"),
        plotOutput("plot_cancer_region", height = "350px"))
    ),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Self-Rated Health: Cancer vs Non-Cancer"),
        plotOutput("plot_self_health", height = "320px")),
      card(card_header("Multimorbidity Distribution"),
        plotOutput("plot_multimorbidity", height = "320px"))
    )
  ),

  # ---- Tab 2: Inequalities ----
  nav_panel("Inequalities", icon = icon("balance-scale"),
    layout_columns(col_widths = c(4, 8),
      card(card_header("Settings"),
        selectInput("ineq_outcome", "Outcome",
          choices = c("Cancer" = "cancer_dx", "Frailty" = "frail_binary",
                      "ADL disability" = "adl_any", "Depression" = "dep_binary",
                      "Hospitalization" = "hospitalized",
                      "Low social participation" = "soc_low",
                      "Multimorbidity (3+)" = "multi3")),
        selectInput("ineq_stratifier", "Stratify by",
          choices = c("Income quartile" = "income_quartile",
                      "Education" = "education_cat",
                      "Race/ethnicity" = "race_label",
                      "Region" = "region_label",
                      "Sex" = "sex_label",
                      "Urban/Rural" = "zona_label")),
        checkboxInput("ineq_adjust", "Adjust for age & sex", value = TRUE)
      ),
      card(card_header("Prevalence by Stratum (survey-weighted)"),
        plotOutput("plot_inequality", height = "450px"))
    ),
    card(card_header("Adjusted Odds Ratios — Forest Plot"),
      plotOutput("plot_forest", height = "400px"))
  ),

  # ---- Tab 3: Social Network ----
  nav_panel("Social Network", icon = icon("users"),
    layout_columns(col_widths = c(3, 3, 3, 3),
      value_box("Score (cancer)", textOutput("vb_soc_cancer"), theme = "danger"),
      value_box("Score (no cancer)", textOutput("vb_soc_no"), theme = "primary"),
      value_box("Difference", textOutput("vb_soc_diff"), theme = "warning"),
      value_box("Social isolation %", textOutput("vb_soc_isolated"), theme = "dark")
    ),
    card(card_header("Social Participation: Cancer vs Non-Cancer (weighted)"),
      plotOutput("plot_social_compare", height = "500px")),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Social Score Distribution"),
        plotOutput("plot_social_dist", height = "320px")),
      card(card_header("Social Score by Age & Cancer Status"),
        plotOutput("plot_social_age", height = "320px"))
    )
  ),

  # ---- Tab 4: Health Profile ----
  nav_panel("Health Profile", icon = icon("heartbeat"),
    layout_columns(col_widths = c(4, 4, 4),
      value_box("Mean Frailty Index", textOutput("vb_fi"), theme = "warning"),
      value_box("Mean QoL Proxy", textOutput("vb_qol"), theme = "success"),
      value_box("Mean Grip (kg)", textOutput("vb_grip"), theme = "info")
    ),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Fried Frailty Phenotype by Cancer Status"),
        plotOutput("plot_frailty_status", height = "380px")),
      card(card_header("Frailty Index Distribution"),
        plotOutput("plot_fi_dist", height = "380px"))
    ),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Comorbidity Prevalence: Cancer vs Non-Cancer"),
        plotOutput("plot_comorbidity", height = "400px")),
      card(card_header("Depression (CES-D ≥ 12) & QoL by Cancer Status"),
        plotOutput("plot_depression_qol", height = "400px"))
    ),
    layout_columns(col_widths = c(6, 6),
      card(card_header("ADL/IADL Disability Count"),
        plotOutput("plot_adl_iadl", height = "350px")),
      card(card_header("BMI Distribution by Cancer Status"),
        plotOutput("plot_bmi", height = "350px"))
    )
  ),

  # ---- Tab 5: Cancer Care ----
  nav_panel("Cancer Care", icon = icon("hospital"),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Treatment Utilization (cancer survivors)"),
        plotOutput("plot_treatment", height = "350px")),
      card(card_header("Cancer Status by Treatment Type"),
        plotOutput("plot_cancer_status", height = "350px"))
    ),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Age at Diagnosis by Income Quartile"),
        plotOutput("plot_dx_age_income", height = "350px")),
      card(card_header("Late Diagnosis (≥70) by Region & Education"),
        plotOutput("plot_late_dx", height = "350px"))
    ),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Healthcare Access: Cancer vs Non-Cancer"),
        plotOutput("plot_healthcare_access", height = "380px")),
      card(card_header("Hospitalization Patterns"),
        plotOutput("plot_hospitalization", height = "380px"))
    )
  ),

  # ---- Tab 6: Statistical Models ----
  nav_panel("Statistical Models", icon = icon("calculator"),
    navset_card_tab(
      nav_panel("Logistic Regression",
        layout_columns(col_widths = c(4, 8),
          card(card_header("Model Setup"),
            selectInput("lr_outcome", "Outcome",
              choices = c("Cancer" = "cancer_dx",
                          "Frailty (Fried)" = "frail_binary",
                          "ADL disability" = "adl_any",
                          "Depression (CES-D≥12)" = "dep_binary",
                          "Hospitalization" = "hospitalized")),
            checkboxGroupInput("lr_predictors", "Predictors",
              choices = c("Age" = "age", "Sex" = "sex",
                          "Education (years)" = "education_yrs",
                          "Income (per capita)" = "income_pc",
                          "Cancer" = "cancer_dx",
                          "N chronic conditions" = "n_chronic",
                          "Social score" = "social_score",
                          "Frailty index" = "frailty_index",
                          "Region" = "region_label",
                          "Race" = "race_label",
                          "Health plan" = "has_health_plan",
                          "BMI" = "bmi"),
              selected = c("age", "sex", "cancer_dx", "education_yrs")),
            actionButton("lr_run", "Run Model", class = "btn-primary w-100")
          ),
          card(card_header("Results"),
            DTOutput("lr_table"),
            plotOutput("lr_forest", height = "400px"))
        )
      ),
      nav_panel("Propensity Score Matching",
        layout_columns(col_widths = c(4, 8),
          card(card_header("PSM Setup"),
            tags$p("Match cancer survivors to controls on demographics."),
            selectInput("psm_outcome", "Compare outcome",
              choices = c("Social score" = "social_score",
                          "Depression score" = "depression_score",
                          "QoL proxy" = "qol_proxy",
                          "Frailty index" = "frailty_index",
                          "ADL disability count" = "adl_disability",
                          "N hospitalizations" = "n_hospitalizations")),
            actionButton("psm_run", "Run PSM", class = "btn-primary w-100")
          ),
          card(card_header("Matched Comparison"),
            verbatimTextOutput("psm_summary"),
            plotOutput("psm_plot", height = "400px"))
        )
      )
    )
  ),

  # ---- Tab 7: ML & Clusters ----
  nav_panel("ML & Clusters", icon = icon("brain"),
    navset_card_tab(
      nav_panel("Clustering",
        layout_columns(col_widths = c(4, 8),
          card(card_header("Settings"),
            selectInput("clust_vars", "Variables", multiple = TRUE,
              choices = c("Social score" = "social_score",
                          "Depression score" = "depression_score",
                          "Frailty index" = "frailty_index",
                          "N chronic conditions" = "n_chronic",
                          "ADL disability" = "adl_disability",
                          "Age" = "age",
                          "Income per capita" = "income_pc",
                          "BMI" = "bmi",
                          "Grip strength" = "grip_strength",
                          "QoL proxy" = "qol_proxy"),
              selected = c("social_score","frailty_index","depression_score","n_chronic")),
            sliderInput("clust_k", "Number of clusters", 2, 6, 3),
            actionButton("clust_run", "Run Clustering", class = "btn-primary w-100")
          ),
          card(card_header("Cluster Visualization"),
            plotOutput("clust_pca", height = "400px"))
        ),
        card(card_header("Cluster Profiles"),
          plotOutput("clust_profiles", height = "350px"))
      ),
      nav_panel("SHAP Analysis",
        layout_columns(col_widths = c(4, 8),
          card(card_header("Random Forest + SHAP"),
            selectInput("shap_outcome", "Predict",
              choices = c("Cancer" = "cancer_dx",
                          "Frailty (binary)" = "frail_binary",
                          "Depression (CES-D≥12)" = "dep_binary",
                          "Hospitalization" = "hospitalized")),
            tags$p(tags$small("Uses ranger RF on 500-sample subset for SHAP.")),
            actionButton("shap_run", "Run SHAP", class = "btn-primary w-100")
          ),
          card(card_header("SHAP Importance"),
            plotOutput("shap_importance", height = "450px"))
        ),
        card(card_header("SHAP Dependence (top variable)"),
          plotOutput("shap_dependence", height = "350px"))
      )
    )
  ),

  # ---- Tab 8: Regional Map ----
  nav_panel("Regional Map", icon = icon("map"),
    layout_columns(col_widths = c(4, 8),
      card(card_header("Indicator"),
        selectInput("map_indicator", "Map indicator",
          choices = c("Cancer prevalence (%)" = "cancer_prev",
                      "Frailty prevalence (%)" = "frail_prev",
                      "Mean social score" = "mean_social",
                      "ADL disability (%)" = "adl_prev",
                      "Hospitalization (%)" = "hosp_prev",
                      "Health plan coverage (%)" = "plan_prev",
                      "Mean frailty index" = "mean_fi",
                      "Depression (%)" = "dep_prev")),
        tags$p(tags$small("Survey-weighted estimates by macro-region."))
      ),
      card(card_header("Brazil — Regional Choropleth"),
        plotOutput("map_choropleth", height = "500px"))
    ),
    card(card_header("Regional Comparison Table"),
      DTOutput("map_table"))
  ),

  # ---- Tab 9: Time to Treatment ----
  nav_panel("Time to Treatment", icon = icon("clock"),
    if (!exists("has_apac") || !has_apac) {
      card(
        card_header("APAC Data Not Found"),
        tags$p("The APAC-Quimioterapia file (apac_qt_sp2022.rds) was not found in the app directory."),
        tags$p("Place the file in the same folder as app.R and restart the app.")
      )
    } else {
      tagList(
        layout_columns(col_widths = c(3, 3, 3, 3),
          value_box("N cases", textOutput("ttt_n"), theme = "primary"),
          value_box("Median TTT (days)", textOutput("ttt_median"), theme = "info"),
          value_box("Compliance ≤60d (%)", textOutput("ttt_comp60"), theme = "success"),
          value_box("Delay >180d (%)", textOutput("ttt_delay180"), theme = "danger")
        ),
        layout_columns(col_widths = c(4, 8),
          card(card_header("Stratify KM curve by"),
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
            ),
            tags$small(tags$em(
              "Data: APAC-Quimioterapia (DATASUS/SIA), all 27 Brazilian states, 2022. Adults \u226550 years. Only SUS chemotherapy."
            ))
          ),
          card(card_header("Kaplan-Meier: Time to First Chemotherapy"),
            plotOutput("ttt_km", height = "380px"))
        ),
        card(card_header("Compliance by Time Band"),
          plotOutput("ttt_bands", height = "280px"))
      )
    }
  ),

  # ---- Tab 10: AI Analyst ----
  nav_panel("AI Analyst", icon = icon("robot"),
    layout_columns(col_widths = c(12),
      card(card_header("Claude-powered Research Assistant"),
        tags$p("Ask questions about the data, request hypothesis generation,
               interpretation of results, or statistical guidance."),
        textAreaInput("ai_question", NULL,
          placeholder = "e.g., What are the main drivers of social isolation among cancer survivors in Brazil?",
          width = "100%", rows = 3),
        layout_columns(col_widths = c(3, 9),
          actionButton("ai_submit", "Ask Claude", class = "btn-primary"),
          tags$small(textOutput("ai_status", inline = TRUE))
        ),
        hr(),
        uiOutput("ai_response")
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # ---- Filtered data ----
  fdata <- reactive({
    d <- df
    if (input$filter_sex != "All") d <- d[d$sex_label == input$filter_sex, ]
    if (input$filter_age != "All") d <- d[d$age_group == input$filter_age, ]
    if (input$filter_region != "All") d <- d[d$region_label == input$filter_region, ]
    if (input$filter_race != "All") d <- d[d$race_label == input$filter_race, ]
    if (input$filter_education != "All") d <- d[d$education_cat == input$filter_education, ]
    d
  })

  # ---- Survey design ----
  svy <- reactive({
    d <- fdata()
    d$frail_binary <- as.integer(!is.na(d$frailty_status) & d$frailty_status == "Frail")
    d$adl_any <- as.integer(d$adl_disability >= 1)
    d$dep_binary <- as.integer(!is.na(d$depression_score) & d$depression_score >= 12)
    d$soc_low <- as.integer(d$social_score <= 1)
    d$multi3 <- as.integer(d$n_chronic >= 3)
    as_survey_design(d, ids = upa, strata = estrato, weights = weight, nest = TRUE)
  })

  # ============================================================
  # TAB 1: OVERVIEW
  # ============================================================
  output$vb_n <- renderText({
    d <- fdata()
    format(round(sum(d$weight, na.rm = TRUE)), big.mark = ",")
  })
  output$vb_cancer <- renderText({
    s <- svy() %>% summarise(p = survey_mean(cancer_dx == 1, na.rm = TRUE))
    paste0(round(s$p * 100, 1), "%")
  })
  output$vb_age <- renderText({
    s <- svy() %>% summarise(m = survey_mean(age, na.rm = TRUE))
    round(s$m, 1)
  })
  output$vb_frail <- renderText({
    s <- svy() %>% summarise(p = survey_mean(frail_binary, na.rm = TRUE))
    paste0(round(s$p * 100, 1), "%")
  })
  output$vb_adl <- renderText({
    s <- svy() %>% summarise(p = survey_mean(adl_any, na.rm = TRUE))
    paste0(round(s$p * 100, 1), "%")
  })
  output$vb_hosp <- renderText({
    s <- svy() %>% summarise(p = survey_mean(hospitalized == 1, na.rm = TRUE))
    paste0(round(s$p * 100, 1), "%")
  })

  output$plot_cancer_age_sex <- renderPlot({
    d <- fdata() %>%
      filter(!is.na(age_group), !is.na(sex_label), !is.na(cancer_label)) %>%
      group_by(age_group, sex_label) %>%
      summarise(prev = weighted.mean(cancer_dx == 1, w = weight, na.rm = TRUE) * 100,
                .groups = "drop")
    ggplot(d, aes(age_group, prev, fill = sex_label)) +
      geom_col(position = "dodge", width = 0.7) +
      scale_fill_manual(values = c("Female" = "#E74C3C", "Male" = "#4A90D9"), name = "Sex") +
      labs(x = "Age group", y = "Cancer prevalence (%)", title = NULL) +
      esmo_theme
  })

  output$plot_cancer_region <- renderPlot({
    d <- fdata() %>%
      filter(!is.na(region_label)) %>%
      group_by(region_label) %>%
      summarise(prev = weighted.mean(cancer_dx == 1, w = weight, na.rm = TRUE) * 100,
                .groups = "drop")
    ggplot(d, aes(reorder(region_label, prev), prev)) +
      geom_col(fill = "#1F3864", width = 0.6) +
      geom_text(aes(label = paste0(round(prev, 1), "%")), hjust = -0.1, size = 4.5) +
      coord_flip() +
      labs(x = NULL, y = "Cancer prevalence (%)") +
      esmo_theme
  })

  output$plot_self_health <- renderPlot({
    d <- fdata() %>%
      filter(!is.na(self_health_label), !is.na(cancer_label)) %>%
      group_by(cancer_label, self_health_label) %>%
      summarise(n = sum(weight, na.rm = TRUE), .groups = "drop_last") %>%
      mutate(pct = n / sum(n) * 100)
    ggplot(d, aes(self_health_label, pct, fill = cancer_label)) +
      geom_col(position = "dodge", width = 0.7) +
      esmo_fill + labs(x = "Self-rated health", y = "%") + esmo_theme
  })

  output$plot_multimorbidity <- renderPlot({
    d <- fdata() %>%
      filter(!is.na(cancer_label), !is.na(multimorbidity))
    ggplot(d, aes(multimorbidity, fill = cancer_label)) +
      geom_bar(aes(weight = weight), position = "dodge", width = 0.7) +
      esmo_fill + labs(x = "Number of chronic conditions", y = "Weighted count") +
      esmo_theme
  })

  # ============================================================
  # TAB 2: INEQUALITIES
  # ============================================================
  output$plot_inequality <- renderPlot({
    s <- svy()
    outcome <- input$ineq_outcome
    strat <- input$ineq_stratifier

    s_data <- s$variables
    s_data$outcome_var <- switch(outcome,
      cancer_dx = s_data$cancer_dx,
      frail_binary = s_data$frail_binary,
      adl_any = s_data$adl_any,
      dep_binary = s_data$dep_binary,
      hospitalized = as.integer(s_data$hospitalized == 1),
      soc_low = s_data$soc_low,
      multi3 = s_data$multi3
    )

    s_data$strat_var <- s_data[[strat]]
    s_upd <- s
    s_upd$variables <- s_data

    res <- s_upd %>%
      filter(!is.na(strat_var), !is.na(outcome_var)) %>%
      group_by(strat_var) %>%
      summarise(
        prev = survey_mean(outcome_var, na.rm = TRUE, vartype = "ci")
      )

    ggplot(res, aes(reorder(strat_var, prev), prev * 100)) +
      geom_col(fill = "#1F3864", width = 0.6) +
      geom_errorbar(aes(ymin = prev_low * 100, ymax = prev_upp * 100),
                    width = 0.2, color = "#E74C3C") +
      coord_flip() +
      labs(x = NULL, y = "Prevalence (%, 95% CI)", title = NULL) +
      esmo_theme
  })

  output$plot_forest <- renderPlot({
    s <- svy()
    outcome <- input$ineq_outcome
    strat <- input$ineq_stratifier

    s_data <- s$variables
    s_data$outcome_var <- switch(outcome,
      cancer_dx = s_data$cancer_dx,
      frail_binary = s_data$frail_binary,
      adl_any = s_data$adl_any,
      dep_binary = s_data$dep_binary,
      hospitalized = as.integer(s_data$hospitalized == 1),
      soc_low = s_data$soc_low,
      multi3 = s_data$multi3
    )
    s_data$strat_var <- factor(s_data[[strat]])
    s_upd <- s
    s_upd$variables <- s_data

    tryCatch({
      fmla <- if (input$ineq_adjust && strat != "sex_label") {
        outcome_var ~ strat_var + age + sex
      } else if (input$ineq_adjust) {
        outcome_var ~ strat_var + age
      } else {
        outcome_var ~ strat_var
      }

      fit <- svyglm(fmla, design = s_upd, family = quasibinomial())
      tidy_fit <- tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(grepl("strat_var", term)) %>%
        mutate(term = gsub("strat_var", "", term))

      ggplot(tidy_fit, aes(estimate, reorder(term, estimate))) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
        geom_point(size = 3, color = "#E74C3C") +
        geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
        labs(x = "Odds Ratio (95% CI)", y = NULL, subtitle = "Reference: first level") +
        esmo_theme
    }, error = function(e) {
      ggplot() + annotate("text", x = 1, y = 1, label = paste("Model error:", e$message), size = 4) +
        theme_void()
    })
  })

  # ============================================================
  # TAB 3: SOCIAL NETWORK
  # ============================================================
  output$vb_soc_cancer <- renderText({
    s <- svy() %>% filter(cancer_dx == 1) %>%
      summarise(m = survey_mean(social_score, na.rm = TRUE))
    round(s$m, 1)
  })
  output$vb_soc_no <- renderText({
    s <- svy() %>% filter(cancer_dx == 0) %>%
      summarise(m = survey_mean(social_score, na.rm = TRUE))
    round(s$m, 1)
  })
  output$vb_soc_diff <- renderText({
    d <- fdata()
    s1 <- weighted.mean(d$social_score[d$cancer_dx == 1], d$weight[d$cancer_dx == 1], na.rm = TRUE)
    s0 <- weighted.mean(d$social_score[d$cancer_dx == 0], d$weight[d$cancer_dx == 0], na.rm = TRUE)
    diff <- s1 - s0
    paste0(ifelse(diff > 0, "+", ""), round(diff, 2))
  })
  output$vb_soc_isolated <- renderText({
    s <- svy() %>% summarise(p = survey_mean(social_score <= 1, na.rm = TRUE))
    paste0(round(s$p * 100, 1), "%")
  })

  output$plot_social_compare <- renderPlot({
    d <- fdata() %>% filter(!is.na(cancer_label))
    soc_long <- d %>%
      select(cancer_label, weight, all_of(names(soc_vars))) %>%
      pivot_longer(-c(cancer_label, weight), names_to = "activity", values_to = "value") %>%
      filter(!is.na(value)) %>%
      group_by(cancer_label, activity) %>%
      summarise(pct = weighted.mean(value == 1, w = weight, na.rm = TRUE) * 100, .groups = "drop") %>%
      mutate(activity_label = soc_vars[activity]) %>%
      mutate(activity_label = factor(activity_label, levels = rev(soc_vars)))
    ggplot(soc_long, aes(pct, activity_label, fill = cancer_label)) +
      geom_col(position = "dodge", width = 0.7) + esmo_fill +
      labs(x = "Participation (%)", y = NULL, title = "Social activities (last 12 months)") +
      esmo_theme + theme(axis.text.y = element_text(size = 12))
  })

  output$plot_social_dist <- renderPlot({
    d <- fdata() %>% filter(!is.na(cancer_label))
    ggplot(d, aes(social_score, fill = cancer_label)) +
      geom_histogram(aes(weight = weight), position = "identity",
                     alpha = 0.6, binwidth = 1, color = "white") +
      esmo_fill + labs(x = "Social score (0-13)", y = "Weighted count") + esmo_theme
  })

  output$plot_social_age <- renderPlot({
    d <- fdata() %>%
      filter(!is.na(age_group), !is.na(cancer_label)) %>%
      group_by(age_group, cancer_label) %>%
      summarise(mean_score = weighted.mean(social_score, w = weight, na.rm = TRUE),
                .groups = "drop")
    ggplot(d, aes(age_group, mean_score, fill = cancer_label)) +
      geom_col(position = "dodge", width = 0.7) + esmo_fill +
      labs(x = "Age group", y = "Mean social score") + esmo_theme
  })

  # ============================================================
  # TAB 4: HEALTH PROFILE
  # ============================================================
  output$vb_fi <- renderText({
    s <- svy() %>% summarise(m = survey_mean(frailty_index, na.rm = TRUE))
    round(s$m, 3)
  })
  output$vb_qol <- renderText({
    s <- svy() %>% summarise(m = survey_mean(qol_proxy, na.rm = TRUE))
    round(s$m, 1)
  })
  output$vb_grip <- renderText({
    s <- svy() %>% summarise(m = survey_mean(grip_strength, na.rm = TRUE))
    round(s$m, 1)
  })

  output$plot_frailty_status <- renderPlot({
    d <- fdata() %>% filter(!is.na(frailty_status), !is.na(cancer_label))
    tab <- d %>%
      group_by(cancer_label, frailty_status) %>%
      summarise(n = sum(weight, na.rm = TRUE), .groups = "drop_last") %>%
      mutate(pct = n / sum(n) * 100)
    ggplot(tab, aes(frailty_status, pct, fill = cancer_label)) +
      geom_col(position = "dodge", width = 0.7) + esmo_fill +
      labs(x = "Frailty Status (Fried)", y = "%") + esmo_theme
  })

  output$plot_fi_dist <- renderPlot({
    d <- fdata() %>% filter(!is.na(frailty_index), !is.na(cancer_label))
    ggplot(d, aes(frailty_index, fill = cancer_label)) +
      geom_density(alpha = 0.5) + esmo_fill +
      geom_vline(xintercept = 0.25, linetype = "dashed", color = "red") +
      annotate("text", x = 0.27, y = 0, label = "Frail threshold", hjust = 0, color = "red") +
      labs(x = "Frailty Index (0-1)", y = "Density") + esmo_theme
  })

  output$plot_comorbidity <- renderPlot({
    d <- fdata() %>% filter(!is.na(cancer_label))
    chronic_vars_ext <- c(
      hypertension = "Hypertension", diabetes = "Diabetes",
      cholesterol = "High cholesterol", heart_failure = "Heart failure",
      arthritis = "Arthritis", copd = "COPD", depression_dx = "Depression",
      osteoporosis = "Osteoporosis", back_problem = "Back problem",
      kidney_disease = "Kidney disease"
    )
    comor_long <- d %>%
      select(cancer_label, weight, all_of(names(chronic_vars_ext))) %>%
      pivot_longer(-c(cancer_label, weight), names_to = "condition", values_to = "value") %>%
      filter(!is.na(value)) %>%
      group_by(cancer_label, condition) %>%
      summarise(pct = weighted.mean(value == 1, w = weight, na.rm = TRUE) * 100, .groups = "drop") %>%
      mutate(condition_label = chronic_vars_ext[condition]) %>%
      mutate(condition_label = factor(condition_label, levels = rev(chronic_vars_ext)))
    ggplot(comor_long, aes(pct, condition_label, fill = cancer_label)) +
      geom_col(position = "dodge", width = 0.7) + esmo_fill +
      labs(x = "Prevalence (%)", y = NULL) +
      esmo_theme + theme(axis.text.y = element_text(size = 12))
  })

  output$plot_depression_qol <- renderPlot({
    d <- fdata() %>% filter(!is.na(cancer_label))
    p1 <- d %>%
      filter(!is.na(depression_score)) %>%
      mutate(dep_cat = ifelse(depression_score >= 12, "CES-D ≥ 12", "CES-D < 12")) %>%
      group_by(cancer_label, dep_cat) %>%
      summarise(n = sum(weight, na.rm = TRUE), .groups = "drop_last") %>%
      mutate(pct = n / sum(n) * 100) %>%
      ggplot(aes(dep_cat, pct, fill = cancer_label)) +
      geom_col(position = "dodge", width = 0.7) + esmo_fill +
      labs(x = NULL, y = "%", title = "Depression Screening") + esmo_theme

    p2 <- d %>%
      filter(!is.na(qol_proxy)) %>%
      ggplot(aes(cancer_label, qol_proxy, fill = cancer_label)) +
      geom_boxplot(width = 0.5, outlier.alpha = 0.3) + esmo_fill +
      labs(x = NULL, y = "QoL Proxy (0-100)", title = "Quality of Life Proxy") +
      esmo_theme + theme(legend.position = "none")

    p1 + p2
  })

  output$plot_adl_iadl <- renderPlot({
    d <- fdata() %>% filter(!is.na(cancer_label))
    p1 <- ggplot(d, aes(factor(pmin(adl_disability, 4)), fill = cancer_label)) +
      geom_bar(aes(weight = weight), position = "dodge") + esmo_fill +
      scale_x_discrete(labels = c("0","1","2","3","4+")) +
      labs(x = "ADL disabilities", y = "Weighted count", title = "ADL") + esmo_theme

    p2 <- ggplot(d, aes(factor(pmin(iadl_disability, 5)), fill = cancer_label)) +
      geom_bar(aes(weight = weight), position = "dodge") + esmo_fill +
      scale_x_discrete(labels = c("0","1","2","3","4","5+")) +
      labs(x = "IADL disabilities", y = "Weighted count", title = "IADL") + esmo_theme

    p1 + p2
  })

  output$plot_bmi <- renderPlot({
    d <- fdata() %>% filter(!is.na(bmi_cat), !is.na(cancer_label))
    tab <- d %>%
      group_by(cancer_label, bmi_cat) %>%
      summarise(n = sum(weight, na.rm = TRUE), .groups = "drop_last") %>%
      mutate(pct = n / sum(n) * 100)
    ggplot(tab, aes(bmi_cat, pct, fill = cancer_label)) +
      geom_col(position = "dodge", width = 0.7) + esmo_fill +
      labs(x = "BMI Category", y = "%") + esmo_theme
  })

  # ============================================================
  # TAB 5: CANCER CARE
  # ============================================================
  tx_vars <- c(
    tx_chemo = "Chemotherapy", tx_surgery = "Surgery",
    tx_radiation = "Radiation", tx_symptom_med = "Symptom medication",
    tx_other = "Other"
  )

  output$plot_treatment <- renderPlot({
    d <- fdata() %>% filter(cancer_dx == 1)
    tx_long <- d %>%
      select(weight, all_of(names(tx_vars))) %>%
      pivot_longer(-weight, names_to = "treatment", values_to = "value") %>%
      filter(!is.na(value), value == 1) %>%
      group_by(treatment) %>%
      summarise(wt_n = sum(weight, na.rm = TRUE), .groups = "drop") %>%
      mutate(treatment_label = tx_vars[treatment],
             pct = wt_n / sum(d$weight, na.rm = TRUE) * 100)
    if (nrow(tx_long) == 0) {
      ggplot() + annotate("text", x = 1, y = 1, label = "No treatment data", size = 6) + theme_void()
    } else {
      ggplot(tx_long, aes(reorder(treatment_label, pct), pct)) +
        geom_col(fill = "#E74C3C", width = 0.6) +
        geom_text(aes(label = paste0(round(pct, 1), "%")), hjust = -0.1, size = 4.5) +
        coord_flip() + labs(x = NULL, y = "% of cancer survivors") + esmo_theme
    }
  })

  output$plot_cancer_status <- renderPlot({
    d <- fdata() %>% filter(cancer_dx == 1, !is.na(cancer_status_label))
    if (nrow(d) == 0) {
      ggplot() + annotate("text", x = 1, y = 1, label = "No data", size = 6) + theme_void()
    } else {
      tab <- d %>%
        group_by(cancer_status_label) %>%
        summarise(n = sum(weight, na.rm = TRUE), .groups = "drop") %>%
        mutate(pct = n / sum(n) * 100)
      ggplot(tab, aes(reorder(cancer_status_label, pct), pct)) +
        geom_col(fill = "#1F3864", width = 0.6) +
        geom_text(aes(label = paste0(round(pct, 1), "%")), hjust = -0.05) +
        coord_flip() + labs(x = NULL, y = "%") + esmo_theme
    }
  })

  output$plot_dx_age_income <- renderPlot({
    d <- fdata() %>%
      filter(cancer_dx == 1, !is.na(cancer_age_dx), cancer_age_dx < 200,
             !is.na(income_quartile))
    if (nrow(d) < 5) {
      ggplot() + annotate("text", x = 1, y = 1, label = "Insufficient data", size = 6) + theme_void()
    } else {
      ggplot(d, aes(income_quartile, cancer_age_dx, fill = income_quartile)) +
        geom_boxplot(width = 0.6) +
        scale_fill_brewer(palette = "Blues") +
        labs(x = "Income Quartile", y = "Age at Diagnosis") +
        esmo_theme + theme(legend.position = "none")
    }
  })

  output$plot_late_dx <- renderPlot({
    d <- fdata() %>%
      filter(cancer_dx == 1, !is.na(late_diagnosis))
    if (nrow(d) < 5) {
      ggplot() + annotate("text", x = 1, y = 1, label = "Insufficient data", size = 6) + theme_void()
    } else {
      tab <- d %>%
        filter(!is.na(region_label)) %>%
        group_by(region_label, late_diagnosis) %>%
        summarise(n = sum(weight, na.rm = TRUE), .groups = "drop_last") %>%
        mutate(pct = n / sum(n) * 100)
      ggplot(tab, aes(region_label, pct, fill = late_diagnosis)) +
        geom_col(position = "stack", width = 0.7) +
        scale_fill_brewer(palette = "YlOrRd", name = "Diagnosis timing") +
        labs(x = NULL, y = "%") + esmo_theme +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
    }
  })

  output$plot_healthcare_access <- renderPlot({
    d <- fdata() %>% filter(!is.na(cancer_label))
    access_vars <- c(
      has_health_plan = "Health plan",
      has_regular_service = "Regular health service",
      easy_appointment = "Easy to get appointment",
      flu_vaccine = "Flu vaccine (12m)"
    )
    acc_long <- d %>%
      select(cancer_label, weight, all_of(names(access_vars))) %>%
      pivot_longer(-c(cancer_label, weight), names_to = "item", values_to = "value") %>%
      filter(!is.na(value)) %>%
      group_by(cancer_label, item) %>%
      summarise(pct = weighted.mean(value == 1, w = weight, na.rm = TRUE) * 100, .groups = "drop") %>%
      mutate(item_label = access_vars[item],
             item_label = factor(item_label, levels = rev(access_vars)))
    ggplot(acc_long, aes(pct, item_label, fill = cancer_label)) +
      geom_col(position = "dodge", width = 0.7) + esmo_fill +
      labs(x = "% Yes", y = NULL) + esmo_theme
  })

  output$plot_hospitalization <- renderPlot({
    d <- fdata() %>% filter(!is.na(cancer_label))
    p1 <- d %>%
      group_by(cancer_label) %>%
      summarise(pct = weighted.mean(hospitalized == 1, w = weight, na.rm = TRUE) * 100,
                .groups = "drop") %>%
      ggplot(aes(cancer_label, pct, fill = cancer_label)) +
      geom_col(width = 0.5) + esmo_fill +
      labs(x = NULL, y = "Hospitalized (%)", title = "12-month hospitalization") +
      esmo_theme + theme(legend.position = "none")

    p2 <- d %>%
      filter(hospitalized == 1, !is.na(hospital_days), hospital_days < 365) %>%
      ggplot(aes(cancer_label, hospital_days, fill = cancer_label)) +
      geom_boxplot(width = 0.5) + esmo_fill +
      labs(x = NULL, y = "Days in hospital", title = "Length of stay") +
      esmo_theme + theme(legend.position = "none")

    p1 + p2
  })

  # ============================================================
  # TAB 6: STATISTICAL MODELS
  # ============================================================

  # Logistic regression
  lr_result <- eventReactive(input$lr_run, {
    s <- svy()
    s_data <- s$variables
    outcome <- input$lr_outcome
    preds <- input$lr_predictors

    s_data$y <- switch(outcome,
      cancer_dx = s_data$cancer_dx,
      frail_binary = s_data$frail_binary,
      adl_any = s_data$adl_any,
      dep_binary = s_data$dep_binary,
      hospitalized = as.integer(s_data$hospitalized == 1)
    )
    s$variables <- s_data

    if (length(preds) == 0) return(NULL)
    fmla <- as.formula(paste("y ~", paste(preds, collapse = " + ")))

    tryCatch({
      fit <- svyglm(fmla, design = s, family = quasibinomial())
      tidy(fit, conf.int = TRUE, exponentiate = TRUE)
    }, error = function(e) NULL)
  })

  output$lr_table <- renderDT({
    res <- lr_result()
    if (is.null(res)) return(NULL)
    res %>%
      mutate(across(where(is.numeric), ~ round(., 3))) %>%
      datatable(options = list(pageLength = 20, dom = "t"), rownames = FALSE)
  })

  output$lr_forest <- renderPlot({
    res <- lr_result()
    if (is.null(res)) return(NULL)
    res_filt <- res %>% filter(term != "(Intercept)")
    ggplot(res_filt, aes(estimate, reorder(term, estimate))) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
      geom_point(size = 3, color = "#1F3864") +
      geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
      labs(x = "Odds Ratio (95% CI)", y = NULL) + esmo_theme
  })

  # Propensity score matching
  psm_result <- eventReactive(input$psm_run, {
    d <- fdata() %>%
      filter(!is.na(cancer_dx), !is.na(age), !is.na(sex),
             !is.na(education_yrs), !is.na(income_pc)) %>%
      mutate(cancer_dx = as.integer(cancer_dx))

    outcome_var <- input$psm_outcome
    d <- d %>% filter(!is.na(.data[[outcome_var]]))

    tryCatch({
      m <- matchit(cancer_dx ~ age + sex + education_yrs + income_pc + region + n_chronic,
                   data = d, method = "nearest", ratio = 3)
      matched <- match.data(m)

      cancer_mean <- mean(matched[[outcome_var]][matched$cancer_dx == 1], na.rm = TRUE)
      control_mean <- mean(matched[[outcome_var]][matched$cancer_dx == 0], na.rm = TRUE)
      tt <- t.test(matched[[outcome_var]] ~ matched$cancer_dx)

      list(
        summary_text = paste0(
          "Matched: ", sum(matched$cancer_dx == 1), " cancer vs ",
          sum(matched$cancer_dx == 0), " controls\n",
          "Cancer mean: ", round(cancer_mean, 2), "\n",
          "Control mean: ", round(control_mean, 2), "\n",
          "Difference: ", round(cancer_mean - control_mean, 2), "\n",
          "p-value: ", format.pval(tt$p.value, digits = 3), "\n",
          "95% CI: [", round(tt$conf.int[1], 2), ", ", round(tt$conf.int[2], 2), "]"
        ),
        matched_data = matched,
        outcome_var = outcome_var
      )
    }, error = function(e) {
      list(summary_text = paste("Error:", e$message), matched_data = NULL, outcome_var = outcome_var)
    })
  })

  output$psm_summary <- renderText({ psm_result()$summary_text })

  output$psm_plot <- renderPlot({
    res <- psm_result()
    if (is.null(res$matched_data)) return(NULL)
    d <- res$matched_data
    d$cancer_label <- ifelse(d$cancer_dx == 1, "Cancer", "No Cancer")
    ggplot(d, aes(cancer_label, .data[[res$outcome_var]], fill = cancer_label)) +
      geom_boxplot(width = 0.5) +
      scale_fill_manual(values = c("No Cancer" = "#4A90D9", "Cancer" = "#E74C3C")) +
      labs(x = NULL, y = res$outcome_var, title = "Matched Comparison") +
      esmo_theme + theme(legend.position = "none")
  })

  # ============================================================
  # TAB 7: ML & CLUSTERS
  # ============================================================

  # Clustering
  clust_result <- eventReactive(input$clust_run, {
    d <- fdata()
    vars <- input$clust_vars
    k <- input$clust_k

    d_clust <- d %>% select(all_of(vars)) %>% na.omit()
    if (nrow(d_clust) < k * 10) return(NULL)

    d_scaled <- scale(d_clust)
    km <- kmeans(d_scaled, centers = k, nstart = 25, iter.max = 100)

    list(km = km, data = d_clust, scaled = d_scaled, vars = vars, k = k)
  })

  output$clust_pca <- renderPlot({
    res <- clust_result()
    if (is.null(res)) return(NULL)
    fviz_cluster(res$km, data = res$scaled, geom = "point", alpha = 0.3,
                 palette = "Set2", ggtheme = esmo_theme)
  })

  output$clust_profiles <- renderPlot({
    res <- clust_result()
    if (is.null(res)) return(NULL)
    d <- res$data
    d$cluster <- factor(res$km$cluster)
    d_long <- d %>%
      pivot_longer(-cluster, names_to = "variable", values_to = "value") %>%
      group_by(cluster, variable) %>%
      summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")
    ggplot(d_long, aes(variable, mean_val, fill = cluster)) +
      geom_col(position = "dodge") +
      scale_fill_brewer(palette = "Set2", name = "Cluster") +
      labs(x = NULL, y = "Mean value") +
      esmo_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  })

  # SHAP
  shap_result <- eventReactive(input$shap_run, {
    withProgress(message = "Running SHAP analysis...", {
      d <- fdata()
      outcome <- input$shap_outcome

      d$y <- switch(outcome,
        cancer_dx = d$cancer_dx,
        frail_binary = as.integer(!is.na(d$frailty_status) & d$frailty_status == "Frail"),
        dep_binary = as.integer(!is.na(d$depression_score) & d$depression_score >= 12),
        hospitalized = as.integer(d$hospitalized == 1)
      )

      feature_cols <- c("age","sex","education_yrs","income_pc","n_chronic",
                        "social_score","frailty_index","adl_disability","iadl_disability",
                        "depression_score","bmi","grip_strength","has_health_plan",
                        "fall_12m","current_smoker")

      d_model <- d %>%
        select(y, all_of(feature_cols)) %>%
        na.omit()

      if (nrow(d_model) < 100) return(NULL)

      # Subsample for speed
      set.seed(42)
      idx <- sample(nrow(d_model), min(500, nrow(d_model)))
      d_sub <- d_model[idx, ]
      bg_idx <- sample(nrow(d_model), min(100, nrow(d_model)))
      bg <- d_model[bg_idx, ]

      incProgress(0.3)

      rf <- ranger(y ~ ., data = d_sub, num.trees = 200, probability = FALSE,
                   importance = "permutation")

      incProgress(0.3)

      X_sub <- d_sub[, feature_cols]
      X_bg <- bg[, feature_cols]

      tryCatch({
        ks <- kernelshap(rf, X = X_sub, bg_X = X_bg)
        sv <- shapviz(ks)
        incProgress(0.3)
        list(sv = sv, rf = rf)
      }, error = function(e) {
        list(sv = NULL, rf = rf, error = e$message)
      })
    })
  })

  output$shap_importance <- renderPlot({
    res <- shap_result()
    if (is.null(res) || is.null(res$sv)) {
      if (!is.null(res$rf)) {
        imp <- data.frame(
          variable = names(res$rf$variable.importance),
          importance = res$rf$variable.importance
        )
        ggplot(imp, aes(importance, reorder(variable, importance))) +
          geom_col(fill = "#1F3864") +
          labs(x = "Permutation Importance", y = NULL,
               title = "Variable Importance (SHAP unavailable, showing permutation)") +
          esmo_theme
      } else {
        ggplot() + annotate("text", x = 1, y = 1, label = "Insufficient data", size = 6) + theme_void()
      }
    } else {
      sv_importance(res$sv, kind = "bee") + esmo_theme
    }
  })

  output$shap_dependence <- renderPlot({
    res <- shap_result()
    if (is.null(res) || is.null(res$sv)) return(NULL)
    top_var <- names(sort(colMeans(abs(res$sv$S)), decreasing = TRUE))[1]
    sv_dependence(res$sv, v = top_var) + esmo_theme
  })

  # ============================================================
  # TAB 8: REGIONAL MAP
  # ============================================================
  output$map_choropleth <- renderPlot({
    if (!has_map) {
      ggplot() + annotate("text", x = 1, y = 1,
        label = "Map data unavailable. Install geobr.", size = 5) + theme_void()
      return()
    }

    d <- fdata()
    indicator <- input$map_indicator

    reg_stats <- d %>%
      group_by(region_label) %>%
      summarise(
        cancer_prev = weighted.mean(cancer_dx == 1, w = weight, na.rm = TRUE) * 100,
        frail_prev = weighted.mean(!is.na(frailty_status) & frailty_status == "Frail",
                                    w = weight, na.rm = TRUE) * 100,
        mean_social = weighted.mean(social_score, w = weight, na.rm = TRUE),
        adl_prev = weighted.mean(adl_disability >= 1, w = weight, na.rm = TRUE) * 100,
        hosp_prev = weighted.mean(hospitalized == 1, w = weight, na.rm = TRUE) * 100,
        plan_prev = weighted.mean(has_health_plan == 1, w = weight, na.rm = TRUE) * 100,
        mean_fi = weighted.mean(frailty_index, w = weight, na.rm = TRUE),
        dep_prev = weighted.mean(!is.na(depression_score) & depression_score >= 12,
                                  w = weight, na.rm = TRUE) * 100,
        .groups = "drop"
      )

    map_data <- br_regions %>% left_join(reg_stats, by = "region_label")

    label_map <- c(
      cancer_prev = "Cancer (%)", frail_prev = "Frailty (%)",
      mean_social = "Mean Social Score", adl_prev = "ADL Disability (%)",
      hosp_prev = "Hospitalization (%)", plan_prev = "Health Plan (%)",
      mean_fi = "Frailty Index", dep_prev = "Depression (%)"
    )

    ggplot(map_data) +
      geom_sf(aes(fill = .data[[indicator]]), color = "white", size = 0.5) +
      scale_fill_distiller(palette = "YlOrRd", direction = 1,
                           name = label_map[indicator]) +
      geom_sf_text(aes(label = paste0(region_label, "\n",
                        round(.data[[indicator]], 1))),
                   size = 3.5, color = "black") +
      theme_void() +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) +
      labs(title = paste("Brazil:", label_map[indicator]))
  })

  output$map_table <- renderDT({
    d <- fdata()
    reg_stats <- d %>%
      group_by(Region = region_label) %>%
      summarise(
        N = n(),
        `Cancer %` = round(weighted.mean(cancer_dx == 1, w = weight, na.rm = TRUE) * 100, 1),
        `Frail %` = round(weighted.mean(!is.na(frailty_status) & frailty_status == "Frail",
                                         w = weight, na.rm = TRUE) * 100, 1),
        `Mean Social` = round(weighted.mean(social_score, w = weight, na.rm = TRUE), 1),
        `ADL Disab %` = round(weighted.mean(adl_disability >= 1, w = weight, na.rm = TRUE) * 100, 1),
        `Hosp %` = round(weighted.mean(hospitalized == 1, w = weight, na.rm = TRUE) * 100, 1),
        `Health Plan %` = round(weighted.mean(has_health_plan == 1, w = weight, na.rm = TRUE) * 100, 1),
        `Frailty Index` = round(weighted.mean(frailty_index, w = weight, na.rm = TRUE), 3),
        .groups = "drop"
      )
    datatable(reg_stats, options = list(dom = "t", pageLength = 5), rownames = FALSE)
  })

  # ============================================================
  # TAB 9: TIME TO TREATMENT
  # ============================================================
  if (has_apac) {

    # National cohort already has ttt, cancer_group, stage_label, etc.
    apac_proc <- local({
      d <- apac
      d$ttt_days <- d$ttt
      d <- d[!is.na(d$ttt_days) & d$ttt_days >= 0 & d$ttt_days <= 730, ]
      d$time_band <- cut(d$ttt_days,
        breaks = c(-Inf, 30, 60, 180, Inf),
        labels = c("0-30d", "31-60d", "61-180d", ">180d"),
        right  = TRUE)
      d
    })

    output$ttt_n <- renderText({
      format(nrow(apac_proc), big.mark = ",")
    })

    output$ttt_median <- renderText({
      round(median(apac_proc$ttt_days, na.rm = TRUE))
    })

    output$ttt_comp60 <- renderText({
      pct <- mean(apac_proc$ttt_days <= 60, na.rm = TRUE) * 100
      paste0(round(pct, 1), "%")
    })

    output$ttt_delay180 <- renderText({
      pct <- mean(apac_proc$ttt_days > 180, na.rm = TRUE) * 100
      paste0(round(pct, 1), "%")
    })

    output$ttt_km <- renderPlot({
      d <- apac_proc
      strat <- input$ttt_strat

      # Cap at 365 days for display
      d$ttt_cap <- pmin(d$ttt_days, 365)
      d$event   <- 1L  # all observed (reached chemo)

      km_df_list <- list()

      if (strat == "overall" || !(strat %in% names(d))) {
        fit <- survfit(Surv(ttt_cap, event) ~ 1, data = d)
        km_df_list[["Overall"]] <- data.frame(
          time  = fit$time,
          surv  = fit$surv,
          group = "Overall"
        )
      } else {
        d[[strat]] <- as.factor(d[[strat]])
        fmla <- as.formula(paste("Surv(ttt_cap, event) ~", strat))
        fit  <- survfit(fmla, data = d)
        strat_levels <- levels(d[[strat]])
        for (lv in strat_levels) {
          idx <- which(names(fit$strata) == paste0(strat, "=", lv))
          if (length(idx) == 0) next
          # Determine row range for this stratum
          starts <- c(1, cumsum(fit$strata) + 1)
          ends   <- cumsum(fit$strata)
          rows   <- starts[idx]:ends[idx]
          km_df_list[[lv]] <- data.frame(
            time  = fit$time[rows],
            surv  = fit$surv[rows],
            group = lv
          )
        }
      }

      km_df <- do.call(rbind, km_df_list)
      km_df$group <- factor(km_df$group)

      n_groups <- length(levels(km_df$group))
      palette  <- if (n_groups <= 2) {
        c("#1F3864", "#E74C3C")
      } else {
        RColorBrewer::brewer.pal(max(3, n_groups), "Set1")[seq_len(n_groups)]
      }

      ggplot(km_df, aes(time, surv, color = group)) +
        geom_step(size = 1) +
        geom_vline(xintercept = 60, linetype = "dashed", color = "grey40", size = 0.7) +
        annotate("text", x = 63, y = 0.98, label = "60d", hjust = 0,
                 color = "grey40", size = 3.5) +
        scale_color_manual(values = palette, name = NULL) +
        scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
        scale_x_continuous(breaks = seq(0, 365, 60)) +
        labs(
          x        = "Days from diagnosis to first chemotherapy",
          y        = "Proportion not yet treated",
          subtitle = "Event = first SUS chemotherapy session; capped at 365 days"
        ) +
        esmo_theme
    })

    output$ttt_bands <- renderPlot({
      d <- apac_proc
      band_tab <- d %>%
        count(time_band) %>%
        mutate(pct = n / sum(n) * 100)

      band_colors <- c(
        "0-30d"   = "#2ECC71",
        "31-60d"  = "#3498DB",
        "61-180d" = "#F39C12",
        ">180d"   = "#E74C3C"
      )

      ggplot(band_tab, aes(time_band, pct, fill = time_band)) +
        geom_col(width = 0.65) +
        geom_text(aes(label = paste0(round(pct, 1), "%")),
                  vjust = -0.4, size = 5, fontface = "bold") +
        scale_fill_manual(values = band_colors, guide = "none") +
        scale_y_continuous(limits = c(0, max(band_tab$pct) * 1.15),
                           labels = scales::percent_format(scale = 1)) +
        labs(x = "Time band", y = "% of patients") +
        esmo_theme
    })

  } else {
    # APAC not loaded — render empty outputs so Shiny doesn't error
    output$ttt_n       <- renderText("—")
    output$ttt_median  <- renderText("—")
    output$ttt_comp60  <- renderText("—")
    output$ttt_delay180 <- renderText("—")
    output$ttt_km      <- renderPlot({
      ggplot() +
        annotate("text", x = 1, y = 1, label = "APAC data not found", size = 6) +
        theme_void()
    })
    output$ttt_bands   <- renderPlot({
      ggplot() +
        annotate("text", x = 1, y = 1, label = "APAC data not found", size = 6) +
        theme_void()
    })
  }

  # ============================================================
  # TAB 10: AI ANALYST
  # ============================================================
  ai_chat <- reactiveValues(response = NULL, status = "")

  observeEvent(input$ai_submit, {
    req(input$ai_question)
    ai_chat$status <- "Thinking..."
    ai_chat$response <- NULL

    d <- fdata()
    context <- paste0(
      "ELSI-Brazil Wave 2 dataset. Nationally representative survey of Brazilians aged 50+.\n",
      "Current sample after filters: n=", nrow(d), "\n",
      "Variables available (186 columns): demographics (age, sex, race, education, income, region, urban/rural), ",
      "cancer (diagnosis, age at dx, treatments, status), ",
      "14 chronic conditions, Fried frailty phenotype (5 criteria), ",
      "deficit-accumulation frailty index, ADL/IADL batteries, ",
      "social participation (13 activities), CES-D 8 depression, ",
      "physical activity, BMI/anthropometry, grip strength, gait speed, ",
      "healthcare utilization (insurance, doctor visits, hospitalization), ",
      "falls, incontinence, sleep quality, smoking, alcohol.\n\n",
      "Key stats (current filter):\n",
      "- Cancer prevalence: ", round(mean(d$cancer_dx == 1, na.rm = TRUE) * 100, 1), "%\n",
      "- Mean age: ", round(mean(d$age, na.rm = TRUE), 1), "\n",
      "- Female: ", round(mean(d$sex == 0, na.rm = TRUE) * 100, 1), "%\n",
      "- Frail (Fried): ", round(mean(d$frailty_status == "Frail", na.rm = TRUE) * 100, 1), "%\n",
      "- Mean social score: ", round(mean(d$social_score, na.rm = TRUE), 1), "/13\n",
      "- ADL disability: ", round(mean(d$adl_disability >= 1, na.rm = TRUE) * 100, 1), "%\n",
      "- Hospitalized 12m: ", round(mean(d$hospitalized == 1, na.rm = TRUE) * 100, 1), "%\n",
      "- Health plan: ", round(mean(d$has_health_plan == 1, na.rm = TRUE) * 100, 1), "%\n",
      "- Mean frailty index: ", round(mean(d$frailty_index, na.rm = TRUE), 3), "\n",
      "- Depression (CES-D>=12): ", round(mean(!is.na(d$depression_score) & d$depression_score >= 12, na.rm = TRUE) * 100, 1), "%\n",
      "\nSurvey weights (peso_calibrado) are available for nationally representative estimates. ",
      "Complex survey design with strata and PSU.\n",
      "The app has: Overview, Inequalities, Social Network, Health Profile (frailty), ",
      "Cancer Care, Statistical Models (survey logistic, PSM), ML/Clusters (SHAP), Regional Maps."
    )

    tryCatch({
      api_key <- Sys.getenv("ANTHROPIC_API_KEY")
      if (nchar(api_key) < 10) {
        ai_chat$response <- "ANTHROPIC_API_KEY not found. Set it in your environment."
        ai_chat$status <- "Error: No API key"
        return()
      }

      resp <- request("https://api.anthropic.com/v1/messages") |>
        req_headers(
          "x-api-key" = api_key,
          "anthropic-version" = "2023-06-01",
          "content-type" = "application/json"
        ) |>
        req_body_json(list(
          model = "claude-sonnet-4-20250514",
          max_tokens = 4096,
          system = paste0(
            "You are an expert epidemiologist and biostatistician analyzing ELSI-Brazil data. ",
            "You help researchers generate hypotheses, interpret results, suggest analyses, ",
            "and provide context from the aging/geriatric/oncology literature. ",
            "Be specific, cite relevant frameworks (Fried frailty, WHO ICOPE, Andersen model), ",
            "and suggest concrete statistical approaches with R code when relevant.\n\n",
            "DATA CONTEXT:\n", context
          ),
          messages = list(list(role = "user", content = input$ai_question))
        )) |>
        req_timeout(60) |>
        req_perform()

      body <- resp_body_json(resp)
      ai_chat$response <- body$content[[1]]$text
      ai_chat$status <- "Done"
    }, error = function(e) {
      ai_chat$response <- paste("API error:", e$message)
      ai_chat$status <- "Error"
    })
  })

  output$ai_status <- renderText(ai_chat$status)

  output$ai_response <- renderUI({
    if (is.null(ai_chat$response)) return(NULL)
    tags$div(
      style = "background: #f8f9fa; border-radius: 8px; padding: 16px; margin-top: 12px;",
      tags$h5(icon("robot"), " Claude's Analysis"),
      tags$div(style = "white-space: pre-wrap; font-family: inherit;",
               HTML(gsub("\n", "<br>", ai_chat$response)))
    )
  })
}

shinyApp(ui, server)
