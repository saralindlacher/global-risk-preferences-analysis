# ==============================================================================
# Cultural Risk Attitudes and Stock Market Volatility
# ------------------------------------------------------------------------------
# Author : Sara Lindlacher Naveira
# Context: Master's thesis (TFM), MSc Financial Risk Management,
#          ICADE - Universidad Pontificia Comillas, 2026
#
# Purpose: Builds national risk-attitude profiles from World Values Survey (WVS)
#          microdata, merges them with the Global Preference Survey (GPS) and
#          World Bank indicators, and estimates their association with stock
#          market volatility (2012 cross-section). Includes diagnostics,
#          robustness checks and an out-of-sample replication on WVS Wave 7.
#
# Inputs : data/raw/   WVS Wave 6, WVS Wave 7, GPS country file
#                      (not included in the repository - see README)
#          World Bank API (downloaded at runtime, internet required)
# Outputs: output/tables/   LaTeX regression and descriptive tables
#          output/figures/  PNG figures
#          output/sessionInfo.txt
#
# Run    : source("cultural_risk_volatility.R") from the repository root.
#          The bootstrap (5,000 resamples) takes a few minutes.
#
# Note   : The longitudinal panel analysis (thesis Section 4.4.2) was run
#          separately on the WVS time-series file and is not included here.
#
# Sections
#   0      Setup and configuration
#   1-3    WVS Wave 6: load, clean, aggregate, explore
#   4-5    GPS data and WVS-GPS merge
#   6-7    World Bank controls and dependent variable
#   8-9    Analysis samples and merged-sample baseline
#   10-15  Main models, diagnostics, log specification, Cyprus check, HC3 SEs
#   16-19  Construct validity, skewness, descriptives, LaTeX export
#   20-22  GDP sensitivity, figures, Cyprus-excluded log table
#   23-29  Economic magnitudes, power, alternative DVs, heterogeneity,
#          bootstrap, delta R-squared, structural controls
#   30-31  WVS Wave 7 replication (2018) and balanced-sample comparison
#   32     Session information
# ==============================================================================


# ---- 0. Setup and configuration ----------------------------------------------

library(tidyverse)
library(moments)      # skewness
library(wbstats)      # World Bank API (standard indicators)
library(httr)         # direct API calls (WGI indicators, structural controls)
library(jsonlite)     # JSON parsing
library(haven)        # read Stata .dta files
library(countrycode)  # ISO country code conversion
library(car)          # variance inflation factors
library(lmtest)       # Breusch-Pagan test, coeftest
library(sandwich)     # HC3 robust covariance matrices
library(stargazer)    # regression tables
library(ggrepel)      # non-overlapping country labels
library(patchwork)    # multi-panel figures
library(boot)         # bootstrap

# All file paths are relative to the repository root, so the script runs
# unchanged on any machine once the raw data are placed in data/raw/.
PATHS <- list(
  wvs_w6  = file.path("data", "raw", "WV6_Data_R_v20201117.rdata"),
  wvs_w7  = file.path("data", "raw", "WVS_Cross-National_Wave_7_Rdata_v6_0.rdata"),
  gps     = file.path("data", "raw", "country_v11.dta"),
  tables  = file.path("output", "tables"),
  figures = file.path("output", "figures")
)

YEAR_MAIN <- 2012       # WVS Wave 6 cross-section (GPS fielded in 2012)
YEAR_REPL <- 2018       # WVS Wave 7 replication (modal fieldwork year)
N_BOOT    <- 5000       # bootstrap resamples
SEED      <- 20260615   # fixed seed so bootstrap results are reproducible

invisible(lapply(c(PATHS$tables, PATHS$figures), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

options(timeout = 120)
httr::set_config(httr::timeout(120))

# Fail early with a clear message if the raw data are missing.
missing_files <- unlist(PATHS[c("wvs_w6", "wvs_w7", "gps")])
missing_files <- missing_files[!file.exists(missing_files)]
if (length(missing_files) > 0) {
  stop("Raw data not found. Place these files in data/raw/ (see README):\n  ",
       paste(missing_files, collapse = "\n  "))
}

# ---- Shared labels for tables -------------------------------------------------
# Defined once and reused by every table. Text and LaTeX versions differ only
# in escaped characters.
MODEL_LABELS <- c("M1: Trust", "M2: Gov Resp", "M3: Competition",
                  "M4: Security", "M5: GPS Risk")

COV_LABELS_TXT <- c("WVS: Social Trust (z)",
                    "WVS: Indiv. Responsibility (z)",
                    "WVS: Pro-Competition (z)",
                    "WVS: Low Security Need (z)",
                    "GPS: Financial Risk",
                    "GDP per Capita (PPP, 2017$)",
                    "Rule of Law (WGI)",
                    "Constant")

COV_LABELS_TEX <- c("WVS: Social Trust (z)",
                    "WVS: Indiv.\\ Responsibility (z)",
                    "WVS: Pro-Competition (z)",
                    "WVS: Low Security Need (z)",
                    "GPS: Financial Risk",
                    "GDP per Capita (PPP, 2017\\$)",
                    "Rule of Law (WGI)",
                    "Constant")

# ---- Helper functions ---------------------------------------------------------

# Load an .rdata file and return its first object, whatever it is named.
load_rdata <- function(path) {
  env <- new.env()
  obj_names <- load(path, envir = env)
  env[[obj_names[1]]]
}

# Fetch one World Bank indicator for one year straight from the v2 API.
# Used for WGI and structural indicators: wbstats forces footnote = "y",
# which makes WGI requests fail. Aggregate regions have no ISO3 code and
# are dropped.
fetch_wb_api <- function(indicator_code, varname, year) {
  resp <- httr::GET(
    paste0("https://api.worldbank.org/v2/en/country/all/indicator/", indicator_code),
    query = list(date = as.character(year), per_page = "20000", format = "json"),
    httr::timeout(60)
  )
  httr::stop_for_status(resp, task = paste("download", indicator_code))
  raw <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"))[[2]]

  out <- tibble(iso3c = raw$countryiso3code, !!varname := raw$value) %>%
    filter(!is.na(.data[[varname]]), iso3c != "")
  cat(sprintf("  %-22s %d: %d countries\n", indicator_code, year, nrow(out)))
  out
}

# HC3 coefficient table. HC3 is the most conservative of the HC family and
# the recommended choice for small cross-country samples.
hc3_table <- function(model) {
  lmtest::coeftest(model, vcov = sandwich::vcovHC(model, type = "HC3"))
}

# HC3 standard errors and p-values as vectors, for stargazer's se = and p =.
# Passing p = explicitly matters: with se = alone, stargazer recomputes
# p-values from a normal approximation, so the stars would not match HC3
# t-tests.
robust_se <- function(model) hc3_table(model)[, "Std. Error"]
robust_p  <- function(model) hc3_table(model)[, "Pr(>|t|)"]

# HC3 estimate, SE and p-value for a single predictor.
get_hc3 <- function(model, var) {
  ct <- hc3_table(model)
  c(beta = ct[var, "Estimate"], se = ct[var, "Std. Error"], p = ct[var, "Pr(>|t|)"])
}

# Significance stars from a p-value (10/5/1% convention).
p_stars <- function(p) {
  dplyr::case_when(p < 0.01 ~ "***", p < 0.05 ~ "**", p < 0.10 ~ "*", TRUE ~ "")
}

# Print a section banner to the console.
banner <- function(title) {
  cat("\n==============================================================\n")
  cat(" ", title, "\n")
  cat("==============================================================\n")
}

# Regression diagnostics: multicollinearity, heteroskedasticity, normality.
run_diagnostics <- function(model, model_name) {
  cat("\n--------------------------------------------------------------\n")
  cat("  Diagnostics:", model_name, "\n")
  cat("--------------------------------------------------------------\n")
  cat("\n[VIF - multicollinearity]\n");          print(round(vif(model), 3))
  cat("\n[Breusch-Pagan - heteroskedasticity]\n"); print(bptest(model))
  cat("\n[Shapiro-Wilk - residual normality]\n");  print(shapiro.test(residuals(model)))
}

# Estimate the five baseline specifications (M1-M5) for a given outcome.
# M1-M4 use the WVS sample and M5 uses the GPS sample. Each sample is kept
# as large as possible rather than restricted to the WVS-GPS overlap.
fit_five_models <- function(dv, wvs_data, gps_data) {
  fit <- function(predictor, data) {
    f <- reformulate(c(predictor, "gdp_per_capita", "rule_of_law"), response = dv)
    # Insert the formula into the call itself, so stargazer reads a clean
    # dependent-variable name from each model.
    eval(bquote(lm(.(f), data = data)))
  }
  list(
    fit("z_risk_trust",       wvs_data),
    fit("z_risk_gov_resp",    wvs_data),
    fit("z_risk_competition", wvs_data),
    fit("z_risk_security",    wvs_data),
    fit("gps_risk",           gps_data)
  )
}

PREDICTORS <- c("z_risk_trust", "z_risk_gov_resp", "z_risk_competition",
                "z_risk_security", "gps_risk")

# Standard five-model stargazer table. If robust = TRUE, both the HC3 standard
# errors and the matching HC3 p-values are passed to stargazer.
five_model_table <- function(models, type, title, dv_label,
                             robust = FALSE, label = "") {
  stargazer(
    models,                                   # a list of lm objects
    type             = type,
    title            = title,
    label            = label,
    dep.var.caption  = "Dependent Variable:",
    dep.var.labels   = dv_label,
    column.labels    = MODEL_LABELS,
    covariate.labels = if (type == "latex") COV_LABELS_TEX else COV_LABELS_TXT,
    se               = if (robust) lapply(models, robust_se) else NULL,
    p                = if (robust) lapply(models, robust_p)  else NULL,
    omit.stat        = c("f", "ser"),
    digits           = 3,
    star.cutoffs     = c(0.10, 0.05, 0.01)
  )
}

# Write a LaTeX table to output/tables/.
save_latex <- function(tex_lines, filename) {
  path <- file.path(PATHS$tables, filename)
  writeLines(tex_lines, path)
  cat("  saved:", path, "\n")
}


# ---- 1. Load WVS Wave 6 microdata --------------------------------------------

wvs_w6 <- load_rdata(PATHS$wvs_w6)
cat(sprintf("WVS Wave 6 loaded: %d respondents x %d variables\n",
            nrow(wvs_w6), ncol(wvs_w6)))


# ---- 2. Select, clean and recode variables -----------------------------------
# Coding follows the WVS-6 codebook. Every risk proxy is rescaled so that
# HIGHER = MORE RISK-TOLERANT, which keeps the signs comparable across models.
#   V24  trust:          1 = most people can be trusted, 2 = need to be careful
#   V98  responsibility: 1 = government should take more responsibility,
#                        10 = people should take more responsibility
#   V99  competition:    1 = competition is good, 10 = competition is harmful
#   V72  security:       "living in secure surroundings is important to this
#                        person"; 1 = very much like me, 6 = not at all like me
#   V258 survey weight;  V242 age;  V248 education (1-9);  V239 income (1-10)

wvs_clean <- wvs_w6 %>%
  select(
    country    = V2,
    weight     = V258,
    v24_raw    = V24,
    v98_raw    = V98,
    v99_raw    = V99,
    v72_raw    = V72,
    age_raw    = V242,
    edu_raw    = V248,
    income_raw = V239
  ) %>%
  # WVS stores non-response (don't know, no answer, not asked, ...) as
  # negative codes. Setting them to NA keeps them out of the averages.
  mutate(across(ends_with("_raw"), ~ ifelse(. < 0, NA, .))) %>%
  mutate(
    risk_trust       = ifelse(v24_raw == 1, 1, 0),  # 1 = trusts others
    risk_gov_resp    = v98_raw,                     # 10 = individual responsibility
    risk_competition = 11 - v99_raw,                # reversed: 10 = pro-competition
    risk_security    = v72_raw,                     # 6 = low need for security
    high_education   = ifelse(edu_raw >= 8, 1, 0),  # some university or above
    high_income      = ifelse(income_raw >= 8, 1, 0)  # top three income steps
  )


# ---- 3. Aggregate to national profiles and standardise -----------------------
# Each country's respondents are collapsed to survey-weighted means (V258), so
# the aggregates represent the national population rather than the raw sample.
# With na.rm = TRUE, each item uses all valid answers to it: a respondent who
# skipped one question still counts towards the others. Nothing is imputed.

national_profiles <- wvs_clean %>%
  group_by(country) %>%
  summarise(
    avg_risk_trust       = weighted.mean(risk_trust,       w = weight, na.rm = TRUE),
    avg_risk_gov_resp    = weighted.mean(risk_gov_resp,    w = weight, na.rm = TRUE),
    avg_risk_competition = weighted.mean(risk_competition, w = weight, na.rm = TRUE),
    avg_risk_security    = weighted.mean(risk_security,    w = weight, na.rm = TRUE),
    control_age          = weighted.mean(age_raw,          w = weight, na.rm = TRUE),
    prop_high_edu        = weighted.mean(high_education,   w = weight, na.rm = TRUE),
    prop_high_income     = weighted.mean(high_income,      w = weight, na.rm = TRUE),
    survey_n             = n(),
    .groups = "drop"
  ) %>%
  # The proxies use different scales (0-1, 1-6, 1-10). Converting to z-scores
  # makes a coefficient read as the effect of a one-SD difference in national
  # attitudes. The z-scores are computed across all Wave 6 countries, not only
  # the later regression sample, so they keep the full cross-country reference.
  mutate(
    z_risk_trust       = as.numeric(scale(avg_risk_trust)),
    z_risk_gov_resp    = as.numeric(scale(avg_risk_gov_resp)),
    z_risk_competition = as.numeric(scale(avg_risk_competition)),
    z_risk_security    = as.numeric(scale(avg_risk_security))
  )

cat(sprintf("National profiles built: %d countries from %d respondents\n",
            nrow(national_profiles), nrow(wvs_clean)))

# Exploratory analysis: distributions, skewness and correlations of the proxies.
summary(national_profiles)

p0_distributions <- national_profiles %>%
  select(avg_risk_trust, avg_risk_gov_resp, avg_risk_competition, avg_risk_security) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  mutate(variable = case_match(variable,
    "avg_risk_trust"       ~ "Trust (0-1)",
    "avg_risk_gov_resp"    ~ "Individual Responsibility (1-10)",
    "avg_risk_competition" ~ "Pro-Competition (1-10)",
    "avg_risk_security"    ~ "Low Security Need (1-6)"
  )) %>%
  ggplot(aes(x = value)) +
  geom_histogram(fill = "steelblue", colour = "white", bins = 15) +
  facet_wrap(~ variable, scales = "free") +
  labs(title = "Distribution of Risk Proxy Variables (National Averages)",
       x = NULL, y = "Count") +
  theme_minimal()
print(p0_distributions)

skew_values <- national_profiles %>%
  summarise(across(c(avg_risk_trust, avg_risk_gov_resp,
                     avg_risk_competition, avg_risk_security), skewness)) %>%
  unlist()

wvs_cor_matrix <- cor(national_profiles %>% select(starts_with("avg_risk")),
                      use = "complete.obs")
print(round(wvs_cor_matrix, 2))


# ---- 4. GPS country-level data -----------------------------------------------
# The GPS risk-taking measure combines incentivised lottery choices with
# survey items (Falk et al., 2018). It is the experimentally validated
# benchmark against which the attitudinal WVS proxies are compared.

gps_clean <- read_dta(PATHS$gps) %>%
  select(
    country_name = country,
    iso3c        = isocode,
    gps_risk     = risktaking,
    gps_trust    = trust
  )


# ---- 5. Merge WVS and GPS on ISO3 codes --------------------------------------
# WVS identifies countries by ISO numeric code and GPS by ISO alpha-3.
# countrycode() translates them to a common key and warns about any code it
# cannot match.

wvs_profiles_mapped <- national_profiles %>%
  mutate(iso3c = countrycode(country, origin = "iso3n", destination = "iso3c",
                             warn = TRUE))

master_risk_dataset <- inner_join(wvs_profiles_mapped, gps_clean, by = "iso3c") %>%
  select(
    iso3c, country_name, gps_risk,
    z_risk_trust, z_risk_gov_resp, z_risk_competition, z_risk_security,
    avg_risk_trust, avg_risk_gov_resp, avg_risk_competition, avg_risk_security,
    control_age, prop_high_edu, prop_high_income, survey_n, gps_trust
  )

# Merge audit: record which countries are lost in the join, and why.
wvs_unmatched <- anti_join(wvs_profiles_mapped, gps_clean, by = "iso3c")
cat(sprintf("WVS-GPS merge: %d matched | %d WVS countries without GPS data: %s\n",
            nrow(master_risk_dataset), nrow(wvs_unmatched),
            paste(sort(na.omit(wvs_unmatched$iso3c)), collapse = ", ")))

validity_matrix <- master_risk_dataset %>%
  select(gps_risk, z_risk_trust, z_risk_gov_resp, z_risk_competition, z_risk_security) %>%
  cor(use = "complete.obs")


# ---- 6. Macroeconomic controls (World Bank) ----------------------------------
# Controls: GDP per capita (PPP, constant 2017 $) for development level, and
# WGI rule of law for institutional quality.

cat("\nDownloading World Bank controls...\n")

wb_gdp_raw <- wb_data(
  indicator  = c("gdp_per_capita" = "NY.GDP.PCAP.PP.KD"),
  start_date = YEAR_MAIN, end_date = YEAR_MAIN
)

rol_clean <- fetch_wb_api("GOV_WGI_RL.EST", "rule_of_law", YEAR_MAIN)

wb_clean <- wb_gdp_raw %>%
  select(iso3c, gdp_per_capita) %>%
  filter(!is.na(iso3c)) %>%
  left_join(rol_clean, by = "iso3c")

master_risk_dataset <- master_risk_dataset %>%
  left_join(wb_clean, by = "iso3c")


# ---- 7. Dependent variable: stock market volatility --------------------------
# GFDD.SM.01: 360-day standard deviation of the national stock market index
# return (World Bank Global Financial Development Database).

volatility_clean <- wb_data(
  indicator  = c("market_volatility" = "GFDD.SM.01"),
  start_date = YEAR_MAIN, end_date = YEAR_MAIN
) %>%
  select(iso3c, market_volatility) %>%
  filter(!is.na(iso3c))

master_risk_dataset <- master_risk_dataset %>%
  left_join(volatility_clean, by = "iso3c")


# ---- 8. Merged-sample baseline (sensitivity check only) -----------------------
# The WVS-GPS overlap is small, so this model is reported only as a
# sensitivity check. The main specifications below use the larger,
# separate samples.

model_volatility <- lm(
  market_volatility ~ gps_risk + z_risk_gov_resp + gdp_per_capita + rule_of_law,
  data = master_risk_dataset
)
banner("Merged-sample baseline")
summary(model_volatility)


# ---- 9. Analysis samples ------------------------------------------------------
# M1-M4 use every WVS country with volatility data and M5 uses every GPS
# country with volatility data. Countries without an observed outcome are
# dropped explicitly, and the resulting sample sizes are printed.

wb_controls <- wb_clean %>%
  left_join(volatility_clean, by = "iso3c")

wvs_reg_data <- wvs_profiles_mapped %>%
  left_join(wb_controls, by = "iso3c") %>%
  filter(!is.na(market_volatility))

gps_reg_data <- gps_clean %>%
  left_join(wb_controls, by = "iso3c") %>%
  filter(!is.na(market_volatility), !is.na(gps_risk))

cat(sprintf("\nAnalysis samples: WVS N = %d | GPS N = %d\n",
            nrow(wvs_reg_data), nrow(gps_reg_data)))

wvs_cor_unmerged <- cor(
  wvs_reg_data %>% select(z_risk_trust, z_risk_gov_resp, z_risk_competition, z_risk_security),
  use = "complete.obs"
)
print(round(wvs_cor_unmerged, 3))


# ---- 10. Baseline OLS models (level dependent variable) ----------------------

models_lev <- fit_five_models("market_volatility", wvs_reg_data, gps_reg_data)
walk2(models_lev, MODEL_LABELS, ~ run_diagnostics(.x, paste(.y, "(level)")))

banner("Regression results - level models (OLS SEs)")
five_model_table(models_lev, "text",
                 "Cultural Risk Attitudes and National Stock Market Volatility (OLS, 2012)",
                 "Market Volatility (360-Day, WB)")


# ---- 11. Outlier diagnosis: Cyprus -------------------------------------------
# The 2012 Cypriot banking crisis makes Cyprus an extreme volatility outlier.
# The table below confirms this before the log transformation and the
# exclusion check.

banner("Outlier diagnosis - WVS sample ranked by volatility")
wvs_reg_data %>%
  select(iso3c, market_volatility, gdp_per_capita, rule_of_law) %>%
  arrange(desc(market_volatility)) %>%
  mutate(across(c(market_volatility, rule_of_law), ~ round(., 3)),
         gdp_per_capita = round(gdp_per_capita, 0)) %>%
  print(n = Inf)

cat(sprintf("\nSample mean volatility: %.2f | Cyprus: %.2f\n",
            mean(wvs_reg_data$market_volatility, na.rm = TRUE),
            filter(wvs_reg_data, iso3c == "CYP")$market_volatility))


# ---- 12. Log-transformed models (main specification) -------------------------
# Volatility is strongly right-skewed. Taking logs compresses the right tail
# and restores residual normality, and coefficients become approximate
# percentage effects. This is the preferred specification.

wvs_reg_data <- wvs_reg_data %>% mutate(log_volatility = log(market_volatility))
gps_reg_data <- gps_reg_data %>% mutate(log_volatility = log(market_volatility))

models_log <- fit_five_models("log_volatility", wvs_reg_data, gps_reg_data)
walk2(models_log, MODEL_LABELS, ~ run_diagnostics(.x, paste(.y, "(log)")))

banner("Table A - log models (OLS SEs)")
five_model_table(models_log, "text",
                 "Cultural Risk Attitudes and Log Stock Market Volatility (OLS, 2012)",
                 "Log Market Volatility")


# ---- 13. Cyprus-excluded robustness check (level dependent variable) ---------
# Re-estimates the level models without Cyprus, to check that residual
# normality is restored and that the coefficients are stable.

wvs_reg_data_ex <- wvs_reg_data %>% filter(iso3c != "CYP")
gps_reg_data_ex <- gps_reg_data %>% filter(iso3c != "CYP")
cat(sprintf("\nCyprus-excluded samples: WVS N = %d | GPS N = %d\n",
            nrow(wvs_reg_data_ex), nrow(gps_reg_data_ex)))

models_ex <- fit_five_models("market_volatility", wvs_reg_data_ex, gps_reg_data_ex)
walk2(models_ex, MODEL_LABELS, ~ run_diagnostics(.x, paste(.y, "(level, excl. Cyprus)")))

banner("Table B - Cyprus-excluded robustness check (level)")
five_model_table(models_ex, "text",
                 "Robustness Check: Excluding Cyprus (OLS, 2012)",
                 "Market Volatility (360-Day, WB)")


# ---- 14. Shapiro-Wilk summary across specifications --------------------------
# Compares residual normality across three specifications: level, log and
# level excluding Cyprus.

sw_stat <- function(m) unname(shapiro.test(residuals(m))$statistic)
sw_p    <- function(m) shapiro.test(residuals(m))$p.value

shapiro_summary <- tibble(
  Model   = MODEL_LABELS,
  Orig_W  = map_dbl(models_lev, sw_stat), Orig_p  = map_dbl(models_lev, sw_p),
  Log_W   = map_dbl(models_log, sw_stat), Log_p   = map_dbl(models_log, sw_p),
  ExCYP_W = map_dbl(models_ex,  sw_stat), ExCYP_p = map_dbl(models_ex,  sw_p)
) %>%
  mutate(
    Orig_pass  = ifelse(Orig_p  > 0.05, "PASS", "FAIL"),
    Log_pass   = ifelse(Log_p   > 0.05, "PASS", "FAIL"),
    ExCYP_pass = ifelse(ExCYP_p > 0.05, "PASS", "FAIL"),
    across(ends_with("_W"), ~ round(., 4)),
    across(ends_with("_p"), ~ formatC(., format = "e", digits = 3))
  )

banner("Shapiro-Wilk: level vs log vs excl. Cyprus (PASS = p > 0.05)")
print(shapiro_summary, n = Inf, width = Inf)


# ---- 15. HC3 heteroskedasticity-robust standard errors -----------------------
# Every reported inference uses HC3. Both the HC3 SEs and the matching HC3
# p-values are passed to stargazer (see robust_p), so the significance
# stars agree with the tests reported in the text.

banner("Table C - log models with HC3 robust SEs (main results)")
five_model_table(models_log, "text",
                 "Cultural Risk Attitudes and Log Stock Market Volatility (OLS, HC3 Robust SEs, 2012)",
                 "Log Market Volatility", robust = TRUE)

banner("Table D - level models with HC3 robust SEs")
five_model_table(models_lev, "text",
                 "Cultural Risk Attitudes and Market Volatility (OLS, HC3 Robust SEs, 2012)",
                 "Market Volatility (360-Day, WB)", robust = TRUE)

# Change in the cultural predictor's SE when moving from OLS to HC3.
# Changes above 20% are flagged for review.
se_change <- function(m, key) {
  ols <- unname(sqrt(diag(vcov(m)))[key])
  hc3 <- unname(robust_se(m)[key])
  c(ols = round(ols, 4), hc3 = round(hc3, 4), pct = round((hc3 - ols) / ols * 100, 1))
}

se_comparison <- pmap_dfr(
  list(MODEL_LABELS, PREDICTORS, models_log, models_lev),
  function(label, pred, m_log, m_lev) {
    lg <- se_change(m_log, pred); lv <- se_change(m_lev, pred)
    tibble(
      Model          = label,
      Log_OLS_SE     = lg[["ols"]], Log_HC3_SE = lg[["hc3"]], `Log_%_change` = lg[["pct"]],
      Log_flag       = ifelse(abs(lg[["pct"]]) > 20, ">20% REVIEW", ""),
      Lev_OLS_SE     = lv[["ols"]], Lev_HC3_SE = lv[["hc3"]], `Lev_%_change` = lv[["pct"]],
      Lev_flag       = ifelse(abs(lv[["pct"]]) > 20, ">20% REVIEW", "")
    )
  }
)

banner("SE comparison: OLS vs HC3 (cultural predictor)")
print(se_comparison, n = Inf, width = Inf)


# ---- 16. Construct validity: GPS vs WVS proxies ------------------------------
# Do the attitudinal WVS proxies converge with the incentivised GPS measure?
# Rule of thumb: |r| >= 0.3 suggests an overlapping construct; |r| < 0.1
# suggests a different one.

banner("Construct validity: GPS risk vs WVS proxies")
print(round(validity_matrix, 2))

gps_cors <- validity_matrix["gps_risk", setdiff(colnames(validity_matrix), "gps_risk")]
cat("\nInterpretation:\n")
walk2(names(gps_cors), gps_cors, function(var, r) {
  flag <- case_when(
    abs(r) >= 0.3 ~ "partial validity (overlapping construct)",
    abs(r) <  0.1 ~ "low convergence (different construct)",
    TRUE          ~ "weak convergence"
  )
  cat(sprintf("  %-20s r = %+.2f  -> %s\n", var, r, flag))
})


# ---- 17. Skewness of the WVS proxies -----------------------------------------

skew_tbl <- tibble(
  Variable = c("Social Trust", "Indiv. Responsibility",
               "Pro-Competition", "Low Security Need"),
  Skewness = round(skew_values, 3),
  Flag     = ifelse(abs(skew_values) > 1, "HIGH SKEW - consider transformation", "")
)
banner("Skewness of WVS risk proxies (national averages)")
print(skew_tbl, n = Inf, width = Inf)


# ---- 18. Merged-sample sensitivity table and descriptive statistics ----------

banner(sprintf("Sensitivity check - merged WVS+GPS sample (N = %d)", nobs(model_volatility)))
stargazer(
  model_volatility,
  type             = "text",
  title            = "Sensitivity Check: Merged WVS+GPS Sample (OLS, 2012)",
  dep.var.caption  = "Dependent Variable:",
  dep.var.labels   = "Market Volatility (360-Day, WB)",
  column.labels    = sprintf("Merged Sample (N = %d)", nobs(model_volatility)),
  covariate.labels = c("GPS: Financial Risk", "WVS: Indiv. Responsibility (z)",
                       "GDP per Capita (PPP, 2017$)", "Rule of Law (WGI)", "Constant"),
  omit.stat        = c("f", "ser"),
  digits           = 3,
  star.cutoffs     = c(0.10, 0.05, 0.01)
)

wvs_desc_data <- wvs_reg_data %>%
  select(market_volatility, log_volatility,
         z_risk_trust, z_risk_gov_resp, z_risk_competition, z_risk_security,
         gdp_per_capita, rule_of_law) %>%
  as.data.frame()

gps_desc_data <- gps_reg_data %>%
  select(market_volatility, log_volatility, gps_risk, gdp_per_capita, rule_of_law) %>%
  as.data.frame()

desc_labels_wvs <- c("Market Volatility", "Log Market Volatility",
                     COV_LABELS_TXT[1:4], COV_LABELS_TXT[6:7])
desc_labels_gps <- c("Market Volatility", "Log Market Volatility",
                     "GPS: Financial Risk Tolerance", COV_LABELS_TXT[6:7])
desc_title_wvs  <- sprintf("Descriptive Statistics: WVS Regression Sample (2012, N = %d)", nrow(wvs_desc_data))
desc_title_gps  <- sprintf("Descriptive Statistics: GPS Regression Sample (2012, N = %d)", nrow(gps_desc_data))
DESC_STATS      <- c("n", "mean", "sd", "min", "median", "max")

banner(desc_title_wvs)
stargazer(wvs_desc_data, type = "text", title = desc_title_wvs,
          covariate.labels = desc_labels_wvs, digits = 3, summary.stat = DESC_STATS)

banner(desc_title_gps)
stargazer(gps_desc_data, type = "text", title = desc_title_gps,
          covariate.labels = desc_labels_gps, digits = 3, summary.stat = DESC_STATS)


# ---- 19. LaTeX export of thesis tables ---------------------------------------

cat("\nWriting LaTeX tables...\n")

save_latex(capture.output(five_model_table(
  models_log, "latex",
  "Cultural Risk Attitudes and Log Stock Market Volatility (OLS, 2012)",
  "Log Market Volatility", label = "tab:log_models"
)), "tableA_log_models.tex")

save_latex(capture.output(five_model_table(
  models_ex, "latex",
  "Robustness Check: Excluding Cyprus (Level Volatility, OLS, 2012)",
  "Market Volatility (360-Day, WB)", label = "tab:cyprus_robustness_level"
)), "tableB_cyprus_robustness_level.tex")

save_latex(capture.output(five_model_table(
  models_log, "latex",
  "Cultural Risk Attitudes and Log Stock Market Volatility (OLS, HC3 Robust SEs, 2012)",
  "Log Market Volatility", robust = TRUE, label = "tab:log_robust"
)), "tableC_log_robust_ses.tex")

save_latex(capture.output(five_model_table(
  models_lev, "latex",
  "Cultural Risk Attitudes and Market Volatility (OLS, HC3 Robust SEs, 2012)",
  "Market Volatility (360-Day, WB)", robust = TRUE, label = "tab:level_robust"
)), "tableD_level_robust_ses.tex")

save_latex(capture.output(stargazer(
  wvs_desc_data, type = "latex", title = desc_title_wvs, label = "tab:desc_wvs",
  covariate.labels = c("Market Volatility", "Log Market Volatility",
                       COV_LABELS_TEX[1:4], COV_LABELS_TEX[6:7]),
  digits = 3, summary.stat = DESC_STATS
)), "desc_stats_wvs.tex")

save_latex(capture.output(stargazer(
  gps_desc_data, type = "latex", title = desc_title_gps, label = "tab:desc_gps",
  covariate.labels = c("Market Volatility", "Log Market Volatility",
                       "GPS: Financial Risk Tolerance", COV_LABELS_TEX[6:7]),
  digits = 3, summary.stat = DESC_STATS
)), "desc_stats_gps.tex")


# ---- 20. GPS sensitivity to the GDP control ----------------------------------
# Is GPS risk simply a proxy for development? If its coefficient changes by
# less than 20% when GDP is dropped, it captures something independent of
# income level.

model5_log    <- models_log[[5]]
model5_no_gdp <- lm(log_volatility ~ gps_risk + rule_of_law, data = gps_reg_data)

beta_full   <- coef(model5_log)["gps_risk"]
beta_no_gdp <- coef(model5_no_gdp)["gps_risk"]
pct_change  <- abs((beta_no_gdp - beta_full) / beta_full) * 100

banner("GPS sensitivity: M5 with vs without GDP control")
stargazer(
  model5_no_gdp, model5_log,
  type             = "text",
  title            = "M5 Sensitivity: GPS Risk Coefficient With and Without GDP Control (OLS, 2012)",
  dep.var.caption  = "Dependent Variable:",
  dep.var.labels   = "Log Market Volatility",
  column.labels    = c("No GDP", "Full M5"),
  covariate.labels = c("GPS: Financial Risk Tolerance", "Rule of Law (WGI)",
                       "GDP per Capita (PPP, 2017$)", "Constant"),
  omit.stat        = c("f", "ser"),
  digits           = 3,
  star.cutoffs     = c(0.10, 0.05, 0.01)
)

cat(sprintf("\nGPS coefficient: full = %.3f | no GDP = %.3f | change = %.1f%%\n",
            beta_full, beta_no_gdp, pct_change))
cat(ifelse(pct_change < 20,
           "Stable: GPS captures attitudes independent of development level.\n",
           "Sensitive: partial overlap with development level, interpret with caution.\n"))


# ---- 21. Figures --------------------------------------------------------------

income_palette <- c(
  "High Income"      = "#2166ac",
  "Upper Middle"     = "#4dac26",
  "Lower Middle/Low" = "#d01c8b"
)

# Approximate World Bank-style income tiers from GDP per capita (PPP).
classify_income <- function(gdp) {
  case_when(gdp > 25000 ~ "High Income",
            gdp > 10000 ~ "Upper Middle",
            TRUE        ~ "Lower Middle/Low")
}

fmt_p <- function(p) {
  case_when(p < 0.001 ~ "p < 0.001",
            p < 0.01  ~ "p < 0.01",
            p < 0.05  ~ "p < 0.05",
            TRUE      ~ paste0("p = ", round(p, 3)))
}

# Scatter of a cultural predictor against log volatility, with the HC3
# estimate shown in the caption.
scatter_plot <- function(data, model, predictor, title, x_label) {
  est <- get_hc3(model, predictor)
  data %>%
    mutate(income_grp = classify_income(gdp_per_capita)) %>%
    ggplot(aes(x = .data[[predictor]], y = log_volatility, colour = income_grp)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "grey40", fill = "grey85", linewidth = 0.8, alpha = 0.4) +
    geom_point(size = 2.8) +
    geom_text_repel(aes(label = iso3c), size = 2.8, max.overlaps = 20,
                    box.padding = 0.35, segment.colour = "grey60") +
    scale_colour_manual(values = income_palette, name = "Income Group") +
    labs(title = title, x = x_label, y = "Log Market Volatility",
         caption = sprintf("HC3 robust SEs  |  β = %.3f, %s  |  N = %d",
                           est[["beta"]], fmt_p(est[["p"]]), nobs(model))) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          plot.caption    = element_text(size = 8, colour = "grey40"))
}

# Figure 1: GPS risk tolerance vs log volatility
p1 <- scatter_plot(gps_reg_data, models_log[[5]], "gps_risk",
                   "GPS Financial Risk Tolerance and Stock Market Volatility (2012)",
                   "GPS Financial Risk Tolerance")
ggsave(file.path(PATHS$figures, "plot1_gps_risk_vs_volatility.png"),
       p1, width = 8, height = 6, dpi = 300)

# Figure 2: pro-competition attitudes vs log volatility
p2 <- scatter_plot(wvs_reg_data, models_log[[3]], "z_risk_competition",
                   "Pro-Competition Attitudes and Stock Market Volatility (2012)",
                   "WVS: Pro-Competition Attitudes (z-score)")
ggsave(file.path(PATHS$figures, "plot2_competition_vs_volatility.png"),
       p2, width = 8, height = 6, dpi = 300)

# Figure 3: coefficient plot, all five log models, 95% CI from HC3 SEs
coef_plot_data <- pmap_dfr(list(MODEL_LABELS, PREDICTORS, models_log),
                           function(label, pred, m) {
  est <- get_hc3(m, pred)
  tibble(model_label = label, beta = est[["beta"]], se_hc3 = est[["se"]], p_val = est[["p"]])
}) %>%
  mutate(ci_lo       = beta - 1.96 * se_hc3,
         ci_hi       = beta + 1.96 * se_hc3,
         significant = p_val < 0.05,
         model_label = factor(model_label, levels = rev(MODEL_LABELS)))

p3 <- ggplot(coef_plot_data, aes(x = beta, y = model_label, colour = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.7) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                width = 0.2, linewidth = 0.8, orientation = "y") +
  geom_point(size = 4) +
  scale_colour_manual(values = c("TRUE" = "#d7191c", "FALSE" = "#999999"),
                      labels = c("TRUE" = "Significant (p < 0.05)", "FALSE" = "Not significant"),
                      name   = NULL) +
  labs(title    = "Effect of Cultural Risk Attitudes on Log Stock Market Volatility",
       subtitle = "Coefficients with 95% CI (HC3 Robust Standard Errors)",
       x = "Coefficient (β)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(PATHS$figures, "plot3_coefficient_plot.png"),
       p3, width = 8, height = 5, dpi = 300)

# Figure 4: correlation heatmap on the WVS-GPS overlap
var_labels <- c(
  "gps_risk"           = "GPS: Financial Risk",
  "z_risk_trust"       = "WVS: Social Trust (z)",
  "z_risk_gov_resp"    = "WVS: Indiv. Resp. (z)",
  "z_risk_competition" = "WVS: Pro-Comp. (z)",
  "z_risk_security"    = "WVS: Low Security (z)",
  "log_volatility"     = "Log Mkt Volatility"
)

heatmap_raw <- master_risk_dataset %>%
  left_join(wvs_reg_data %>% select(iso3c, log_volatility), by = "iso3c") %>%
  select(all_of(names(var_labels))) %>%
  filter(complete.cases(.))
cat(sprintf("Heatmap N = %d (WVS-GPS overlap)\n", nrow(heatmap_raw)))

p4 <- cor(heatmap_raw) %>%
  as.data.frame() %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "correlation") %>%
  mutate(var1 = factor(var_labels[var1], levels = unname(var_labels)),
         var2 = factor(var_labels[var2], levels = rev(unname(var_labels)))) %>%
  ggplot(aes(x = var1, y = var2, fill = correlation)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = round(correlation, 2), colour = abs(correlation) > 0.6), size = 3.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d7191c",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
  labs(title = "Correlation Matrix: Risk Proxies and Market Volatility", x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
ggsave(file.path(PATHS$figures, "plot4_correlation_heatmap.png"),
       p4, width = 8, height = 6, dpi = 300)

# Combined 2x2 panel
combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title   = "Cultural Risk Attitudes and Stock Market Volatility (2012)",
    caption = "Data: WVS Wave 6, GPS, World Bank. OLS with HC3 robust SEs.",
    theme   = theme(plot.title   = element_text(size = 13, face = "bold", hjust = 0.5),
                    plot.caption = element_text(size = 8, colour = "grey40"))
  )
ggsave(file.path(PATHS$figures, "thesis_plots_combined.png"),
       combined, width = 16, height = 12, dpi = 300)
cat("Figures saved to", PATHS$figures, "\n")


# ---- 22. Cyprus-excluded log models: HC3 LaTeX table -------------------------
# This table is built by hand so that the main thesis table shows HC3 SEs and
# stars on the log scale, in the same layout as the stargazer output.

models_ex_log <- fit_five_models("log_volatility", wvs_reg_data_ex, gps_reg_data_ex)

# One LaTeX cell pair (estimate with stars, SE in brackets) using HC3.
fmt_coef_tex <- function(model, var) {
  ct <- hc3_table(model)
  if (!var %in% rownames(ct)) return(c("", ""))
  est   <- ct[var, "Estimate"]
  stars <- p_stars(ct[var, "Pr(>|t|)"])
  est_txt <- formatC(abs(est), format = "f", digits = 3)
  c(paste0(ifelse(est < 0, "$-$", ""), est_txt, ifelse(stars == "", "", paste0("$^{", stars, "}$"))),
    paste0("(", formatC(ct[var, "Std. Error"], format = "f", digits = 3), ")"))
}

# Estimate row, SE row and spacer row for one variable across all models.
tex_rows <- function(var, label, spacer = TRUE) {
  cells <- map(models_ex_log, fmt_coef_tex, var = var)
  c(paste0(label, " & ", paste(map_chr(cells, 1), collapse = " & "), " \\\\"),
    paste0(" & ",        paste(map_chr(cells, 2), collapse = " & "), " \\\\"),
    if (spacer) "& & & & & \\\\")
}

row_vars   <- c(PREDICTORS, "gdp_per_capita", "rule_of_law", "(Intercept)")
body_rows  <- unlist(map2(row_vars, COV_LABELS_TEX,
                          ~ tex_rows(.x, .y, spacer = .x != "(Intercept)")))

tex_b_log <- c(
  "\\begin{table}[!htbp] \\centering",
  "  \\caption{Robustness Check: Excluding Cyprus (Log Volatility, HC3 Robust SEs, 2012)}",
  "  \\label{tab:cyprus_robustness}",
  "\\begin{tabular}{@{\\extracolsep{5pt}}lccccc}",
  "\\\\[-1.8ex]\\hline",
  "\\hline \\\\[-1.8ex]",
  " & \\multicolumn{5}{c}{Dependent Variable:} \\\\",
  "\\cline{2-6}",
  "\\\\[-1.8ex] & \\multicolumn{5}{c}{Log Market Volatility} \\\\",
  paste0(" & ", paste(MODEL_LABELS, collapse = " & "), " \\\\"),
  "\\\\[-1.8ex] & (1) & (2) & (3) & (4) & (5) \\\\",
  "\\hline \\\\[-1.8ex]",
  body_rows,
  "\\hline \\\\[-1.8ex]",
  paste0("Observations & ", paste(map_int(models_ex_log, nobs), collapse = " & "), " \\\\"),
  paste0("R$^{2}$ & ", paste(map_chr(models_ex_log, ~ formatC(summary(.x)$r.squared,
                                                              format = "f", digits = 3)),
                             collapse = " & "), " \\\\"),
  "\\hline",
  "\\hline \\\\[-1.8ex]",
  "\\textit{Note:}  & \\multicolumn{5}{r}{HC3 robust SEs. $^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01} \\\\",
  "\\end{tabular}",
  "\\end{table}"
)
save_latex(tex_b_log, "tableB_cyprus_robustness.tex")


# ---- 23. Economic magnitudes --------------------------------------------------
# Converts log-scale coefficients into percentage changes in volatility:
# %change = (exp(beta * delta_x) - 1) * 100. Every coefficient is read from
# the fitted models, so these numbers update automatically with the data.

magnitude_report <- function(data, model, predictor, label, delta_one_sd) {
  beta <- get_hc3(model, predictor)[["beta"]]
  ext  <- data %>% arrange(.data[[predictor]]) %>%
    select(iso3c, all_of(predictor), market_volatility) %>% slice(c(1, n()))

  one_sd_pct <- (exp(beta * delta_one_sd) - 1) * 100
  range_pct  <- (exp(beta * diff(range(data[[predictor]], na.rm = TRUE))) - 1) * 100

  banner(sprintf("Economic magnitude - %s (beta = %+.3f)", label, beta))
  print(ext)
  cat(sprintf("One-SD increase:            %+.1f%% change in volatility\n", one_sd_pct))
  cat(sprintf("Full range (%s -> %s):     %+.1f%% change in volatility\n",
              ext$iso3c[1], ext$iso3c[2], range_pct))
  invisible(c(one_sd_pct = one_sd_pct, range_pct = range_pct))
}

# WVS proxies are z-scores, so one SD is a change of 1 unit. GPS risk is in
# its original units, so one SD is its standard deviation in the GPS sample.
m3_magnitude <- magnitude_report(wvs_reg_data, models_log[[3]], "z_risk_competition",
                                 "M3: Pro-Competition", delta_one_sd = 1)
m4_magnitude <- magnitude_report(wvs_reg_data, models_log[[4]], "z_risk_security",
                                 "M4: Low Security Need", delta_one_sd = 1)
m5_magnitude <- magnitude_report(gps_reg_data, models_log[[5]], "gps_risk",
                                 "M5: GPS Risk", delta_one_sd = sd(gps_reg_data$gps_risk, na.rm = TRUE))

# Benchmark: the volatility gap between developed and emerging markets.
benchmark <- gps_reg_data %>%
  mutate(developed = gdp_per_capita > 25000) %>%
  group_by(developed) %>%
  summarise(mean_vol = mean(market_volatility, na.rm = TRUE),
            mean_log_vol = mean(log_volatility, na.rm = TRUE),
            n = n(), .groups = "drop")
banner("Benchmark: developed vs emerging volatility gap")
print(benchmark)
cat(sprintf("Mean log-volatility gap implies a %.1f%% volatility difference\n",
            abs((exp(diff(benchmark$mean_log_vol)) - 1) * 100)))


# ---- 24. Power analysis: minimum detectable effect ---------------------------
# Smallest |beta| detectable at alpha = 0.05 (two-tailed) with 80% power,
# given the observed HC3 SE:  MDE = SE_HC3 * (t_{0.975, df} + t_{0.80, df}).
# This separates "no effect" from "an effect too small for N to detect".

power_analysis <- function(model, predictor, label) {
  est <- get_hc3(model, predictor)
  df  <- model$df.residual
  mde <- est[["se"]] * (qt(0.975, df) + qt(0.80, df))
  ratio <- abs(est[["beta"]]) / mde
  tibble(
    Model  = label,
    Beta   = round(est[["beta"]], 3),
    SE_HC3 = round(est[["se"]], 3),
    df     = df,
    MDE_80 = round(mde, 3),
    Ratio  = round(ratio, 2),
    Status = case_when(
      ratio >= 1.0 ~ "Detected",
      ratio >= 0.5 ~ "Underpowered (between half and full MDE)",
      ratio >= 0.2 ~ "Far below MDE (<50% of detectable threshold)",
      TRUE         ~ "Negligible relative to MDE (likely true null)"
    )
  )
}

results_power <- pmap_dfr(list(models_log, PREDICTORS, MODEL_LABELS), power_analysis)
banner("Power analysis - minimum detectable effect (alpha = 0.05, power = 0.80)")
print(results_power, n = Inf, width = Inf)


# ---- 25. Alternative dependent variables -------------------------------------
# Do M3 and M5 hold for other measures of national financial instability?
#   Alt 1: bank Z-score, inverted so that higher = riskier
#   Alt 2: absolute CPI inflation (nominal instability)
#   Alt 3: SD of real GDP growth over 2010-2014 (macro volatility)
# Every outcome is logged, as in the main specification.

cat("\nDownloading alternative dependent variables...\n")

zscore_raw <- wb_data(indicator = c("bank_zscore" = "GFDD.SI.05"),
                      start_date = YEAR_MAIN, end_date = YEAR_MAIN) %>%
  select(iso3c, bank_zscore) %>%
  filter(!is.na(iso3c), !is.na(bank_zscore)) %>%
  # A higher Z-score means a more stable banking sector. The model logs the
  # outcome, so storing 1/Z gives log(1/Z) = -log(Z): higher = riskier.
  mutate(bank_risk_level = 1 / bank_zscore)

cpi_raw <- wb_data(indicator = c("cpi_inflation" = "FP.CPI.TOTL.ZG"),
                   start_date = YEAR_MAIN, end_date = YEAR_MAIN) %>%
  select(iso3c, cpi_inflation) %>%
  filter(!is.na(iso3c), !is.na(cpi_inflation)) %>%
  mutate(cpi_abs = abs(cpi_inflation))

gdp_growth_raw <- wb_data(indicator = c("gdp_growth" = "NY.GDP.MKTP.KD.ZG"),
                          start_date = YEAR_MAIN - 2, end_date = YEAR_MAIN + 2) %>%
  select(iso3c, date, gdp_growth) %>%
  filter(!is.na(iso3c), !is.na(gdp_growth)) %>%
  group_by(iso3c) %>%
  summarise(gdp_volatility = sd(gdp_growth), n_years = n(), .groups = "drop") %>%
  filter(n_years >= 4)   # require at least 4 of the 5 years

cat(sprintf("  bank Z-score: %d | CPI: %d | GDP growth SD: %d countries\n",
            nrow(zscore_raw), nrow(cpi_raw), nrow(gdp_growth_raw)))

add_alt_dvs <- function(data) {
  data %>%
    left_join(zscore_raw,     by = "iso3c") %>%
    left_join(cpi_raw,        by = "iso3c") %>%
    left_join(gdp_growth_raw, by = "iso3c")
}
wvs_reg_alt <- add_alt_dvs(wvs_reg_data)
gps_reg_alt <- add_alt_dvs(gps_reg_data)

# Regress log(outcome) on one predictor plus the standard controls (HC3).
# Countries with a missing or non-positive outcome are dropped, and the model
# is skipped if fewer than 15 countries remain.
run_alt_regression <- function(data, outcome, predictor, model_label, dv_label) {
  d <- data %>% filter(!is.na(.data[[outcome]]), .data[[outcome]] > 0)
  if (nrow(d) < 15) {
    cat(sprintf("  [%s on %s] skipped: only %d countries\n", model_label, dv_label, nrow(d)))
    return(NULL)
  }
  d   <- d %>% mutate(log_dv = log(.data[[outcome]]))
  fit <- lm(reformulate(c(predictor, "gdp_per_capita", "rule_of_law"), "log_dv"), data = d)
  est <- get_hc3(fit, predictor)
  tibble(Model = model_label, DV = dv_label, N = nobs(fit),
         Beta = round(est[["beta"]], 3), SE_HC3 = round(est[["se"]], 3),
         p_value = round(est[["p"]], 4), Stars = p_stars(est[["p"]]),
         R2 = round(summary(fit)$r.squared, 3))
}

alt_outcomes <- c(
  "Main: Stock Vol (GFDD.SM.01)"     = "market_volatility",
  "Alt 1: Bank Risk (-log Z-score)"  = "bank_risk_level",
  "Alt 2: |CPI Inflation| (2012)"    = "cpi_abs",
  "Alt 3: GDP Growth SD (2010-14)"   = "gdp_volatility"
)

results_alt <- imap_dfr(alt_outcomes, function(outcome, dv_label) {
  bind_rows(
    run_alt_regression(wvs_reg_alt, outcome, "z_risk_competition", "M3: Competition", dv_label),
    run_alt_regression(gps_reg_alt, outcome, "gps_risk",           "M5: GPS Risk",    dv_label)
  )
})

banner("Alternative dependent variables - M3 and M5 (HC3)")
print(results_alt, n = Inf, width = Inf)

robust_summary <- results_alt %>%
  group_by(Model) %>%
  summarise(DVs_tested       = n(),
            DVs_same_sign    = sum(sign(Beta) == sign(first(Beta))),
            DVs_signif_5pct  = sum(p_value < 0.05),
            DVs_signif_10pct = sum(p_value < 0.10),
            .groups = "drop")
print(robust_summary, width = Inf)


# ---- 26. Heterogeneity: OECD vs emerging markets ------------------------------
# Are M3 and M5 driven by one type of market? Sub-sample estimates, plus a
# Chow-style F-test of an interaction between culture and OECD membership.

# OECD members as of 2012 (34 countries)
OECD_2012 <- c(
  "AUS", "AUT", "BEL", "CAN", "CHL", "CZE", "DNK", "EST", "FIN", "FRA",
  "DEU", "GRC", "HUN", "ISL", "IRL", "ISR", "ITA", "JPN", "KOR", "LUX",
  "MEX", "NLD", "NZL", "NOR", "POL", "PRT", "SVK", "SVN", "ESP", "SWE",
  "CHE", "TUR", "GBR", "USA"
)

wvs_reg_oecd <- wvs_reg_data %>% mutate(oecd = iso3c %in% OECD_2012)
gps_reg_oecd <- gps_reg_data %>% mutate(oecd = iso3c %in% OECD_2012)

banner("Heterogeneity: OECD vs emerging markets")
cat(sprintf("WVS sample: %d OECD, %d emerging | GPS sample: %d OECD, %d emerging\n",
            sum(wvs_reg_oecd$oecd), sum(!wvs_reg_oecd$oecd),
            sum(gps_reg_oecd$oecd), sum(!gps_reg_oecd$oecd)))

run_subsample <- function(data, predictor, is_oecd, model_label) {
  d <- data %>% filter(oecd == is_oecd)
  group_label <- ifelse(is_oecd, "OECD", "Emerging")
  if (nrow(d) < 12) {
    cat(sprintf("  [%s, %s] skipped: only %d countries\n", model_label, group_label, nrow(d)))
    return(NULL)
  }
  fit <- lm(reformulate(c(predictor, "gdp_per_capita", "rule_of_law"), "log_volatility"), data = d)
  est <- get_hc3(fit, predictor)
  tibble(Model = model_label, Group = group_label, N = nobs(fit),
         Beta = round(est[["beta"]], 3), SE_HC3 = round(est[["se"]], 3),
         p_value = round(est[["p"]], 4), Stars = p_stars(est[["p"]]),
         R2 = round(summary(fit)$r.squared, 3))
}

heterogeneity <- bind_rows(
  run_subsample(wvs_reg_oecd, "z_risk_competition", TRUE,  "M3: Competition"),
  run_subsample(wvs_reg_oecd, "z_risk_competition", FALSE, "M3: Competition"),
  run_subsample(gps_reg_oecd, "gps_risk",           TRUE,  "M5: GPS Risk"),
  run_subsample(gps_reg_oecd, "gps_risk",           FALSE, "M5: GPS Risk")
)
print(heterogeneity, n = Inf, width = Inf)

# H0: the cultural effect is the same in OECD and emerging markets.
interaction_test <- function(data, predictor, label) {
  base <- reformulate(c(predictor, "oecd", "gdp_per_capita", "rule_of_law"), "log_volatility")
  full <- reformulate(c(paste0(predictor, " * oecd"), "gdp_per_capita", "rule_of_law"), "log_volatility")
  res  <- anova(lm(base, data = data), lm(full, data = data))
  cat(sprintf("\n%s interaction F-test: F = %.2f, p = %.4f\n", label, res$F[2], res$`Pr(>F)`[2]))
  invisible(res)
}
interaction_test(wvs_reg_oecd, "z_risk_competition", "M3 (competition)")
interaction_test(gps_reg_oecd, "gps_risk",           "M5 (GPS)")


# ---- 27. Bootstrap percentile confidence intervals ---------------------------
# A non-parametric check on the HC3 inference: case (row) resampling with
# N_BOOT replications and percentile 95% CIs. The fixed seed makes the
# intervals reproducible.

set.seed(SEED)

bootstrap_ci <- function(data, formula_text, predictor, label, n_boot = N_BOOT) {
  f <- as.formula(formula_text)
  stat_fn <- function(d, i) {
    cf <- coef(lm(f, data = d[i, ]))
    if (predictor %in% names(cf)) cf[[predictor]] else NA_real_
  }
  boot_draws <- boot(data, stat_fn, R = n_boot)$t
  boot_draws <- boot_draws[!is.na(boot_draws)]

  fit    <- lm(f, data = data)
  est    <- get_hc3(fit, predictor)
  t_crit <- qt(0.975, fit$df.residual)
  boot_q <- quantile(boot_draws, c(0.025, 0.975))

  tibble(
    Model       = label,
    Beta        = round(est[["beta"]], 3),
    HC3_lo      = round(est[["beta"]] - t_crit * est[["se"]], 3),
    HC3_hi      = round(est[["beta"]] + t_crit * est[["se"]], 3),
    Boot_lo     = round(unname(boot_q[1]), 3),
    Boot_hi     = round(unname(boot_q[2]), 3),
    Boot_excl_0 = (boot_q[1] > 0) == (boot_q[2] > 0),
    Boot_bias   = round(mean(boot_draws) - est[["beta"]], 3),
    Boot_SD     = round(sd(boot_draws), 3),
    HC3_SE      = round(est[["se"]], 3),
    SD_ratio    = round(sd(boot_draws) / est[["se"]], 2)
  )
}

results_boot <- bind_rows(
  bootstrap_ci(wvs_reg_data, "log_volatility ~ z_risk_competition + gdp_per_capita + rule_of_law",
               "z_risk_competition", "M3: Pro-Competition"),
  bootstrap_ci(wvs_reg_data, "log_volatility ~ z_risk_security + gdp_per_capita + rule_of_law",
               "z_risk_security", "M4: Low Security Need"),
  bootstrap_ci(gps_reg_data, "log_volatility ~ gps_risk + gdp_per_capita + rule_of_law",
               "gps_risk", "M5: GPS Risk")
)

banner(sprintf("Bootstrap percentile 95%% CIs (%d case resamples) vs HC3", N_BOOT))
print(results_boot, n = Inf, width = Inf)


# ---- 28. Marginal explanatory power of culture (delta R-squared) -------------
# Baseline model: log_vol ~ GDP + rule of law.
# Full model:     baseline + one cultural variable.
# Reports the change in R-squared and a nested F-test.

compute_delta_r2 <- function(data, cultural_var, label) {
  m_base <- lm(log_volatility ~ gdp_per_capita + rule_of_law, data = data)
  m_full <- lm(reformulate(c(cultural_var, "gdp_per_capita", "rule_of_law"), "log_volatility"),
               data = data)
  s_base <- summary(m_base); s_full <- summary(m_full)
  ftest  <- anova(m_base, m_full)
  delta  <- s_full$r.squared - s_base$r.squared
  tibble(
    Model        = label,
    R2_baseline  = round(s_base$r.squared, 3),
    R2_full      = round(s_full$r.squared, 3),
    Delta_R2     = round(delta, 3),
    Delta_R2_pct = round(delta / s_full$r.squared * 100, 1),
    Delta_AdjR2  = round(s_full$adj.r.squared - s_base$adj.r.squared, 3),
    F_stat       = round(ftest$F[2], 2),
    F_pvalue     = round(ftest$`Pr(>F)`[2], 4),
    Verdict      = ifelse(ftest$`Pr(>F)`[2] < 0.05, "Improves fit", "No significant improvement")
  )
}

results_dr2 <- pmap_dfr(
  list(list(wvs_reg_data, wvs_reg_data, wvs_reg_data, wvs_reg_data, gps_reg_data),
       PREDICTORS, MODEL_LABELS),
  compute_delta_r2
)
banner("Marginal contribution of cultural variables")
print(results_dr2, n = Inf, width = Inf)


# ---- 29. Robustness to omitted structural variables ---------------------------
# Adds standard cross-country controls that plausibly come before culture in
# the causal chain:
#   private credit / GDP (financial depth), trade openness (exposure to
#   external shocks), CPI inflation (nominal instability).
# Market capitalisation and turnover are left out on purpose. They plausibly
# mediate the path culture -> market depth -> volatility, which would make
# them "bad controls" (Angrist & Pischke, 2008).

cat("\nDownloading structural controls...\n")
credit_raw    <- fetch_wb_api("FS.AST.PRVT.GD.ZS", "private_credit", YEAR_MAIN)
trade_raw     <- fetch_wb_api("NE.TRD.GNFS.ZS",    "trade_openness", YEAR_MAIN)
inflation_raw <- fetch_wb_api("FP.CPI.TOTL.ZG",    "inflation",      YEAR_MAIN)

add_structural <- function(data) {
  data %>%
    left_join(credit_raw,    by = "iso3c") %>%
    left_join(trade_raw,     by = "iso3c") %>%
    left_join(inflation_raw, by = "iso3c")
}
wvs_reg_str <- add_structural(wvs_reg_data)
gps_reg_str <- add_structural(gps_reg_data)

# Re-estimate with each structural control added on its own, then all three
# together. Reports the HC3 coefficient on culture and its % change from the
# baseline.
compare_structural <- function(data, cultural_var, label) {
  base_terms <- c(cultural_var, "gdp_per_capita", "rule_of_law")
  ok_credit  <- with(data, !is.na(private_credit) & private_credit > 0)
  ok_trade   <- with(data, !is.na(trade_openness) & trade_openness > 0)
  ok_infl    <- !is.na(data$inflation)

  specs <- list(
    "Baseline"               = list(terms = base_terms,                                         rows = TRUE),
    "+ private credit / GDP" = list(terms = c(base_terms, "log(private_credit)"),               rows = ok_credit),
    "+ trade openness"       = list(terms = c(base_terms, "log(trade_openness)"),               rows = ok_trade),
    "+ inflation"            = list(terms = c(base_terms, "inflation"),                         rows = ok_infl),
    "+ all three"            = list(terms = c(base_terms, "log(private_credit)",
                                              "log(trade_openness)", "inflation"),              rows = ok_credit & ok_trade & ok_infl)
  )

  out <- imap_dfr(specs, function(s, spec_name) {
    fit <- lm(reformulate(s$terms, "log_volatility"), data = data[s$rows, ])
    est <- get_hc3(fit, cultural_var)
    tibble(Model = label, Specification = spec_name, N = nobs(fit),
           Beta = est[["beta"]], SE_HC3 = est[["se"]], p_value = est[["p"]])
  }) %>%
    mutate(Pct_change = abs((Beta - Beta[1]) / Beta[1]) * 100,
           Stars      = p_stars(p_value))
  out
}

results_structural <- bind_rows(
  compare_structural(wvs_reg_str, "z_risk_competition", "M3: Pro-Competition"),
  compare_structural(wvs_reg_str, "z_risk_security",    "M4: Low Security Need"),
  compare_structural(gps_reg_str, "gps_risk",           "M5: GPS Risk")
)

banner("Robustness to omitted structural variables (HC3)")
results_structural %>%
  mutate(across(c(Beta, SE_HC3), ~ round(., 3)),
         p_value = round(p_value, 3), Pct_change = round(Pct_change, 1)) %>%
  print(n = Inf, width = Inf)

# Interpretation for M5, generated from the results rather than typed by hand.
m5_all <- results_structural %>% filter(Model == "M5: GPS Risk", Specification == "+ all three")
cat(sprintf(paste0("\nM5 with all three structural controls: beta = %+.3f, p = %.3f, ",
                   "%.1f%% change from baseline (N = %d).\n"),
            m5_all$Beta, m5_all$p_value, m5_all$Pct_change, m5_all$N))


# ---- 30. Out-of-sample replication: WVS Wave 7 (2018) ------------------------
# Re-estimates M1-M4 on WVS Wave 7 (fieldwork 2017-2022, modal year 2018)
# with World Bank 2018 data. M5 cannot be replicated because the GPS exists
# only for 2012.
#
# Item mapping, Wave 6 -> Wave 7 (checked against the WVS-7 codebook):
#   V24 -> Q57   trust, same wording and coding
#   V98 -> Q108  responsibility, same 1-10 scale (10 = people responsible)
#   V99 -> Q109  competition, same 1-10 scale (1 = competition is good)
#   V72 -> Q150  NOT the same item. Wave 7 does not repeat the 6-point
#                Schwartz "secure surroundings" item. Q150 is a forced choice:
#                1 = freedom, 2 = security. It is recoded to 1 = freedom (low
#                security need), 0 = security, so that higher still means more
#                risk-tolerant. The M4 replication therefore uses an
#                approximate proxy and should be read with caution.

banner("WVS Wave 7 replication (2018 cross-section)")

wvs_w7 <- load_rdata(PATHS$wvs_w7)
cat(sprintf("WVS Wave 7 loaded: %d respondents x %d variables\n", nrow(wvs_w7), ncol(wvs_w7)))

wvs_w7_clean <- wvs_w7 %>%
  select(country  = B_COUNTRY,
         weight   = W_WEIGHT,
         q57_raw  = Q57,
         q108_raw = Q108,
         q109_raw = Q109,
         q150_raw = Q150) %>%
  mutate(across(ends_with("_raw"), ~ ifelse(. < 0, NA, .))) %>%
  mutate(
    risk_trust       = ifelse(q57_raw == 1, 1, 0),
    risk_gov_resp    = q108_raw,
    risk_competition = 11 - q109_raw,
    risk_security    = case_when(q150_raw == 1 ~ 1,   # freedom  = low security need
                                 q150_raw == 2 ~ 0,   # security = high security need
                                 TRUE          ~ NA_real_)
  )

national_profiles_w7 <- wvs_w7_clean %>%
  group_by(country) %>%
  summarise(
    avg_risk_trust       = weighted.mean(risk_trust,       w = weight, na.rm = TRUE),
    avg_risk_gov_resp    = weighted.mean(risk_gov_resp,    w = weight, na.rm = TRUE),
    avg_risk_competition = weighted.mean(risk_competition, w = weight, na.rm = TRUE),
    avg_risk_security    = weighted.mean(risk_security,    w = weight, na.rm = TRUE),
    survey_n             = n(),
    .groups = "drop"
  ) %>%
  mutate(
    z_risk_trust       = as.numeric(scale(avg_risk_trust)),
    z_risk_gov_resp    = as.numeric(scale(avg_risk_gov_resp)),
    z_risk_competition = as.numeric(scale(avg_risk_competition)),
    z_risk_security    = as.numeric(scale(avg_risk_security)),
    iso3c = countrycode(country, origin = "iso3n", destination = "iso3c", warn = TRUE)
  )

cat(sprintf("Wave 7 national profiles: %d countries from %d respondents\n",
            nrow(national_profiles_w7), nrow(wvs_w7_clean)))

cat("\nDownloading World Bank data for", YEAR_REPL, "...\n")
vol_repl <- wb_data(indicator = c("market_volatility" = "GFDD.SM.01"),
                    start_date = YEAR_REPL, end_date = YEAR_REPL) %>%
  select(iso3c, market_volatility) %>% filter(!is.na(market_volatility))
gdp_repl <- wb_data(indicator = c("gdp_per_capita" = "NY.GDP.PCAP.PP.KD"),
                    start_date = YEAR_REPL, end_date = YEAR_REPL) %>%
  select(iso3c, gdp_per_capita) %>% filter(!is.na(gdp_per_capita))
rol_repl <- fetch_wb_api("GOV_WGI_RL.EST", "rule_of_law", YEAR_REPL)

wvs_reg_w7 <- national_profiles_w7 %>%
  left_join(vol_repl, by = "iso3c") %>%
  left_join(gdp_repl, by = "iso3c") %>%
  left_join(rol_repl, by = "iso3c") %>%
  filter(!is.na(market_volatility), market_volatility > 0,
         !is.na(gdp_per_capita), !is.na(rule_of_law)) %>%
  mutate(log_volatility = log(market_volatility))

cat(sprintf("Wave 7 regression sample: N = %d\n", nrow(wvs_reg_w7)))

# Compare each Wave 7 estimate with its Wave 6 counterpart, which is taken
# directly from the fitted Wave 6 models.
replicate_w7 <- function(w6_model, predictor, label) {
  fit <- lm(reformulate(c(predictor, "gdp_per_capita", "rule_of_law"), "log_volatility"),
            data = wvs_reg_w7)
  w7 <- get_hc3(fit, predictor)
  w6 <- get_hc3(w6_model, predictor)
  tibble(Model = label, N_W7 = nobs(fit),
         W7_Beta = round(w7[["beta"]], 3), W7_SE = round(w7[["se"]], 3),
         W7_p = round(w7[["p"]], 4), W7_Stars = p_stars(w7[["p"]]),
         W6_Beta = round(w6[["beta"]], 3), W6_p = round(w6[["p"]], 4),
         Sign_match = ifelse(sign(w7[["beta"]]) == sign(w6[["beta"]]), "Same", "Opposite"))
}

results_w7 <- pmap_dfr(
  list(models_log[1:4], PREDICTORS[1:4],
       c(MODEL_LABELS[1:3], "M4: Security (W7: freedom vs security)")),
  replicate_w7
)
print(results_w7, n = Inf, width = Inf)


# ---- 31. Balanced-sample comparison: same countries in both waves -------------
# Restricting both waves to the countries they share separates genuine change
# over time from changes in sample composition.

common_iso <- intersect(wvs_reg_data$iso3c, wvs_reg_w7$iso3c)
banner(sprintf("Balanced sample: Wave 6 vs Wave 7 (N = %d common countries)", length(common_iso)))

waves <- list(`Wave 6 (2012)` = wvs_reg_data, `Wave 7 (2018)` = wvs_reg_w7)

balanced <- expand_grid(Predictor = c("z_risk_competition", "z_risk_security"),
                        Wave      = names(waves)) %>%
  pmap_dfr(function(Predictor, Wave) {
    fit <- lm(reformulate(c(Predictor, "gdp_per_capita", "rule_of_law"), "log_volatility"),
              data = filter(waves[[Wave]], iso3c %in% common_iso))
    est <- get_hc3(fit, Predictor)
    tibble(Predictor = Predictor, Wave = Wave, N = nobs(fit),
           Beta = round(est[["beta"]], 3), SE_HC3 = round(est[["se"]], 3),
           p_value = round(est[["p"]], 3))
  })
print(balanced, n = Inf, width = Inf)


# ---- 32. Session information --------------------------------------------------
# Records the R and package versions used, for reproducibility.

writeLines(capture.output(sessionInfo()), file.path("output", "sessionInfo.txt"))

banner("Pipeline complete")
cat("Tables:  ", PATHS$tables, "\n")
cat("Figures: ", PATHS$figures, "\n")
