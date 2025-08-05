# Trends in and Factors Associated with Long COVID and Recovery among US Adults
# Data source is the National Health Interview Survey (NHIS), 2022-2024
# Rishi Shah
# July 20, 2025

# load necessary libraries
library(data.table)
library(foreign)
library(dplyr)
library(survey)
library(srvyr)
library(data.table)
library(epiDisplay)
library(ggplot2)
library(gtsummary)
library(gt)
library(tidyverse)
library(stringr)

# read in 2022 NHIS survey data
df_2022 <- fread('adult22.csv') %>%
  rename_all(tolower) %>%
  mutate(
    year = 2022,
    ever_covid = case_when(
      cvddiag_a == 1 | postest_a == 1 ~ 1,
      cvddiag_a == 2 & postest_a == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    ever_longcovid = case_when(
      ever_covid == 1 & longcvd_a == 1 ~ 1,
      ever_covid == 1 & longcvd_a == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    recovered_longcovid = case_when(
      longcvd_a == 1 & sympnow_a == 1 ~ 0,   # still has symptoms
      longcvd_a == 1 & sympnow_a == 2 ~ 1,   # recovered
      TRUE ~ NA_real_
    )
  )

# read in 2023 NHIS survey data
df_2023 <- fread('adult23.csv') %>%
  rename_all(tolower) %>%
  mutate(
    year = 2023,
    ever_covid = case_when(
      evercovd_a == 1 ~ 1,
      evercovd_a == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    ever_longcovid = case_when(
      ever_covid == 1 & longcovd1_a == 1 ~ 1,
      ever_covid == 1 & longcovd1_a == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    recovered_longcovid = case_when(
      longcovd1_a == 1 & sympnow1_a == 1 ~ 0,   # still has symptoms
      longcovd1_a == 1 & sympnow1_a == 2 ~ 1,   # recovered
      TRUE ~ NA_real_
    )
  )

# read in 2024 NHIS survey data
df_2024 <- fread('adult24.csv') %>%
  rename_all(tolower) %>%
  mutate(
    year = 2024,
    ever_covid = case_when(
      evercovd_a == 1 ~ 1,
      evercovd_a == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    ever_longcovid = case_when(
      ever_covid == 1 & longcovd2_a == 1 ~ 1,
      ever_covid == 1 & longcovd2_a == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    recovered_longcovid = case_when(
      longcovd2_a == 1 & sympnow1_a == 1 ~ 0,   # still has symptoms
      longcovd2_a == 1 & sympnow1_a == 2 ~ 1,   # recovered
      TRUE ~ NA_real_
    )
  )

# bind data
df <- bind_rows(df_2022, df_2023, df_2024)

# harmonize urbanicity
df$urban <- ifelse(df$year == 2024, df$urbrrl23, df$urbrrl)

# plot trends in rates of COVID-19, long COVID, and recovery
covid_color <- "#374e55"
longcovid_color <- "#df8f44" 
recovery_color <- "#b24745"

calculate_covid_estimates <- function(data, year_val) {
  
  year_data <- data %>% 
    filter(year == year_val, !is.na(wtfa_a), wtfa_a > 0)
  
  svy_design <- svydesign(
    ids = ~ppsu,
    strata = ~pstrat, 
    weights = ~wtfa_a,
    data = year_data,
    nest = TRUE
  )
  
  covid_est <- svymean(~ever_covid, svy_design, na.rm = TRUE)
  longcovid_est <- svymean(~ever_longcovid, svy_design, na.rm = TRUE)
  
  year_data$recovered_from_longcovid_pop <- ifelse(
    year_data$ever_longcovid == 1 & year_data$recovered_longcovid == 1, 1, 0
  )
  
  svy_design <- svydesign(
    ids = ~ppsu,
    strata = ~pstrat, 
    weights = ~wtfa_a,
    data = year_data,
    nest = TRUE
  )
  
  recovery_pop_est <- svymean(~recovered_from_longcovid_pop, svy_design, na.rm = TRUE)
  
  covid_total <- svytotal(~ever_covid, svy_design, na.rm = TRUE)
  longcovid_total <- svytotal(~ever_longcovid, svy_design, na.rm = TRUE)
  recovery_total <- svytotal(~recovered_from_longcovid_pop, svy_design, na.rm = TRUE)
  
  results <- data.frame(
    year = year_val,
    covid_pct = as.numeric(covid_est[1]) * 100,
    covid_se = as.numeric(SE(covid_est)[1]) * 100,
    longcovid_pct = as.numeric(longcovid_est[1]) * 100,
    longcovid_se = as.numeric(SE(longcovid_est)[1]) * 100,
    recovery_pct = as.numeric(recovery_pop_est[1]) * 100,
    recovery_se = as.numeric(SE(recovery_pop_est)[1]) * 100,
    covid_total = as.numeric(covid_total[1]) / 1000000,
    covid_total_se = as.numeric(SE(covid_total)[1]) / 1000000,
    longcovid_total = as.numeric(longcovid_total[1]) / 1000000,
    longcovid_total_se = as.numeric(SE(longcovid_total)[1]) / 1000000,
    recovery_total = as.numeric(recovery_total[1]) / 1000000,
    recovery_total_se = as.numeric(SE(recovery_total)[1]) / 1000000
  )
  
  results$covid_ci_low <- results$covid_pct - 1.96 * results$covid_se
  results$covid_ci_upp <- results$covid_pct + 1.96 * results$covid_se
  results$longcovid_ci_low <- results$longcovid_pct - 1.96 * results$longcovid_se
  results$longcovid_ci_upp <- results$longcovid_pct + 1.96 * results$longcovid_se
  results$recovery_ci_low <- results$recovery_pct - 1.96 * results$recovery_se
  results$recovery_ci_upp <- results$recovery_pct + 1.96 * results$recovery_se
  
  results$covid_total_ci_low <- results$covid_total - 1.96 * results$covid_total_se
  results$covid_total_ci_upp <- results$covid_total + 1.96 * results$covid_total_se
  results$longcovid_total_ci_low <- results$longcovid_total - 1.96 * results$longcovid_total_se
  results$longcovid_total_ci_upp <- results$longcovid_total + 1.96 * results$longcovid_total_se
  results$recovery_total_ci_low <- results$recovery_total - 1.96 * results$recovery_total_se
  results$recovery_total_ci_upp <- results$recovery_total + 1.96 * results$recovery_total_se
  
  return(results)
}

estimates_2022 <- calculate_covid_estimates(df, 2022)
estimates_2023 <- calculate_covid_estimates(df, 2023)
estimates_2024 <- calculate_covid_estimates(df, 2024)

covid_trends <- bind_rows(estimates_2022, estimates_2023, estimates_2024)

totals_table <- covid_trends %>%
  select(year, covid_total, longcovid_total, recovery_total) %>%
  mutate(across(c(covid_total, longcovid_total, recovery_total), ~round(.x, 1)))

plot_data_totals <- covid_trends %>%
  select(year, 
         covid_total, covid_total_ci_low, covid_total_ci_upp,
         longcovid_total, longcovid_total_ci_low, longcovid_total_ci_upp,
         recovery_total, recovery_total_ci_low, recovery_total_ci_upp) %>%
  tidyr::pivot_longer(
    cols = -year,
    names_to = c("measure", "stat"),
    names_pattern = "(.+)_total(_ci_low|_ci_upp|$)",
    values_to = "value"
  ) %>%
  mutate(stat = ifelse(stat == "", "total", gsub("_", "", stat))) %>%
  tidyr::pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    measure_label = case_when(
      measure == "covid" ~ "Ever had COVID-19",
      measure == "longcovid" ~ "Ever had Long COVID",
      measure == "recovery" ~ "Recovered from Long COVID"
    ),
    color = case_when(
      measure == "covid" ~ covid_color,
      measure == "longcovid" ~ longcovid_color,
      measure == "recovery" ~ recovery_color
    )
  )

plot_covid_totals <- ggplot(plot_data_totals, aes(x = year, y = total, color = measure_label, fill = measure_label)) +
  geom_ribbon(aes(ymin = cilow, ymax = ciupp), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("%.1fM", total), y = total + max(plot_data_totals$total) * 0.03),
            color = "black", size = 4, vjust = 0.2) +
  scale_color_manual(values = c("Ever had COVID-19" = covid_color,
                                "Ever had Long COVID" = longcovid_color,
                                "Recovered from Long COVID" = recovery_color)) +
  scale_fill_manual(values = c("Ever had COVID-19" = covid_color,
                               "Ever had Long COVID" = longcovid_color, 
                               "Recovered from Long COVID" = recovery_color)) +
  labs(
    title = "Number of US adults with COVID-19, long COVID,\nand recovery from long COVID, 2022–2024",
    x = "Year",
    y = "Number of US adults (millions)",
    color = "",
    fill = ""
  ) +
  scale_x_continuous(breaks = 2022:2024) +
  scale_y_continuous(limits = c(0, max(plot_data_totals$ciupp) * 1.1), 
                     breaks = seq(0, ceiling(max(plot_data_totals$ciupp)/10)*10, by = 20)) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 14),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.text = element_text(size = 11),
    legend.margin = margin(6, 6, 6, 6),
    legend.key.size = unit(0.8, "cm")
  )


################################
#### 1. DATA PRE-PROCESSING ####
################################

# set the seed for reproducibility
set.seed(07202025)

# covariates 

# race/ethnicity
df$race_simple <- NA
df$race_simple[df$hispallp_a == 1] <- 1
df$race_simple[df$hispallp_a == 2] <- 2
df$race_simple[df$hispallp_a == 3] <- 3
df$race_simple[df$hispallp_a == 4] <- 4
df$race_simple[df$hispallp_a >= 5 & df$hispallp_a <= 7] <- 5

levels <- c("Hispanic", "NH White", "NH Black", "NH Asian", "NH Other")
labels <- c(1, 2, 3, 4, 5)
df$race_simple <- factor(df$race_simple, levels = labels, labels = levels)

# binarize sex
# 1 = female
df$female <- 0
df$female[df$sex_a == 2] <- 1

levels <- c("Male", "Female")
labels <- c(0, 1)
df$female <- factor(df$female, levels = labels, labels = levels)

# binarize less than or greater than high school education
# 1 = some college or higher
df$educ_binary <- NA
df$educ_binary[df$educp_a >= 1 & df$educp_a < 3] <- 0
df$educ_binary[df$educp_a >= 3 & df$educp_a < 99] <- 1

levels <- c("Less HS", "More HS")
labels <- c(0, 1)
df$educ_binary <- factor(df$educ_binary, levels = labels, labels = levels)

# binarize low and middle/high family income
# low = < 200% of the federal poverty limit = 0
# middle/high = >= 200% of the federal poverty limit = 1
df$inc_binary <- NA
df$inc_binary[df$ratcat_a >= 1 & df$ratcat_a <= 7] <- 0
df$inc_binary[df$ratcat_a >= 8 & df$ratcat_a <= 14] <- 1

levels <- c("Low income", "High income")
labels <- c(0, 1)
df$inc_binary <- factor(df$inc_binary, levels = labels, labels = levels)

# make three age categories
df$agecat <- NA
df$agecat[df$agep_a < 40] <- 1
df$agecat[df$agep_a >= 40 & df$agep_a < 65] <- 2
df$agecat[df$agep_a >= 65 & df$agep_a < 99] <- 3

levels <- c("< 40 years", "40-64 years", "≥ 65 years")
labels <- c(1, 2, 3)
df$agecat <- factor(df$agecat, levels = labels, labels = levels)

# make binary variables for each region of the US
df$neast <- ifelse(df$region == 1, 1, 0)
df$midwest <- ifelse(df$region == 2, 1, 0)
df$south <- ifelse(df$region == 3, 1, 0)
df$west <- ifelse(df$region == 4, 1, 0)

levels <- c("Northeast", "Midwest", "South", "West")
labels <- c(1, 2, 3, 4)
df$region <- factor(df$region, levels = labels, labels = levels)

# make binary variables for each urban/rural classification
df$lrge_cent_metr <- ifelse(df$urban == 1, 1, 0)
df$lrge_frin_metr <- ifelse(df$urban == 2, 1, 0)
df$med_sml_metr <- ifelse(df$urban == 3, 1, 0)
df$non_metr <- ifelse(df$urban == 4, 1, 0)

df$urbrur <- ifelse(df$urban %in% c(1, 2, 3), 1, ifelse(df$urban == 4, 0, NA))
levels <- c("Nonmetropolitan", "Metropolitan")
labels <- c(0, 1)
df$urbrur <- factor(df$urbrur, levels = labels, labels = levels)

# binarize employment status
df$empstat <- ifelse(df$empwrklsw1_a == 1, 1, ifelse(df$empwrklsw1_a == 2, 0, NA))

levels <- c("Unemployed", "Employed")
labels <- c(0, 1)
df$empstat <- factor(df$empstat, levels = labels, labels = levels)

# covid severity (only available for 2022)
df$cvdsev_a <- ifelse(df$cvdsev_a %in% c(7, 9), NA, df$cvdsev_a)
levels <- c("Mild symptoms", "Moderate symptoms", "Severe symptoms")
labels <- c(2, 3, 4)
df$cvdsev_a <- factor(df$cvdsev_a, levels = labels, labels = levels)

# LC-associated activity limitations (only available for 2023 and 2024)
df$lcvdact_a <- ifelse(df$lcvdact_a %in% c(7, 8, 9), NA, df$lcvdact_a)
levels <- c("Not at all", "A little", "A lot")
labels <- c(1, 2, 3)
df$lcvdact_a <- factor(df$lcvdact_a, levels = labels, labels = levels)

# insurance coverage
# make three insurance status categories (public, private, uninsured)
# 1 = public, 2 = private, 3 = uninsured
df$insurancecat <- NA
df$insurancecat[df$medicare_a == 1 | df$medicare_a == 2 | df$medicaid_a == 1 | df$medicaid_a == 2 | df$chip_a == 1 | df$chip_a == 2 | df$othpub_a == 1 | df$othpub_a == 2 | df$othgov_a == 1 | df$othgov_a == 2 | df$military_a == 1 | df$military_a == 2] <- 1
df$insurancecat[df$private_a == 1 | df$private_a == 2] <- 2
df$insurancecat[df$notcov_a == 1] <- 3

levels <- c("Public", "Private", "Uninsured")
labels <- c(1, 2, 3)
df$insurancecat <- factor(df$insurancecat, levels = labels, labels = levels)

################################
####    2. DATA ANALYSIS    ####
################################

options(survey.lonely.psu = 'adjust')
data_2022 <- as_survey(subset(df, year == 2022), id = ppsu, weight = wtfa_a, strata = pstrat, nest = TRUE)
data_2023 <- as_survey(subset(df, year == 2023), id = ppsu, weight = wtfa_a, strata = pstrat, nest = TRUE)
data_2024 <- as_survey(subset(df, year == 2024), id = ppsu, weight = wtfa_a, strata = pstrat, nest = TRUE)

# only include those who ever had COVID-19 and a yes/no answer to the LC question
subset_data_2022 <- data_2022 %>%
  filter(ever_covid == 1 & !is.na(ever_longcovid))

subset_data_2023 <- data_2023 %>%
  filter(ever_covid == 1 & !is.na(ever_longcovid))

subset_data_2024 <- data_2024 %>%
  filter(ever_covid == 1 & !is.na(ever_longcovid))

# what are the characteristics/covariates associated with LC and recovery each year from 2022-2024?

# make summary table with national estimates for LC
svy_table_lc_2022 <- tbl_svysummary(subset_data_2022, by = ever_longcovid, 
               include = c(agecat, female, race_simple, empstat, inc_binary,
                           insurancecat, educ_binary, urbrur, region, cvdsev_a)) %>% add_p()

svy_table_lc_2023 <- tbl_svysummary(subset_data_2023, by = ever_longcovid, 
               include = c(agecat, female, race_simple, empstat, inc_binary,
                           insurancecat, educ_binary, urbrur, region)) %>% add_p()

svy_table_lc_2024 <- tbl_svysummary(subset_data_2024, by = ever_longcovid, 
               include = c(agecat, female, race_simple, empstat, inc_binary,
                           insurancecat, educ_binary, urbrur, region)) %>% add_p()

# normal raw counts table without national estimates for LC
raw_counts_lc_2022 <- df %>%
  filter(year == 2022, ever_covid == 1 & !is.na(ever_longcovid)) %>%
  select(agecat, female, race_simple, empstat, inc_binary, insurancecat,
         educ_binary, urbrur, region, cvdsev_a, ever_longcovid) %>%
  tbl_summary(by = ever_longcovid, missing = "no") %>%
  add_p()

raw_counts_lc_2023 <- df %>%
  filter(year == 2023, ever_covid == 1 & !is.na(ever_longcovid)) %>%
  select(agecat, female, race_simple, empstat, inc_binary, insurancecat,
         educ_binary, urbrur, region, ever_longcovid) %>%
  tbl_summary(by = ever_longcovid, missing = "no") %>%
  add_p()

raw_counts_lc_2024 <- df %>%
  filter(year == 2024, ever_covid == 1 & !is.na(ever_longcovid)) %>%
  select(agecat, female, race_simple, empstat, inc_binary, insurancecat,
         educ_binary, urbrur, region, ever_longcovid) %>%
  tbl_summary(by = ever_longcovid, missing = "no") %>%
  add_p()

# only include those who ever had LC and a yes/no answer to the symptoms now question
subset_data_2022_recovery <- subset_data_2022 %>%
  filter(ever_longcovid == 1 & !is.na(recovered_longcovid))

subset_data_2023_recovery <- subset_data_2023 %>%
  filter(ever_longcovid == 1 & !is.na(recovered_longcovid))

subset_data_2024_recovery <- subset_data_2024 %>%
  filter(ever_longcovid == 1 & !is.na(recovered_longcovid))

# make summary table with national estimates for LC recovery
svy_table_lc_recovery_2022 <- tbl_svysummary(subset_data_2022_recovery, by = recovered_longcovid, 
               include = c(agecat, female, race_simple, empstat, inc_binary,
                           insurancecat, educ_binary, urbrur, region, cvdsev_a)) %>% add_p()

svy_table_lc_recovery_2023 <- tbl_svysummary(subset_data_2023_recovery, by = recovered_longcovid, 
               include = c(agecat, female, race_simple, empstat, inc_binary,
                           insurancecat, educ_binary, urbrur, region, lcvdact_a)) %>% add_p()

svy_table_lc_recovery_2024 <- tbl_svysummary(subset_data_2024_recovery, by = recovered_longcovid, 
               include = c(agecat, female, race_simple, empstat, inc_binary,
                           insurancecat, educ_binary, urbrur, region, lcvdact_a)) %>% add_p()

# normal raw counts table without national estimates for LC recovery
raw_counts_lc_recovery_2022 <- df %>%
  filter(year == 2022, ever_longcovid == 1 & !is.na(recovered_longcovid)) %>%
  select(agecat, female, race_simple, empstat, inc_binary, insurancecat,
         educ_binary, urbrur, region, cvdsev_a, recovered_longcovid) %>%
  tbl_summary(by = recovered_longcovid, missing = "no") %>%
  add_p()

raw_counts_lc_recovery_2023 <- df %>%
  filter(year == 2023, ever_longcovid == 1 & !is.na(recovered_longcovid)) %>%
  select(agecat, female, race_simple, empstat, inc_binary, insurancecat,
         educ_binary, urbrur, region, lcvdact_a, recovered_longcovid) %>%
  tbl_summary(by = recovered_longcovid, missing = "no") %>%
  add_p()

raw_counts_lc_recovery_2024 <- df %>%
  filter(year == 2024, ever_longcovid == 1 & !is.na(recovered_longcovid)) %>%
  select(agecat, female, race_simple, empstat, inc_binary, insurancecat,
         educ_binary, urbrur, region, lcvdact_a, recovered_longcovid) %>%
  tbl_summary(by = recovered_longcovid, missing = "no") %>%
  add_p()


# multivariate logistic regression model to examine which factors are associated with LC and recovery

# fit multivariate logistic regression model for factors associated with LC by year
fit_lc_2022 <- svyglm(ever_longcovid ~ 
                       relevel(agecat, ref = "< 40 years") +
                       relevel(female, ref = "Male") +
                       relevel(race_simple, ref = "NH White") + 
                       relevel(urbrur, ref = "Metropolitan") +
                       relevel(educ_binary, ref = "More HS"),
                       # relevel(cvdsev_a, ref = "Mild symptoms"),
                       design = subset_data_2022,
                     family = quasibinomial(link = 'logit'))

fit_lc_2023 <- svyglm(ever_longcovid ~ 
                   relevel(agecat, ref = "< 40 years") +
                   relevel(female, ref = "Male") +
                   relevel(race_simple, ref = "NH White") + 
                   relevel(urbrur, ref = "Metropolitan") +
                   relevel(educ_binary, ref = "More HS"),
                   design = subset_data_2023,
                 family = quasibinomial(link = 'logit'))

fit_lc_2024 <- svyglm(ever_longcovid ~ 
                        relevel(agecat, ref = "< 40 years") +
                        relevel(female, ref = "Male") +
                        relevel(race_simple, ref = "NH White") + 
                        relevel(urbrur, ref = "Metropolitan") +
                        relevel(educ_binary, ref = "More HS"),
                        design = subset_data_2024,
                      family = quasibinomial(link = 'logit'))

# fit multivariate logistic regression model for factors associated with LC recovery by year
fit_lc_recovery_2022 <- svyglm(recovered_longcovid ~ 
                                relevel(agecat, ref = "< 40 years") +
                                relevel(female, ref = "Male") +
                                relevel(race_simple, ref = "NH White") + 
                                relevel(urbrur, ref = "Metropolitan") +
                                relevel(educ_binary, ref = "More HS"),
                                # relevel(cvdsev_a, ref = "Mild symptoms"),
                                design = subset_data_2022_recovery,
                              family = quasibinomial(link = 'logit'))

fit_lc_recovery_2023 <- svyglm(recovered_longcovid ~ 
                                 relevel(agecat, ref = "< 40 years") +
                                 relevel(female, ref = "Male") +
                                 relevel(race_simple, ref = "NH White") + 
                                 relevel(urbrur, ref = "Metropolitan") +
                                 relevel(educ_binary, ref = "More HS"),
                                 design = subset_data_2023_recovery,
                               family = quasibinomial(link = 'logit'))


fit_lc_recovery_2024 <- svyglm(recovered_longcovid ~ 
                                relevel(agecat, ref = "< 40 years") +
                                relevel(female, ref = "Male") +
                                relevel(race_simple, ref = "NH White") + 
                                relevel(urbrur, ref = "Metropolitan") +
                                relevel(educ_binary, ref = "More HS"),
                                design = subset_data_2024_recovery,
                              family = quasibinomial(link = 'logit'))

# get odds ratios and confidence intervals
extract_aors <- function(model) {
  broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    select(term, estimate, conf.low, conf.high) %>%
    mutate(aor = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high)) %>%
    select(term, aor)
}

aor_lc_2022 <- extract_aors(fit_lc_2022)
aor_recovery_2022 <- extract_aors(fit_lc_recovery_2022)

aor_lc_2023 <- extract_aors(fit_lc_2023)
aor_recovery_2023 <- extract_aors(fit_lc_recovery_2023)

aor_lc_2024 <- extract_aors(fit_lc_2024)
aor_recovery_2024 <- extract_aors(fit_lc_recovery_2024)

# export all data
df_counts_2022 <- as_tibble(raw_counts_lc_2022$table_body)
df_weighted_2022 <- as_tibble(svy_table_lc_2022$table_body)
df_aor_2022 <- aor_lc_2022

df_counts_2023 <- as_tibble(raw_counts_lc_2023$table_body)
df_weighted_2023 <- as_tibble(svy_table_lc_2023$table_body)
df_aor_2023 <- aor_lc_2023

df_counts_2024 <- as_tibble(raw_counts_lc_2024$table_body)
df_weighted_2024 <- as_tibble(svy_table_lc_2024$table_body)
df_aor_2024 <- aor_lc_2024


df_counts_recovery_2022 <- as_tibble(raw_counts_lc_recovery_2022$table_body)
df_weighted_recovery_2022 <- as_tibble(svy_table_lc_recovery_2022$table_body)
df_aor_recovery_2022 <- aor_recovery_2022

df_counts_recovery_2023 <- as_tibble(raw_counts_lc_recovery_2023$table_body)
df_weighted_recovery_2023 <- as_tibble(svy_table_lc_recovery_2023$table_body)
df_aor_recovery_2023 <- aor_recovery_2023

df_counts_recovery_2024 <- as_tibble(raw_counts_lc_recovery_2024$table_body)
df_weighted_recovery_2024 <- as_tibble(svy_table_lc_recovery_2024$table_body)
df_aor_recovery_2024 <- aor_recovery_2024


clean_aor_terms <- function(aor_df) {
  aor_df %>%
    filter(term != "(Intercept)") %>%
    mutate(
      variable = case_when(
        str_detect(term, "agecat") ~ "agecat",
        str_detect(term, "female") ~ "female",
        str_detect(term, "race_simple") ~ "race_simple",
        str_detect(term, "urbrur") ~ "urbrur",
        str_detect(term, "educ_binary") ~ "educ_binary",
        str_detect(term, "cvdsev_a") ~ "cvdsev_a",
        TRUE ~ NA_character_
      ),
      label = term %>%
        str_replace_all('relevel\\([^\\)]+\\)', '') %>%
        str_replace_all('\"', '') %>%
        str_trim()
    ) %>%
    select(variable, label, aor)
}

aor_clean_2022 <- clean_aor_terms(df_aor_2022)
aor_clean_2023 <- clean_aor_terms(df_aor_2023)
aor_clean_2024 <- clean_aor_terms(df_aor_2024)

merge_lc_tables <- function(counts, weighted, aor_clean) {
  counts_slim <- counts %>%
    filter(row_type == "level") %>%
    select(variable, label, raw_stat_1 = stat_1, raw_stat_2 = stat_2)
  
  weighted_slim <- weighted %>%
    filter(row_type == "level") %>%
    select(variable, label, weighted_stat_1 = stat_1, weighted_stat_2 = stat_2)
  
  p_values <- weighted %>%
    filter(row_type == "label") %>%
    select(variable, p.value)
  
  merged <- counts_slim %>%
    left_join(weighted_slim, by = c("variable", "label")) %>%
    left_join(p_values, by = "variable") %>%
    left_join(aor_clean, by = c("variable", "label"))
  
  return(merged)
}

final_lc_2022 <- merge_lc_tables(df_counts_2022, df_weighted_2022, aor_clean_2022)
final_lc_2023 <- merge_lc_tables(df_counts_2023, df_weighted_2023, aor_clean_2023)
final_lc_2024 <- merge_lc_tables(df_counts_2024, df_weighted_2024, aor_clean_2024)

clean_aor_terms_recovery <- function(aor_df) {
  aor_df %>%
    filter(term != "(Intercept)") %>%
    mutate(
      variable = case_when(
        str_detect(term, "agecat") ~ "agecat",
        str_detect(term, "female") ~ "female",
        str_detect(term, "race_simple") ~ "race_simple",
        str_detect(term, "urbrur") ~ "urbrur",
        str_detect(term, "educ_binary") ~ "educ_binary",
        str_detect(term, "cvdsev_a") ~ "cvdsev_a",
        str_detect(term, "lcvdact_a") ~ "lcvdact_a",
        TRUE ~ NA_character_
      ),
      label = term %>%
        str_replace_all('relevel\\([^\\)]+\\)', '') %>%
        str_replace_all('\"', '') %>%
        str_trim()
    ) %>%
    select(variable, label, aor)
}

aor_clean_recovery_2022 <- clean_aor_terms_recovery(df_aor_recovery_2022)
aor_clean_recovery_2023 <- clean_aor_terms_recovery(df_aor_recovery_2023)
aor_clean_recovery_2024 <- clean_aor_terms_recovery(df_aor_recovery_2024)

merge_recovery_tables <- function(counts, weighted, aor_clean) {
  counts_slim <- counts %>%
    filter(row_type == "level") %>%
    select(variable, label, raw_stat_1 = stat_1, raw_stat_2 = stat_2)
  
  weighted_slim <- weighted %>%
    filter(row_type == "level") %>%
    select(variable, label, weighted_stat_1 = stat_1, weighted_stat_2 = stat_2)
  
  p_values <- weighted %>%
    filter(row_type == "label") %>%
    select(variable, p.value)
  
  merged <- counts_slim %>%
    left_join(weighted_slim, by = c("variable", "label")) %>%
    left_join(p_values, by = "variable") %>%
    left_join(aor_clean, by = c("variable", "label"))
  
  return(merged)
}

final_recovery_2022 <- merge_recovery_tables(df_counts_recovery_2022, df_weighted_recovery_2022, aor_clean_recovery_2022)
final_recovery_2023 <- merge_recovery_tables(df_counts_recovery_2023, df_weighted_recovery_2023, aor_clean_recovery_2023)
final_recovery_2024 <- merge_recovery_tables(df_counts_recovery_2024, df_weighted_recovery_2024, aor_clean_recovery_2024)


write_csv(final_lc_2022, "longcovid_table_2022.csv")
write_csv(final_lc_2023, "longcovid_table_2023.csv")
write_csv(final_lc_2024, "longcovid_table_2024.csv")

write_csv(final_recovery_2022, "recovery_table_2022.csv")
write_csv(final_recovery_2023, "recovery_table_2023.csv")
write_csv(final_recovery_2024, "recovery_table_2024.csv")