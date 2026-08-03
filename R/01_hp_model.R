# ---------------------------------------------------------------------------
# 01 - Subregional H. pylori seroprevalence model
#
# Fits a negative-binomial model to age-specific seroprevalence estimates
# pooled from published studies, then predicts a complete subregion x period x
# age-group grid under a reference study profile.
#
# Run from the project root:  source("R/01_hp_model.R")
# ---------------------------------------------------------------------------

library(data.table)
library(readxl)
library(dplyr)
library(INLA)
library(ggplot2)
library(tidyr)
library(ggpubr)

source("R/utils_model.R")
source("R/utils_plot.R")

dir.create("output", showWarnings = FALSE)


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

# pooled study estimates
hpPrev <- read_excel("data/age_spec_serology.xlsx", 1) |>
  select(
    Study_design_3, Estimated_period, Setting, age_spec_group,
    Detection_method_2, Area1, HDI_2019, ROB, Total, HP, HP_prevalence
  ) |>
  mutate(
    # collapse the one study period that straddles two bands
    Estimated_period = case_when(
      Estimated_period == "2015~2021" ~ "2015~2019",
      TRUE ~ Estimated_period
    )
  ) |>
  rename(population = Total)

# country -> subregion lookup. The right-hand labels match the country strings
# as they appear in the extraction sheet, including "Omen", which is a
# misspelling carried in the source data and must be preserved for the join.
subregionMapping <- fread(
  "data/subregionMapping.csv"
) |>
  select(country, Subregion) |>
  mutate(
    country = case_when(
      country == "United States of America" ~ "USA",
      country == "Myanmar"                  ~ "Burma OR Myanmar",
      country == "Viet Nam"                 ~ "Vietnam",
      country == "Oman"                     ~ "Omen",
      country == "Brunei Darussalam"        ~ "Brunei",
      TRUE ~ country
    )
  ) |>
  distinct()

yearMapping <- data.frame(
  Estimated_period = c(
    "1990_before", "1990~1994", "1995~1999", "2000~2004", "2005~2009",
    "2010~2014", "2015~2019", "2020~2024", "2025~2029", "2030~2035",
    "2035~2039"
  ),
  year_idx = 0:10
)

ageMapping <- data.frame(
  age_spec_group = c("a_20~50", "b_50~70", "c_70"),
  agp_idx        = 0:2
)


# ---------------------------------------------------------------------------
# HDI: observed values, then linear extrapolation on the logit scale
# ---------------------------------------------------------------------------

# working on the logit scale keeps extrapolated HDI inside (0, 1)
hdi <- fread("data/HDI_subregion.csv") |>
  mutate(hdi_logit = qlogis(hdi_mean)) |>
  rename(year_idx = Year_dummy)

model_hdi <- lm(hdi_logit ~ Subregion * year_idx, data = hdi)

newdata_hdi <- expand.grid(
  Subregion = levels(as.factor(hdi$Subregion)),
  year_idx  = 5:10
)

hdi_pred <- cbind(
  newdata_hdi,
  pred_hdi = plogis(predict(model_hdi, newdata_hdi))
) |>
  mutate(
    HDI_rank_dummy = case_when(
      pred_hdi < 0.550 ~ 0,
      pred_hdi < 0.700 ~ 1,
      pred_hdi < 0.800 ~ 2,
      TRUE             ~ 3
    )
  )

hdiMapping <- bind_rows(hdi, hdi_pred) |>
  select(year_idx, Subregion, HDI_rank_dummy) |>
  unique()


# ---------------------------------------------------------------------------
# Population denominators, collapsed to the three modelled age groups
# ---------------------------------------------------------------------------

popYearMapping <- data.frame(
  year     = paste0(seq(1985, 2050, 5), "-", seq(1989, 2054, 5)),
  year_idx = 0:13
)

pop <- readRDS("data/pop_subregion.RDS") |>
  ungroup() |>
  na.omit() |>
  select(Subregion, age, sex, year, population) |>
  mutate(
    agp_idx = case_when(
      age %in% c("20-24", "25-29", "30-34", "35-39", "40-44", "45-49") ~ 0,
      age %in% c("50-54", "55-59", "60-64", "65-69")                   ~ 1,
      age %in% c("70-74", "75-79", "80-84", "85+")                     ~ 2,
      TRUE ~ NA_real_
    )
  ) |>
  left_join(popYearMapping, by = "year") |>
  summarise(
    population = sum(population, na.rm = TRUE),
    .by = c(Subregion, agp_idx, year_idx)
  ) |>
  na.omit()

rm(hdi, hdi_pred, model_hdi, newdata_hdi, popYearMapping)


# ---------------------------------------------------------------------------
# Assemble the modelling frame
#
# Observed study rows and the prediction grid are stacked into one frame so
# INLA fits and predicts in a single pass. `row_type` records which is which;
# every downstream step selects on that flag rather than on row position.
# ---------------------------------------------------------------------------

input <- hpPrev |>
  left_join(subregionMapping, by = c("Area1" = "country")) |>
  left_join(yearMapping, by = "Estimated_period") |>
  left_join(ageMapping, by = "age_spec_group") |>
  left_join(hdiMapping, by = c("year_idx", "Subregion")) |>
  select(
    Study_design_3, Setting, Detection_method_2, Subregion, year_idx,
    agp_idx, HDI_rank_dummy, HP, population, HP_prevalence
  ) |>
  mutate(row_type = "observed")

# reference profile for prediction: population-based cross-sectional
# serology studies
pred_grid <- expand.grid(
  Subregion          = unique(na.omit(input$Subregion)),
  year_idx           = as.integer(0:10),
  agp_idx            = as.integer(0:2),
  Study_design_3     = "cs_cn",
  Setting            = "a_Population",
  Detection_method_2 = "a_Serology",
  stringsAsFactors   = FALSE
) |>
  mutate(HP = NA_integer_, row_type = "predicted") |>
  left_join(hdiMapping, by = c("year_idx", "Subregion")) |>
  left_join(pop, by = c("Subregion", "agp_idx", "year_idx"))

input2 <- bind_rows(input, pred_grid) |>
  mutate(
    HP                 = as.integer(round(HP)),
    Subregion_idx      = as.integer(as.factor(Subregion)),
    Study_design_3     = as.factor(Study_design_3),
    HDI_rank_dummy     = as.factor(HDI_rank_dummy),
    Setting            = as.factor(Setting),
    Detection_method_2 = as.factor(Detection_method_2),
    agp_idx            = as.factor(agp_idx),
    year_linear        = year_idx,
    age_sub_idx        = as.integer(as.factor(paste0(agp_idx, "_", Subregion_idx)))
  )


# ---------------------------------------------------------------------------
# Fit
# ---------------------------------------------------------------------------

fit <- fit_hp_prevalence(input2)
summary(fit)

output <- input2 |>
  mutate(
    HP_Prev_smoothed = pmin(fit$summary.fitted.values$mean, 1),
    HP_Prev_lower    = pmin(fit$summary.fitted.values$`0.025quant`, 1),
    HP_Prev_upper    = pmin(fit$summary.fitted.values$`0.975quant`, 1)
  )

write.csv(output, "output/hpPrev_output.csv", row.names = FALSE)

# prediction grid on its own, for downstream burden calculations
output_pred <- output |> filter(row_type == "predicted")
write.csv(output_pred, "output/hpPrev_prediction_grid.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# Fit checks
# ---------------------------------------------------------------------------

fitcheck <- hp_pearson_residuals(input2, fit)

plot_pearson_resid(fitcheck, year_idx,  "Year")
plot_pearson_resid(fitcheck, agp_idx,   "Age group")
plot_pearson_resid(fitcheck, Subregion, "Subregion", rotate_x = TRUE)


# ---------------------------------------------------------------------------
# Trend figures
# ---------------------------------------------------------------------------

plot_prev_trend(output, 0)
plot_prev_trend(output, 1)
plot_prev_trend(output, 2)
