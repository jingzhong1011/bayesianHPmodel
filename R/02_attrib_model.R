# ---------------------------------------------------------------------------
# 02 - Non-attributable gastric cancer burden model
#
# For each burden-draw iteration, anatomical subsite and sex, fits a Poisson
# space-age-period model with a Knorr-Held type IV interaction and retains
# posterior draws of the non-attributable rate.
#
# 2026-07-17: type IV interaction corrected via generic0 + Kronecker structure.
#
# Run from the project root:  source("R/02_attrib_model.R")
# ---------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(tidyr)
library(INLA)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(future.apply)

source("R/utils_model.R")
source("R/utils_plot.R")

# tune to the machine; each INLA call is pinned to one thread via num.threads
N_WORKERS <- 4
N_DRAWS   <- 100
N_ITERS   <- 100

TEMP_DIR <- "output/temp_results"
LOG_DIR  <- "output/logs"

plan(multisession, workers = N_WORKERS)

dir.create(TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR,  recursive = TRUE, showWarnings = FALSE)


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

nonattrGC_subregion <- readRDS(
  "data/nonattrGC_input.RDS"
) |>
  ungroup() |>
  na.omit()

pop <- readRDS("data/pop_subregion.RDS") |>
  ungroup() |>
  na.omit() |>
  select(Subregion, age, sex, year, population)

yearMapping <- data.frame(
  year     = paste0(seq(2000, 2035, 5), "-", seq(2004, 2039, 5)),
  year_idx = 1:8
)

ageMapping <- data.frame(
  age     = c(paste0(seq(20, 80, 5), "-", seq(24, 84, 5)), "85+"),
  age_idx = 1:14
)

# age_idx2 and age_idx3 are deliberate copies: INLA needs a distinct index
# variable for each f() term that uses the same underlying covariate.
add_model_indices <- function(df) {
  df |>
    mutate(
      Subregion_idx = as.integer(as.factor(Subregion)),
      age_idx2      = as.integer(age_idx),
      age_idx3      = as.integer(age_idx),
      year_idx2     = as.integer(year_idx)
    )
}

obs_all_iters <- nonattrGC_subregion |>
  group_by(iter, year, Subregion, sex, age, V2) |>
  summarise(
    non_attrib_idx = sum(round(pmax(0, non_attrib))),
    population     = sum(pop_total),
    .groups = "drop"
  ) |>
  left_join(yearMapping, by = "year") |>
  left_join(ageMapping, by = "age") |>
  add_model_indices()

pred_combo <- expand_grid(
  # character, not factor: joining a factor column against the character
  # `year` in pop/yearMapping relies on silent coercion
  year      = paste0(seq(2000, 2035, 5), "-", seq(2004, 2039, 5)),
  Subregion = levels(as.factor(nonattrGC_subregion$Subregion)),
  sex       = c("Female", "Male"),
  age       = levels(as.factor(nonattrGC_subregion$age)),
  V2        = levels(as.factor(nonattrGC_subregion$V2))
)

pred_data <- pred_combo |>
  left_join(pop, by = c("year", "Subregion", "sex", "age")) |>
  filter(!is.na(population)) |>
  mutate(non_attrib_idx = NA_real_) |>
  left_join(yearMapping, by = "year") |>
  left_join(ageMapping, by = "age") |>
  add_model_indices()


# ---------------------------------------------------------------------------
# Type IV interaction structure
#
# Depends only on the grid dimensions, which are identical across all 400
# tasks, so it is built once here rather than inside the worker.
# ---------------------------------------------------------------------------

n_age_global  <- nrow(ageMapping)   # 14
n_year_global <- nrow(yearMapping)  # 8

typeIV <- build_typeIV_components(n_main = n_age_global, n_group = n_year_global)

cat(sprintf(
  "Type IV interaction structure pre-computed: %d x %d, rankdef = %d\n",
  nrow(typeIV$Cmatrix), ncol(typeIV$Cmatrix), typeIV$rankdef
))


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------

task_grid <- expand.grid(
  iter_idx = 1:N_ITERS,
  v2_cat   = c("Cardia", "Noncardia"),
  sex_cat  = c("Male", "Female"),
  stringsAsFactors = FALSE
)

# results are written per task so an interrupted run can be resumed
invisible(future_lapply(seq_len(nrow(task_grid)), function(task_id) {

  current_iter <- task_grid$iter_idx[task_id]
  current_v2   <- task_grid$v2_cat[task_id]
  current_sex  <- task_grid$sex_cat[task_id]

  out_file <- file.path(
    TEMP_DIR,
    sprintf("res_iter%d_%s_%s.RDS", current_iter, current_v2, current_sex)
  )

  if (file.exists(out_file)) return(NULL)

  res <- run_attrib_inla(
    iter_idx   = current_iter,
    target_v2  = current_v2,
    target_sex = current_sex,
    obs_df     = obs_all_iters,
    pred_df    = pred_data,
    typeIV     = typeIV,
    n_draws    = N_DRAWS,
    log_dir    = LOG_DIR
  )

  if (!is.null(res)) saveRDS(res, out_file)

  rm(res)
  gc()

  NULL
}, future.seed = TRUE))


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

read_logs <- function(kind) {
  files <- list.files(LOG_DIR, pattern = paste0("^", kind, "_"), full.names = TRUE)
  if (length(files) == 0) return(NULL)
  rbindlist(lapply(files, fread), fill = TRUE)
}

fail_log <- read_logs("fail")
if (!is.null(fail_log)) {
  cat(sprintf("%d / %d iteration-strata combinations failed to fit.\n",
              nrow(fail_log), nrow(task_grid)))
}

conv_log <- read_logs("convergence")
if (!is.null(conv_log)) {
  cat(sprintf("%d / %d iteration-strata combinations flagged non-convergence.\n",
              nrow(conv_log), nrow(task_grid)))
}

completed <- length(list.files(TEMP_DIR, pattern = "^res_"))
cat(sprintf("%d / %d tasks produced output.\n", completed, nrow(task_grid)))


# ---------------------------------------------------------------------------
# Combine and summarise
# ---------------------------------------------------------------------------

files <- list.files(TEMP_DIR, pattern = "^res_", full.names = TRUE)
nonattrGC_output <- rbindlist(lapply(files, readRDS))

draw_cols <- grep("^draw_", names(nonattrGC_output), value = TRUE)
id_cols   <- c("year", "Subregion", "sex", "age", "V2", "iter",
               "non_attrib_rate_obs")

nonattrGC_output_long <- melt(
  nonattrGC_output,
  id.vars       = id_cols,
  measure.vars  = draw_cols,
  variable.name = "draw_id",
  value.name    = "draw_value"
)

# credible intervals pool across both posterior draws and burden iterations
nonattrGC_result <- nonattrGC_output_long[, .(
  pred_mean  = mean(draw_value),
  pred_lower = quantile(draw_value, probs = 0.025),
  pred_upper = quantile(draw_value, probs = 0.975),
  obs_mean   = mean(non_attrib_rate_obs, na.rm = TRUE),
  obs_lower  = quantile(non_attrib_rate_obs, probs = 0.025, na.rm = TRUE),
  obs_upper  = quantile(non_attrib_rate_obs, probs = 0.975, na.rm = TRUE)
), by = .(year, Subregion, sex, age, V2)]

saveRDS(nonattrGC_output, "output/nonattrGC_output.RDS")
write.csv(nonattrGC_result, "output/nonattrGC_result.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# Trend figures
# ---------------------------------------------------------------------------

plot_nonattr_trend(nonattrGC_result, "Cardia",    "Male",   "25-29")
# plot_nonattr_trend(nonattrGC_result, "Noncardia", "Male",   "25-29")
# plot_nonattr_trend(nonattrGC_result, "Noncardia", "Male",   "40-44")
# plot_nonattr_trend(nonattrGC_result, "Noncardia", "Female", "40-44")
