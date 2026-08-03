# ---------------------------------------------------------------------------
# Model utilities shared by R/01_hp_model.R and R/02_attrib_model.R
#
# Sourced, not run directly. Requires: INLA, Matrix, data.table.
# ---------------------------------------------------------------------------


#' Penalised-complexity prior for a precision hyperparameter
#'
#' Encodes the statement P(sigma > u) = a on the marginal standard deviation
#' of the random effect.
#'
#' @param u Scale of the standard deviation.
#' @param a Tail probability.
#' @return A list suitable for the `hyper` argument of `INLA::f()`.
pc_hpr <- function(u, a) {
  list(prec = list(prior = "pc.prec", param = c(u, a)))
}


#' Structure matrix of a first-order random walk
#'
#' @param m Number of ordered levels.
#' @return An `m x m` symmetric structure matrix.
build_rw1_structure <- function(m) {
  stopifnot(m >= 2)

  R <- matrix(0, m, m)
  diag(R) <- c(1, rep(2, m - 2), 1)

  for (i in seq_len(m - 1)) {
    R[i, i + 1] <- -1
    R[i + 1, i] <- -1
  }

  R
}


#' Sum-to-zero constraints for a Knorr-Held type IV interaction
#'
#' A type IV interaction (RW1 x RW1) has a rank deficiency of
#' `n_main + n_group - 1`, which INLA cannot infer on its own. This builds the
#' matching set of linear constraints: the interaction must sum to zero within
#' every level of the main index and within every level of the grouping index
#' (one group constraint is dropped, since it is implied by the others).
#'
#' Index convention: position = (group - 1) * n_main + main, i.e. the main
#' index varies fastest. `build_typeIV_components()` builds the structure
#' matrix under the same convention.
#'
#' @param n_main Number of levels of the fast-varying index (e.g. age).
#' @param n_group Number of levels of the slow-varying index (e.g. year).
#' @return A list with constraint matrix `A` and right-hand side `e`.
build_typeIV_constraint <- function(n_main, n_group) {
  m  <- n_main
  n  <- n_group
  mn <- m * n

  A <- matrix(0, nrow = n + m - 1, ncol = mn)

  # one constraint per level of the main index
  for (i in seq_len(m)) {
    A[i, seq(1, n * m, by = m) + (i - 1)] <- 1
  }

  # one constraint per level of the grouping index, minus one
  for (i in seq_len(n - 1)) {
    A[i + m, (1:m) + (i - 1) * m] <- 1
  }

  list(A = A, e = rep(0, n + m - 1))
}


#' Pre-compute the components of a type IV interaction term
#'
#' The structure matrix, constraints and rank deficiency depend only on the
#' dimensions of the age-year grid, which are fixed across all strata. Building
#' them once and passing them into the fitting function avoids recomputing an
#' identical Kronecker product for every task.
#'
#' @param n_main Number of levels of the fast-varying index (e.g. age).
#' @param n_group Number of levels of the slow-varying index (e.g. year).
#' @return A list with `Cmatrix`, `constr` and `rankdef`, ready to be handed to
#'   `INLA::f(..., model = "generic0")`.
build_typeIV_components <- function(n_main, n_group) {
  R_main  <- build_rw1_structure(n_main)
  R_group <- build_rw1_structure(n_group)

  # kronecker(group, main) puts the main index in the fast-varying position
  R_full <- Matrix::Matrix(kronecker(R_group, R_main), sparse = TRUE)

  constr  <- build_typeIV_constraint(n_main, n_group)
  rankdef <- n_group + n_main - 1

  stopifnot(ncol(constr$A) == ncol(R_full))
  stopifnot(rankdef == nrow(constr$A))

  list(
    Cmatrix = INLA::inla.scale.model(R_full, constr = constr),
    constr  = constr,
    rankdef = rankdef,
    n_main  = n_main,
    n_group = n_group
  )
}


#' Append a one-row diagnostic log for a single fitting task
#'
#' Each task writes its own file rather than appending to a shared one. Under
#' `future_lapply()` several workers run concurrently, and concurrent appends to
#' a single CSV can interleave mid-line and corrupt the file.
#'
#' @param log_dir Directory to write into; created if absent.
#' @param kind Log category, used as the filename prefix (e.g. "fail").
#' @param iter_idx,target_v2,target_sex Task identifiers.
#' @param ... Further columns to record.
#' @return Invisibly, the path written.
write_task_log <- function(log_dir, kind, iter_idx, target_v2, target_sex, ...) {
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  path <- file.path(
    log_dir,
    sprintf("%s_iter%d_%s_%s.csv", kind, iter_idx, target_v2, target_sex)
  )

  utils::write.csv(
    data.frame(
      iter = iter_idx,
      v2   = target_v2,
      sex  = target_sex,
      time = Sys.time(),
      ...,
      stringsAsFactors = FALSE
    ),
    path,
    row.names = FALSE
  )

  invisible(path)
}


#' Fit the H. pylori prevalence model
#'
#' Negative-binomial model for seropositive counts with population as an
#' offset. Fixed effects for study design, setting, HDI rank, detection method
#' and age group; a linear year trend; IID subregion and age-by-subregion
#' effects; and a subregion-replicated RW1 over period.
#'
#' @param data Stacked observed and prediction-grid rows.
#' @return An `inla` object.
fit_hp_prevalence <- function(data) {
  INLA::inla(
    HP ~ 1 +
      Study_design_3 + Setting + HDI_rank_dummy + Detection_method_2 +
      agp_idx +
      f(year_linear,   model = "linear", mean.linear = 0, prec.linear = 1) +
      f(Subregion_idx, model = "iid", hyper = pc_hpr(0.50, 0.50)) +
      f(age_sub_idx,   model = "iid", hyper = pc_hpr(0.50, 0.50)) +
      f(year_idx,      model = "rw1", hyper = pc_hpr(0.50, 0.50),
        scale.model = TRUE, replicate = Subregion_idx, constr = TRUE),
    family = "nbinomial",
    E      = population,
    data   = data,
    control.predictor = list(compute = TRUE, link = 1),
    control.compute   = list(
      config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE,
      return.marginals.predictor = FALSE
    )
  )
}


#' Pearson residuals for a fitted negative-binomial prevalence model
#'
#' @param data The data frame passed to `fit_hp_prevalence()`.
#' @param fit The resulting `inla` object.
#' @return `data` with a `pearson_resid` column, rows with no residual dropped.
hp_pearson_residuals <- function(data, fit) {
  prev_fitted <- fit$summary.fitted.values$mean
  mu_fitted   <- data$population * prev_fitted

  phi_fitted <- fit$summary.hyperpar[
    "size for the nbinomial observations (1/overdispersion)", "mean"
  ]

  nb_sd <- sqrt(mu_fitted + (mu_fitted^2 / phi_fitted))

  out <- data
  out$pearson_resid <- (out$HP - mu_fitted) / nb_sd
  out[!is.na(out$pearson_resid), ]
}


#' Fit the non-attributable gastric cancer model for one iteration-stratum
#'
#' Poisson model for non-attributable case counts with population as an offset,
#' fitted separately for each combination of burden-draw iteration, anatomical
#' subsite and sex. Returns posterior draws of the rate on the prediction grid.
#'
#' @param iter_idx Burden-draw iteration to fit.
#' @param target_v2 Anatomical subsite ("Cardia" or "Noncardia").
#' @param target_sex "Male" or "Female".
#' @param obs_df Observed counts across all iterations.
#' @param pred_df Prediction grid with population attached.
#' @param typeIV Output of `build_typeIV_components()`.
#' @param n_draws Number of posterior draws to retain.
#' @param log_dir Directory for failure and convergence logs.
#' @return A data frame of one row per prediction cell with `n_draws` draw
#'   columns, or `NULL` if the stratum was skipped or the fit failed.
run_attrib_inla <- function(iter_idx,
                            target_v2,
                            target_sex,
                            obs_df,
                            pred_df,
                            typeIV,
                            n_draws = 100,
                            log_dir = "output/logs") {

  n_age  <- typeIV$n_main
  n_year <- typeIV$n_group

  data.table::setDT(obs_df)
  data.table::setDT(pred_df)

  this_obs  <- obs_df[iter == iter_idx & V2 == target_v2 & sex == target_sex]
  this_pred <- pred_df[V2 == target_v2 & sex == target_sex]

  this_input <- merge(
    this_pred[, -c("non_attrib_idx")],
    this_obs[, .(year, Subregion, sex, age, V2, non_attrib_idx)],
    by    = c("year", "Subregion", "sex", "age", "V2"),
    all.x = TRUE
  )
  this_input[, iter := iter_idx]

  if (nrow(this_input) == 0 || all(is.na(this_input$non_attrib_idx))) {
    warning(sprintf(
      "Task skipped: iter = %d, v2 = %s, sex = %s has no valid data.",
      iter_idx, target_v2, target_sex
    ))
    return(NULL)
  }

  # combined age-year index for the type IV interaction; age varies fastest,
  # matching the Kronecker order used in build_typeIV_components()
  this_input[, spacetime_idx := (year_idx - 1) * n_age + age_idx3]

  # the pre-computed structure matrix assumes a complete age-year grid, so a
  # stratum with holes in it would be fitted against the wrong constraints
  grid_check <- data.table::uniqueN(this_input[, .(age_idx3, year_idx)])
  if (grid_check != n_age * n_year) {
    warning(sprintf(
      paste0("iter = %d, v2 = %s, sex = %s: age-year grid incomplete ",
             "(%d combos, expected %d). Skipping."),
      iter_idx, target_v2, target_sex, grid_check, n_age * n_year
    ))
    return(NULL)
  }

  # age_idx / age_idx2 / age_idx3 hold identical values. INLA requires a
  # distinct index variable per f() term, so the main effect, the
  # subregion-replicated effect and the interaction each get their own copy.
  formula <- non_attrib_idx ~ 1 +
    f(Subregion_idx, model = "iid", hyper = pc_hpr(0.5, 0.01)) +
    f(age_idx,  model = "rw1", hyper = pc_hpr(0.5, 0.01),
      scale.model = TRUE, constr = TRUE) +
    f(year_idx, model = "rw1", hyper = pc_hpr(0.5, 0.01),
      scale.model = TRUE, constr = TRUE) +
    f(age_idx2,  model = "rw1", replicate = Subregion_idx,
      hyper = pc_hpr(1, 0.01), scale.model = TRUE, constr = TRUE) +
    f(year_idx2, model = "rw1", replicate = Subregion_idx,
      hyper = pc_hpr(1, 0.01), scale.model = TRUE, constr = TRUE) +
    # generic0 with an explicit rank deficiency and extra constraints: INLA's
    # own constr = TRUE cannot express the type IV structure correctly
    f(spacetime_idx, model = "generic0",
      Cmatrix = typeIV$Cmatrix, rankdef = typeIV$rankdef,
      hyper = pc_hpr(1, 0.01), constr = FALSE, extraconstr = typeIV$constr,
      diagonal = 1e-6)

  fit <- tryCatch(
    INLA::inla(
      formula,
      family = "poisson",
      E      = population,
      data   = this_input,
      control.predictor = list(compute = TRUE, link = 1),
      control.compute   = list(config = TRUE, return.marginals.predictor = FALSE),
      num.threads = "1:1"
    ),
    error = function(e) {
      write_task_log(log_dir, "fail", iter_idx, target_v2, target_sex,
                     error = e$message)
      warning(sprintf(
        "INLA failed at iter = %d, v2 = %s, sex = %s. Error: %s",
        iter_idx, target_v2, target_sex, e$message
      ))
      NULL
    }
  )

  if (is.null(fit)) return(NULL)

  if (fit$mode$mode.status != 0) {
    write_task_log(log_dir, "convergence", iter_idx, target_v2, target_sex,
                   mode_status = fit$mode$mode.status)
    warning(sprintf(
      "Non-convergence flag at iter = %d, v2 = %s, sex = %s",
      iter_idx, target_v2, target_sex
    ))
  }

  samples      <- INLA::inla.posterior.sample(n = n_draws, result = fit)
  pred_indices <- grep("^Predictor:", rownames(samples[[1]]$latent))

  draws_matrix <- sapply(samples, function(x) exp(x$latent[pred_indices, 1]))
  colnames(draws_matrix) <- paste0("draw_", seq_len(n_draws))

  result_df <- this_input |>
    dplyr::mutate(
      non_attrib_rate_obs = non_attrib_idx / population,
      iter = iter_idx
    ) |>
    dplyr::select(
      year, Subregion, sex, age, V2, iter,
      population, non_attrib_idx, non_attrib_rate_obs
    ) |>
    dplyr::bind_cols(as.data.frame(draws_matrix))

  rm(fit, samples)
  gc()

  result_df
}
