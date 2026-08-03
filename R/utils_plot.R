# ---------------------------------------------------------------------------
# Plotting utilities shared by R/01_hp_model.R and R/02_attrib_model.R
#
# Sourced, not run directly. Requires: ggplot2, dplyr, data.table.
# ---------------------------------------------------------------------------


#' House theme for every figure in this project
#'
#' @param base_size Base font size.
#' @return A ggplot2 theme.
theme_hp <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
      strip.text       = ggplot2::element_text(face = "bold", size = 8),
      panel.grid.minor = ggplot2::element_blank()
    )
}


#' Boxplot of Pearson residuals against a grouping variable
#'
#' Medians close to zero indicate the model is tracking the data. Complex
#' random effects make the posterior non-Gaussian, which inflates p_W in WAIC;
#' the residuals are the more reliable check of fit here.
#'
#' @param data Output of `hp_pearson_residuals()`.
#' @param group_var Bare column name to group by.
#' @param label Axis label for the grouping variable.
#' @param rotate_x Rotate x-axis labels 45 degrees.
#' @return A ggplot object.
plot_pearson_resid <- function(data, group_var, label, rotate_x = FALSE) {
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = as.factor({{ group_var }}), y = pearson_resid)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
    ggplot2::labs(
      title = paste0("Pearson residuals vs ", label),
      x     = label,
      y     = "Pearson residual"
    ) +
    theme_hp()

  if (rotate_x) {
    p <- p + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  }

  p
}


#' Fitted H. pylori prevalence trend by subregion, for one age group
#'
#' Observed study estimates are drawn as points, the smoothed prediction grid
#' as a line with a 95% credible ribbon. Rows are separated using the
#' `row_type` flag set in `R/01_hp_model.R` rather than by position, so
#' adding or dropping studies cannot silently shift which rows are plotted.
#'
#' @param data Model output containing `row_type`, `agp_idx`, `year_idx` and
#'   the `HP_Prev_*` columns.
#' @param target_agp Age group index: 0 (20-50), 1 (50-70) or 2 (70+).
#' @return A ggplot object.
plot_prev_trend <- function(data, target_agp) {
  age_text <- switch(
    as.character(target_agp),
    "0" = "20-50",
    "1" = "50-70",
    "2" = "70+",
    stop("target_agp must be 0, 1 or 2.")
  )

  in_group <- data |>
    dplyr::filter(as.character(agp_idx) == as.character(target_agp))

  obs_data <- in_group |>
    dplyr::filter(row_type == "observed", !is.na(HP_prevalence)) |>
    dplyr::arrange(Subregion, year_idx)

  pred_data <- in_group |>
    dplyr::filter(row_type == "predicted") |>
    dplyr::arrange(Subregion, year_idx)

  if (nrow(pred_data) == 0) {
    stop("No prediction rows found for this age group.")
  }

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = pred_data,
      ggplot2::aes(
        x = year_idx, ymin = HP_Prev_lower, ymax = HP_Prev_upper,
        fill = Subregion
      ),
      alpha = 0.2, colour = NA
    ) +
    ggplot2::geom_line(
      data = pred_data,
      ggplot2::aes(x = year_idx, y = HP_Prev_smoothed, colour = Subregion),
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = obs_data,
      ggplot2::aes(x = year_idx, y = HP_prevalence, colour = Subregion),
      size = 1, alpha = 0.7
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(min(pred_data$year_idx), max(pred_data$year_idx), by = 2)
    ) +
    ggplot2::labs(
      title = paste0("HP prevalence trend by subregion (age ", age_text, ")"),
      x     = "Year index",
      y     = "Prevalence"
    ) +
    theme_hp() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::facet_wrap(~ Subregion, ncol = 5)
}


#' Non-attributable gastric cancer rate trend by subregion
#'
#' @param data Summarised output of `R/02_attrib_model.R`, one row per
#'   year-subregion-sex-age-subsite with `pred_*` and `obs_*` columns.
#' @param target_v2 Anatomical subsite ("Cardia" or "Noncardia").
#' @param target_sex "Male" or "Female".
#' @param target_age Five-year age band label, e.g. "25-29".
#' @return A ggplot object.
plot_nonattr_trend <- function(data, target_v2, target_sex, target_age) {
  data.table::setDT(data)

  plot_data <- data[V2 == target_v2 & sex == target_sex & age == target_age]

  if (nrow(plot_data) == 0) {
    stop("No rows match the requested subsite, sex and age band.")
  }

  # rates are modelled per person; report per 100,000
  rate_cols <- c("pred_mean", "pred_lower", "pred_upper",
                 "obs_lower", "obs_upper")
  plot_data[, (rate_cols) := lapply(.SD, function(x) x * 1e5), .SDcols = rate_cols]

  plot_data <- plot_data[order(Subregion, year)]

  if (is.character(plot_data$year) || is.factor(plot_data$year)) {
    plot_data[, year_idx := as.numeric(as.factor(year))]
    x_breaks <- unique(plot_data$year_idx)
    x_labels <- unique(as.character(plot_data$year))
  } else {
    plot_data[, year_idx := year]
    x_breaks <- ggplot2::waiver()
    x_labels <- ggplot2::waiver()
  }

  ggplot2::ggplot(plot_data, ggplot2::aes(x = year_idx)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = pred_lower, ymax = pred_upper, fill = Subregion),
      alpha = 0.2, colour = NA
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = pred_mean, colour = Subregion), linewidth = 1
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = obs_lower, ymax = obs_upper, colour = Subregion),
      alpha = 0.6, linewidth = 0.8
    ) +
    ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    ggplot2::labs(
      title = sprintf(
        "Non-attributable rate trend by subregion\n(%s | %s | Age: %s)",
        target_v2, target_sex, target_age
      ),
      x = "Year",
      y = "Rate (per 100,000 population)"
    ) +
    theme_hp() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::facet_wrap(~ Subregion, ncol = 5, scales = "free_y")
}
