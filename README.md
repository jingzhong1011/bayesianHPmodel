# Bayesian estimation of *H. pylori* prevalence and non-attributable gastric cancer burden

Subnational Bayesian models for two linked quantities: 
- Trend in
*Helicobacter pylori* prevalence across world subregions
- Share of gastric cancer burden that is *not* attributable to *H. pylori* infection.  
  
Both are fitted with Integrated Nested Laplace Approximation (INLA).

## Models

### Prevalence (`R/01_hp_model.R`)

Seropositive counts are modelled as negative-binomial with population as an
offset, pooling age-specific estimates from published serological studies.

- Fixed effects for study design, setting, detection method, HDI rank and age group, so that heterogeneous study protocols can be projected onto a common reference profile (population-based cross-sectional serology)
- A linear period trend, plus a subregion-replicated first-order random walk over period
- IID subregion and age-by-subregion effects
- Penalised-complexity priors throughout
- HDI is extrapolated forward on the logit scale, which keeps projected values inside the unit interval

Observed study rows and the prediction grid are stacked into a single frame
and fitted in one pass, so predictions carry the full posterior uncertainty. A
`row_type` flag distinguishes the two.

### Non-attributable burden (`R/02_attrib_model.R`)

Non-attributable case counts are modelled as Poisson with population as an
offset, fitted separately for each combination of burden-draw iteration,
anatomical subsite (cardia / noncardia) and sex — 400 fits in total.

- Main effects: IID subregion, RW1 over age, RW1 over period
- Subregion-replicated RW1 over age and over period
- A Knorr-Held **type IV** age-period interaction

The type IV term is the substantive methodological piece. INLA's built-in
`constr = TRUE` cannot express its constraint structure, so the interaction is
specified as `generic0` with an explicitly supplied structure matrix
(the Kronecker product of two RW1 structure matrices), an explicit rank
deficiency of `n_age + n_period - 1`, and a matching set of extra sum-to-zero
constraints. That structure depends only on the grid dimensions, so it is
built once and reused across all 400 fits rather than recomputed each time.

Uncertainty is propagated by drawing from the posterior of each fit and
pooling across both posterior draws and burden iterations.

## Repository layout

```
R/
  utils_model.R   PC priors, RW1 and type IV structures, INLA fitting functions
  utils_plot.R    shared theme and figure functions
  01_hp_model.R   seroprevalence model, end to end
  02_attrib_model.R  non-attributable burden model, end to end
data/
output/
```

Scripts are run from the project root and source the `R/` files themselves.

## Running it

The burden model runs 400 INLA fits in parallel via `future.apply`. Each fit is
pinned to a single thread, so `N_WORKERS` at the top of
`R/02_attrib_model.R` should be set to the number of cores available.
Results are written per task under `output/temp_results/`, so an interrupted
run resumes where it stopped rather than starting over. Failures and
non-convergence flags are logged per task under `output/logs/` and summarised
at the end of the script.

```r
source("R/01_hp_model.R")
source("R/02_attrib_model.R")
```

## Dependencies

`INLA`, `Matrix`, `data.table`, `dplyr`, `tidyr`, `readxl`, `ggplot2`,
`ggpubr`, `future.apply`.

INLA is not on CRAN:

```r
install.packages(
  "INLA",
  repos = c(INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE
)
```

## Data

All inputs derive from published literature, UN population estimates and UNDP HDI estimates. No individual-level records are used. The input files
themselves are not distributed with this repository; `README.md` documents
the file names, required columns and coding conventions the scripts expect.

## License

This work is licensed under [CC BY 4.0](LICENSE).
