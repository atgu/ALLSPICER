###############################################################################
## figureS22.R
##
## Figure S22 -- relationship between the AF-weighted variant effect correlation
## and the ALLSPICE p-value across the simulation settings (effect-generating
## mode x effect-size sigma x number of variants). Shows that, under the
## proportional null, estimation noise alone does not produce significance, and
## that the departure from proportionality needed to reach significance shrinks
## with the number of variants.
##
## Extracted from `R/simulation_analysis.R` (Parts 1-3). Self-contained: reads the
## per-replicate simulation outputs produced by figures/figure1-2_figureS1-S7.R
## (<stub>_test_results.csv and <stub>_beta.csv in result_path).
##
## OUTPUTS (to result_path / figure_path)
##   sim_effect_corr_vs_pvalue.csv   per-cell median p (+ IQR/range) vs effect corr
##   sim_detectable_effect_corr.csv  effect corr at which median p hits a target
##                                   (pooled across the observed sim modes)
##   figureS22.png                   Figure S22 (the per-setting summary)
###############################################################################

## ---------------------------------------------------------------------------
## Setup
## ---------------------------------------------------------------------------
source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/constants.R')
source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/simulations.R')
## figure_path and result_path come from constants.R (the main figure directory).
library(ggpubr)

## Re-attach the tidyverse LAST so its verbs win on the search path. The sourced
## files pull in plyr/MASS, which mask dplyr::summarize/select; in particular
## plyr::summarise silently drops the grouping columns.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(tibble)
  library(ggplot2)
})

PFLOOR <- 1e-300   # floor for p-values that underflow to 0, so -log10 stays finite

#' AF-weighted, uncentered ("cosine") correlation of two effect-size vectors.
#'
#' This is the quantity ALLSPICE effectively tests: cos angle between the
#' AF-weighted effect vectors, == +-1 iff the effects are exactly proportional
#' (b1 = c*b2). Unlike Pearson correlation it does NOT center the vectors
#' (ALLSPICE's null is a line through the origin), and it weights each variant by
#' its allele frequency (ALLSPICE uses A = 2*diag(AF)).
#' @param x,y numeric effect-size vectors (one entry per variant).
#' @param w   non-negative weights (allele frequency); non-finite entries dropped.
#' @return scalar in [-1, 1]; NA if a vector has zero weighted norm.
w_cos <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w)
  x <- x[ok]; y <- y[ok]; w <- w[ok]
  d <- sqrt(sum(w * x^2) * sum(w * y^2))
  if (d == 0) NA_real_ else sum(w * x * y) / d
}

#' NA-safe min / max. Cells whose p-values are all NA (e.g. r = +-1, where
#' ALLSPICE is degenerate) return NA instead of +-Inf with a warning.
safe_min <- function(x) { x <- x[is.finite(x)]; if (length(x)) min(x) else NA_real_ }
safe_max <- function(x) { x <- x[is.finite(x)]; if (length(x)) max(x) else NA_real_ }

## ---------------------------------------------------------------------------
## Part 1: summarize each simulation setting (effect correlation vs ALLSPICE p)
## ---------------------------------------------------------------------------

## The 12 settings used in figures/figure1-2_figureS1-S7.R: 4 effect-generating
## modes x 3 (sigma) parameter combos. `stub` is the file prefix on disk; `sigma`
## is parsed from the combo name for labelling.
##   null_mle   : effects exactly proportional (the null is TRUE) -> calibration
##   alt        : effects independent  (cor ~ 0)
##   correlated : effects correlated with a fixed rho
##   nonlinear  : effects in a non-linear (b2 + b2^2) relationship
settings <- tidyr::expand_grid(
  sim_type = c('null_mle', 'alt', 'correlated', 'nonlinear'),
  combo    = c('par_combo1_sigma1', 'par_combo2_sigma0.1', 'par_combo3_sigma0.01')
) %>%
  mutate(stub  = paste0(sim_type, '_sim100_', combo),
         sigma = sub('.*_sigma', '', combo))

#' Summarize one simulation setting: variant effect correlation vs ALLSPICE p.
#'
#' Reads the setting's beta + test-results files, computes a per-replicate
#' variant effect correlation from the beta file (both ordinary Pearson and the
#' ALLSPICE-aligned AF-weighted cosine; both observed-sumstats and true-effect
#' versions), attaches each replicate's ALLSPICE p-value, then aggregates to one
#' row per parameter cell (par_set = unique n_var x r x pi [x c]).
#'
#' @param stub file prefix, e.g. "correlated_sim100_par_combo1_sigma1".
#' @return tibble with one row per par_set: medians of the four effect-correlation
#'   measures; n_rep; n_valid_p (finite p-values); median / IQR / min / max of the
#'   p-value; frac_p_lt_0.05; and `stub`. Returns NULL if either file is missing.
summarize_setting <- function(stub) {
  res_f  <- paste0(result_path, stub, '_test_results.csv')
  beta_f <- paste0(result_path, stub, '_beta.csv')
  if (!file.exists(res_f) || !file.exists(beta_f)) {
    message('skip (missing files): ', stub); return(NULL)
  }
  res  <- readr::read_csv(res_f,  show_col_types = FALSE)   # one row per replicate
  beta <- readr::read_csv(beta_f, show_col_types = FALSE)   # one row per variant

  ## (1) Per-replicate variant effect correlation across the gene's variants.
  ##     *_obs  use the estimated summary-stat effects (what ALLSPICE sees);
  ##     *_true use the simulated true effects (the underlying setting);
  ##     *_afw  is the AF-weighted cosine (ALLSPICE-aligned, == +-1 iff proportional);
  ##     the plain versions are ordinary Pearson correlation.
  eff <- beta %>%
    dplyr::group_by(par_set, simulation, n_var, r, pi) %>%
    dplyr::summarize(
      eff_corr_obs      = suppressWarnings(cor(b1_sum_stats,  b2_sum_stats)),
      eff_corr_true     = suppressWarnings(cor(b1_true_value, b2_true_value)),
      eff_corr_obs_afw  = w_cos(b1_sum_stats,  b2_sum_stats,  AF),
      eff_corr_true_afw = w_cos(b1_true_value, b2_true_value, AF),
      .groups = 'drop')

  ## (2) Attach the matching ALLSPICE p-value, then aggregate per parameter cell.
  ##     IQR (q25/q75) and min/max describe the p-value *range*; n_valid_p counts
  ##     replicates with a finite p (0 at the r = +-1 degenerate cells).
  eff %>%
    dplyr::inner_join(res %>% dplyr::select(par_set, simulation, pvalue),
                      by = c('par_set', 'simulation')) %>%
    dplyr::group_by(par_set, n_var, r, pi) %>%
    dplyr::summarize(
      n_rep                    = dplyr::n(),
      median_eff_corr_obs      = median(eff_corr_obs,      na.rm = TRUE),
      median_eff_corr_true     = median(eff_corr_true,     na.rm = TRUE),
      median_eff_corr_obs_afw  = median(eff_corr_obs_afw,  na.rm = TRUE),
      median_eff_corr_true_afw = median(eff_corr_true_afw, na.rm = TRUE),
      n_valid_p                = sum(is.finite(pvalue)),
      median_pvalue            = median(pvalue, na.rm = TRUE),
      pvalue_q25               = quantile(pvalue, 0.25, na.rm = TRUE),
      pvalue_q75               = quantile(pvalue, 0.75, na.rm = TRUE),
      pvalue_min               = safe_min(pvalue),
      pvalue_max               = safe_max(pvalue),
      frac_p_lt_0.05           = mean(pvalue < 0.05, na.rm = TRUE),
      .groups = 'drop') %>%
    mutate(stub = stub)
}

## Run over all settings and assemble the master table (degenerate cells kept;
## flagged by n_valid_p == 0).
summary_tab <- purrr::map_dfr(settings$stub, summarize_setting) %>%
  dplyr::left_join(settings, by = 'stub') %>%
  dplyr::select(sim_type, sigma, par_set, n_var, r, pi, n_rep, n_valid_p,
                median_eff_corr_obs, median_eff_corr_true,
                median_eff_corr_obs_afw, median_eff_corr_true_afw,
                median_pvalue, pvalue_q25, pvalue_q75, pvalue_min, pvalue_max,
                frac_p_lt_0.05) %>%
  dplyr::arrange(sim_type, sigma, n_var, r, pi)

readr::write_csv(summary_tab, paste0(result_path, 'sim_effect_corr_vs_pvalue.csv'))
cat('settings summarized:', length(unique(paste(summary_tab$sim_type, summary_tab$sigma))), '\n')
print(summary_tab, n = 60)

## ---------------------------------------------------------------------------
## Part 2: invert median-p -> required effect correlation
## ---------------------------------------------------------------------------

#' Effect correlation at which the median p-value reaches a target.
#'
#' Within a group, the rows supply (median_eff_corr_obs_afw, median_pvalue)
#' points; this monotone-interpolates the AF-weighted effect correlation at which
#' the median p equals `target_p`. Non-finite rows (degenerate cells) are dropped;
#' returns NA if there are < 2 usable points or the target is outside the observed
#' p-range (no extrapolation).
#' @param d        data frame for one group (needs median_pvalue, median_eff_corr_obs_afw).
#' @param target_p target p-value (e.g. 0.05, or 4.23e-6 for genome-wide).
#' @return scalar AF-weighted effect correlation, or NA.
invert_p <- function(d, target_p) {
  d <- d %>% dplyr::arrange(median_pvalue) %>%
    dplyr::filter(is.finite(median_pvalue), is.finite(median_eff_corr_obs_afw))
  if (nrow(d) < 2 || target_p < min(d$median_pvalue) || target_p > max(d$median_pvalue)) return(NA_real_)
  stats::approx(d$median_pvalue, d$median_eff_corr_obs_afw, xout = target_p, ties = mean)$y
}

## Pool ACROSS sim modes within each (sigma, n_var): the modes place points at
## different effect correlations (null ~ +-1, correlated ~ rho, nonlinear ~
## intermediate, independent ~ 0), which together form the correlation axis to
## invert. Caveat: this assumes p depends on effect correlation regardless of the
## *shape* of the departure -- the rho-sweep (Figure S21) avoids that assumption.
detect_tab <- summary_tab %>%
  group_by(sigma, n_var) %>%
  group_modify(~ tibble::tibble(corr_afw_at_p_0.05 = invert_p(.x, 0.05),
                                corr_afw_at_p_gw   = invert_p(.x, 4.23e-6))) %>%
  ungroup() %>% arrange(sigma, n_var)
readr::write_csv(detect_tab, paste0(result_path, 'sim_detectable_effect_corr.csv'))
cat('\n-- detectable AF-weighted effect correlation at p thresholds (pooled across modes) --\n')
print(detect_tab, n = 40)

## ---------------------------------------------------------------------------
## Part 3: figure -- median ALLSPICE p-value vs variant effect correlation
## ---------------------------------------------------------------------------
## ===========================================================================
## FIGURE S22 CAPTION (for the Supplementary Information)
## ---------------------------------------------------------------------------
## Figure S22. Relationship between the AF-weighted variant effect correlation and
## the ALLSPICE p-value across simulation settings. Each point is one parameter
## cell (a unique combination of number of variants, phenotypic correlation r, and
## mixture weight pi), summarized over 100 replicate simulations: the x-axis is the
## median AF-weighted, uncentered correlation between the two phenotypes'
## per-variant effect vectors (= 1 iff exactly proportional) and the y-axis is the
## median ALLSPICE -log10(p), with vertical bars spanning the p-value interquartile
## range across replicates. Colour denotes the number of variants per gene
## (5, 10, 20, 100). Panels are arranged by effect-size standard deviation (rows:
## sigma = 1, 0.1, 0.01) and effect-generating mechanism (columns: Null /
## proportional, Independent, Correlated, Non-linear). The dashed line marks
## nominal significance (p = 0.05). Under the proportional null (left column) the
## median p stays near 0.5 across the entire correlation axis -- including where
## estimation noise pushes the observed correlation away from 1 -- so estimation
## noise around a truly proportional relationship does not by itself produce
## significance; in the other mechanisms the p-value falls as the effects depart
## from proportionality, and the departure needed to reach significance shrinks as
## the number of variants grows. Degenerate cells (r = +-1, where ALLSPICE's
## c^2 - 2cr + 1 denominator collapses and p is undefined) are omitted.
## Simulations use n_ind = 1,000 individuals.
## ===========================================================================

## warn about / plot missing values.
plot_dat <- summary_tab %>%
  dplyr::filter(is.finite(median_pvalue), is.finite(median_eff_corr_obs_afw)) %>%
  mutate(sim_type = factor(sim_type, levels = c('null_mle', 'alt', 'correlated', 'nonlinear'),
                           labels = c('Null (proportional)', 'Independent',
                                      'Correlated', 'Non-linear')),
         sigma = factor(paste0('sigma = ', sigma),
                        levels = c('sigma = 1', 'sigma = 0.1', 'sigma = 0.01')),
         n_var = factor(n_var))

## x = AF-weighted effect correlation (=1 iff exactly proportional);
## y = -log10(median p); error bars span the p-value IQR; colour = n_var;
## facets = sigma (rows) x mode (cols). Dashed line = nominal p = 0.05.
p <- ggplot(plot_dat, aes(x = median_eff_corr_obs_afw, y = -log10(pmax(median_pvalue, PFLOOR)),
                          color = n_var)) +
  geom_errorbar(aes(ymin = -log10(pmax(pvalue_q75, PFLOOR)),
                    ymax = -log10(pmax(pvalue_q25, PFLOOR))), width = 0.03, alpha = 0.5) +
  geom_point(size = 2) +
  geom_hline(yintercept = -log10(0.05), lty = 2, colour = 'grey40') +
  facet_grid(sigma ~ sim_type, scale = 'free_y') +
  scale_color_viridis_d(name = 'Number of variants', option = 'C', end = 0.9) +
  labs(x = 'Median AF-weighted variant effect correlation',
       y = expression(-log[10](median~p_ALLSPICE)~"(IQR error bars)")) +
  theme_classic(base_size = 12) +
  theme(legend.position = 'top',
        strip.background = element_blank(), strip.text = element_text(face = 'bold'))

png(paste0(figure_path, 'figureS22.png'), height = 7, width = 9, units = 'in', res = 300)
print(p)
dev.off()
cat('\nwrote:', paste0(figure_path, 'figureS22.png'), '\n')
