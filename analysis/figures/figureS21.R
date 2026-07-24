###############################################################################
## figureS21.R
##
## Figure S21 -- ALLSPICE detection surface: how the smallest detectable
## departure from proportionality (1 - rho) depends on the number of variants,
## the effect size (sigma), and the phenotypic correlation (r).
##
## Extracted from `R/simulation_analysis.R` (Part 4: the dedicated rho-sweep).
## Self-contained: sourcing constants.R + simulations.R supplies the ALLSPICE
## helpers (get_ac_mat/get_af_mat/get_geno_mat/get_true_beta/get_beta_hat/
## get_c_hat/get_mle_beta/get_likelihood_test_stats) and `all_n_var`.
##
## OUTPUTS (to result_path / figure_path)
##   sim_rho_sweep_raw.csv             raw per-replicate output of the rho sweep
##   sim_rho_sweep_summary.csv         per-(rho, n_var, sigma, r) medians + power
##   sim_rho_sweep_detectable_corr.csv detectable rho / AF-weighted corr per cell
##   sim_rho_sweep_detection_curve.png Figure S21 (the detection surface)
###############################################################################

## ---------------------------------------------------------------------------
## Setup
## ---------------------------------------------------------------------------
source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/constants.R')
source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/simulations.R')  # get_*, all_n_var
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
#'
#' @param x,y numeric effect-size vectors (one entry per variant).
#' @param w   non-negative weights (allele frequency). Non-finite entries in
#'            x/y/w are dropped pairwise.
#' @return scalar in [-1, 1]; NA if a vector has zero weighted norm.
w_cos <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w)
  x <- x[ok]; y <- y[ok]; w <- w[ok]
  d <- sqrt(sum(w * x^2) * sum(w * y^2))
  if (d == 0) NA_real_ else sum(w * x * y) / d
}

## ---------------------------------------------------------------------------
## Part 4: dedicated rho-sweep -> smooth correlation->p detection curve
## ---------------------------------------------------------------------------
## The existing settings give only a few points on the effect-correlation axis.
## To map the relationship continuously we vary the TRUE effect correlation
## directly using `mode = 'correlated'`, where
##     b1 = g2 * (rho * z2 + sqrt(1 - rho^2) * z1)   =>   cor(b1, b2) = rho,
## sweep rho on a fine grid x n_var, and read off the smallest departure from
## proportionality (1 - rho) that ALLSPICE detects at each threshold.
## The sweep crosses rho x n_var x sigma x r to give a full detection SURFACE
## (how the detectable departure depends on effect size sigma and phenotypic
## correlation r), not just one slice. pi and n_ind are held fixed.
set.seed(1)                       # reproducible draws
rho_grid    <- c(0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.97, 0.99, 0.995, 0.999, 0.9995)
n_var_grid  <- all_n_var          # c(5, 10, 20, 100)
sigma_grid  <- c(1, 0.1, 0.01)    # effect-size SD (matches the par_combo settings)
r_grid      <- c(0, 0.5)          # phenotypic correlation
SWEEP_TIMES <- 200                # replicates per (rho, n_var, sigma, r); raise to smooth
SWEEP_PI <- 0.5; SWEEP_NIND <- 1000   # fixed nuisance params

#' Fast drop-in for simulations.R::get_pheno_pair().
#'
#' The original draws the per-individual residual with one mvtnorm::rmvnorm call
#' PER INDIVIDUAL (apply over n_ind columns) -- the dominant cost of the sweep.
#' Since the 2x2 residual covariance R is identical for every individual, we draw
#' all n_ind residuals in a single call. Statistically identical; ~100-1000x faster.
#' @param b 2 x n_var true-effect matrix; @param X n_var x n_ind genotypes;
#' @param r phenotypic correlation. @return 2 x n_ind phenotype matrix.
get_pheno_pair_fast <- function(b, X, r) {
  R  <- matrix(c(1, r, r, 1), nrow = 2)
  MU <- b %*% X                                   # 2 x n_ind expected phenotypes
  MU + t(mvtnorm::rmvnorm(n = ncol(MU), sigma = R))  # add all residuals at once
}

#' Fast drop-in for simulations.R::pheno_corr_test().
#'
#' Identical logic/return to pheno_corr_test() (reuses get_ac_mat, get_af_mat,
#' get_geno_mat, get_true_beta, get_beta_hat, get_c_hat, get_mle_beta,
#' get_likelihood_test_stats) but generates phenotypes with get_pheno_pair_fast().
#' @return list(results = c(c_hat, lambda, pvalue), beta = per-variant data.frame).
pheno_corr_test_fast <- function(n_ind, n_var, c, r, pi, sigma, mle = TRUE, null = TRUE,
                                 mode = 'linear', rho = 0.5, sigma_eps = sigma) {
  AC <- get_ac_mat(n_var)
  A  <- get_af_mat(AC, n_ind)
  X  <- get_geno_mat(AC, n_ind)
  b  <- get_true_beta(n_var, c, pi, sigma, null = null, mode = mode, rho = rho, sigma_eps = sigma_eps)
  Y  <- get_pheno_pair_fast(b, X, r)
  b_hat  <- get_beta_hat(Y, X, A, n_ind)
  b1_hat <- matrix(b_hat[1, ], nrow = 1); b2_hat <- matrix(b_hat[2, ], nrow = 1)
  c_hat  <- get_c_hat(b1_hat, b2_hat, A, r)
  b_mle  <- get_mle_beta(b1_hat, b2_hat, c, r, null = null)
  df <- n_var - 1
  if (null) { c_hat <- if_else(mle, c_hat, c); df <- if_else(mle, df, n_var) }
  lambda <- get_likelihood_test_stats(n_ind, r, b1_hat, b2_hat, c_hat, A)
  pvalue <- pchisq(lambda, df, lower.tail = FALSE)
  beta <- data.frame(b1_true_value = b[1, ], b1_sum_stats = b_hat[1, ], b1_mle = b_mle[1, ],
                     b2_true_value = b[2, ], b2_sum_stats = b_hat[2, ], b2_mle = b_mle[2, ],
                     AF = c(diag(A)))
  list(results = c(c_hat = c_hat, lambda = lambda, pvalue = pvalue), beta = beta)
}

#' Run the rho-sweep simulation over rho x n_var x sigma x r.
#'
#' For each cell, simulate `times` ALLSPICE tests under the correlated-effects
#' alternative and record, per replicate, the p-value and the AF-weighted effect
#' correlation (observed and true).
#' @param rho_grid,n_var_grid,sigma_grid,r_grid grids to cross.
#' @param times replicates per cell.
#' @param n_ind,pi fixed simulation parameters (sample size, causal-variant proportion).
#' @return long tibble: one row per replicate (rho, n_var, sigma, r, pvalue, eff_corr_*_afw).
run_rho_sweep <- function(rho_grid, n_var_grid, sigma_grid, r_grid, times, n_ind, pi) {
  grid <- tidyr::expand_grid(rho = rho_grid, n_var = n_var_grid, sigma = sigma_grid, r = r_grid)
  purrr::pmap_dfr(grid, function(rho, n_var, sigma, r) {
    reps <- replicate(times,
      pheno_corr_test_fast(n_ind = n_ind, n_var = n_var, c = 0, r = r, pi = pi, sigma = sigma,
                           mle = TRUE, null = FALSE, mode = 'correlated', rho = rho),  # c unused when null=FALSE
      simplify = FALSE)
    tibble::tibble(
      rho = rho, n_var = n_var, sigma = sigma, r = r,
      pvalue            = vapply(reps, function(z) z$results[['pvalue']], numeric(1)),
      eff_corr_obs_afw  = vapply(reps, function(z) w_cos(z$beta$b1_sum_stats,  z$beta$b2_sum_stats,  z$beta$AF), numeric(1)),
      eff_corr_true_afw = vapply(reps, function(z) w_cos(z$beta$b1_true_value, z$beta$b2_true_value, z$beta$AF), numeric(1)))
  })
}

sweep <- run_rho_sweep(rho_grid, n_var_grid, sigma_grid, r_grid, SWEEP_TIMES, SWEEP_NIND, SWEEP_PI)
readr::write_csv(sweep, paste0(result_path, 'sim_rho_sweep_raw.csv'))

## Per-(rho, n_var, sigma, r) medians + power (fraction significant) at nominal
## and genome-wide thresholds.
sweep_summ <- sweep %>%
  dplyr::group_by(rho, n_var, sigma, r) %>%
  dplyr::summarize(
    n_rep = dplyr::n(),
    median_eff_corr_obs_afw = median(eff_corr_obs_afw, na.rm = TRUE),
    median_pvalue = median(pvalue, na.rm = TRUE),
    pvalue_q25 = quantile(pvalue, 0.25, na.rm = TRUE),
    pvalue_q75 = quantile(pvalue, 0.75, na.rm = TRUE),
    power_0.05 = mean(pvalue < 0.05,    na.rm = TRUE),
    power_gw   = mean(pvalue < 4.23e-6, na.rm = TRUE),
    .groups = 'drop')
readr::write_csv(sweep_summ, paste0(result_path, 'sim_rho_sweep_summary.csv'))

#' Interpolate the column `xcol` at which the median p hits `target_p`.
#' Same idea as `invert_p` but for the sweep, where the x-axis can be either the
#' true rho or the AF-weighted observed correlation.
#' @param d group data frame (needs median_pvalue and `xcol`).
#' @param xcol name of the x column to read off ("rho" or "median_eff_corr_obs_afw").
#' @param target_p target p-value.
#' @return scalar value of `xcol` at the crossing, or NA (no extrapolation).
invert_sweep <- function(d, xcol, target_p) {
  d <- d %>% dplyr::arrange(median_pvalue) %>%
    dplyr::filter(is.finite(median_pvalue), is.finite(.data[[xcol]]))
  if (nrow(d) < 2 || target_p < min(d$median_pvalue) || target_p > max(d$median_pvalue)) return(NA_real_)
  stats::approx(d$median_pvalue, d[[xcol]], xout = target_p, ties = mean)$y
}

#' Explain why a detectable-correlation crossing is or isn't found, so that NA
#' values in the table are interpretable rather than ambiguous.
#'   "detectable"          -> the median p crosses `target_p` within the swept
#'                            rho grid (the interpolated value is reported);
#'   "no_power"            -> the median p never drops to `target_p` even at
#'                            rho = 0 (maximal heterogeneity) -> no power in this
#'                            (sigma, r, n_var) regime;
#'   "always_significant" -> the median p is below `target_p` across the whole
#'                            grid (significant even at rho ~ 1) -> the crossing
#'                            lies beyond the swept rho range.
#' @return one of the three status strings.
sweep_status <- function(d, target_p) {
  mp <- d$median_pvalue[is.finite(d$median_pvalue)]
  if (length(mp) < 2) return("undefined")
  if (target_p < min(mp)) return("no_power")
  if (target_p > max(mp)) return("always_significant")
  "detectable"
}

## For each (sigma, r, n_var), the detectable true rho (and AF-weighted obs corr)
## at nominal and genome-wide thresholds. Larger rho = ALLSPICE detects subtler
## departures from proportionality (smaller 1 - rho). The status_* columns label
## why a rho_at_* entry is NA (no_power vs always_significant; see sweep_status).
detect_sweep <- sweep_summ %>% group_by(sigma, r, n_var) %>%
  group_modify(~ tibble::tibble(
    rho_at_p_0.05      = invert_sweep(.x, 'rho', 0.05),
    rho_at_p_gw        = invert_sweep(.x, 'rho', 4.23e-6),
    corr_afw_at_p_0.05 = invert_sweep(.x, 'median_eff_corr_obs_afw', 0.05),
    corr_afw_at_p_gw   = invert_sweep(.x, 'median_eff_corr_obs_afw', 4.23e-6),
    status_p_0.05      = sweep_status(.x, 0.05),
    status_p_gw        = sweep_status(.x, 4.23e-6))) %>%
  ungroup() %>% arrange(sigma, r, n_var)
readr::write_csv(detect_sweep, paste0(result_path, 'sim_rho_sweep_detectable_corr.csv'))
cat('\n-- rho-sweep: detectable true rho / AF-weighted corr at thresholds (per sigma, r, n_var) --\n')
print(detect_sweep, n = 100)

## ===========================================================================
## FIGURE S21 CAPTION (for the Supplementary Information)
## ---------------------------------------------------------------------------
## Figure S21. Simulation-based calibration of the heterogeneity magnitude that
## ALLSPICE can detect. Each panel plots the median ALLSPICE -log10(p) (over 200
## simulated replicates per cell) against the true cross-trait effect correlation
## rho, where the departure from proportionality is 1 - rho; one line per number
## of variants per gene (5, 10, 20, 100). Panels are arranged by effect-size
## standard deviation (rows: sigma = 1, 0.1, 0.01) and phenotypic correlation
## (columns: r = 0, 0.5). Dashed and dotted horizontal lines mark the nominal
## (p = 0.05) and study-wide (p = 4.23 x 10^-6) significance thresholds; the rho
## at which a curve crosses a threshold is the smallest departure from
## proportionality detectable in that regime. The y-axis is capped at 50 (top row
## only) for legibility, and triangles denote points beyond the cap (i.e. more
## significant). Simulations use n_ind = 1,000 individuals and causal proportion
## pi = 0.5, with effects drawn under the correlated model (cor(beta1, beta2) =
## rho exactly). The detectable departure shrinks as the number of variants grows;
## at sigma <= 0.1 the test has essentially no power at this sample size; and with
## 5 variants study-wide significance is not reached even for fully independent
## effects (rho = 0). These thresholds are specific to the simulated effect size
## and sample size and will shift with both.
## ===========================================================================

## Detection SURFACE: median -log10(p) vs true rho, one line per n_var, faceted by
## sigma (rows) x phenotypic r (cols). The rho where a line crosses a threshold is
## the minimum detectable correlation in that (sigma, r) regime.
## - y is capped at YCAP so the steep threshold-crossing region is legible rather
##   than dominated by p-values that underflow to 0 (which carry no real magnitude);
## - points whose true -log10(p) exceeds the cap are drawn as triangles ("off-scale,
##   more significant"), and their (meaningless) IQR error bars are suppressed.
YCAP <- 50
cap  <- function(p) pmin(-log10(pmax(p, PFLOOR)), YCAP)
sweep_plot <- sweep_summ %>%
  mutate(n_var     = factor(n_var),
         floored   = -log10(pmax(median_pvalue, PFLOOR)) > YCAP,   # off the cap
         y         = cap(median_pvalue),
         ymin      = cap(pvalue_q75),
         ymax      = cap(pvalue_q25),
         sigma_lab = factor(paste0('sigma = ', sigma), levels = paste0('sigma = ', sigma_grid)),
         r_lab     = paste0('phenotypic r = ', r))

## Reference lines drawn per-facet ONLY where the panel's data reaches them, so the
## genome-wide line (5.37) doesn't balloon the empty low-signal (low-sigma) panels.
panel_max <- sweep_plot %>% dplyr::group_by(sigma_lab, r_lab) %>%
  dplyr::summarize(ymax = max(y, na.rm = TRUE), .groups = 'drop')
ref_lines <- tidyr::expand_grid(
    panel_max,
    tibble::tibble(lvl  = c('nominal (0.05)', 'genome-wide (4.23e-6)'),
                   yint = c(-log10(0.05), -log10(4.23e-6)))) %>%
  dplyr::filter(yint <= ymax + 0.3)

p_sweep <- ggplot(sweep_plot, aes(rho, y, color = n_var, group = n_var)) +
  geom_errorbar(data = function(d) dplyr::filter(d, !floored),          # drop floored bars
                aes(ymin = ymin, ymax = ymax), width = 0.012, alpha = 0.4) +
  geom_line() +
  geom_point(aes(shape = floored), size = 1.9) +
  geom_hline(data = ref_lines, aes(yintercept = yint, linetype = lvl), colour = 'grey40') +
  facet_grid(sigma_lab ~ r_lab, scales = 'free_y') +
  scale_color_viridis_d(name = 'Number of variants', option = 'C', end = 0.9) +
  scale_linetype_manual(values = c(`nominal (0.05)` = 'dashed',
                                   `genome-wide (4.23e-6)` = 'dotted'), name = 'Threshold') +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = 'none') +
  labs(x = expression('True effect correlation '*rho*''),
       y = expression(-log[10](median~p_ALLSPICE))) +
  theme_classic(base_size = 12) + themes +
  theme(legend.position = 'top', strip.background = element_blank(),
        strip.text = element_text(face = 'bold'), plot.caption = element_text(hjust = 0))
png(paste0(figure_path, 'figureS21.png'), height = 7.5, width = 7.5, units = 'in', res = 300)
print(p_sweep)
dev.off()
cat('\nwrote:', paste0(figure_path, 'figureS21.png'), '\n')
