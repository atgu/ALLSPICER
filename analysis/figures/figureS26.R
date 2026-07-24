###############################################################################
## figureS26.R
##
## Figure S26 -- ALLSPICE null-calibration robustness (R2): type-I error under
## the exact proportional null for ultra-rare variants (singletons/doubletons ->
## strongly non-normal effect estimates) and heavy-tailed (t) residual noise,
## with the phenotypic correlation correctly specified. Shows ultra-rare variants
## are conservative and heavy-tailed noise gives only mild inflation at high
## variant counts.
##
## Extracted from `R/simulation_analysis.R` (Part 5, R2). Self-contained:
## sourcing constants.R + simulations.R supplies the ALLSPICE helpers
## (get_ac_mat/get_af_mat/get_geno_mat/get_true_beta/get_beta_hat/get_c_hat/
## get_likelihood_test_stats) and `all_n_var`.
##
## OUTPUTS (to result_path / figure_path)
##   sim_null_calibration_ultrarare.csv  type-I error (0.05, 0.01) + lambda_GC per cell
##   figureS26.png                       Figure S26 (type-I error, ultra-rare / non-normal)
###############################################################################

## ---------------------------------------------------------------------------
## Setup
## ---------------------------------------------------------------------------
source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/constants.R')
source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/simulations.R')  # get_*, all_n_var
## figure_path and result_path come from constants.R (the main figure directory).
library(ggpubr)

## Re-attach the tidyverse LAST so its verbs win on the search path (the sourced
## files pull in plyr/MASS, which mask dplyr::summarize/select).
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(tibble)
  library(ggplot2)
})

## ---------------------------------------------------------------------------
## Part 5 (shared): proportional-null calibration helpers
## ---------------------------------------------------------------------------
## Data are generated under the exact proportional null (b1 = c*b2, zero
## heterogeneity); we measure ALLSPICE type-I error and the genomic-inflation
## factor lambda_GC when the model assumptions are violated. If type-I error
## stays nominal and lambda_GC ~ 1, a significant ALLSPICE p cannot be an
## artifact of those violations.
set.seed(2)
CALIB_TIMES <- 1000               # replicates per cell (type-I at 0.05 + lambda_GC)
CALIB_C  <- 0.5                   # proportionality constant under the null
CALIB_SIGMA <- 1; CALIB_PI <- 0.5; CALIB_NIND <- 1000

#' One proportional-null ALLSPICE test, with controllable misspecification.
#'
#' Data are generated under the exact null b1 = c*b2 (no heterogeneity). The test
#' may be given a phenotypic correlation r_test != the generating r_gen, the
#' allele-count ceiling can be lowered to force ultra-rare variants, and the
#' residual noise can be heavy-tailed (multivariate t, scaled to unit variance).
#' @return c(lambda, pvalue, df).
null_calib_test <- function(n_ind, n_var, c, r_gen, r_test, sigma,
                            min_cnt = 1, max_cnt = 100, noise = c('normal', 't'), df_t = 3, pi = 0.5) {
  noise <- match.arg(noise)
  ## allele counts drawn uniformly in [min_cnt, max_cnt) (cf. get_ac_mat, which
  ## fixes min_cnt = 1); min_cnt lets us carve disjoint AC bins (e.g. AC 3-99).
  AC <- diag(floor(runif(n = n_var, min = min_cnt, max = max_cnt)), ncol = n_var, nrow = n_var)
  A  <- get_af_mat(AC, n_ind)
  X  <- get_geno_mat(AC, n_ind)
  b  <- get_true_beta(n_var, c, pi, sigma, null = TRUE)     # exact proportional null
  R  <- matrix(c(1, r_gen, r_gen, 1), nrow = 2)
  MU <- b %*% X
  resid <- if (noise == 'normal') {
    t(mvtnorm::rmvnorm(n = ncol(MU), sigma = R))
  } else {
    t(mvtnorm::rmvt(n = ncol(MU), sigma = R * (df_t - 2) / df_t, df = df_t))  # cov == R
  }
  Y <- MU + resid
  b_hat  <- get_beta_hat(Y, X, A, n_ind)
  b1_hat <- matrix(b_hat[1, ], nrow = 1); b2_hat <- matrix(b_hat[2, ], nrow = 1)
  c_hat  <- get_c_hat(b1_hat, b2_hat, A, r_test)            # test uses r_test (maybe wrong)
  lambda <- get_likelihood_test_stats(n_ind, r_test, b1_hat, b2_hat, c_hat, A)
  df <- n_var - 1                                           # MLE c_hat under the null
  c(lambda = lambda, pvalue = pchisq(lambda, df, lower.tail = FALSE), df = df)
}

#' Replicate null_calib_test and summarize calibration for one cell.
#' @return tibble: type-I error at 0.05 / 0.01 and genomic inflation lambda_GC
#'   (= median(lambda) / qchisq(0.5, df); ~1 if well-calibrated).
calib_summary <- function(times, ...) {
  args <- list(...)
  out  <- vapply(seq_len(times), function(i) do.call(null_calib_test, args), numeric(3))
  pv   <- out['pvalue', ]; lam <- out['lambda', ]; df <- out['df', 1]
  ok  <- is.finite(pv)
  tibble::tibble(
    n_valid    = sum(ok),
    type1_0.05 = mean(pv[ok] < 0.05),
    type1_0.01 = mean(pv[ok] < 0.01),
    lambda_gc  = median(lam[ok], na.rm = TRUE) / qchisq(0.5, df))
}

## --- R2: ultra-rare variants + heavy-tailed (non-normal) noise -------------
## Three allele-count regimes:
##   AC 1-2    ultra-rare only (singletons/doubletons; non-normal effect estimates)
##   AC 3-99   the complementary, NON-ultra-rare rare/uncommon range (disjoint)
##   AC 1-99   the full baseline spectrum used in every other simulation / the real
##             AF < 1e-4 analysis (includes the ultra-rare end)
ac_levels <- c('ultra-rare (AC 1-2)', 'moderately rare (AC 3-99)', 'baseline (AC 1-99)')
r2_grid <- tidyr::expand_grid(
  ac     = ac_levels,
  noise  = c('normal', 't (df=3)'),
  n_var  = all_n_var) %>%
  mutate(min_cnt   = dplyr::case_when(ac == 'moderately rare (AC 3-99)' ~ 3, TRUE ~ 1),
         max_cnt   = dplyr::case_when(ac == 'ultra-rare (AC 1-2)'       ~ 3, TRUE ~ 100),
         noise_arg = ifelse(noise == 'normal', 'normal', 't'))

calib_rare <- purrr::pmap_dfr(r2_grid, function(ac, noise, n_var, min_cnt, max_cnt, noise_arg) {
  dplyr::bind_cols(
    tibble::tibble(ac = ac, noise = noise, n_var = n_var),
    calib_summary(CALIB_TIMES, n_ind = CALIB_NIND, n_var = n_var, c = CALIB_C,
                  r_gen = 0.5, r_test = 0.5, sigma = CALIB_SIGMA,    # r correctly specified
                  min_cnt = min_cnt, max_cnt = max_cnt, noise = noise_arg, df_t = 3, pi = CALIB_PI))
})
readr::write_csv(calib_rare, paste0(result_path, 'sim_null_calibration_ultrarare.csv'))
cat('\n-- R2 null calibration: type-I error / lambda_GC for ultra-rare / non-normal --\n')
print(calib_rare, n = 100)

## ===========================================================================
## FIGURE S26 CAPTION (for the Supplementary Information)
## ---------------------------------------------------------------------------
## Figure S26. ALLSPICE type-I error under the proportional null across allele-count
## regimes and noise distributions, with the phenotypic correlation correctly
## specified. Data were simulated under the exact proportional null (β1 = c·β2 with
## c = 0.5; residual correlation both generated and supplied to the test at r = 0.5)
## at n_ind = 1,000 individuals, effect-size SD σ = 1, and mixture weight π = 0.5,
## with 1,000 replicates per cell. Bars show the empirical type-I error at α = 0.05
## (dashed line) by number of variants per gene (m = 5, 10, 20, 100). Panels contrast
## three allele-count regimes: ultra-rare only (AC 1–2, singletons/doubletons, for
## which the per-variant effect estimate is driven by ~one carrier and is strongly
## non-normal), the complementary non-ultra-rare range (AC 3–99), and the full
## baseline spectrum (AC 1–99) used in every other simulation and matching the real
## AF < 10⁻⁴ analysis. Colour contrasts Gaussian residual noise with heavy-tailed
## multivariate-t (df = 3) noise scaled to unit variance. Both the ultra-rare and the
## moderately-rare regimes are conservative (type-I at or below nominal) rather than
## anti-conservative; heavy-tailed noise produces only mild inflation that grows with
## variant count and is attenuated, not amplified, in the ultra-rare regime.
## ===========================================================================

p_calib_rare <- ggplot(calib_rare %>% mutate(n_var = factor(n_var),
                                             ac = factor(ac, levels = ac_levels)),
                       aes(n_var, type1_0.05, fill = noise)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.9) +
  geom_hline(yintercept = 0.05, lty = 2, colour = 'grey40') +
  facet_wrap(~ ac) +
  scale_fill_brewer(palette = 'Set2', name = 'Residual noise') +
  labs(x = 'Number of variants', y = 'Type-I error at alpha = 0.05') +
  theme_classic(base_size = 12) +
  theme(legend.position = 'top', strip.background = element_blank(),
        strip.text = element_text(face = 'bold'), plot.caption = element_text(hjust = 0))
png(paste0(figure_path, 'figureS26.png'), height = 5, width = 7.5, units = 'in', res = 300)
print(p_calib_rare)
dev.off()
cat('\nwrote:', paste0(figure_path, 'figureS26.png'), '\n')
