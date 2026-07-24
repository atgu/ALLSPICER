###############################################################################
## figureS25.R
##
## Figure S25 -- ALLSPICE null-calibration robustness (R1): type-I error under
## the exact proportional null as the phenotypic correlation r approaches the
## singular boundary (c^2 - 2cr + 1 -> 0) and is mis-specified (test given the
## wrong r). Shows the test stays calibrated when r is correct -- even at
## r = 0.9/0.99 -- and inflates only when r is over-estimated toward 1.
##
## Extracted from `R/simulation_analysis.R` (Part 5, R1). Self-contained:
## sourcing constants.R + simulations.R supplies the ALLSPICE helpers
## (get_ac_mat/get_af_mat/get_geno_mat/get_true_beta/get_beta_hat/get_c_hat/
## get_likelihood_test_stats) and `all_n_var`.
##
## OUTPUTS (to result_path / figure_path)
##   sim_null_calibration_r.csv  type-I error (0.05, 0.01) + lambda_GC per cell
##   figureS25.png               Figure S25 (type-I error vs r and mis-spec)
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
                            max_cnt = 100, noise = c('normal', 't'), df_t = 3, pi = 0.5) {
  noise <- match.arg(noise)
  AC <- get_ac_mat(n_var, max_cnt = max_cnt)
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

## --- R1: r near the boundary + r mis-specification -------------------------
r1_grid <- tidyr::expand_grid(
  r_gen = c(0, 0.5, 0.9, 0.99),
  delta = c(-0.1, -0.05, 0, 0.05, 0.1),     # r_test = r_gen + delta (mis-specification)
  n_var = all_n_var) %>%
  mutate(r_test = pmin(pmax(r_gen + delta, -0.99), 0.99))

calib_r <- purrr::pmap_dfr(r1_grid, function(r_gen, delta, n_var, r_test) {
  dplyr::bind_cols(
    tibble::tibble(r_gen = r_gen, delta = delta, n_var = n_var, r_test = r_test),
    calib_summary(CALIB_TIMES, n_ind = CALIB_NIND, n_var = n_var, c = CALIB_C,
                  r_gen = r_gen, r_test = r_test, sigma = CALIB_SIGMA,
                  max_cnt = 100, noise = 'normal', pi = CALIB_PI))
})
readr::write_csv(calib_r, paste0(result_path, 'sim_null_calibration_r.csv'))
cat('\n-- R1 null calibration: type-I error / lambda_GC vs r and mis-specification --\n')
print(calib_r, n = 100)

## ===========================================================================
## FIGURE S25 CAPTION (for the Supplementary Information)
## ---------------------------------------------------------------------------
## Figure S25. ALLSPICE type-I error under the proportional null as the phenotypic
## correlation approaches the singular boundary and is mis-specified. Data were
## simulated under the exact proportional null (β1 = c·β2 with c = 0.5; zero
## heterogeneity) at n_ind = 1,000 individuals, effect-size SD σ = 1, and mixture
## weight π = 0.5, with 1,000 replicates per cell. Each panel is a number of
## variants per gene (m = 5, 10, 20, 100). The x-axis is the generating phenotypic
## correlation r ∈ {0, 0.5, 0.9, 0.99} (approaching the boundary at which ALLSPICE's
## denominator c² − 2cr + 1 becomes singular); colour is the mis-specification
## δ = r_test − r_gen ∈ {−0.1, −0.05, 0, 0.05, 0.1} supplied to the test; the
## y-axis is the empirical type-I error at α = 0.05 (dashed line). With r correctly
## specified (δ = 0) the type-I error stays near nominal across all r, including the
## near-boundary values r = 0.9 and 0.99 — the near-singular denominator does not by
## itself break calibration. The test is anti-conservative only when r is
## over-estimated (δ > 0) and conservative when under-estimated (δ < 0), with the
## inflation increasing toward the boundary (largest at true r = 0.9 with δ = +0.1).
## ===========================================================================

p_calib_r <- ggplot(calib_r %>% mutate(n_var = factor(n_var), delta = factor(delta)),
                    aes(r_gen, type1_0.05, color = delta, group = delta)) +
  geom_line() + geom_point(size = 1.8) +
  geom_hline(yintercept = 0.05, lty = 2, colour = 'grey40') +
  facet_wrap(~ n_var, labeller = labeller(n_var = function(x) paste0('N variants = ', x))) +
  scale_color_brewer(palette = 'RdBu') +
  labs(x = 'Generating phenotypic correlation r', y = 'Type-I error at alpha = 0.05', color = expression(''*delta*'')) +
  theme_classic(base_size = 12) +
  theme(legend.position = 'right', strip.background = element_blank(),
        strip.text = element_text(face = 'bold'), plot.caption = element_text(hjust = 0))
png(paste0(figure_path, 'figureS25.png'), height = 6, width = 9, units = 'in', res = 300)
print(p_calib_r)
dev.off()
cat('\nwrote:', paste0(figure_path, 'figureS25.png'), '\n')
