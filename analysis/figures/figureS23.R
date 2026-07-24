###############################################################################
## figureS23.R
##
## Figure S23 -- AF-weighted effect-size correlation (rho_w) vs ALLSPICE
## significance for the real study-wide-significant gene-trait pairs. Each
## significant missense|LC pair (phenotypic correlation r < 0.8) is plotted by its
## rho_w (magnitude of departure from proportionality; 1 - |rho_w|) against its
## ALLSPICE -log10(p), showing the significant hits span a wide range of rho_w
## (they are not near-proportional).
##
## Extracted from `R/real_gene_rho.R` (the Fig S{real} portion). Self-contained:
## sourcing constants.R supplies read_pleiotropy_results/modify_results_table/
## annotation_types and the paths; read_gcs_fast is inlined below (needs gsutil).
##
## INPUT  gs://ukb-diverse-pops/wlu/genebass_notebooks/sig_gene_variant_betas_500k.txt.bgz
##        (per-variant betas for the significant genes; see subset_sig_gene_variants
##        / subset_pleiotropy_data.py for how it was produced)
## OUTPUTS
##   ajhg_revision/data/real_gene_effect_corr_magnitude.csv  (Table S{realgene})
##   figureS23.png  (Figure S23; main figure_path)
###############################################################################

## ---------------------------------------------------------------------------
## Setup
## ---------------------------------------------------------------------------
source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/constants.R')
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(tibble); library(ggplot2)
})

PFLOOR <- 1e-300   # floor for p-values that underflow to 0, so -log10 stays finite

#' Fast GCS reader: gsutil -> tempfile -> read -> cleanup (inlined from the
#' aou_gwas constants so this script does not depend on that file).
read_gcs_fast <- function(gcs_uri, format = NULL, ..., gsutil = NULL) {
  if (is.null(gsutil)) {
    candidates <- c("/Users/wlu/google-cloud-sdk/bin/gsutil", Sys.which("gsutil"))
    candidates <- unique(candidates[nzchar(candidates)])
    gsutil <- candidates[file.exists(candidates)][1]
  }
  if (is.na(gsutil) || length(gsutil) == 0 || !nzchar(gsutil) || !file.exists(gsutil)) {
    stop("No usable gsutil found.")
  }
  uri_path <- sub("^gs://[^/]+/", "", gcs_uri)
  ext <- tools::file_ext(uri_path)
  base_ext <- if (tolower(ext) == "gz") tools::file_ext(tools::file_path_sans_ext(uri_path)) else ext
  fmt <- tolower(if (!is.null(format)) format else if (tolower(ext) == "gz") base_ext else ext)
  tmp <- tempfile(fileext = if (nzchar(ext)) paste0(".", ext) else "")
  on.exit(unlink(tmp), add = TRUE)
  res <- tryCatch(system2(gsutil, c("cp", gcs_uri, tmp), stdout = TRUE, stderr = TRUE),
                  error = function(e) e)
  status <- attr(res, "status")
  if (inherits(res, "error") || (!is.null(status) && status != 0) || !file.exists(tmp)) {
    stop(paste(c("gsutil copy failed.", if (length(res)) res else NULL), collapse = "\n"))
  }
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Need data.table")
  data.table::fread(tmp, ...)
}

w_cos <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w)
  x <- x[ok]; y <- y[ok]; w <- w[ok]
  d <- sqrt(sum(w * x^2) * sum(w * y^2))
  if (d == 0) NA_real_ else sum(w * x * y) / d
}

RMAX <- 0.8   # restrict to pairs with phenotypic correlation corr < RMAX.
              # Rationale (a priori, not calibration-driven): (1) ALLSPICE is meant
              # to compare two DIFFERENT phenotypes -- near-perfectly correlated
              # traits (corr -> 1) are effectively the same phenotype; (2) the
              # likelihood denominator c^2 - 2*r*c + 1 degenerates as r -> 1.
              # One-sided (signed) bound: high positive correlation only.

raw_results_500k <- read_pleiotropy_results('burden', '500k')
results_500k <- modify_results_table(raw_results_500k, 'burden', '500k') %>%
  mutate(annotation = factor(annotation, levels = annotation_types)) %>%
  filter(pvalue < 4.23e-6 & annotation == 'missense|LC' & corr < RMAX)

## Variant-level effect sizes that actually fed ALLSPICE, stream-filtered to the
## significant genes (long format: one row per (variant, phenotype); pheno
## identity rebuilt as phenoname below).
top_triplets <- read_gcs_fast('gs://ukb-diverse-pops/wlu/genebass_notebooks/sig_gene_variant_betas_500k.txt.bgz') %>%
  mutate(coding = if_else(is.na(coding), '', coding),
         modifier = if_else(is.na(modifier), '', modifier),
         phenoname = paste0(trait_type, '_', phenocode, '_', pheno_sex, '_', coding, '_', modifier))

sig_pairs <- results_500k %>% distinct(gene, pheno1, pheno2, description1, description2)
cat('strictly-significant missense|LC pairs:', nrow(sig_pairs), '\n')

real_corr <- purrr::pmap_dfr(sig_pairs, function(gene, pheno1, pheno2, description1, description2) {
  v <- top_triplets %>%
    filter(gene == !!gene, AF < 1e-4, phenoname %in% c(pheno1, pheno2)) %>%
    mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
    filter(annotation == 'missense|LC') %>%
    pivot_wider(names_from = phenoname, values_from = BETA,
                id_cols = c(locus, alleles, AF), values_fn = mean)
  ok <- all(c(pheno1, pheno2) %in% names(v))
  if (ok) v <- v %>% filter(!is.na(.data[[pheno1]]), !is.na(.data[[pheno2]]))
  tibble::tibble(gene = gene, pheno1 = pheno1, pheno2 = pheno2,
                 description1 = description1, description2 = description2,
                 n_var_used = if (ok) nrow(v) else 0L,
                 rho_w = if (ok && nrow(v) >= 2) w_cos(v[[pheno1]], v[[pheno2]], v$AF) else NA_real_)
}) %>%
  mutate(departure_from_proportionality = 1 - abs(rho_w)) %>%
  left_join(results_500k %>%
              distinct(gene, pheno1, pheno2, allspice_pvalue = pvalue, c_hat, n_var_test = n_var),
            by = c('gene', 'pheno1', 'pheno2')) %>%
  arrange(allspice_pvalue)

out <- '~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/ajhg_revision/data/real_gene_effect_corr_magnitude.csv'
readr::write_csv(real_corr, out)
cat('wrote:', out, '\n\n')
print(as.data.frame(real_corr %>%
  select(gene, description1, description2, n_var_used, rho_w,
         departure_from_proportionality, allspice_pvalue)), row.names = FALSE)

## --- Magnitude distribution across the significant pairs (numbers cited in the
##     Supplementary Note / response: median departure, fraction > 0.5, smallest |rho_w|) ---
a <- abs(real_corr$rho_w); a <- a[is.finite(a)]; dep <- 1 - a
cat(sprintf('\n-- significant-pair magnitude summary (N = %d) --\n', length(a)))
cat(sprintf('  |rho_w|: min %.3f, median %.3f, max %.3f\n', min(a), median(a), max(a)))
cat(sprintf('  departure 1-|rho_w|: median %.3f (IQR %.3f-%.3f)\n',
            median(dep), quantile(dep, .25), quantile(dep, .75)))
cat(sprintf('  departure > 0.5: %d/%d (%.0f%%); |rho_w| < 0.1: %d\n',
            sum(dep > 0.5), length(a), 100 * mean(dep > 0.5), sum(a < 0.1)))
cat(sprintf('  smallest |rho_w|: %s\n',
            paste(round(sort(a)[seq_len(min(6, length(a)))], 3), collapse = ', ')))

## ===========================================================================
## FIGURE S23 CAPTION (for the Supplementary Information)
## ---------------------------------------------------------------------------
## Figure S23. AF-weighted effect-size correlation (rho_w) versus ALLSPICE
## significance for the study-wide-significant gene-trait pairs. Each point is one
## study-wide-significant (p < 4.23 x 10^-6) burden missense|LC gene-trait pair
## with phenotypic correlation r < 0.8 (N = 49 pairs across 12 genes). rho_w is the
## allele-frequency-weighted, uncentered correlation between the two phenotypes'
## rare-variant (AF < 10^-4) effect estimates (= +-1 iff the effects are exactly
## proportional; the departure from proportionality is 1 - |rho_w|); the y-axis is
## -log10(ALLSPICE p) and point colour is the number of variants (log scale). The
## y-axis is capped at -log10(p) = 30; points exceeding it (p < 10^-30) are drawn as
## triangles. The dashed line marks the study-wide significance threshold
## (4.23 x 10^-6). The significant pairs span a wide range of rho_w (|rho_w| from
## 0.003 to 0.81, with
## 43 of 49 having a departure 1 - |rho_w| > 0.5), i.e. they are not
## near-proportional. Because rho_w is computed from estimated effects, sampling
## noise attenuates it toward 0, so it is a conservative measure of proportionality.
## ===========================================================================

## --- Supporting panel: significance is set by variant count, not by departure ---
YCAP <- 30
pdat <- real_corr %>% filter(!is.na(rho_w)) %>%
  mutate(floored   = -log10(pmax(allspice_pvalue, PFLOOR)) > YCAP)
sp <- suppressWarnings(c(
  p_dep   = cor(-log10(pmax(pdat$allspice_pvalue, 1e-300)), 1 - abs(pdat$rho_w), method = 'spearman'),
  p_nvar  = cor(-log10(pmax(pdat$allspice_pvalue, 1e-300)), pdat$n_var_test,     method = 'spearman')))
cat(sprintf('\nSpearman(-log10 p, departure) = %.2f ; Spearman(-log10 p, n_var) = %.2f\n',
            sp['p_dep'], sp['p_nvar']))

p_real <- ggplot(pdat, aes(rho_w, -log10(pmax(allspice_pvalue, 1e-30)), color = n_var_test)) +
  geom_point(aes(shape = floored), size = 1.9) +
  geom_hline(yintercept = -log10(4.23e-6), lty = 2, colour = 'grey40') +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = 'none') +
  scale_color_continuous(name = 'Number of variants', trans = 'log10') +
  labs(x = expression('AF-weighted effect correlation '*rho[w]*'  ( = 1 iff exactly proportional )'),
       y = expression(-log[10]('p_ALLSPICE'))) +
  theme_classic(base_size = 12) + theme(plot.caption = element_text(hjust = 0), legend.position = 'top')

fig_out <- paste0(figure_path, 'figureS23.png')
png(path.expand(fig_out), height = 5, width = 7.5, units = 'in', res = 300)
print(p_real)
dev.off()
cat('wrote:', fig_out, '\n')
