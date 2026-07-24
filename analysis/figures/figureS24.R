###############################################################################
## figureS24.R
##
## Figure S24 -- magnitude (rho_w) of study-wide-significant vs non-significant
## gene-trait pairs WITHIN the same genes (a within-gene control). Significant
## pairs do not have systematically larger departures from proportionality than
## non-significant pairs of the same genes, i.e. significance is not enriched for
## trivial near-proportional signals.
##
## Extracted from `R/real_gene_rho.R` (the sig-vs-non-significant comparison).
## Self-contained: sourcing constants.R supplies read_pleiotropy_results/
## modify_results_table/annotation_types and the paths; read_gcs_fast is inlined
## below (needs gsutil).
##
## INPUT  gs://ukb-diverse-pops/wlu/genebass_notebooks/sig_gene_variant_betas_500k.txt.bgz
## OUTPUTS
##   ajhg_revision/data/real_gene_sig_vs_nonsig_magnitude.csv
##   ajhg_revision/data/real_gene_sig_vs_nonsig_by_gene.csv
##   figureS24.png  (Figure S24; per-gene faceted; main figure_path)
##   real_gene_sig_vs_nonsig_overview.png  (supporting; single-panel pooled overview, main figure_path)
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

RMAX <- 0.8   # restrict to phenotypically distinct pairs (corr < RMAX); see figureS23.R.

raw_results_500k <- read_pleiotropy_results('burden', '500k')

## Variant-level effect sizes that fed ALLSPICE, stream-filtered to the significant
## genes (long format: one row per (variant, phenotype)).
top_triplets <- read_gcs_fast('gs://ukb-diverse-pops/wlu/genebass_notebooks/sig_gene_variant_betas_500k.txt.bgz') %>%
  mutate(coding = if_else(is.na(coding), '', coding),
         modifier = if_else(is.na(modifier), '', modifier),
         phenoname = paste0(trait_type, '_', phenocode, '_', pheno_sex, '_', coding, '_', modifier))

## =====================================================================
## Comparison group: magnitude for NON-significant pairs of the SAME genes.
## Variant-level betas exist only for the genes in top_triplets (the significant
## set), so the non-significant comparison is necessarily WITHIN those genes
## (other phenotype-pairs that did not reach study-wide significance). This is a
## within-gene control (variant count comparable); it cannot sample genome-wide
## non-significant pairs.
## =====================================================================
genes_with_data <- unique(top_triplets$gene)

compute_rho_one <- function(gene, pheno1, pheno2) {
  v <- top_triplets %>%
    filter(gene == !!gene, AF < 1e-4, phenoname %in% c(pheno1, pheno2)) %>%
    mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
    filter(annotation == 'missense|LC') %>%
    pivot_wider(names_from = phenoname, values_from = BETA,
                id_cols = c(locus, alleles, AF), values_fn = mean)
  ok <- all(c(pheno1, pheno2) %in% names(v))
  if (ok) v <- v %>% filter(!is.na(.data[[pheno1]]), !is.na(.data[[pheno2]]))
  c(n = if (ok) nrow(v) else 0, rho = if (ok && nrow(v) >= 2) w_cos(v[[pheno1]], v[[pheno2]], v$AF) else NA_real_)
}

cand <- modify_results_table(raw_results_500k, 'burden', '500k') %>%
  mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
  filter(annotation == 'missense|LC', gene %in% genes_with_data, corr < RMAX) %>%
  distinct(gene, pheno1, pheno2, allspice_pvalue = pvalue, n_var_test = n_var) %>%
  mutate(significant = allspice_pvalue < 4.23e-6)

rw <- purrr::pmap_dfr(cand %>% select(gene, pheno1, pheno2), function(gene, pheno1, pheno2) {
  o <- compute_rho_one(gene, pheno1, pheno2)
  tibble::tibble(gene = gene, pheno1 = pheno1, pheno2 = pheno2,
                 n_var_used = as.integer(o[['n']]), rho_w = o[['rho']])
})
cmp <- cand %>% left_join(rw, by = c('gene', 'pheno1', 'pheno2')) %>%
  filter(!is.na(rho_w)) %>%
  mutate(departure = 1 - abs(rho_w),
         group = if_else(significant, 'study-wide significant', 'non-significant'))

cmp_out <- '~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/ajhg_revision/data/real_gene_sig_vs_nonsig_magnitude.csv'
readr::write_csv(cmp, cmp_out)
cat('\nwrote:', cmp_out, '\n')

cat('\n-- magnitude: study-wide significant vs non-significant pairs (within the', length(genes_with_data), 'genes with variant data) --\n')
print(cmp %>% group_by(group) %>%
        summarize(n = dplyr::n(),
                  median_departure = round(median(departure), 3),
                  median_nvar = median(n_var_test), .groups = 'drop'))
wt <- suppressWarnings(wilcox.test(departure ~ group, data = cmp))
cat(sprintf('Wilcoxon departure (sig vs non-sig): p = %.3g\n', wt$p.value))
cat('Spearman(-log10 p, departure) across ALL', nrow(cmp), 'pairs =',
    round(suppressWarnings(cor(-log10(pmax(cmp$allspice_pvalue, 1e-300)), cmp$departure, method = 'spearman')), 3), '\n')

## ===========================================================================
## FIGURE S24 CAPTION (for the Supplementary Information)
## ---------------------------------------------------------------------------
## Figure S24. Magnitude of departure from proportionality (rho_w) for study-wide
## significant versus non-significant gene-trait pairs, shown per gene. Each panel is
## one of the 12 genes that carry a study-wide-significant missense|LC pair (the genes
## with variant-level data); within a panel, every missense|LC phenotype pair with
## phenotypic correlation r < 0.8 is plotted by its AF-weighted effect-size correlation
## rho_w (x) and -log10(ALLSPICE p) (y), coloured by whether it reached study-wide
## significance (p < 4.23 x 10^-6, red) or not (grey) and sized by the number of
## variants (log scale). The y-axis is capped at -log10(p) = 30; points exceeding it
## (p < 10^-30) are drawn as triangles, and the dashed line marks the study-wide
## threshold. This is a within-gene control: within each gene the significant and
## non-significant pairs span a comparable range of rho_w, so significance is not
## enriched for larger departures from proportionality (overall, their departures
## 1 - |rho_w| do not differ significantly; Wilcoxon test in the text) — consistent
## with rho_w (magnitude) and the p-value (a power-dependent detection statistic)
## being complementary. The comparison is necessarily within these genes because
## variant-level betas are available only for them.
## ===========================================================================

p_cmp <- ggplot(cmp, aes(rho_w, -log10(pmax(allspice_pvalue, 1e-300)),
                         colour = group, size = n_var_test)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(4.23e-6), lty = 2, colour = 'grey40') +
  scale_size_continuous(name = 'Number of\nvariants', trans = 'log10') +
  scale_colour_manual(values = c('study-wide significant' = '#D62728',
                                 'non-significant' = 'grey55'), name = NULL) +
  labs(x = expression('AF-weighted effect correlation '*rho[w]*'  ( = 1 iff exactly proportional )'),
       y = expression(-log[10]('ALLSPICE p-value')),
       title = 'Significant vs non-significant pairs (within the same genes)',
       caption = sprintf('n = %d significant, %d non-significant; dashed = genome-wide 4.23e-6.',
                         sum(cmp$significant), sum(!cmp$significant))) +
  theme_classic(base_size = 12) + theme(plot.caption = element_text(hjust = 0))

## single-panel overview (pooled across genes), kept as a supporting figure;
## figureS24.png itself is the per-gene faceted version written below.
cmp_fig <- paste0(figure_path, 'real_gene_sig_vs_nonsig_overview.png')
png(path.expand(cmp_fig), height = 5, width = 8, units = 'in', res = 300)
print(p_cmp)
dev.off()
cat('wrote:', cmp_fig, '\n')

## --- Per-gene breakdown (within-gene control; variant count comparable within a gene) ---
YCAP <- 30
per_gene <- cmp %>%
  group_by(gene, group) %>%
  summarize(n_pairs = dplyr::n(),
            median_rho_w = round(median(rho_w), 3),
            median_departure = round(median(departure), 3),
            median_nvar = median(n_var_test), .groups = 'drop') %>%
  arrange(gene, group)
pg_out <- '~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/ajhg_revision/data/real_gene_sig_vs_nonsig_by_gene.csv'
readr::write_csv(per_gene, pg_out)
cat('\nwrote:', pg_out, '\n\n-- per-gene: significant vs non-significant magnitude --\n')
print(as.data.frame(per_gene), row.names = FALSE)

p_pergene <- ggplot(cmp %>%
                      mutate(floored   = -log10(pmax(allspice_pvalue, PFLOOR)) > YCAP), aes(rho_w, -log10(pmax(allspice_pvalue, 1e-30)),
                             colour = group, size = n_var_test)) +
  geom_point(aes(shape = floored), size = 1.9) +
  geom_hline(yintercept = -log10(4.23e-6), lty = 2, colour = 'grey40') +
  facet_wrap(~ gene, scales = 'free', nrow = 3) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = 'none') +
  scale_colour_manual(values = c('study-wide significant' = '#D62728',
                                 'non-significant' = 'grey55'), name = NULL) +
  scale_size_continuous(name = 'Number of variants', trans = 'log10') +
  labs(x = expression('AF-weighted effect correlation '*rho[w]),
       y = expression(-log[10]('p_ALLSPICE'))) +
  theme_classic(base_size = 10) +
  theme(legend.position = 'top', strip.background = element_blank(),
        strip.text = element_text(face = 'bold'))
fig_out <- paste0(figure_path, 'figureS24.png')
png(path.expand(fig_out), height = 5, width = 7.5, units = 'in', res = 300)
print(p_pergene)
dev.off()
cat('wrote:', fig_out, '\n')
