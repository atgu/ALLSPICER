###############################################################################
## figureS27.R
##
## Figure S27 + Supplementary Table S27 -- ALLSPICE synonymous negative-control /
## single-variant inflation (reviewer Comment 3). (A) Trait-specific synonymous
## single-variant genomic-control lambda_GC (median [IQR] across traits) and the
## ALLSPICE-level synonymous lambda_GC (diagnostic) vs AF ceiling. (B) Real (missense|LC,
## pLoF) triplet -log10(p) at the primary AF<1e-4 ceiling, before vs after applying the
## trait-specific synonymous lambda_GC correction (n_ind -> n_ind / lambda_GC).
##
## SELF-CONTAINED and the SOLE owner of this analysis (supersedes the former
## R/figureS27.R). It computes the synonymous inflation + re-runs every real triplet
## across AF ceilings (cached to results/modified_test/), writes Supplementary
## Table S27, and renders Figure S27. Re-running reuses the cached per-triplet scan
## (figureS27_scan_real_triplets.csv); delete it to force a full refresh.
##
## INPUTS  : hybrid sweep per-variant betas (GCS) + published results CSV (n_ind
##           reconstruction) + the ALLSPICER package.
## OUTPUTS : figureS27.{png,pdf} -> figure_path ;
##           figureS27_table.csv (= Table S27, inflation + variants-remaining + strict
##           counts), figureS27_scan_real_triplets.csv -> results/modified_test/
##
## ---------------------------------------------------------------------------
## MANUSCRIPT TEXT (reviewer Comment 3) -- incorporated from the former
## figureS27_text.md; source data figureS27_{table,scan_real_triplets}.csv.
## ---------------------------------------------------------------------------
## ANALYSIS. The ALLSPICE statistic assumes a per-variant sampling variance
## var(beta_hat_i) = 1/(n_ind * 2*AF_i) -- the "1/n" term. Under-estimation of this
## variance (e.g. from residual population structure) would inflate the test. We estimate
## the empirical single-variant inflation from SYNONYMOUS variants as a TRAIT-SPECIFIC
## genomic-control factor lambda_GC (median synonymous single-variant chi^2, df = 1, per
## trait, / qchisq(0.5,1)), substitute it for the theoretical scale (n_ind -> n_ind /
## lambda_GC, floored at 1 so genomic control only deflates), and re-run ALLSPICE on the
## real (missense|LC, pLoF) triplets to ask whether the heterogeneity signals survive.
## Repeated across AF ceilings.
##
## FIGURE S27 (legend). Heterogeneity signals re-tested after a trait-specific synonymous
## genomic-control correction. (A) Negative-control genomic control versus allele-frequency
## (AF) ceiling. The trait-specific SINGLE-VARIANT lambda_GC (median-based genomic control of
## the synonymous single-variant chi^2, df = 1, computed per trait) is mild and rises with
## the ceiling (median across traits = 1.02 / 1.05 / 1.07 / 1.08 at AF < 1e-5 / 1e-4 / 1e-3 /
## 1e-2; error bars = IQR across traits) -- common variants carry more inflation, as expected
## for residual structure. As a direct test-level check, the ALLSPICE statistic run on
## synonymous gene-trait PAIRS (genomic-controlled on the p-value scale) is CONSERVATIVE
## (lambda_gc_allspice = 0.03 / 0.06 / 0.07 / 0.38, all < 1; median synonymous ALLSPICE p ~
## 0.9), i.e. no test-level inflation. Dashed line: no inflation (lambda_GC = 1); values
## below it indicate conservatism. (B) ALLSPICE -log10(p) for the real (missense|LC, pLoF)
## gene-trait pairs at the primary AF < 1e-4 ceiling, before (x) versus after (y) the
## trait-specific synonymous correction (for a pair the two traits' lambda_GC are averaged
## and floored at 1; n_ind -> n_ind / lambda_GC). Dashed line: y = x; dotted red lines:
## strict threshold p = 4.23e-6. The strongest signals remain strictly significant; a
## fraction of marginal calls fall below threshold, reflecting the test's sensitivity to
## small per-variant inflation when aggregating over many variants.
##
## TABLE S27 (caption). Per AF ceiling: the number of traits (n_traits) and synonymous
## single-variant tests (n_syn_tests) used for the trait-specific genomic control, and the
## distribution of the trait-specific synonymous lambda_GC across traits (median
## lambda_trait_median; IQR lambda_trait_q25-lambda_trait_q75); the ALLSPICE-level
## synonymous genomic control (lambda_gc_allspice), which is < 1 (the test is conservative
## on the matched negative control) and is a calibration diagnostic, not applied as a
## correction; the number of variants entering ALLSPICE per gene-trait pair (n_var_median
## with IQR n_var_q25-n_var_q75); and, for missense|LC and pLoF gene-trait pairs, the number
## of strictly significant findings (p < 4.23e-6) before (n_strict_raw) and after
## (n_strict_traitspec_corr) the trait-specific synonymous correction.
##
## METHODS (add to Supplementary Methods). To assess sensitivity to residual single-variant
## inflation, including residual population structure, we used synonymous variants as a
## negative control. The ALLSPICE likelihood assumes each effect-size estimate has sampling
## variance 1/(n_ind * 2*AF_i); under-estimation of this variance inflates the test. For
## each phenotype we computed a trait-specific synonymous genomic-control factor
## lambda_GC = median(synonymous single-variant chi^2, df = 1) / qchisq(0.5, 1) over
## synonymous variants below each AF ceiling (1e-5, 1e-4 [primary], 1e-3, 1e-2; traits with
## >= 10 synonymous tests). For each real gene-trait pair we averaged the two traits'
## lambda_GC, floored it at 1 (so the correction only deflates), substituted it for the
## theoretical scale (n_ind -> n_ind / lambda_GC), and re-ran ALLSPICE, comparing the number
## of strictly significant findings (p < 4.23e-6) before and after. As an additional, more
## direct check we also ran ALLSPICE on the synonymous gene-trait pairs themselves and
## genomic-controlled the ALLSPICE statistic on the p-value scale (each synonymous ALLSPICE
## p -> 1-df chi^2; lambda_gc_allspice = median / qchisq(0.5, 1)); this came out < 1 (the
## test is conservative on the matched null), so it is reported as a diagnostic rather than
## applied as a correction.
##
## RESULTS (add to main text or Supplement). Synonymous inflation was mild at rare-variant
## thresholds. At the primary AF < 1e-4 threshold the median trait-specific synonymous
## genomic-control factor was 1.05 (IQR 0.98-1.12), based on 48,359 synonymous single-variant
## association tests across 85 traits (Fig S27A), rising modestly with the ceiling (median
## lambda_GC = 1.02 / 1.05 / 1.07 / 1.08 at AF < 1e-5 / 1e-4 / 1e-3 / 1e-2) -- common variants
## carry more inflation, as expected for residual structure, reinforcing AF < 1e-4 as the
## primary threshold. After applying the trait-specific synonymous correction, strictly
## significant pLoF pairs decreased from 53 to 37 and missense pairs from 70 to 40 at
## AF < 1e-4 (Table S27; at the other ceilings missense 32->26, 125->69, 178->94 and pLoF
## 45->39, 54->35, 54->35 for AF < 1e-5 / 1e-3 / 1e-2). This indicates that marginal calls
## are sensitive to residual variance calibration, as expected for an aggregate test near a
## study-wide threshold. Importantly, the strongest findings -- including the ALB and ALPL
## examples and the APOB/TET2 pLoF signals -- remained strictly significant after this
## correction. A direct test-level check reinforced this: ALLSPICE run on the synonymous
## pairs themselves was conservative, not inflated (lambda_gc_allspice = 0.03-0.38, all < 1;
## median synonymous ALLSPICE p ~ 0.9), so a matched negative control shows no inflation of
## the actual statistic; because lambda_gc_allspice < 1 it is reported as a diagnostic, not
## applied. Population stratification severe enough to manufacture the heterogeneity signal
## would inflate the rare-variant null far more, so the rare-variant findings are unlikely
## to be artifacts of residual structure. This is a worst-case bound, not an estimate of
## false positives: lambda_GC is computed from synonymous variants in the same significant
## genes, which themselves tag real signal, so it over-states the true structure-only
## inflation.
###############################################################################

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggpubr) })

result_path <- path.expand("~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/results/")
figure_path <- path.expand("~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/figure_2024/")
mt_dir      <- file.path(result_path, "modified_test"); dir.create(mt_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_path, showWarnings = FALSE, recursive = TRUE)
STRICT_P <- 4.23e-6; ORIG_AF <- 1e-4; AF_GRID <- c(1e-5, 1e-4, 1e-3, 1e-2)
SCAN_CSV     <- file.path(mt_dir, "figureS27_scan_real_triplets.csv")  # per-triplet x AF (cache)
TABLE_CSV    <- file.path(mt_dir, "figureS27_table.csv")               # = Supplementary Table S27
# hybrid per-variant betas, built by python/pull_sig_variant_betas_full_af.py --mode hybrid
# (AF<1e-4 from make_sig_variant_betas_all_annotations.py; AF>=1e-4 from --mode full_af).
VARIANTS_GCS <- "gs://ukb-diverse-pops/wlu/genebass_notebooks/sig_gene_variant_betas_500k_all_annotations_sweep.txt.gz"

# ---- compute the synonymous inflation (infl) + per-triplet AF scan (used when cache absent) ----
compute_scan <- function() {
  suppressMessages({ library(magrittr); library(mvtnorm) })
  if (!suppressMessages(suppressWarnings(require(ALLSPICER, quietly = TRUE)))) {
    pkg_dir <- Filter(dir.exists, c(
      path.expand("~/Partners HealthCare Dropbox/Wenhan Lu/github_repo/ALLSPICER/R"),
      path.expand("~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/R")))[1]
    source(file.path(pkg_dir, "utils.R")); source(file.path(pkg_dir, "ALLSPICE.R"))
  }
  read_gcs_fast <- function(gcs_uri) {
    cand <- unique(c("/Users/wlu/google-cloud-sdk/bin/gsutil", Sys.which("gsutil")))
    gsutil <- cand[nzchar(cand) & file.exists(cand)][1]
    if (is.na(gsutil)) stop("No usable gsutil found.")
    tmp <- tempfile(fileext = ".bgz"); on.exit(unlink(tmp), add = TRUE)
    system2(gsutil, c("cp", gcs_uri, tmp), stdout = TRUE, stderr = TRUE)
    if (!file.exists(tmp)) stop("gsutil cp failed: ", gcs_uri); fread(tmp)
  }
  phenoname_of <- function(d) paste0(d$trait_type, "_", d$phenocode, "_", d$pheno_sex, "_",
                                     ifelse(is.na(d$coding), "", d$coding), "_", d$modifier)
  res <- fread(file.path(result_path, "continuous_final_AF_1e_4_burden_syn_var_corr_500k_2023_results.csv"))
  res[, annotation := as.character(annotation)]; res <- res[sig_gene == 2 & n_var > 1 & !is.na(pvalue)]
  tri <- read_gcs_fast(VARIANTS_GCS)
  tri[, phenoname := phenoname_of(tri)]
  tri[, annotation := ifelse(annotation %in% c("missense", "LC"), "missense|LC", annotation)]

  triplet_variants <- function(g, anno, p1, p2) {
    sub <- tri[gene == g & annotation == anno & phenoname %in% c(p1, p2)]; if (nrow(sub) == 0) return(NULL)
    sub <- sub[, .(BETA = BETA[1]), by = .(locus, alleles, AF, phenoname)]
    w <- dcast(sub, locus + alleles + AF ~ phenoname, value.var = "BETA")
    if (!all(c(p1, p2) %in% names(w))) return(NULL)
    setnames(w, c(p1, p2), c("beta1", "beta2"))
    w <- w[complete.cases(w[, .(AF, beta1, beta2)])]; if (nrow(w) == 0) NULL else w
  }
  reconstruct_n_ind <- function(w, r, c_hat, lambda) {
    v <- w[AF < ORIG_AF]; if (nrow(v) < 2) return(NA_real_)
    A <- 2 * diag(v$AF); b1 <- t(as.matrix(v$beta1)); b2 <- t(as.matrix(v$beta2))
    Q <- c((b1 - c_hat * b2) %*% A %*% t(b1 - c_hat * b2)) / (c_hat^2 - 2 * c_hat * r + 1)
    lambda / Q
  }
  run_triplet <- function(w, r, n_ind, g, p1, p2, af_max) {
    v <- w[AF < af_max]; if (nrow(v) < 2 || is.na(n_ind)) return(NULL)
    ALLSPICE(data.frame(beta1 = v$beta1, beta2 = v$beta2, AF = v$AF),
             pheno_corr = r, n_ind = n_ind, gene = g, pheno1 = p1, pheno2 = p2,
             beta1_field = "beta1", beta2_field = "beta2", af_field = "AF")
  }
  prep <- function(rows) Filter(Negate(is.null), lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i]; w <- triplet_variants(r$gene, r$annotation, r$pheno1, r$pheno2); if (is.null(w)) return(NULL)
    n <- reconstruct_n_ind(w, r$corr, r$c_hat, r$lambda); if (is.na(n)) return(NULL)
    list(row = r, w = w, n_ind = n)
  }))

  # (1) TRAIT-SPECIFIC synonymous single-variant genomic control. For each phenotype
  #     (trait), lambda_GC = median(synonymous single-variant chi^2, df=1) / qchisq(0.5,1)
  #     over synonymous variants below each AF ceiling (require >= MIN_SYN tests for a
  #     stable per-trait median). Summarized per AF ceiling as median [IQR] across traits.
  MIN_SYN <- 10
  syn <- tri[annotation == "synonymous" & !is.na(Pvalue) & Pvalue > 0 & Pvalue <= 1 & !is.na(AF),
             .(phenoname, AF, chi2 = qchisq(Pvalue, df = 1, lower.tail = FALSE))]
  trait_lambda <- rbindlist(lapply(AF_GRID, function(af) {
    syn[AF < af, .(n_tests = .N, lambda_trait = median(chi2) / qchisq(0.5, 1)), by = phenoname][
        n_tests >= MIN_SYN][, af_threshold := af][]
  }))
  infl <- trait_lambda[, .(n_traits = .N, n_syn_tests = sum(n_tests),
                           lambda_trait_median = median(lambda_trait),
                           lambda_trait_q25 = quantile(lambda_trait, 0.25),
                           lambda_trait_q75 = quantile(lambda_trait, 0.75)),
                       by = af_threshold][order(af_threshold)]
  cat("== trait-specific synonymous single-variant GC (median [IQR] across traits) ==\n"); print(infl)

  # (1b) Diagnostic: run ALLSPICE ITSELF on the synonymous gene-trait pairs and
  #      genomic-control the ALLSPICE statistic (p-value scale, since df = n_var-1 varies).
  #      This matched test-level negative control comes out < 1 (conservative), so it is
  #      reported as a diagnostic, not applied as a correction.
  syn_prep <- prep(res[annotation == "synonymous"])
  cat(sprintf("covered synonymous ALLSPICE pairs: %d\n", length(syn_prep)))
  syn_allspice <- rbindlist(lapply(AF_GRID, function(af) {
    ps <- vapply(syn_prep, function(t) {
      o <- run_triplet(t$w, t$row$corr, t$n_ind, t$row$gene, t$row$pheno1, t$row$pheno2, af_max = af)
      if (is.null(o)) NA_real_ else o$pvalue
    }, numeric(1))
    ps <- ps[is.finite(ps) & ps > 0 & ps <= 1]
    data.table(af_threshold = af, lambda_gc_allspice = median(qchisq(ps, df = 1, lower.tail = FALSE)) / qchisq(0.5, 1))
  }))
  infl <- merge(infl, syn_allspice, by = "af_threshold")
  cat("== + synonymous ALLSPICE-level genomic control (diagnostic) ==\n"); print(infl)

  # (2) re-run real triplets and apply the TRAIT-SPECIFIC synonymous correction: for a
  #     pair (t1,t2), lambda_pair = mean of the two traits' lambda_GC, floored at 1 so the
  #     genomic control only deflates (never inflates); n_ind -> n_ind / lambda_pair.
  lam_key <- trait_lambda[, .(phenoname, af_threshold, lambda_trait)]
  get_lam <- function(ph, af) { v <- lam_key[phenoname == ph & af_threshold == af, lambda_trait]; if (length(v)) v[1] else NA_real_ }
  real_prep <- prep(res[annotation %in% c("missense|LC", "pLoF")])
  cat(sprintf("covered real triplets: %d\n", length(real_prep)))
  scan <- rbindlist(lapply(AF_GRID, function(af) {
    la <- infl[af_threshold == af, lambda_gc_allspice]; im <- infl[af_threshold == af]
    rbindlist(lapply(real_prep, function(t) {
      raw <- run_triplet(t$w, t$row$corr, t$n_ind, t$row$gene, t$row$pheno1, t$row$pheno2, af_max = af)
      if (is.null(raw)) return(NULL)
      df <- raw$n_var - 1  # correction only rescales n_ind -> recompute p from the same lambda/df
      l1 <- get_lam(t$row$pheno1, af); l2 <- get_lam(t$row$pheno2, af)
      lam_pair <- mean(c(l1, l2), na.rm = TRUE); if (!is.finite(lam_pair)) lam_pair <- 1
      lam_used <- max(lam_pair, 1)  # genomic control only deflates, never inflates
      data.table(af_threshold = af, lambda_gc_allspice = la,
                 lambda_trait_median = im$lambda_trait_median, lambda_trait_q25 = im$lambda_trait_q25,
                 lambda_trait_q75 = im$lambda_trait_q75, n_traits = im$n_traits, n_syn_tests = im$n_syn_tests,
                 gene = t$row$gene, annotation = t$row$annotation, pheno1 = t$row$pheno1, pheno2 = t$row$pheno2,
                 n_var = raw$n_var, lambda_pair = lam_pair, lambda_used = lam_used, p_raw = raw$pvalue,
                 p_traitspec_corr = 1 - pchisq(raw$lambda / lam_used, df))
    }))
  }))
  fwrite(scan, SCAN_CSV)
  scan
}

# ---- load per-triplet scan (cache) or compute ----
scan <- if (file.exists(SCAN_CSV)) fread(SCAN_CSV) else compute_scan()
if (!"p_traitspec_corr" %in% names(scan)) scan <- compute_scan()  # cache predates the trait-specific correction
infl <- unique(scan[, .(af_threshold, n_traits, n_syn_tests, lambda_trait_median,
                        lambda_trait_q25, lambda_trait_q75, lambda_gc_allspice)])[order(af_threshold)]

# ---- Supplementary Table S27: inflation + variants-remaining + strict counts ----
# For AF ceilings that act as the "strata" here, "variants remaining" = number of
# variants entering ALLSPICE per gene-trait pair at that ceiling (median [IQR]).
supp_table <- merge(
  infl,
  scan[, .(n_pairs = .N,
           n_var_median = as.double(median(n_var)),
           n_var_q25 = as.double(quantile(n_var, 0.25)),
           n_var_q75 = as.double(quantile(n_var, 0.75)),
           n_strict_raw            = sum(p_raw            < STRICT_P),
           n_strict_traitspec_corr = sum(p_traitspec_corr < STRICT_P)),
       by = .(af_threshold, annotation)],
  by = "af_threshold")[order(annotation, af_threshold)]
fwrite(supp_table, TABLE_CSV)   # = Supplementary Table S27
cat("== Table S27: inflation, variants remaining per pair, and strict counts ==\n")
print(supp_table)

## ---- Panel A: trait-specific synonymous single-variant lambda_GC (median [IQR] across
## traits; the correction applied) vs the ALLSPICE statistic run on synonymous pairs
## (lambda_gc_allspice, < 1 => the actual test is conservative on a matched null). ----
pA_dat <- infl[order(af_threshold)]
pA <- ggplot(pA_dat, aes(factor(af_threshold))) +
  geom_hline(yintercept = 1, lty = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = lambda_trait_q25, ymax = lambda_trait_q75,
                    colour = "trait-specific single-variant lambda_GC (correction)"), width = 0.15) +
  geom_line(aes(y = lambda_trait_median, colour = "trait-specific single-variant lambda_GC (correction)", group = 1)) +
  geom_point(aes(y = lambda_trait_median, colour = "trait-specific single-variant lambda_GC (correction)"), size = 2.6) +
  geom_line(aes(y = lambda_gc_allspice, colour = "ALLSPICE-level lambda_GC (synonymous pairs)", group = 1)) +
  geom_point(aes(y = lambda_gc_allspice, colour = "ALLSPICE-level lambda_GC (synonymous pairs)"), size = 2.6, shape = 17) +
  scale_x_discrete(labels = c("1e-5", "1e-4", "1e-3", "1e-2")) +
  scale_colour_manual(values = c("trait-specific single-variant lambda_GC (correction)" = "#2C7FB8",
                                 "ALLSPICE-level lambda_GC (synonymous pairs)" = "#31A354"), name = NULL) +
  ylim(0, max(1.2, max(pA_dat$lambda_trait_q75, na.rm = TRUE) * 1.05)) +
  labs(x = "AF ceiling", y = expression("Genomic-control "*lambda[GC]*" (median [IQR] across traits)")) +
  theme_classic(base_size = 12) +
  theme(legend.position = "inside", legend.position.inside = c(0.5, 0.55))

## ---- Panel B: real triplets re-tested at the primary AF<1e-4 ceiling, before vs after
## the trait-specific synonymous correction. ----
d4 <- scan[af_threshold == ORIG_AF]
pB <- ggplot(d4, aes(-log10(p_raw), -log10(p_traitspec_corr), colour = annotation)) +
  geom_abline(slope = 1, lty = 2, colour = "grey60") +
  geom_hline(yintercept = -log10(STRICT_P), lty = 3, colour = "red") +
  geom_vline(xintercept = -log10(STRICT_P), lty = 3, colour = "red") +
  geom_point(size = 1.9, alpha = 0.8) +
  annotation_color_scale2+
  labs(x = expression("Raw "*-log[10](p)),
       y = expression(atop("After trait-specific synonymous", lambda[GC]*" correction  "*-log[10](p)))) +
  theme_classic(base_size = 12) +
  theme(legend.position = "inside", legend.position.inside = c(0.78, 0.16))

fig <- ggpubr::ggarrange(pA + theme(plot.margin = margin(0.8, 0.1, 0.1, 0.1, "cm")),
                         pB + theme(plot.margin = margin(0.8, 0.1, 0.1, 0.1, "cm")), nrow = 1, labels = c("(A) Negative control (dashed = no inflation; < 1 = conservative)",
                                                      "(B) Real triplets re-tested (AF < 1e-4)"),
                 font.label = list(size = 12, face = "bold"), hjust = 0, vjust =1)
# Draw via device + print(): ggsave() calls ggplot_build() in ggplot2 >= 3.5, which has
# no method for a ggarrange (list-based) object -> "no applicable method for 'ggplot_build'".
png(file.path(figure_path, "figureS27.png"), width = 10, height = 4.6, units = "in", res = 300); print(fig); dev.off()
pdf(file.path(figure_path, "figureS27.pdf"), width = 10, height = 4.6); print(fig); dev.off()
cat("wrote figureS27.{png,pdf} to", figure_path, "\n")

