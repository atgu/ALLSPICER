###############################################################################
## figureS28.R
##
## Figure S28 (two panels) -- pLoF heterogeneity in APOB and TET2 (reviewer
## Comment 2 + follow-up):
##   (A) Annotation-stratified analysis (Task 2 follow-up): re-runs the strictly-
##       significant APOB/TET2 pLoF triplets after restricting variants by consequence
##       class, gnomAD match, PASS filter, LCR/segdup overlap, LOFTEE confidence/flags,
##       and NMD targeting, to localize which variant class carries the surviving signal.
##   (B) COUNT-MATCHED SUBSAMPLING CONTROL: because ALLSPICE aggregates a per-variant
##       deviation, its power grows with the number of variants tested, so a stratum
##       that keeps fewer variants loses strict pairs mechanically. We build a
##       rarefaction null -- strict pairs vs number of variants, from random subsets of
##       the FULL pLoF set -- and, per stratum, resample the full set down to that
##       stratum's per-pair size and compare the observed strict count to the size-
##       matched null. On the curve: ON = count-driven; BELOW = the class carries less
##       heterogeneity per variant; ABOVE = more.
##
## SELF-CONTAINED and the SOLE owner of this analysis (supersedes the former
## R/task2_stratified_apob_tet2.R and figureS28_S29.R). Caches to results/modified_test/;
## delete a cache CSV to force that step to refresh.
##
## INPUTS  : hybrid sweep per-variant betas (GCS) + gnomAD pLoF annotation TSV
##           (python/pull_gnomad_plof_annotations.py) + the ALLSPICER package +
##           the published results CSV (for n_ind reconstruction).
## OUTPUTS : figureS28.{png,pdf} (2-panel A|B) -> figure_path ;
##           task2_stratified_apob_tet2.csv, task2_stratified_summary.csv (=Table S28),
##           task2_stratified_heatmap.png, figureS28_rarefaction_null.csv,
##           tableS29_countmatched.csv (=Table S29) -> results/modified_test/
##
## ---------------------------------------------------------------------------
## FIGURE CAPTION (manuscript)
## ---------------------------------------------------------------------------
## Figure S28. Annotation-stratified analysis of pLoF heterogeneity in APOB and TET2,
## with a count-matched control.
## (A) Bars show the number of strictly significant pLoF association pairs (of 53)
## retained after restricting variants by consequence class (frameshift only, stop-gain/
## splice only, SNVs only), gnomAD v4.1.1 match status, PASS-filter status, low-complexity
## or segmental-duplication overlap, LOFTEE confidence and flag status, and predicted NMD
## targeting. The signals persist after filters targeting variant quality (PASS only, 50
## of 53), repeat-region mapping artifacts (outside LCR/segdup, 53), hard LOFTEE flags
## (high-confidence with hard flags removed, 51), and NMD escape (NMD-targeted only, 51).
## The strata that retain few pairs are those that keep few variants (SNVs only and stop-
## gain/splice only, ~65-68 variants, 2 each; frameshift only, 230 variants, 35) or that
## additionally remove the soft PHYLOCSF_WEAK conservation flag ("LOFTEE HC, all flags
## removed", 13) rather than reflecting low variant quality. Dashed line, total strictly
## significant pLoF pairs (53).
## (B) Count-matched control. Number of strictly significant pairs (of 53) as a function
## of the number of variants entering the test. Grey line and ribbon, median and 5-95%
## interval of the strict count from random subsets of the full pLoF set at each size
## (rarefaction null over 500 resamples) -- the count-only expectation. Points, the
## annotation-defined strata at their median per-pair variant count and observed strict
## count, coloured by filter category. Strata lying on the grey curve are explained by the
## number of variants alone; strata falling below the size-matched null (Table S29) would
## carry less heterogeneity per variant than a random pLoF subset of equal size -- none do,
## and frameshift-only sits modestly above, indicating the heterogeneity is distributed
## across pLoF variants roughly in proportion to their number and is not concentrated in
## (or an artifact of) any single variant class.
###############################################################################

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggpubr) })

result_path <- path.expand("~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/results/")
figure_path <- path.expand("~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/figure_2024/")
mt_dir      <- file.path(result_path, "modified_test"); dir.create(mt_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_path, showWarnings = FALSE, recursive = TRUE)
STRICT_P <- 4.23e-6; ORIG_AF <- 1e-4
RESULTS_CSV  <- file.path(mt_dir, "task2_stratified_apob_tet2.csv")  # per-triplet x stratum (cache)
SUMMARY_CSV  <- file.path(mt_dir, "task2_stratified_summary.csv")    # = Supplementary Table S28
RAREF_CSV    <- file.path(mt_dir, "figureS28_rarefaction_null.csv")  # panel B rarefaction null (cache)
MATCH_CSV    <- file.path(mt_dir, "tableS29_countmatched.csv")       # = Supplementary Table S29 (cache)
# VARIANTS_GCS: hybrid per-variant betas, built by python/pull_sig_variant_betas_full_af.py --mode hybrid
#   (AF<1e-4 half = published set from make_sig_variant_betas_all_annotations.py; AF>=1e-4 half = --mode full_af).
# ANNOT_GCS: gnomAD v4.1.1 pLoF annotations for APOB/TET2, built by python/pull_gnomad_plof_annotations.py.
VARIANTS_GCS <- "gs://ukb-diverse-pops/wlu/genebass_notebooks/sig_gene_variant_betas_500k_all_annotations_sweep.txt.gz"
ANNOT_GCS    <- "gs://aou_wlu/allspice/apob_tet2_plof_gnomad_v4.1.1.tsv.bgz"
B_REP  <- 500L                                                       # resamples per size / stratum
K_GRID <- c(8, 15, 25, 40, 55, 70, 90, 120, 160, 200, 240, 290)      # rarefaction variant counts

# Stratum predicates on the joined variant table (applied at AF < 1e-4). Variants absent
# from gnomAD (matched = NA) are excluded from any annotation-based stratum.
fl <- function(x) !is.na(x) & x
STRATA <- list(
  all = function(w) rep(TRUE, nrow(w)),
  snv_only = function(w) w$is_snv,                                    # = Task 2 (indels removed)
  frameshift_only = function(w) fl(w$is_frameshift),
  ptv_non_frameshift = function(w) fl(w$matched) & !is.na(w$lof) & !fl(w$is_frameshift),
  gnomad_matched = function(w) fl(w$matched),
  loftee_HC = function(w) w$lof == "HC" & fl(w$matched),
  loftee_HC_noflag = function(w) w$lof == "HC" & (is.na(w$lof_flags) | w$lof_flags == "") & fl(w$matched),
  # only HARD flags removed -- KEEPS the soft PHYLOCSF_WEAK conservation warning (which
  # dominates the flagged set); the contrast with loftee_HC_noflag isolates its effect.
  loftee_HC_no_hardflag = function(w) w$lof == "HC" & fl(w$matched) &
    (is.na(w$lof_flags) | w$lof_flags == "" |
     !grepl("NAGNAG|END_TRUNC|DONOR|SPLICE|ANC_ALLELE|NON_CAN|SINGLE_EXON", w$lof_flags)),
  pass_only = function(w) fl(w$pass_qc),
  outside_lcr_segdup = function(w) fl(w$matched) & !fl(w$lcr) & !fl(w$segdup),
  nmd_targeted_only = function(w) fl(w$matched) & !fl(w$in_last_exon)
)
STRATA_LEVELS <- names(STRATA)

# ---------------------------------------------------------------------------
# Shared helpers (used by both the stratified compute and the rarefaction control)
# ---------------------------------------------------------------------------
load_allspice <- function() {
  suppressMessages({ library(magrittr); library(mvtnorm) })
  if (!suppressMessages(suppressWarnings(require(ALLSPICER, quietly = TRUE)))) {
    pkg_dir <- Filter(dir.exists, c(
      path.expand("~/Partners HealthCare Dropbox/Wenhan Lu/github_repo/ALLSPICER/R"),
      path.expand("~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/R")))[1]
    source(file.path(pkg_dir, "utils.R")); source(file.path(pkg_dir, "ALLSPICE.R"))
  }
}
read_gcs_fast <- function(gcs_uri) {
  gsutil <- c("/Users/wlu/google-cloud-sdk/bin/gsutil", Sys.which("gsutil"))
  gsutil <- gsutil[nzchar(gsutil) & file.exists(gsutil)][1]
  if (is.na(gsutil)) stop("No usable gsutil found.")
  tmp <- tempfile(fileext = ".bgz"); on.exit(unlink(tmp), add = TRUE)
  system2(gsutil, c("cp", gcs_uri, tmp), stdout = TRUE, stderr = TRUE)
  if (!file.exists(tmp)) stop("gsutil cp failed: ", gcs_uri)
  fread(tmp)
}
phenoname_of <- function(d) paste0(d$trait_type, "_", d$phenocode, "_", d$pheno_sex, "_",
                                   ifelse(is.na(d$coding), "", d$coding), "_", d$modifier)
reconstruct_n_ind <- function(v, r, c_hat, lambda) {   # v = variant table at AF < ORIG_AF
  if (nrow(v) < 2) return(NA_real_)
  A <- 2 * diag(v$AF); b1 <- t(as.matrix(v$beta1)); b2 <- t(as.matrix(v$beta2))
  Q <- c((b1 - c_hat * b2) %*% A %*% t(b1 - c_hat * b2)) / (c_hat^2 - 2 * c_hat * r + 1)
  lambda / Q
}
# Closed-form ALLSPICE p-value (A = 2*diag(af); quadratic forms reduce to weighted sums).
# Matches ALLSPICE() exactly but avoids data.frame overhead for the resampling loops.
allspice_p <- function(beta1, beta2, af, r, n_ind) {
  n <- length(beta1); if (n < 2) return(NA_real_)
  u <- sum(af * beta2 * (beta1 - r * beta2))
  v <- sum(af * (beta2^2 - beta1^2))
  w <- sum(af * beta1 * (r * beta1 - beta2))
  disc <- v^2 - 4 * u * w; if (u == 0 || disc < 0) return(NA_real_)
  c1 <- (-v + sqrt(disc)) / (2 * u); c2 <- (-v - sqrt(disc)) / (2 * u)
  chat <- if (u > 0) max(c1, c2) else min(c1, c2)
  lambda <- (n_ind / (chat^2 - 2 * chat * r + 1)) * (2 * sum(af * (beta1 - chat * beta2)^2))
  1 - pchisq(lambda, n - 1)
}

# Build one entry per strictly-significant APOB/TET2 pLoF pair: the full pLoF variant
# table (AF<1e-4, gnomAD annotations joined) + reconstructed n_ind + phenotype corr.
build_triplets <- function() {
  load_allspice()
  res <- fread(file.path(result_path, "continuous_final_AF_1e_4_burden_syn_var_corr_500k_2023_results.csv"))
  res[, annotation := as.character(annotation)]; res <- res[sig_gene == 2 & n_var > 1 & !is.na(pvalue)]
  tri <- read_gcs_fast(VARIANTS_GCS)
  tri[, phenoname := phenoname_of(tri)]
  tri[, annotation := ifelse(annotation %in% c("missense", "LC"), "missense|LC", annotation)]
  ann <- read_gcs_fast(ANNOT_GCS); ann[, matched := TRUE]
  if ("is_snv" %in% names(ann)) ann[, is_snv := NULL]
  setkey(ann, variant_id)

  triplet_variants <- function(g, p1, p2) {
    sub <- tri[gene == g & annotation == "pLoF" & phenoname %in% c(p1, p2)]
    if (nrow(sub) == 0) return(NULL)
    sub <- sub[, .(BETA = BETA[1]), by = .(locus, alleles, AF, phenoname)]
    w <- dcast(sub, locus + alleles + AF ~ phenoname, value.var = "BETA")
    if (!all(c(p1, p2) %in% names(w))) return(NULL)
    setnames(w, c(p1, p2), c("beta1", "beta2"))
    w <- w[complete.cases(w[, .(AF, beta1, beta2)])]; if (nrow(w) == 0) return(NULL)
    a <- gsub('^\\[|\\]$|"', "", w$alleles)
    w[, ref := tstrsplit(a, ",", fixed = TRUE)[[1]]]; w[, alt := tstrsplit(a, ",", fixed = TRUE)[[2]]]
    w[, is_snv := nchar(ref) == 1L & nchar(alt) == 1L]
    w[, variant_id := paste(locus, ref, alt, sep = ":")]
    merge(w, ann, by = "variant_id", all.x = TRUE)
  }
  rows <- res[gene %in% c("APOB", "TET2") & annotation == "pLoF" & pvalue < STRICT_P]
  cat(sprintf("strictly-significant APOB/TET2 pLoF triplets: %d\n", nrow(rows)))
  tlist <- lapply(seq_len(nrow(rows)), function(i) {
    row <- rows[i]
    w <- triplet_variants(row$gene, row$pheno1, row$pheno2); if (is.null(w)) return(NULL)
    n_ind <- reconstruct_n_ind(w[AF < ORIG_AF], row$corr, row$c_hat, row$lambda)
    if (is.na(n_ind)) return(NULL)
    list(gene = row$gene, p1 = row$pheno1, p2 = row$pheno2, r = row$corr, n_ind = n_ind, w = w)
  })
  Filter(Negate(is.null), tlist)
}

# Per-triplet x stratum ALLSPICE results (panel A / Table S28).
compute_strata <- function(tlist) {
  results <- rbindlist(lapply(tlist, function(t) {
    w <- t$w
    rbindlist(lapply(STRATA_LEVELS, function(st) {
      v <- w[STRATA[[st]](w) & AF < ORIG_AF]
      p <- if (nrow(v) < 2) NA_real_ else allspice_p(v$beta1, v$beta2, v$AF, t$r, t$n_ind)
      data.table(gene = t$gene, pheno1 = t$p1, pheno2 = t$p2, stratum = st,
                 n_var = nrow(v), pvalue = p)
    }))
  }))
  results[, strict := pvalue < STRICT_P]
  fwrite(results, RESULTS_CSV)
  results
}

# Rarefaction null + per-stratum count-matched test (panel B / Table S29).
compute_rarefaction <- function(tlist, results) {
  full <- lapply(tlist, function(t) t$w[AF < ORIG_AF, .(beta1, beta2, AF)])  # full pLoF set / pair
  r    <- vapply(tlist, `[[`, numeric(1), "r")
  nind <- vapply(tlist, `[[`, numeric(1), "n_ind")
  tkey <- vapply(tlist, function(t) paste(t$gene, t$p1, t$p2), character(1))

  # one resample: for each pair draw `sizes[i]` variants (capped at its pool), count strict pairs
  draw_strict <- function(sizes) {
    sum(vapply(seq_along(full), function(i) {
      m <- nrow(full[[i]]); k <- min(sizes[i], m); if (k < 2) return(FALSE)
      idx <- sample.int(m, k)
      p <- allspice_p(full[[i]]$beta1[idx], full[[i]]$beta2[idx], full[[i]]$AF[idx], r[i], nind[i])
      !is.na(p) && p < STRICT_P
    }, logical(1)))
  }
  pool <- vapply(full, nrow, numeric(1))

  # (i) rarefaction: same target size k for every pair (capped at its pool)
  raref <- rbindlist(lapply(K_GRID, function(k) {
    d <- vapply(seq_len(B_REP), function(b) draw_strict(rep(k, length(full))), integer(1))
    data.table(k = k, mean_n_var = mean(pmin(k, pool)),
               n_strict_mean = mean(d), q05 = quantile(d, .05), q50 = median(d), q95 = quantile(d, .95))
  }))

  # (ii) count-matched: per stratum, match each pair's own retained size, resample from its full pool
  matched <- rbindlist(lapply(setdiff(STRATA_LEVELS, "all"), function(st) {
    sz <- results[stratum == st]; sz[, kk := paste(gene, pheno1, pheno2)]
    tgt <- sz$n_var[match(tkey, sz$kk)]; tgt[is.na(tgt)] <- 0
    obs <- sum(results[stratum == st, strict], na.rm = TRUE)
    d <- vapply(seq_len(B_REP), function(b) draw_strict(tgt), integer(1))
    data.table(stratum = st, median_size = as.double(median(sz$n_var)), observed = obs,
               null_mean = mean(d), null_q05 = quantile(d, .05), null_q50 = median(d), null_q95 = quantile(d, .95),
               p_depleted = mean(d <= obs), p_enriched = mean(d >= obs))
  }))
  matched[, verdict := fifelse(observed < null_q05, "below null (class < random of same size)",
                        fifelse(observed > null_q95, "above null (class > random of same size)",
                                "within null (count-driven)"))]
  fwrite(raref, RAREF_CSV); fwrite(matched, MATCH_CSV)
  list(raref = raref, matched = matched)
}

# ---------------------------------------------------------------------------
# Orchestration: build once if any step needs recomputing; otherwise read caches.
# ---------------------------------------------------------------------------
need_strata <- !file.exists(RESULTS_CSV)
need_raref  <- !(file.exists(RAREF_CSV) && file.exists(MATCH_CSV))
if (need_strata || need_raref) tlist <- build_triplets()

results <- if (need_strata) compute_strata(tlist) else fread(RESULTS_CSV)
results[, stratum := factor(stratum, levels = STRATA_LEVELS)]

# One human-readable line per stratum: what it restricts to, # pairs tested, # strict
# retained, and the number of variants remaining per pair (median [IQR]). = Table S28.
STRAT_DESC <- c(
  all = "All pLoF variants (AF<1e-4)",
  snv_only = "SNVs only (indels removed)",
  frameshift_only = "Frameshift indels only",
  ptv_non_frameshift = "Stop-gain / splice only (non-frameshift PTV)",
  gnomad_matched = "Present in gnomAD v4.1.1",
  loftee_HC = "LOFTEE high-confidence (HC)",
  loftee_HC_noflag = "LOFTEE HC, all flags removed (drops soft PHYLOCSF_WEAK)",
  loftee_HC_no_hardflag = "LOFTEE HC, hard flags removed (keeps PHYLOCSF_WEAK)",
  pass_only = "gnomAD PASS filter",
  outside_lcr_segdup = "Outside low-complexity / segmental-duplication regions",
  nmd_targeted_only = "NMD-targeted (not last-exon)")
summ <- results[, .(
  description  = STRAT_DESC[as.character(stratum[1])],
  n_pairs      = .N,
  n_testable   = sum(!is.na(pvalue)),
  n_strict     = sum(strict, na.rm = TRUE),
  n_var_median = as.double(median(n_var)),
  n_var_q25    = as.double(quantile(n_var, 0.25)),
  n_var_q75    = as.double(quantile(n_var, 0.75)),
  n_var_min    = min(n_var),
  n_var_max    = max(n_var)
), by = stratum][order(stratum)]
fwrite(summ, SUMMARY_CSV)   # = Supplementary Table S28
cat("== Table S28: APOB/TET2 pLoF strict pairs + variants remaining per stratum ==\n")
print(summ)

# ============================================================================
# Panel A: strict pLoF pairs retained per stratum
# ============================================================================
TOTAL <- summ[stratum == "all", n_strict]
lab <- c(all = "All pLoF variants", gnomad_matched = "In gnomAD v4.1.1",
         frameshift_only = "Frameshift only", ptv_non_frameshift = "Stop-gain / splice only",
         snv_only = "SNVs only (indels removed)", pass_only = "gnomAD PASS only",
         outside_lcr_segdup = "Outside LCR / segdup", loftee_HC = "LOFTEE high-confidence",
         loftee_HC_no_hardflag = "LOFTEE HC, hard flags removed",
         loftee_HC_noflag = "LOFTEE HC, all flags removed", nmd_targeted_only = "NMD-targeted only")
cat_map <- c(all = "Baseline", gnomad_matched = "Baseline",
             frameshift_only = "Consequence class", ptv_non_frameshift = "Consequence class",
             snv_only = "Consequence class", pass_only = "Quality / region",
             outside_lcr_segdup = "Quality / region", loftee_HC = "LOFTEE confidence",
             loftee_HC_no_hardflag = "LOFTEE confidence", loftee_HC_noflag = "LOFTEE confidence",
             nmd_targeted_only = "NMD")
CAT_LEVELS <- c("Baseline", "Consequence class", "Quality / region", "LOFTEE confidence", "NMD")
CAT_COLS <- c("Baseline" = "grey55", "Consequence class" = "#D55E00",
              "Quality / region" = "#0072B2", "LOFTEE confidence" = "#009E73", "NMD" = "#CC79A7")
b <- summ[match(STRATA_LEVELS, stratum)]
b[, label := factor(lab[as.character(stratum)], levels = rev(unname(lab[STRATA_LEVELS])))]
b[, category := factor(cat_map[as.character(stratum)], levels = CAT_LEVELS)]

pA <- ggplot(b, aes(label, n_strict, fill = category)) +
  geom_col(width = 0.72, alpha = 0.92) +
  geom_hline(yintercept = TOTAL, lty = 2, colour = "grey40") +
  geom_text(aes(label = n_strict), hjust = -0.25, size = 3.4) + coord_flip() +
  scale_fill_manual(values = CAT_COLS, name = NULL) +
  scale_y_continuous(limits = c(0, TOTAL * 1.08), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = sprintf("Strictly significant pLoF pairs retained (of %d)", TOTAL)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 10)) +
  guides(fill = guide_legend(nrow = 2))

# ---- secondary QC view (not a panel): -log10(p) heatmap (triplet x stratum) ----
d <- copy(results)
d[, pair := paste0(gene, ": ", sub("continuous_", "", pheno1), " x ", sub("continuous_", "", pheno2))]
d[, mlogp := -log10(pvalue)]
ph <- ggplot(d, aes(stratum, pair, fill = mlogp)) +
  geom_tile(colour = "white") +
  geom_tile(data = d[strict == TRUE], colour = "black", linewidth = 0.7, fill = NA) +
  geom_text(aes(label = ifelse(is.na(mlogp), "NA", sprintf("%.0f\n(%d)", mlogp, n_var))), size = 2.6) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey85", name = expression(-log[10](p))) +
  labs(x = NULL, y = NULL, title = "APOB/TET2 pLoF heterogeneity by variant stratum",
       subtitle = "black border = strictly significant (p < 4.23e-6); cell = -log10(p) (n_var)") +
  theme_classic(base_size = 11) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(mt_dir, "task2_stratified_heatmap.png"), ph,
       width = 10, height = 1.2 + 0.42 * length(unique(d$pair)), dpi = 300, limitsize = FALSE)
cat("wrote task2_stratified_heatmap.png to", mt_dir, "\n")

# ============================================================================
# Panel B: count-matched control (rarefaction null + per-stratum size-matched test)
# ============================================================================
# Validate the fast closed-form allspice_p against the cached ALLSPICE 'all' p-values
# (also confirms build_triplets reproduces the full pLoF set used for the barplot).
if (exists("tlist")) {
  chk <- vapply(tlist, function(t) { v <- t$w[AF < ORIG_AF]; allspice_p(v$beta1, v$beta2, v$AF, t$r, t$n_ind) }, numeric(1))
  key <- vapply(tlist, function(t) paste(t$gene, t$p1, t$p2), character(1))
  ar  <- results[stratum == "all"]; ar[, kk := paste(gene, pheno1, pheno2)]
  ref <- ar$pvalue[match(key, ar$kk)]
  cat(sprintf("check: allspice_p vs cached 'all' pvalue -- max |delta(-log10 p)| = %.2g\n",
              max(abs(-log10(chk) + log10(ref)), na.rm = TRUE)))
}
if (need_raref) {
  rr <- compute_rarefaction(tlist, results); raref <- rr$raref; matched <- rr$matched
} else { raref <- fread(RAREF_CSV); matched <- fread(MATCH_CSV) }
cat("\n== Table S29: count-matched control (observed strict vs size-matched null) ==\n")
print(matched[, .(stratum, median_size, observed, null_mean, null_q05, null_q95, p_depleted, verdict)])

# Strata plotted at (median variants, observed strict), over the rarefaction null band.
pts <- summ[stratum != "all", .(stratum, n_var_median, n_strict)]
pts[, category := factor(cat_map[as.character(stratum)], levels = CAT_LEVELS)]
pts[, label := lab[as.character(stratum)]]
pts <- merge(pts, matched[, .(stratum, verdict)], by = "stratum", all.x = TRUE)
pts[, off_null := grepl("below|above", verdict)]

pB <- ggplot() +
  geom_ribbon(data = raref, aes(mean_n_var, ymin = q05, ymax = q95), fill = "grey80", alpha = 0.6) +
  geom_line(data = raref, aes(mean_n_var, q50), colour = "grey45") +
  geom_point(data = pts, aes(n_var_median, n_strict, colour = category, shape = off_null), size = 3) +
  scale_colour_manual(values = CAT_COLS, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
  labs(x = "Variants entering the test per pair (stratum median)",
       y = sprintf("Strictly significant pLoF pairs (of %d)", TOTAL)) +
  theme_classic(base_size = 12) + theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2))
pB <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  pB + ggrepel::geom_text_repel(data = pts, aes(n_var_median, n_strict, label = label, colour = category),
                                size = 3, show.legend = FALSE, max.overlaps = 20, seed = 1)
} else {
  pB + geom_text(data = pts, aes(n_var_median, n_strict, label = label, colour = category),
                 vjust = -0.8, size = 2.8, show.legend = FALSE)
}

# ============================================================================
# Figure S28 = two panels (A: stratified barplot, B: count-matched control)
# ============================================================================
fig <- ggpubr::ggarrange(pA, pB, ncol = 2, widths = c(1, 1), labels = c("A", "B"),
                 font.label = list(size = 14, face = "bold"))
# Draw via device + print(): ggsave() calls ggplot_build() in ggplot2 >= 3.5, which has
# no method for a ggarrange (list-based) object -> "no applicable method for 'ggplot_build'".
png(file.path(figure_path, "figureS28.png"), width = 14, height = 6, units = "in", res = 300); print(fig); dev.off()
pdf(file.path(figure_path, "figureS28.pdf"), width = 14, height = 6); print(fig); dev.off()
cat("wrote figureS28.{png,pdf} (2-panel) to", figure_path, "\n")
