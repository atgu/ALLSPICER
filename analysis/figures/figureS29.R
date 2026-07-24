# =====================================================================
# Figure S29. Functional-domain and regional missense constraint (RMC)
# annotation of missense-driven ALLSPICE signals.
#
# Rare missense variants (AF < 1e-4) from strictly-significant ALLSPICE
# gene-trait pairs are annotated by (i) protein-domain membership and (ii)
# regional missense O/E (from the exploded public gnomAD RMC; see
# explode_rmc_to_variants.py / gnomad_methods PR #789). Panels:
#   A  gene-level constraint (missense Z, RMC min O/E): sig vs rest
#   B  enrichment of nominally trait-associated variants in domains / constrained regions
#   C  regional O/E vs variant effect magnitude (|beta|max)
#   D  beta-beta scatter coloured by protein domain (ALB, single-region example)
#   E  pair-level cross-trait heterogeneity (1-|rho_w|) vs gene constraint (missense Z, RMC min O/E)
#   F  lollipop of variant effect along the protein with RMC context (HECTD4,
#      multi-region; ALB/ALPL are single-region RMC so are not shown as lollipops)
#
# Domain and constraint annotations enrich for phenotypically associated
# variants but do not explain the magnitude of non-proportional cross-trait
# heterogeneity.
# =====================================================================
#
# ---------------------------------------------------------------------
# FIGURE LEGEND (for the supplement)
# ---------------------------------------------------------------------
# Figure S29. Functional-domain and regional missense constraint (RMC) annotation
# of missense-driven ALLSPICE signals.
# (A) Gene-level constraint (missense Z-score; minimum regional missense O/E) for
#     strictly-significant genes versus other tested genes; boxplots over violins,
#     with per-group gene counts and Wilcoxon rank-sum p-values.
# (B) Fraction of rare missense variants in an annotated protein domain or in a
#     significantly constrained RMC region (chi-square p < 1e-3), split by whether
#     the variant is nominally trait-associated (p < 0.05 for either trait, used
#     descriptively); Fisher's exact odds ratios and p-values and per-bar variant
#     counts are shown.
# (C) Variant level: regional missense O/E versus variant effect magnitude
#     (maximum |beta| across the two traits); each point is one rare missense
#     variant, coloured by nominal association; Spearman correlation shown.
# (D) Variant-level beta-beta scatter for ALB (Albumin vs Calcium); each point is a
#     rare missense variant coloured by protein domain and shaped by nominal
#     significance in each trait. ALB is single-region for RMC, so protein domains
#     rather than regional O/E are shown here.
# (E) Pair level: cross-trait effect heterogeneity (1 - |rho_w|, where rho_w is the
#     allele-frequency-weighted cosine correlation ALLSPICE tests) versus gene-level
#     constraint (missense Z; minimum regional O/E); each point is a significant
#     gene-trait pair. Dashed lines are visual trends only; Spearman rho and p are
#     shown per subpanel.
# (F) Variant effects along the protein for the multi-region example HECTD4
#     (Reticulocyte count vs Mean sphered cell volume), against genomic position;
#     background shaded by regional O/E, stems/points coloured by nominal
#     significance.
# Panels operate at different levels: C is variant-level effect magnitude versus
# regional O/E, whereas E is pair-level heterogeneity versus gene-level constraint.
# Domain and constraint annotations enrich for phenotypically associated variants (B)
# but do not explain variant effect magnitude (C) or the pair-level effect
# heterogeneity underlying the ALLSPICE signal (E). Related to SuppTableS32.
# ---------------------------------------------------------------------
source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/constants.R')
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2)
                                 library(readr); library(purrr); library(stringr); library(ggpubr) })
set.seed(42)

data_dir <- '~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/ajhg_revision/data/'
fig_out  <- '~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/figure_2024/'
RMC_VAR_URI <- 'gs://ukb-diverse-pops/wlu/genebass_notebooks/allspice/variants_with_rmc.tsv.bgz'

NOMINAL <- 0.05; RMC_SIG <- 1e-3; STRICT_P <- 4.23e-6
accent  <- '#C0392B'; neutral <- 'grey72'
fill_vals <- c('Rest' = neutral, 'Sig' = accent, 'Not assoc.' = neutral, 'Assoc.' = accent)

# shared theme for a consistent, clean look across panels
base_theme <- theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = 'bold', size = 12.5),
        plot.title.position = 'plot',
        plot.subtitle = element_text(size = 9, colour = 'grey35'),
        plot.margin = margin(10, 12, 10, 12),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold', size = 10),
        axis.title = element_text(size = 11),
        legend.title = element_text(size = 10), legend.text = element_text(size = 9),
        legend.key.size = unit(.4, 'cm'))

# ---- helpers ---------------------------------------------------------
read_gcs_fast <- function(gcs_uri, ..., gsutil = NULL) {
  if (is.null(gsutil)) {
    cand <- c("/Users/wlu/google-cloud-sdk/bin/gsutil", Sys.which("gsutil"))
    cand <- unique(cand[nzchar(cand)]); gsutil <- cand[file.exists(cand)][1]
  }
  if (is.na(gsutil) || !nzchar(gsutil) || !file.exists(gsutil)) stop("No usable gsutil found.")
  ext <- tools::file_ext(sub("^gs://[^/]+/", "", gcs_uri))
  tmp <- tempfile(fileext = if (nzchar(ext)) paste0(".", ext) else ""); on.exit(unlink(tmp), add = TRUE)
  res <- tryCatch(system2(gsutil, c("cp", gcs_uri, tmp), stdout = TRUE, stderr = TRUE), error = function(e) e)
  if (inherits(res, "error") || (!is.null(attr(res, "status")) && attr(res, "status") != 0) || !file.exists(tmp))
    stop(paste(c("gsutil copy failed.", res), collapse = "\n"))
  data.table::fread(tmp, ...)
}

split_variant_fields <- function(df, locus_col = "locus", alleles_col = "alleles") {
  lp <- strsplit(as.character(df[[locus_col]]), ":", fixed = TRUE)
  df$chrom <- vapply(lp, `[`, character(1), 1); df$pos <- as.integer(vapply(lp, `[`, character(1), 2))
  ac <- gsub("^\\[|\\]$", "", as.character(df[[alleles_col]])); ac <- gsub('"', "", ac)
  ap <- strsplit(ac, ",", fixed = TRUE)
  df$ref <- vapply(ap, `[`, character(1), 1); df$alt <- vapply(ap, `[`, character(1), 2); df
}

annotate_domain_rmc <- function(var_info, domain_map, domain_anno, var_rmc) {
  v <- var_info %>% mutate(variant_id = paste(chrom, pos, ref, alt, sep = ':')) %>%
    left_join(domain_map %>% distinct(chrom, pos, ref, alt, aa_pos, enst, uniprot_isoform),
              by = c('chrom', 'pos', 'ref', 'alt'))
  dom_ref <- domain_anno %>% transmute(uniprot_isoform, aa_start, aa_end,
    family = dplyr::coalesce(dplyr::na_if(domain_name, ''), dplyr::na_if(interpro_id, ''),
                             dplyr::na_if(gene3d_id, ''), 'domain'),
    domain_label = paste0(family, ' (', aa_start, '-', aa_end, ')'))
  dom <- v %>% distinct(variant_id, uniprot_isoform, aa_pos) %>%
    filter(!is.na(uniprot_isoform), !is.na(aa_pos)) %>%
    inner_join(dom_ref, by = 'uniprot_isoform', relationship = 'many-to-many') %>%
    filter(aa_pos >= aa_start, aa_pos <= aa_end) %>% arrange(variant_id, aa_start) %>%
    group_by(variant_id) %>% summarize(domain_label = first(domain_label), .groups = 'drop')
  v %>% left_join(dom, by = 'variant_id') %>% left_join(var_rmc, by = c('chrom', 'pos', 'ref', 'alt')) %>%
    mutate(domain_label = if_else(is.na(domain_label), 'No domain', domain_label))
}

# =====================================================================
# Load inputs
# =====================================================================
message('Loading inputs ...')
domain_anno         <- read_tsv(paste0(data_dir, 'uniprot_superfamily_clean.tsv'), show_col_types = FALSE)
locus_domain_subset <- read_tsv(paste0(data_dir, 'missense_subset_grch38.tsv'), show_col_types = FALSE)
gene_constrain      <- read_tsv('~/Dropbox (Partners HealthCare)/0_very_often_used/gnomad.v4.1.constraint_metrics.tsv', show_col_types = FALSE)
rmc_metric          <- read_tsv(paste0(data_dir, 'gnomad_v4.1.1_all_mcrs.tsv'), show_col_types = FALSE)

var_rmc <- read_gcs_fast(RMC_VAR_URI) %>% tibble::as_tibble() %>% split_variant_fields() %>%
  transmute(chrom, pos, ref, alt, region_oe = rmc_region_oe, region_oe_p = rmc_region_oe_chisq_p,
            region_start = rmc_region_start, region_end = rmc_region_end,
            rmc_region = if_else(!is.na(rmc_region_start), paste0(rmc_region_start, '-', rmc_region_end), NA_character_))

raw <- read_pleiotropy_results('burden', '500k')
res_all <- modify_results_table(raw, 'burden', '500k')
results_500k <- res_all %>% mutate(annotation = factor(annotation, levels = annotation_types)) %>%
  filter(pvalue < STRICT_P, annotation == 'missense|LC')
sig_pairs <- results_500k %>% distinct(gene, pheno1, pheno2, description1, description2)

top_triplets <- read_delim(paste0(data_path, 'top_significant_triplets.txt.bgz'), delim = '\t',
                           col_types = cols(phenocode = col_character())) %>%
  mutate(coding = if_else(is.na(coding), '', coding),
         phenoname = paste0(trait_type, '_', phenocode, '_', pheno_sex, '_', coding, '_', modifier)) %>%
  split_variant_fields()

triplets_annot <- annotate_domain_rmc(
  top_triplets %>% filter(AF < 1e-4) %>%
    mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)),
  locus_domain_subset, domain_anno, var_rmc) %>%
  mutate(in_domain = domain_label != 'No domain',
         in_rmc = !is.na(region_oe_p) & region_oe_p < RMC_SIG)

# =====================================================================
# Panel A: gene-level constraint, sig vs rest (missense Z, RMC min O/E)
# =====================================================================
gene_sig <- res_all %>% filter(annotation == 'missense|LC') %>% group_by(gene) %>%
  summarize(minp = min(pvalue, na.rm = TRUE), .groups = 'drop') %>%
  mutate(group = factor(if_else(minp < STRICT_P, 'Sig', 'Rest'), levels = c('Rest', 'Sig')))
fna <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[[1]] else NA }
constr <- gene_constrain %>%
  mutate(across(c(canonical, mane_select), ~ as.character(.) %in% c('TRUE', 'true', 'True'))) %>%
  filter(mane_select | canonical) %>% group_by(gene) %>%
  summarize(mis_z = fna(`mis.z_score`), .groups = 'drop')
gene_enst <- gene_constrain %>% filter(grepl('^ENST', transcript)) %>% distinct(gene, transcript)
rmc_min <- rmc_metric %>% inner_join(gene_enst, by = 'transcript', relationship = 'many-to-many') %>%
  group_by(gene) %>% summarize(rmc_min_oe = min(region_oe, na.rm = TRUE), .groups = 'drop')
gc_dat <- gene_sig %>% left_join(constr, by = 'gene') %>% left_join(rmc_min, by = 'gene')

pa_dat <- bind_rows(
  gc_dat %>% transmute(group, score = mis_z, metric = 'Missense Z'),
  gc_dat %>% transmute(group, score = rmc_min_oe, metric = 'RMC min O/E')) %>%
  filter(!is.na(score))
n_lab_A <- pa_dat %>% group_by(metric, group) %>% summarize(y = min(score), n = n(), .groups = 'drop')
pA <- ggplot(pa_dat, aes(group, score, fill = group)) +
  geom_violin(alpha = .35, colour = NA, scale = 'width') +
  geom_boxplot(width = .16, outlier.shape = NA, alpha = .9) +
  geom_text(data = n_lab_A, aes(group, y, label = paste0(n, ' genes')), vjust = 1, size = 3.2, inherit.aes = FALSE) +
  facet_wrap(~ metric, scales = 'free_y') +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', comparisons = list(c('Rest', 'Sig')), size = 4, vjust = 0.2) +
  scale_fill_manual(values = fill_vals, guide = 'none') +
  labs(x = NULL, y = 'Gene constraint score') +
  base_theme + theme(axis.text.x = element_text(face = 'bold'))

# =====================================================================
# Panel B: enrichment of associated variants in domains / constrained regions
# =====================================================================
pair_var <- function(g, p1, p2) {
  d <- triplets_annot %>% filter(gene == g, phenoname %in% c(p1, p2), annotation == 'missense|LC') %>%
    pivot_wider(names_from = phenoname, values_from = c(BETA, Pvalue),
                id_cols = c(variant_id, gene, pos, aa_pos, domain_label, in_domain, in_rmc,
                            region_oe, region_oe_p, region_start, region_end))
  b1 <- paste0('BETA_', p1); b2 <- paste0('BETA_', p2); q1 <- paste0('Pvalue_', p1); q2 <- paste0('Pvalue_', p2)
  if (!all(c(b1, b2, q1, q2) %in% names(d))) return(NULL)
  d$beta1 <- d[[b1]]; d$beta2 <- d[[b2]]; d$p1 <- d[[q1]]; d$p2 <- d[[q2]]
  d %>% filter(!is.na(beta1), !is.na(beta2)) %>%
    mutate(assoc = (p1 <= NOMINAL) | (p2 <= NOMINAL), beta_max = pmax(abs(beta1), abs(beta2)), gene = g)
}
pt <- pmap(sig_pairs, function(gene, pheno1, pheno2, description1, description2) pair_var(gene, pheno1, pheno2))
pt <- pt[!vapply(pt, is.null, logical(1))]; allv <- bind_rows(pt)
gv <- allv %>% group_by(gene, variant_id) %>%
  summarize(region_oe = first(region_oe), in_rmc = any(in_rmc), in_domain = any(in_domain),
            assoc = any(assoc), beta_max = max(beta_max, na.rm = TRUE), .groups = 'drop')

or_lab <- function(flag) { ft <- fisher.test(table(gv$assoc, gv[[flag]]))
  sprintf('OR=%.2f, p=%.1g', unname(ft$estimate), ft$p.value) }
enr <- gv %>% group_by(assoc) %>%
  summarize(`In annotated domain` = mean(in_domain), `In constrained RMC region` = mean(in_rmc), .groups = 'drop') %>%
  pivot_longer(-assoc, names_to = 'feature', values_to = 'frac') %>%
  mutate(assoc = factor(if_else(assoc, 'Assoc.', 'Not assoc.'), levels = c('Not assoc.', 'Assoc.')))
or_ann <- tibble(feature = c('In annotated domain', 'In constrained RMC region'),
                 lab = c(or_lab('in_domain'), or_lab('in_rmc')))
n_lab_B <- tidyr::crossing(
  feature = c('In annotated domain', 'In constrained RMC region'),
  gv %>% count(assoc, name = 'n') %>%
    mutate(assoc = factor(if_else(assoc, 'Assoc.', 'Not assoc.'), levels = c('Not assoc.', 'Assoc.'))))
pB <- ggplot(enr, aes(feature, frac, fill = assoc)) +
  geom_col(position = position_dodge(.7), width = .65, alpha = .9) +
  geom_text(data = n_lab_B, aes(feature, y = 0.015, label = paste0('n=', n), group = assoc),
            position = position_dodge(.7), vjust = 0, size = 2.9, inherit.aes = FALSE) +
  geom_text(data = or_ann, aes(feature, y = 0.98, label = lab), inherit.aes = FALSE, size = 3.8, vjust = 1) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_fill_manual(name = NULL, values = fill_vals) +
  labs(x = NULL, y = 'Fraction of variants') +
  base_theme + theme(legend.position = 'top', axis.text.x = element_text(size = 10))

# =====================================================================
# Panel C: regional O/E vs variant effect magnitude
# =====================================================================
# Each point = one rare missense variant (deduped to gene x variant); coloured by whether
# it is nominally trait-associated (p<0.05 for either trait of any significant pair).
gvo <- gv %>% filter(!is.na(region_oe)) %>%
  mutate(assoc_lab = factor(if_else(assoc, 'Nominally assoc.', 'Not assoc.'),
                            levels = c('Not assoc.', 'Nominally assoc.')))
sp <- suppressWarnings(cor.test(gvo$region_oe, gvo$beta_max, method = 'spearman'))
pC <- ggplot(gvo, aes(region_oe, beta_max, colour = assoc_lab)) +
  geom_point(alpha = .3, size = .7) +
  geom_smooth(aes(group = 1), method = 'lm', se = FALSE, lty = 2, colour = 'black') +
  scale_colour_manual(name = NULL, values = c('Not assoc.' = neutral, 'Nominally assoc.' = accent)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
  labs(x = 'Regional missense O/E (lower = more constrained)', y = 'Variant |effect| (max across traits)',) +
  base_theme + theme(legend.position = 'top')

# =====================================================================
# Panel D: beta-beta scatter coloured by protein domain (variant-level)
# =====================================================================
# pick a specific gene-trait pair by trait descriptions (either order); else the first.
pick_pair <- function(gene_name, d1 = NULL, d2 = NULL) {
  prs <- sig_pairs %>% filter(gene == gene_name)
  if (!is.null(d1)) prs <- prs %>%
    filter((description1 == d1 & description2 == d2) | (description1 == d2 & description2 == d1))
  prs %>% slice(1)
}
beta_scatter <- function(gene_name, color_by, d1 = NULL, d2 = NULL) {
  pr <- pick_pair(gene_name, d1, d2)
  d <- pair_var(pr$gene, pr$pheno1, pr$pheno2)
  if (is.null(d) || nrow(d) == 0) return(ggplot() + theme_void())
  d <- d %>% mutate(sig = factor(case_when(
    p1 <= NOMINAL & p2 <= NOMINAL ~ 'Both', p1 > NOMINAL & p2 <= NOMINAL ~ pr$description2,
    p1 <= NOMINAL & p2 > NOMINAL ~ pr$description1, TRUE ~ 'None'),
    levels = c('Both', pr$description1, pr$description2, 'None')))
  base <- ggplot(d, aes(beta1, beta2)) +
    geom_vline(xintercept = 0, lty = 2, lwd = .25) + geom_hline(yintercept = 0, lty = 2, lwd = .25) +
    scale_shape_manual(name = paste0('p<', NOMINAL), values = c(16, 17, 15, 1), drop = FALSE) +
    labs(x = paste0('Effect on ', pr$description1, ' (β)'),
         y = paste0('Effect on ', pr$description2, ' (β)')) + base_theme +
    theme(legend.position = 'right', legend.box = 'vertical')
  if (color_by == 'domain') {
    lv <- sort(setdiff(unique(d$domain_label), 'No domain'))
    cols <- c(setNames(scales::hue_pal()(length(lv)), lv), 'No domain' = 'grey80')
    base + geom_point(aes(shape = sig, colour = domain_label), size = 2.2) +
      scale_colour_manual(name = 'Protein domain', values = cols, breaks = c(lv, 'No domain'), na.value = 'grey80') +
      guides(colour = guide_legend(ncol = 1, order = 1), shape = guide_legend(ncol = 1, order = 2))
  } else {
    base + geom_point(aes(shape = sig, colour = region_oe), size = 2.2) +
      scale_colour_viridis_c(name = 'Regional\nmissense O/E', direction = -1, na.value = 'grey80') +
      guides(colour = guide_colourbar(order = 1), shape = guide_legend(ncol = 1, order = 2)) +
      labs(title = paste0(gene_name, ' — regional missense O/E'))
  }
}
pD <- beta_scatter('ALB', 'domain')      # single-region example, domain colouring

# =====================================================================
# Panel E: pair-level effect heterogeneity (1 - |rho_w|) vs gene constraint
#   rho_w = AF-weighted cosine correlation ALLSPICE tests; departure = 1 - |rho_w|.
#   Tests whether missense tolerance explains cross-trait effect heterogeneity.
# =====================================================================
w_cos <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w); x <- x[ok]; y <- y[ok]; w <- w[ok]
  d <- sqrt(sum(w * x^2) * sum(w * y^2)); if (d == 0) NA_real_ else sum(w * x * y) / d
}
real_corr <- pmap_dfr(sig_pairs, function(gene, pheno1, pheno2, description1, description2) {
  v <- top_triplets %>% filter(gene == !!gene, AF < 1e-4, phenoname %in% c(pheno1, pheno2)) %>%
    mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
    filter(annotation == 'missense|LC') %>%
    pivot_wider(names_from = phenoname, values_from = BETA, id_cols = c(locus, alleles, AF), values_fn = mean)
  ok <- all(c(pheno1, pheno2) %in% names(v))
  if (ok) v <- v %>% filter(!is.na(.data[[pheno1]]), !is.na(.data[[pheno2]]))
  tibble(gene = gene, rho_w = if (ok && nrow(v) >= 2) w_cos(v[[pheno1]], v[[pheno2]], v$AF) else NA_real_)
}) %>% mutate(departure = 1 - abs(rho_w)) %>%
  left_join(constr, by = 'gene') %>% left_join(rmc_min, by = 'gene') %>% filter(!is.na(departure))
sp_z <- suppressWarnings(cor.test(real_corr$departure, real_corr$mis_z, method = 'spearman'))
sp_o <- suppressWarnings(cor.test(real_corr$departure, real_corr$rmc_min_oe, method = 'spearman'))
he_long <- bind_rows(
  real_corr %>% transmute(departure, x = mis_z,
    metric = sprintf('Gene missense Z (higher = constrained)\nSpearman rho = %.2f, p = %.2f', unname(sp_z$estimate), sp_z$p.value)),
  real_corr %>% transmute(departure, x = rmc_min_oe,
    metric = sprintf('Gene RMC min O/E (lower = constrained)\nSpearman rho = %.2f, p = %.2f', unname(sp_o$estimate), sp_o$p.value))) %>%
  filter(!is.na(x))
pE <- ggplot(he_long, aes(x, departure)) +
  geom_smooth(method = 'lm', se = FALSE, lty = 2, colour = 'grey55', linewidth = .6) +   # visual trend only
  geom_point(alpha = .8, size = 1.9, colour = accent) +
  facet_wrap(~ metric, scales = 'free_x', strip.position = 'bottom') +
  labs(x = NULL, y = 'Cross-trait heterogeneity\n(1 - |rho_w|)') +
  base_theme +
  theme(strip.placement = 'outside', strip.text = element_text(face = 'bold', size = 9.5))

# =====================================================================
# Panels F-H: lollipops (variant effect along protein; domain + RMC context)
# =====================================================================
lollipop <- function(gene_name, d1 = NULL, d2 = NULL) {
  pr <- pick_pair(gene_name, d1, d2)
  d <- pair_var(pr$gene, pr$pheno1, pr$pheno2)
  if (is.null(d) || nrow(d) == 0) return(ggplot() + theme_void() + labs(title = paste0(gene_name, ' (no data)')))
  xvar <- if (sum(!is.na(d$aa_pos)) >= 0.5 * nrow(d)) 'aa_pos' else 'pos'
  long <- bind_rows(
    d %>% transmute(x = .data[[xvar]], beta = beta1, p = p1, region_oe, trait = pr$description1),
    d %>% transmute(x = .data[[xvar]], beta = beta2, p = p2, region_oe, trait = pr$description2)) %>% filter(!is.na(x))
  if (xvar == 'aa_pos') {
    reg <- d %>% filter(!is.na(region_oe), !is.na(aa_pos)) %>%
      group_by(region_start, region_end, region_oe) %>% summarize(xmin = min(aa_pos), xmax = max(aa_pos), .groups = 'drop')
  } else {
    reg <- d %>% filter(!is.na(region_oe)) %>% distinct(region_start, region_end, region_oe) %>%
      mutate(xmin = region_start, xmax = region_end)
  }
  ggplot(long, aes(x, beta)) +
    { if (nrow(reg)) geom_rect(data = reg, inherit.aes = FALSE,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = region_oe), alpha = .5) } +
    scale_fill_viridis_c(name = 'Region O/E', direction = -1, na.value = 'grey90') +
    geom_hline(yintercept = 0, lty = 2, lwd = .25) +
    geom_segment(aes(xend = x, yend = 0, colour = p <= NOMINAL), lwd = .35) +
    geom_point(aes(colour = p <= NOMINAL), size = 1.4) +
    scale_x_continuous(label = comma) +
    scale_colour_manual(name = paste0('p<', NOMINAL), values = c('grey60', '#D62728')) +
    facet_wrap(~ trait, ncol = 1, scales = 'free_y') +
    labs(x = if (xvar == 'aa_pos') 'Amino-acid position' else 'Genomic position', y = 'Variant effect (β)') +
    base_theme + theme(strip.text = element_text(face = 'bold', size = 11), legend.position = 'right')
}
# ALB/ALPL are single-region RMC (uniform O/E) -> not informative as lollipops; show only
# the multi-region example (HECTD4), where regional O/E varies within the gene.
pF <- lollipop('HECTD4', 'Reticulocyte count', 'Mean sphered cell volume')

# =====================================================================
# Supplementary table: per gene-trait pair summary behind Figure S29
#   ALLSPICE stats + rare-missense counts + domain/RMC enrichment (B) +
#   effect heterogeneity 1-|rho_w| (E) + gene-level constraint (A/E).
# =====================================================================
rho_tab <- pmap_dfr(sig_pairs, function(gene, pheno1, pheno2, description1, description2) {
  v <- top_triplets %>% filter(gene == !!gene, AF < 1e-4, phenoname %in% c(pheno1, pheno2)) %>%
    mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
    filter(annotation == 'missense|LC') %>%
    pivot_wider(names_from = phenoname, values_from = BETA, id_cols = c(locus, alleles, AF), values_fn = mean)
  ok <- all(c(pheno1, pheno2) %in% names(v))
  if (ok) v <- v %>% filter(!is.na(.data[[pheno1]]), !is.na(.data[[pheno2]]))
  rw <- if (ok && nrow(v) >= 2) w_cos(v[[pheno1]], v[[pheno2]], v$AF) else NA_real_
  tibble(gene, pheno1, pheno2, rho_w = rw, departure = 1 - abs(rw))
})
enr_tab <- pmap_dfr(sig_pairs, function(gene, pheno1, pheno2, description1, description2) {
  d <- pair_var(gene, pheno1, pheno2); n <- if (is.null(d)) 0L else nrow(d)
  tibble(gene, pheno1, pheno2, n_missense_var = n,
         n_nominally_assoc = if (n) sum(d$assoc) else 0L,
         pct_in_domain = if (n) 100 * mean(d$in_domain) else NA_real_,
         pct_in_constrained_region = if (n) 100 * mean(d$in_rmc) else NA_real_)
})
supp_tab <- sig_pairs %>%
  left_join(results_500k %>% distinct(gene, pheno1, pheno2, allspice_pvalue = pvalue, c_hat, n_var_test = n_var),
            by = c('gene', 'pheno1', 'pheno2')) %>%
  left_join(rho_tab, by = c('gene', 'pheno1', 'pheno2')) %>%
  left_join(enr_tab, by = c('gene', 'pheno1', 'pheno2')) %>%
  left_join(constr, by = 'gene') %>% left_join(rmc_min, by = 'gene') %>%
  transmute(gene, trait1 = description1, trait2 = description2,
            allspice_pvalue, c_hat, n_var_test,
            n_missense_var, n_nominally_assoc,
            pct_in_domain = round(pct_in_domain, 1), pct_in_constrained_region = round(pct_in_constrained_region, 1),
            rho_w = round(rho_w, 3), departure_1_minus_abs_rho_w = round(departure, 3),
            gene_missense_z = round(mis_z, 2), gene_rmc_min_oe = round(rmc_min_oe, 3)) %>%
  arrange(allspice_pvalue)
write_csv(supp_tab, paste0(data_dir, 'tableS_figureS29_per_pair.csv'))
cat('wrote', paste0(data_dir, 'tableS_figureS29_per_pair.csv'), '(', nrow(supp_tab), 'pairs )\n')

# =====================================================================
# Assemble
# =====================================================================
row1 <- ggpubr::ggarrange(pA, pB, pC, ncol = 3, labels = c('(A) Gene level — constraint of significant vs other genes',
                                                           '(B) Variant level — enrichment of associated variants',
                                                           '(C) Variant level — effect magnitude vs regional O/E'), font.label = list(size = 13), hjust = 0, vjust =1)
row2 <- ggpubr::ggarrange(pD, pE, ncol = 2, labels = c('(D) Variant level — ALB effects by protein domain',
                                                       '(E) Pair level — heterogeneity vs gene constraint'), font.label = list(size = 13), hjust = 0, vjust =1)
row3 <- ggpubr::ggarrange(pF, ncol = 1, labels = '(F) Variant level — HECTD4 effects over multi-region RMC (genomic position)', font.label = list(size = 13), hjust = 0, vjust =1)
figS29 <- ggpubr::ggarrange(row1, row2, row3, nrow = 3, heights = c(1, 1.3, 1.15))

ggsave(paste0(fig_out, 'figureS29.png'), figS29, height = 12, width = 15, dpi = 300, limitsize = FALSE, bg = 'white')
ggsave(paste0(fig_out, 'figureS29.pdf'), figS29, height = 12, width = 15, limitsize = FALSE)
cat('wrote', paste0(fig_out, 'figureS29.png / .pdf'), '\n')

