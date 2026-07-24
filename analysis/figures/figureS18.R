source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/constants.R')

raw_results_500k <- read_csv('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/ALLSPICE/continuous_ALL_AC_5_burden_syn_var_corr_500k_2024_results_corr.csv')
results_500k <- modify_results_table(raw_results_500k, 'burden', '500k')
results_500k %>%
  group_by(annotation) %>%
  dplyr::summarize(sig05 = sum(pvalue < 0.05))
sig_results_500k <- results_500k %>% filter(pvalue < 4.24e-6)

## --- Genome-wide strict pairs (corr < 0.8): ultra-rare AC<=5 vs the AF<1e-4 primary.
## This is the "52 of 98 strictly-significant pairs persisted at AC<=5" statement
## (SI / response: restricting to ultra-rare variants weakened but did not abolish
## the genome-wide signal). Both counts use the study-wide threshold p < 4.23e-6.
ac5_strict_corr08  <- nrow(results_500k %>% filter(pvalue < 4.23e-6, corr < 0.8))
primary_results    <- modify_results_table(
  read_csv(paste0(result_path, 'continuous_final_AF_1e_4_burden_syn_var_corr_500k_2023_results.csv')),
  'burden', '500k')
primary_strict_corr08 <- nrow(primary_results %>% filter(pvalue < 4.23e-6, corr < 0.8))
cat(sprintf('Genome-wide strict pairs (corr<0.8): AC<=5 = %d ; AF<1e-4 primary = %d  (%d of %d persist at AC<=5)\n',
            ac5_strict_corr08, primary_strict_corr08, ac5_strict_corr08, primary_strict_corr08))


figure_real_data_qq <- function(results, name=NULL, save=TRUE){
  results <- results %>%
    group_by(annotation) %>%
    arrange(pvalue) %>%
    add_count() %>%
    mutate(observed = -log10(pvalue),
           rank = order(pvalue),
           expected = -(log10(rank / (n+1))),
           annotation = factor(annotation, levels = annotation_types))
  mx <- as.numeric(max(max(results[results$pvalue>0, 'observed']),
                       max(results[results$pvalue>0, 'expected'])))

  figure <- results %>%
    ggplot+ aes(y=observed,x=expected, color = annotation, label = gene) +
    geom_point(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1) +
    geom_hline(yintercept = -log10(4.23e-6), lty=2) +
    labs(x=expression(Expected -log[10](p)), y=expression(Observed -log[10](p)), color = 'Annotation') +
    annotation_color_scale + annotation_fill_scale +
    xlim(0, mx) +
    ylim(0, mx) +
    themes +
    facet_grid(~annotation, labeller = label_type)
  # geom_text_repel(max.overlaps = 2)

  if(save){
    png(paste0(figure_path, name, "_qqplot.png"), width=5, height=3.5, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

p1 <- figure_real_data_qq(results_500k %>% filter(sig_gene == 2 & corr < 0.9), save = F)
p2 <- results_500k %>% filter(sig_gene == 2 & corr < 0.9)%>%
  mutate(annotation = factor(annotation, levels = annotation_types)) %>%
  ggplot +
  aes(x = corr, y = -log10(pvalue), color = annotation, size = n_var) +
  labs(x = 'Phenotypic correlation', y = expression(-log[10](p)), size = 'Number of variants', color = 'Annotation') +
  geom_point(alpha = 0.5) +
  annotation_color_scale +
  # scale_y_log10() +
  geom_hline(yintercept = -log10(0.05/nrow(results_500k)), lty = 2) +
  # geom_vline(xintercept = c(-0.1, 0.1, 0.8, 1), lty = 2) +
  scale_size_continuous(range = c(0.01, 4)) +
  facet_grid(~annotation, labeller = labeller(annotation = annotation_names)) + themes + theme(legend.position = 'top')

p3 <- results_500k %>%
  # filter(annotation != 'synonymous') %>%
  mutate(sig = if_else(pvalue < 4.24e-6, 'Strictly significant', if_else(pvalue < 0.05, 'Nominally significant', 'Not significant')),
         annotation = factor(annotation, levels = annotation_types)) %>%
  mutate(sig = factor(sig, levels = c('Not significant', 'Nominally significant', 'Strictly significant'))) %>%
  ggplot + aes(x = sig, y = n_var, color = annotation) +
  scale_y_log10(label = comma) + labs(y = 'Number of variants', x = NULL) +
  geom_boxplot(width = 0.2, size = 0.75) + annotation_color_scale +
  facet_wrap(~annotation, labeller=label_type, nrow=1, scale = 'free') +
  guides(color = 'none') + themes

pLoF_sig <- format_sig_result_matrix(sig_results_500k %>% filter(corr < 0.8), annot='pLoF')
p41 <- plot_sig_result_matrix(pLoF_sig, annot='pLoF')

mis_sig <- format_sig_result_matrix(sig_results_500k %>% filter(corr < 0.8), annot='missense|LC')
p42 <- plot_sig_result_matrix(mis_sig, annot='missense|LC')

p4 <- ggpubr::ggarrange(p41, p42, nrow=1, hjust = 0, align = 'h', widths = c(1.5, 1.2),
                             font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')
)

figure <- ggpubr::ggarrange(p1 +
                              theme(axis.title = element_text(face = 'plain', size = 11),
                                    plot.margin = unit(c(1,0,0,0.5), "cm"), legend.position = 'none'),
                            p2 + guides(color = "none") +
                              theme(axis.title = element_text(face = 'plain', size = 11),
                                    plot.margin = unit(c(0.7,0,0,0.5), "cm")),
                            p3 +
                              theme(axis.title = element_text(face = 'plain', size = 11),
                                    plot.margin = unit(c(1,0,0,0.5), "cm"), legend.position = 'none'),
                            p4 +
                              theme(axis.title = element_text(face = 'plain', size = 11),
                                    plot.margin = unit(c(1,0,0,0.5), "cm"), legend.position = 'none'),
                            labels = c('(A) QQ plots of ALLSPICE test results across high-quality phenotypes',
                                       '(B) Relationship between phenotypic correlation and ALLSPICE p-value',
                                       '(C) Number of variants in triplets across significance levels',
                                       '(D) Triplets strictly significant in ALLSPICE test (p-value < 4.24e-6)'
                                       ),
                            ncol=1, vjust = 2, hjust = 0, font.label = list(size = 10, color = "black", face = "bold", family = NULL), heights = c(0.18, 0.2, 0.18, 0.4))

png(paste0(figure_path,'figureS18.png'), height = 15, width = 12, units = 'in', res = 300)
print(figure)
dev.off()
