## Make main figure4
source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/constants.R')
library(magick)

gene_data <- read.csv(paste0(data_path, 'gene_phewas_burden_sig_count_239.csv'), sep = '\t') %>%
  mutate(annotation = factor(annotation, levels = annotation_types2[c(1,2,4,3)]),
         interval = get_freq_interval(CAF))

gene_name_label <- gene_data %>%
  filter(n_phewas_sig > 10) %>%
  filter(annotation != 'pLoF|missense|LC') %>%
  mutate(annotation = factor(annotation, levels = annotation_types2[c(1,2,4,3)]))

label_type = labeller(annotation = annotation_names2)
annotation_fill_scale2 = scale_fill_manual(name = 'Annotation', values = colors2, breaks = annotation_types2[c(1,2,4,3)], labels = annotation_names2[c(1,2,4,3)])

figureA <- gene_data %>%
  filter(n_phewas_sig > 1) %>%
  ggplot + aes(x = n_phewas_sig, fill=annotation) +
  geom_histogram(binwidth = 1, position = 'dodge', stat ='count', color='white')  +
  annotation_color_scale2 +
  annotation_fill_scale2 +
  themes + theme_classic() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16)) +
  theme(legend.position = 'top',
        legend.direction = 'horizontal') +
    geom_text(data = gene_name_label, aes(x = n_phewas_sig, y = 3, label = gene_symbol, color = annotation), size = 3, show.legend = FALSE)  +
  geom_text(data = gene_data %>% filter(n_phewas_sig >= 1) %>% group_by(annotation) %>% dplyr::summarize(n_pleiotropy = sum(n_phewas_sig>1), n = n(), p = sum(n_phewas_sig>1)/n()),
                  aes(x = Inf, y = 70, color = annotation,label = paste0('Pleiotropic genes (%):\n', n_pleiotropy, '/', n, '=', round(p*100, 2), '%')), hjust = 1, show.legend = FALSE) +
  labs(x = 'N associations', y = 'N genes', color = NULL, fill= NULL) +
  facet_grid(~annotation, labeller = label_type)+
  theme(plot.margin = unit(c(0.5,1,0.5,0.2), "cm"))

raw_results_500k <- read_pleiotropy_results('burden', '500k')
results_500k <- modify_results_table(raw_results_500k, 'burden', '500k') %>%
  mutate(sig = if_else(pvalue < 4.23e-6, 'Strictly significant', if_else(pvalue < 0.05, 'Nominally significant', 'Not significant')),
         annotation = factor(annotation, levels = annotation_types)) %>%
  mutate(sig = factor(sig, levels = c('Not significant', 'Nominally significant', 'Strictly significant')))
sig_results_500k <- results_500k %>% filter(pvalue < 4.23e-6)

# group_sum <- results_500k %>%
#   mutate(sig = pvalue < 4.23e-6,
#          annotation = factor(annotation, levels = annotation_types)) %>%
#   group_by(annotation) %>%
#   dplyr::summarize(
#     nvar = mean(n_var),
#     sig_total = sum(sig),
#     non_sig = sum(!sig),
#     prop = sum(sig)/n()
#   )
#
# pairwise.t.test(results_500k[, 'n_var'], results_500k[, 'sig'],
#                 p.adjust.method = "BH",
#                 pool.sd = FALSE)
#
# pairwise.t.test(results_500k[results_500k$annotation == 'pLoF', 'n_var'], results_500k[results_500k$annotation == 'pLoF', 'sig'],
#                 p.adjust.method = "BH",
#                 pool.sd = FALSE)
#
# pairwise.t.test(results_500k[results_500k$annotation == 'missense|LC', 'n_var'], results_500k[results_500k$annotation == 'missense|LC', 'sig'],
#                 p.adjust.method = "BH",
#                 pool.sd = FALSE)
cor_test <- cor.test(results_500k$n_var, results_500k$pvalue)
cor_test$p.value


# prop.test(
#   x = group_sum$sig_total,  # number of "successes"
#   n = group_sum$total       # total per group
# )

# pairwise.prop.test(
#   x = group_sum$sig_total,
#   n = group_sum$total,
#   p.adjust.method = "bonferroni"
# )
#
# counts <- cbind(
#   Significant = group_sum$sig_total,
#   NotSignificant = group_sum$total - group_sum$sig_total
# )
#
# rownames(counts) <- group_sum$annotation
#
# chisq.test(counts)

figureB <- results_500k %>%
  filter(annotation != 'synonymous')  %>%
  ggplot + aes(x = sig, y = n_var, color = annotation) +
  scale_y_log10(label = comma) + labs(y = 'Number of variants', x = NULL) +
  geom_boxplot(width = 0.2, size = 0.75) + annotation_color_scale +
  facet_wrap(~annotation, labeller=label_type, ncol=1, scale = 'free') +
  guides(color = 'none')

pLoF_sig <- format_sig_result_matrix(sig_results_500k %>% filter(corr < 0.8), annot='pLoF')
figureC1 <- plot_sig_result_matrix(pLoF_sig, annot='pLoF')

mis_sig <- format_sig_result_matrix(sig_results_500k %>% filter(corr < 0.8), annot='missense|LC')
figureC2 <- plot_sig_result_matrix(mis_sig, annot='missense|LC')

figureC <- ggpubr::ggarrange(figureC1, figureC2, ncol=1, hjust = 0, align = 'v', heights = c(1, 1),
                             font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')
)
# figureB
# png(paste0(figure_path,'figure4/figure4B.png'), height = 5, width =10, units = 'in', res = 300)
# print(figureB)
# dev.off()
#

# figure3 <- ggpubr::ggarrange(ggpubr::ggarrange(figureA + theme(plot.margin = unit(c(0.5, 1, 0, 0.5), 'cm')),
#                                                figureB + theme(plot.margin = unit(c(0.8, 0.3, 0, 0), 'cm')),
#                                                nrow=1, hjust = 0, widths = c(2, 1), vjust = 1.2,
#                                                labels = c('(A) Schematic overview of triplet, gene-level horizontal and vertical pleiotropy',
#                                                           '(B) Number of variants in triplets across significance levels'),
#                                                font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')) ,
#                              figureC + theme(plot.margin = unit(c(0.5, 0, 0, 0), 'cm')), ncol=1, hjust = 0, heights = c(1, 1.8), vjust = 1,
#                              labels = c('',
#                                         '(C) Triplets strictly significant in ALLSPICE test (p-value < 4.23e-6)'),
#                              font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')
# )
# figure3
# png(paste0(figure_path,'figure4.png'), height = 10, width =12, units = 'in', res = 300)
# print(figure3)
# dev.off()


bottom_row <- ggpubr::ggarrange(figureC1 + theme(plot.margin = unit(c(0.5, 1, 0, 1), 'cm')),
                                figureC2 + theme(plot.margin = unit(c(0.5, 1, 0, 1), 'cm')),
                                nrow = 1, hjust = 0, vjust = 1.2,
                                labels = c('(B) Pairs of associations strictly significant in ALLSPICE test (pLoF; p-value < 4.23e-6)',
                                           '(C) Pairs of associations strictly significant in ALLSPICE test (Missense; p-value < 4.23e-6)'),
                                font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'))

figure4 <- ggpubr::ggarrange(figureA + theme(plot.margin = unit(c(0.5, 0, 0, 0), 'cm')),
                              bottom_row,
                              ncol = 1, hjust = 0, vjust = 1.2, heights = c(1, 2.2),
                              labels = c('(A) Number of associations per gene', ''),
                              font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'))
figure4
png(paste0(figure_path,'figure4.png'), height = 10, width =13, units = 'in', res = 300)
print(figure4)
dev.off()
