library(ggplot2)

set.seed(1234)

n <- 50

## Core cluster near zero
x_core <- rnorm(n, mean = 0, sd = 0.15)
y_core <- rnorm(n, mean = 0, sd = 0.15)

## Add a few non-null / outlier effects
k <- 6
x_out <- rnorm(k, mean = 0.2, sd = 0.3)
y_out <- rnorm(k, mean = 0.4, sd = 1)

## Combine
df <- data.frame(
  trait1 = c(x_core, x_out),
  trait2 = c(y_core, y_out)
)

## Plot
p <- ggplot(df, aes(trait1, trait2)) +
  geom_point(size = 2.5, color = color_lof) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  # coord_cartesian(xlim = c(-0.6, 0.8), ylim = c(-0.6, 1.2)) +
  labs(
    x = "Variant effect size (Trait A)",
    y = "Variant effect size (Trait B)",
    title = NULL
  ) +
  theme_classic(base_size = 14)
png(paste0(figure_path,'figure2/horizontal_example.png'), height =3, width = 6, units = 'in', res = 300)
print(p)
dev.off()


set.seed(456)
n <- 50

## Trait 1 effects: mostly small, a few larger
x_core <- rnorm(n - 5, mean = 0, sd = 0.15)
x_out  <- rnorm(5, mean = 0.4, sd = 0.5)
x <- c(x_core, x_out)

## Linear relationship + noise (null: proportional effects)
beta <- 1.2
y <- beta * x + rnorm(n, mean = 0, sd = 0.1)

df <- data.frame(
  trait1 = x,
  trait2 = y
)

p <- ggplot(df, aes(trait1, trait2)) +
  geom_point(size = 2.5, color = color_lof) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  # coord_cartesian(xlim = c(-0.4, 0.9), ylim = c(-0.4, 1.1)) +
  labs(
    x = "Variant effect size (Trait A)",
    y = "Variant effect size (Trait B)",
    title = NULL
  ) +
  theme_classic(base_size = 14)
png(paste0(figure_path,'figure2/vertical_example.png'), height =3, width = 6, units = 'in', res = 300)
print(p)
dev.off()



## Make main figure2
source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/constants.R')
library(magick)

figureA <- ggdraw() +
  draw_image(image_read(paste0('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/figure_2024/figure2/figure2_scheme.png')))+
  theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm"))

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

pLoF_sig <- format_sig_result_matrix(sig_results_500k, annot='pLoF')
figureC1 <- plot_sig_result_matrix(pLoF_sig, annot='pLoF')

mis_sig <- format_sig_result_matrix(sig_results_500k, annot='missense|LC')
figureC2 <- plot_sig_result_matrix(mis_sig, annot='missense|LC')

figureC <- ggpubr::ggarrange(figureC1, figureC2, ncol=1, hjust = 0, align = 'v', heights = c(1, 1),
                             font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')
)
# figureB
# png(paste0(figure_path,'figure2/figure2B.png'), height = 5, width =10, units = 'in', res = 300)
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
# png(paste0(figure_path,'figure2.png'), height = 10, width =12, units = 'in', res = 300)
# print(figure3)
# dev.off()


figure2 <-ggpubr::ggarrange(ggpubr::ggarrange(figureA + theme(plot.margin = unit(c(0.5, 0, 0, 0), 'cm')),
                             figureC1 + theme(plot.margin = unit(c(0.5, 1, 0, 1), 'cm')), ncol=1, hjust = 0,
                             heights = c(1, 1.8), vjust = 1.2,
                             labels = c('(A) Schematic overview of triplet, gene-level horizontal and vertical pleiotropy',
                                        '(B) Triplets strictly significant in ALLSPICE test (pLoF; p-value < 4.23e-6)'),
                             font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')
                             ), figureC2 + theme(plot.margin = unit(c(0.5, 0, 0, 0), 'cm')),
                            nrow=1, hjust = 0,  widths = c(1, 1), vjust = 1.2,
                            labels = c('', '(C) Triplets strictly significant in ALLSPICE test (Missense; p-value < 4.23e-6)'),
                            font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'))
figure2
png(paste0(figure_path,'figure2.png'), height = 7, width =13, units = 'in', res = 300)
print(figure2)
dev.off()
