source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/constants.R')
library(magick)

# Figure 3A: Gene-annotation pair and pleiotropy scheme
figure3A <- ggdraw() +
  draw_image(image_read(paste0('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/figure_2024/figure3/figure3_scheme_1.png'))) +
  theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm"))

# Figure 3B: ALLSPICE triplet scheme
figure3B <- ggdraw() +
  draw_image(image_read(paste0('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/figure_2024/figure3/figure3_scheme_2.png'))) +
  theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm"))

# Combine
figure <- ggpubr::ggarrange(figure1A, figure2A,
                            labels = c('(A) Definitions of gene-annotation pair',
                                       '(B) Schematic overview of association pairs, gene-level horizontal and vertical pleiotropy'),
                            ncol = 1, hjust = 0, vjust = 2,
                            font.label = list(size = 10, color = "black", face = "bold", family = NULL),
                            heights = c(1, 1))
png(paste0(figure_path, 'figure3.png'), height = 6, width = 8, units = 'in', res = 300)
print(figure)
dev.off()
