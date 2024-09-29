source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')
library(magick)

raw_results_500k <- read_pleiotropy_results('burden', '500k') 
results_500k <- modify_results_table(raw_results_500k, 'burden', '500k')

top_triplets <- read_delim(paste0(data_path, "top_significant_triplets.txt.bgz"), delim='\t', col_types = cols(phenocode = col_character())) %>% 
  mutate(coding = if_else(is.na(coding), '', coding)) %>%
  mutate(phenoname = paste0(trait_type, '_', phenocode, '_',  pheno_sex, '_',  coding, '_',  modifier)) %>%
  mutate(BETA_adjusted = sqrt(2*AF*(1-AF))*BETA)

## Protein figures 
library(cowplot)
# p2 <- ggplot_pdf(image_read(paste0(figure_path, 'figure3_S17/ukb_pleiotropy_ALB_variants.1AO6.by_ca.colored.png'))) +
#   theme(plot.margin = unit(c(0.7,0,0,0), "cm")) 

p1 <- ggdraw() + 
  draw_image(image_read(paste0(figure_path, 'figure3_S17/ukb_pleiotropy_ALB_variants.1AO6.by_ca.colored.png'))) + 
  theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm")) + 
  draw_label("Calcium\nbinding\nsites", x = 0.15, y = 0.53, hjust = 0.5, vjust = 0.5, color = 'orange', size = 8) +
  draw_line(
    x = c(0.2, 0.3), y = c(0.48, 0.43),  # Adjust these coordinates as needed
    color = 'orange',
    arrow = arrow(length = unit(0.03, "npc"))
  )
p1

p2 <- ggplot_pdf(image_read_pdf(paste0(figure_path, 'figure3_S17/ukb_pleiotropy_ALB_variants.P02768.by_phe_p0.05_0exl1.3d_dist_ca.pdf')))+
  theme(plot.margin = unit(c(0.7,0,0,0), "cm")) 
# p4 <- ggplot_pdf(image_read(paste0(figure_path, 'figure3_S17/ukb_pleiotropy_ALB_variants.1AO6.by_ca_betasign.colored.png')))+
#   theme(plot.margin = unit(c(0.7,0,0,0), "cm"))

p3 <- ggdraw() + 
  draw_image(image_read(paste0(figure_path, 'figure3_S17/ukb_pleiotropy_ALB_variants.1AO6.by_ca_betasign.colored.png'))) + 
  theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm")) + 
  draw_label("Calcium\nbinding\nsites", x = 0.15, y = 0.53, hjust = 0.5, vjust = 0.5, color = 'orange', size = 8) +
  draw_line(
    x = c(0.2, 0.3), y = c(0.48, 0.43),  # Adjust these coordinates as needed
    color = 'orange',
    arrow = arrow(length = unit(0.03, "npc"))
  )
p3
p4 <- ggplot_pdf(image_read_pdf(paste0(figure_path, 'figure3_S17/ukb_pleiotropy_ALB_variants.P02768.by_phe30680_betasign_pmin0.05_exl.3d_dist_ca.legend.pdf')))+
  theme(plot.margin = unit(c(0.7,0,0,0), "cm"))

figure = ggpubr::ggarrange(ggpubr::ggarrange(p1, p2, 
                                             labels = c('(A) ALB missense variants\n      on protein structure 1AO6\n      (by pvalue)', 
                                                        '(B) Distance between ALB missense variants\n      and calcium binding sites (by pvalue)'), 
                                             nrow=1, hjust=0, widths = c(1,1.8), vjust = 1,
                                             font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')), 
                           ggpubr::ggarrange(p3, p4, 
                                             labels = c('(C) ALB missense variants\n      on protein structure 1AO6\n      (by effect direction)', 
                                                        '(D) Distance between ALB missense variants\n      and calcium binding sites (by effect direction)'), 
                                             nrow=1, hjust=0, widths = c(1,1.8), vjust = 1,
                                             font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')), 
                           nrow=2, heights = c( 0.16, 0.15), labels = c('', ''), hjust = 0, 
                           font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')
                           )
png(paste0(figure_path,'figure3.png'), height = 5.5, width = 7.5, units = 'in', res = 300)
print(figure)
dev.off()



## Protein figures 
library(cowplot)
# p2 <- ggplot_pdf(image_read(paste0(figure_path, 'figure3_S17/ukb_pleiotropy_ALB_variants.1AO6.by_ca.colored.png'))) +
#   theme(plot.margin = unit(c(0.7,0,0,0), "cm")) 

p1 <- ggdraw() + 
  draw_image(image_read(paste0(figure_path, 'figure3_S17/ukb_pleiotropy_ALB_variants_ac5.1AO6.by_phe_p0.05_0exl1.colored.png'))) + 
  theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm")) + 
  draw_label("Calcium\nbinding\nsites", x = 0.15, y = 0.53, hjust = 0.5, vjust = 0.5, color = 'orange', size = 8) +
  draw_line(
    x = c(0.2, 0.3), y = c(0.48, 0.43),  # Adjust these coordinates as needed
    color = 'orange',
    arrow = arrow(length = unit(0.03, "npc"))
  )
p1

p2 <- ggplot_pdf(image_read_pdf(paste0(figure_path, 'figure3_S17/ukb_pleiotropy_ALB_variants_ac5.P02768.by_phe_p0.05_0exl1.3d_dist_ca.pdf')))+
  theme(plot.margin = unit(c(0.7,0,0,0), "cm")) 
# p4 <- ggplot_pdf(image_read(paste0(figure_path, 'figure3_S17/ukb_pleiotropy_ALB_variants.1AO6.by_ca_betasign.colored.png')))+
#   theme(plot.margin = unit(c(0.7,0,0,0), "cm"))

p3 <- ggdraw() + 
  draw_image(image_read(paste0(figure_path, 'figure3_S17/ukb_pleiotropy_ALB_variants_ac5.1AO6.by_phe_betasign_pmin0.05_exl.colored.png'))) + 
  theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm")) + 
  draw_label("Calcium\nbinding\nsites", x = 0.15, y = 0.53, hjust = 0.5, vjust = 0.5, color = 'orange', size = 8) +
  draw_line(
    x = c(0.2, 0.3), y = c(0.48, 0.43),  # Adjust these coordinates as needed
    color = 'orange',
    arrow = arrow(length = unit(0.03, "npc"))
  )
p3
p4 <- ggplot_pdf(image_read_pdf(paste0(figure_path, 'figure3_S17/ukb_pleiotropy_ALB_variants_ac5.P02768.by_phe30680_betasign_pmin0.05_exl.3d_dist_ca.legend.pdf')))+
  theme(plot.margin = unit(c(0.7,0,0,0), "cm"))

figure = ggpubr::ggarrange(ggpubr::ggarrange(p1, p2, 
                                             labels = c('(A) ALB missense variants\n      on protein structure 1AO6\n      (by pvalue)', 
                                                        '(B) Distance between ALB missense variants\n      and calcium binding sites (by pvalue)'), 
                                             nrow=1, hjust=0, widths = c(1,1.8), vjust = 1,
                                             font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')), 
                           ggpubr::ggarrange(p3, p4, 
                                             labels = c('(C) ALB missense variants\n      on protein structure 1AO6\n      (by effect direction)', 
                                                        '(D) Distance between ALB missense variants\n      and calcium binding sites (by effect direction)'), 
                                             nrow=1, hjust=0, widths = c(1,1.8), vjust = 1,
                                             font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')), 
                           nrow=2, heights = c( 0.16, 0.15), labels = c('', ''), hjust = 0, 
                           font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')
)
png(paste0(figure_path,'figureS17.png'), height = 5.5, width = 7.5, units = 'in', res = 300)
print(figure)
dev.off()
