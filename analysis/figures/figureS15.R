source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')
library(magick)
library(emojifont)

raw_results_500k <- read_pleiotropy_results('burden', '500k') 
results_500k <- modify_results_table(raw_results_500k, 'burden', '500k')

p1 <- results_500k %>%
  filter(gene == 'ALPL') %>%
  select(description1=description2, description2=description1, corr) %>%
  distinct() %>%
  rbind(., results_500k %>%
          filter(gene == 'ALPL') %>%
          select(description1, description2, corr) %>%
          distinct()) %>%
  
  ggplot() +
  geom_tile(aes(x = description2, y = description1, fill = corr)) +
  geom_text(aes(x = description2, y = description1, label = round(corr,2)), size=4) +
  coord_equal() +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-0.3,0.3)) +
  labs(x =NULL, y = NULL) + theme_classic() + 
  theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") )+
  themes + 
  theme_classic()  + 
  theme(legend.title = element_text(face = 'plain', size = 11), 
        axis.title = element_text(face = 'plain', size = 11), 
        plot.margin = unit(c(0.5,0,0,0), "cm"),
        axis.text.x = element_text(angle = 15, hjust =1))

p2 <- ggplot_pdf(image_read(paste0(figure_path, 'figureS14/ukb_pleiotropy_ALPL_variants.7YIV.by_phe_p0.05_exl.colored.png')))+
  theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
p3 <- ggplot_pdf(image_read(paste0(figure_path, 'figureS14/ukb_pleiotropy_ALPL_variants.7YIV.by_phe_betasign_exl.colored.png')))+
  theme(plot.margin = unit(c(0.5,0,0,0), "cm"))



figure = ggpubr::ggarrange(p1, 
                           ggpubr::ggarrange(p2, p3, 
                                             labels = c('(B) ALPL missense variants on protein structure 7YIV\n      (by pvalue)', 
                                                        '(C) ALPL missense variants on protein structure 7YIV\n      (by effect direction)'), 
                                             nrow=1, hjust=0, widths = c(1, 1), 
                                             font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')), 
                           nrow=2, heights = c(0.1, 0.12), labels = c('(A) ALPL pleiotropy domain', ''), hjust = 0, 
                           font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')
)
png(paste0(figure_path,'figureS15.png'), height = 6, width = 7.5, units = 'in', res = 300)
print(figure)
dev.off()