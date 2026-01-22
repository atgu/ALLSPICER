source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/constants.R')
library(magick)
library(ggvenn)
library(ggVennDiagram)

gene_data <- read.csv(paste0(data_path, 'gene_phewas_burden_sig_count_239.csv'), sep = '\t') %>%
  # filter(annotation != 'pLoF|missense|LC') %>%
  mutate(annotation = factor(annotation, levels = annotation_types2[c(1,2,4,3)]),
         interval = get_freq_interval(CAF))


gene_name_label <- gene_data %>%
  filter(n_phewas_sig > 10) %>%
  filter(annotation != 'pLoF|missense|LC') %>%
  mutate(annotation = factor(annotation, levels = annotation_types2[c(1,2,4,3)]))

figureA <- ggdraw() +
  draw_image(image_read(paste0('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/figure_2024/figure1/figure1_scheme.png')))+
  theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm"))

label_type = labeller(annotation = annotation_names2)
annotation_fill_scale2 = scale_fill_manual(name = 'Annotation', values = colors2, breaks = annotation_types2[c(1,2,4,3)], labels = annotation_names2[c(1,2,4,3)])

figureB <- gene_data %>%
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

figureB


# pLoF = gene_data %>% filter(n_phewas_sig > 1 & annotation == 'pLoF') %$% gene_symbol
# Missense = gene_data %>% filter(n_phewas_sig > 1 & annotation == 'missense|LC') %$% gene_symbol
# Synonymous = gene_data %>% filter(n_phewas_sig > 1 & annotation == 'synonymous') %$% gene_symbol
# pLoFMis = gene_data %>% filter(n_phewas_sig > 1 & annotation == 'pLoF|missense|LC') %$% gene_symbol
#
# Reduce(intersect, list(pLoF, Missense))
# Reduce(intersect, list(pLoF, Synonymous)) # 0
# Reduce(intersect, list(pLoF, pLoFMis))
# Reduce(intersect, list(Missense, Synonymous))
# Reduce(intersect, list(Missense, pLoFMis))
# Reduce(intersect, list(Synonymous, pLoFMis))
# Reduce(intersect, list(pLoF, Missense, Synonymous)) # 0
# Reduce(intersect, list(pLoF, Missense, pLoFMis))
# Reduce(intersect, list(Missense, Synonymous, pLoFMis))
# Reduce(intersect, list(pLoF, Missense, Synonymous, pLoFMis)) # 0
#
# gene_sets = list(pLoF=pLoF, Missense = Missense, Synonymous=Synonymous, pLoFMis=pLoFMis)
# ggVennDiagram(gene_sets, label = "count")
#
# data_599 <- read_csv(paste0(data_path, 'pleiotropy_2024_gene_burden_sig_cnt_summary_annotated.csv'))
# pLoF = data_599 %>% filter(n_sig_gene > 1 & annotation == 'pLoF') %$% gene_symbol
# Missense = data_599 %>% filter(n_sig_gene > 1 & annotation == 'missense|LC') %$% gene_symbol
# Synonymous = data_599 %>% filter(n_sig_gene > 1 & annotation == 'synonymous') %$% gene_symbol
# pLoFMis = data_599 %>% filter(n_sig_gene > 1 & annotation == 'pLoF|missense|LC') %$% gene_symbol
#
# Reduce(intersect, list(pLoF, Missense))
# Reduce(intersect, list(pLoF, Synonymous)) # 0
# Reduce(intersect, list(pLoF, pLoFMis))
# Reduce(intersect, list(Missense, Synonymous))
# Reduce(intersect, list(Missense, pLoFMis))
# Reduce(intersect, list(Synonymous, pLoFMis))
# Reduce(intersect, list(pLoF, Missense, Synonymous)) # 0
# Reduce(intersect, list(pLoF, Missense, pLoFMis))
# Reduce(intersect, list(Missense, Synonymous, pLoFMis))
# Reduce(intersect, list(pLoF, Missense, Synonymous, pLoFMis)) # 0
#
# gene_sets = list(pLoF=pLoF, Missense = Missense, pLoFMis=pLoFMis)
# ggVennDiagram(gene_sets, label = "count")

# figureC <- ggdraw() +
#   draw_image(image_read(paste0('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/figure_2024/figure1/figure1_venn.png')))+
#   theme(plot.margin = unit(c(0.3, 0, 0.7, 0), "cm"))

data_239 <- read.csv(paste0(data_path, 'gene_phewas_burden_sig_count_239.csv'), sep = '\t') %>%
  mutate(interval = get_freq_interval(freq=CAF))  %>%
  mutate(interval = factor(interval, levels = c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'),
                           labels = c('(0.01%, 0.1%]', '(0.1%, 1%]', '(1%, 10%]', paste0('(10%, ', bquote("\U221E"), ' )') )))%>%
  filter(n_phewas_sig >= 1)

pLoF_mis_239 <- data_239 %>% filter(annotation %in% c('pLoF', 'missense|LC'))
no_combine_239 <- data_239 %>% filter(annotation != 'pLoF|missense|LC')

figureC <- gene_list_pleiotropy_figure(pLoF_mis_239, test='burden_239', panel = F, filter_cat=T, overwrite = T) +
  theme(legend.title = element_text(face = 'plain', size = 11),
        axis.title = element_text(face = 'plain', size = 11),
        plot.margin = unit(c(0.5,1,0.5,0.2), "cm"),
        axis.title.y = element_blank())


figure = ggpubr::ggarrange(figureA, figureB, figureC,
                           labels = c('(A) Definitions of gene-annotation pair and pleiotropy',
                                      '(B) Number of associations per gene',
                                      '(C) Proportion of pleiotropic genes (association >1/>=1) across gene categories'),
                           ncol=1, common.legend=TRUE, vjust = 0.5, hjust = 0,
                           font.label = list(size = 10, color = "black", face = "bold", family = NULL),
                           heights = c(0.18, 0.18, 0.23))
png(paste0(figure_path,'figure1.png'), height = 8, width = 8, units = 'in', res = 300)
print(figure)
dev.off()
