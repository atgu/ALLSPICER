source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')

gene_data <- read.csv(paste0(data_path, 'gene_phewas_burden_sig_count_239.csv'), sep = '\t') %>%
  filter(annotation != 'pLoF|missense|LC') %>%
  mutate(annotation = factor(annotation, levels = annotation_types),
         interval = get_freq_interval(CAF))


gene_name_label <- gene_data %>% 
  filter(n_phewas_sig > 10) %>%
  filter(annotation != 'pLoF|missense|LC') %>%
  mutate(annotation = factor(annotation, levels = annotation_types))

figureA <- gene_data %>%
  filter(n_phewas_sig > 1) %>%
  ggplot + aes(x = n_phewas_sig, fill=annotation) +
  geom_histogram(binwidth = 1, alpha = 0.8, position = 'dodge', stat ='count', color='white')  +
  annotation_color_scale2 + 
  annotation_fill_scale2 + 
  themes + theme_classic() + 
  theme(legend.position = 'top',
        legend.direction = 'horizontal') +
    geom_text(data = gene_name_label, aes(x = n_phewas_sig, y = 3, label = gene_symbol, color = annotation), size = 3, show.legend = FALSE)  +
  geom_text(data = gene_data %>% filter(n_phewas_sig >= 1) %>% group_by(annotation) %>% dplyr::summarize(n_pleiotropy = sum(n_phewas_sig>1), n = n(), p = sum(n_phewas_sig>1)/n()),
                  aes(x = Inf, y = 60, color = annotation,label = paste0('Pleiotropic genes (%):\n', n_pleiotropy, '/', n, '=', round(p*100, 2), '%')), hjust = 1) + 
  labs(x = 'N associations', y = 'N genes', color = 'Annotation', fill= 'Annotation') +
  facet_grid(~annotation, labeller = label_type, scale = 'free')+
  theme(plot.margin = unit(c(0.5,1,0.5,0.2), "cm"))
# figureA

data_239 <- read.csv(paste0(data_path, 'gene_phewas_burden_sig_count_239.csv'), sep = '\t') %>%
  mutate(interval = get_freq_interval(freq=CAF))  %>%
  mutate(interval = factor(interval, levels = c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'),
                           labels = c('(0.01%, 0.1%]', '(0.1%, 1%]', '(1%, 10%]', paste0('(10%, ', bquote("\U221E"), ' )') )))%>% 
  filter(n_phewas_sig >= 1)

pLoF_mis_239 <- data_239 %>% filter(annotation %in% c('pLoF', 'missense|LC'))
no_combine_239 <- data_239 %>% filter(annotation != 'pLoF|missense|LC')

figureB <- gene_list_pleiotropy_figure(pLoF_mis_239, test='burden_239', panel = F, filter_cat=T, overwrite = T) +
  theme(legend.title = element_text(face = 'plain', size = 11), 
        axis.title = element_text(face = 'plain', size = 11), 
        plot.margin = unit(c(0.5,1,0.5,0.2), "cm"),
        axis.title.y = element_blank())



figure = ggpubr::ggarrange(figureA, figureB,  
                           labels = c('(A) Number of associations per gene', 
                                      '(B) Proportion of pleiotropic genes (association >1/>=1) across gene categories'), 
                           nrow=2, common.legend=TRUE, vjust = 0, hjust = 0, 
                           font.label = list(size = 10, color = "black", face = "bold", family = NULL),
                           heights = c(0.18, 0.23))
png(paste0(figure_path,'figure1.png'), height = 6, width = 8, units = 'in', res = 300)
print(figure)
dev.off()