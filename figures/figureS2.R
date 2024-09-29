source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')

gene_data <- read.csv(paste0(data_path, 'gene_phewas_burden_sig_count_239.csv'), sep = '\t') %>%
  filter(annotation != 'pLoF|missense|LC') %>%
  mutate(annotation = factor(annotation, levels = annotation_types),
         interval = get_freq_interval(CAF)) 

gene_prop_summary <- gene_data %>%
  group_by(interval, annotation) %>%
  dplyr::summarize(multi_trait_prop = sum(n_phewas_sig >1)/n(),
                   single_trait_prop = sum(n_phewas_sig == 1)/n(),) %>%
  pivot_longer(cols = c('single_trait_prop', 'multi_trait_prop'))  %>%
  mutate(name = factor(name, levels = c('single_trait_prop', 'multi_trait_prop'), labels = c('Trait-specific', 'Pleiotropic')))
gene_n_summary <- gene_data %>%
  group_by(interval, annotation) %>%
  dplyr::summarize(multi_trait_n = sum(n_phewas_sig >1),
                   single_trait_n = sum(n_phewas_sig == 1),)%>%
  pivot_longer(cols = c('single_trait_n', 'multi_trait_n'))  %>%
  mutate(name = factor(name, levels = c('single_trait_n', 'multi_trait_n'), labels = c('Trait-specific', 'Pleiotropic')))
max_n <- max(gene_n_summary$value)*15
max_p <- max(gene_prop_summary$value)
gene_summary = merge(gene_prop_summary, gene_n_summary, by = c('interval', 'annotation', 'name'))

figure <- gene_summary  %>%
  mutate(interval = factor(interval, levels = c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'),
                           labels = c('(0.01%, 0.1%]', '(0.1%, 1%]', '(1%, 10%]', paste0('(10%, ', bquote("\U221E"), ' )') ))) %>%
  ggplot +
  labs(x = 'CAF interval', y = 'Gene (%)', lty = NULL) + 
  geom_point(aes(x = interval, y = value.x, color = annotation, group=interaction(annotation))) + 
  geom_line(aes(x = interval, y = value.x, color = annotation, group=interaction(annotation))) +
  geom_col(aes(x = interval, y = value.y/max_n, color=annotation, fill=annotation), alpha=0.2) + 
  annotation_color_scale + 
  annotation_fill_scale + 
  scale_y_continuous(label=percent, sec.axis = sec_axis(~.*max_n, name="N genes")) +
  theme(legend.position = 'top', legend.direction = 'horizontal') +
  facet_grid(name~annotation, scale = 'free', labeller = label_type)+
  theme(axis.text.x = element_text(angle=15, hjust=1, size = 6))
figure

png(paste0(figure_path,'figureS2.png'), height = 4, width = 6, units = 'in', res = 300)
print(figure)
dev.off()