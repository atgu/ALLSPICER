source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')

data <- read_csv(paste0(data_path, 'pleiotropy_2024_gene_burden_sig_cnt_summary_annotated.csv'))
gene_info <- read_delim('~/gene_lists/lists/data/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', delim = '\t') %>%
  dplyr::select(gene_symbol = gene, gene_id, pLI, oe_lof_upper_bin, cds_length, oe_lof_upper)

data <- data %>%
  merge(., gene_info, by = c('gene_symbol', 'gene_id'), all.x = T) %>%
  mutate(pleiotropy = if_else(n_sig_gene > 1, 'Pleiotropic', 'Non-pleiotropic'),
         domain_category = case_when(
           n_sig_gene == 0 ~ 'Non-associated',
           n_sig_gene == 1 ~ 'Trait-specific',
           n_sig_gene > 1 & n_domain_associated == 1 ~ 'Domain-specific',
           n_sig_gene > 1 & n_domain_associated > 1 ~ 'Multi-domain'
         ),
         disease_domain_category = case_when(
           n_disease_associated == 0 ~ 'Non-associated',
           n_disease_associated == 1 ~ 'Trait-specific',
           n_disease_associated > 1 & n_disease_domain_associated == 1 ~ 'Domain-specific',
           n_disease_associated > 1 & n_disease_domain_associated > 1 ~ 'Multi-domain'
         ),
         
         )

p1 <- data %>% 
  mutate(domain_category = factor(domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain'))) %>%
  filter(domain_category != 'Non-associated' & annotation != 'pLoF|missense|LC') %>%
  mutate(annotation = factor(annotation, levels= annotation_types)) %>%
  ggplot + aes(x = annotation, color = annotation, fill=annotation, alpha = domain_category) +
  geom_bar(position = "fill") +
  scale_alpha_discrete(name = NULL, range=c(0.25, 1)) + 
  scale_x_discrete(labels = annotation_names) +
  annotation_color_scale + annotation_fill_scale + themes +
  labs(x = NULL, y=NULL, color = NULL, fill=NULL) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 8),
        plot.margin = unit(c(1,0,0,0), "cm")) + 
  guides(color = 'none', fill='none')
p1

png(paste0('~/Desktop/figure1_proportion_domain_specific.png'), width=4, height=4, units = 'in', res = 300)
print(p1)
dev.off()


annotation_types3 = c('Watanabe et al.2019', annotation_types)
annotation_names3 = c('Watanabe et al.2019', annotation_names)
names(annotation_names3) = annotation_types3
colors3 = c(colors, 'Watanabe et al.2019' = '#1F77B4')
annotation_fill_scale3 = scale_fill_manual(name = 'Annotation', values = colors3, breaks = annotation_types3, labels = annotation_names3)
annotation_color_scale3 = scale_color_manual(name = 'Annotation', values = colors3, breaks = annotation_types3, labels = annotation_names3)

p12 <- data %>% 
  filter(domain_category != 'Non-associated' & annotation != 'pLoF|missense|LC') %>%
  group_by(annotation) %>%
  mutate(n = n()) %>%
  group_by(domain_category, annotation, n) %>%
  dplyr::summarize(cnt = n()) %>%
  mutate(prop = cnt/n) %>%
  dplyr::select(-n, -cnt) %>%
  rbind(., data.frame(domain_category = c('Trait-specific', 'Domain-specific', 'Multi-domain'),
                      annotation = rep('Watanabe et al.2019', 3),
                      prop = c(0.1880, 0.1403, 0.6717))) %>%
  mutate(domain_category = factor(domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain'))) %>%
  
  mutate(annotation = factor(annotation, levels= annotation_types3)) %>%
  ggplot + aes(x = annotation, y = prop, color = annotation, fill=annotation, alpha = domain_category) +
  geom_bar(stat = "identity") +
  scale_alpha_discrete(name = NULL, range=c(0.25, 1)) + 
  scale_x_discrete(labels = annotation_names3) +
  annotation_color_scale3 + annotation_fill_scale3 + themes +
  labs(x = NULL, y=NULL, color = NULL, fill=NULL) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 8),
        plot.margin = unit(c(1,0,0,0), "cm")) + 
  guides(color = 'none', fill='none')
p12

png(paste0('~/Desktop/figure1_proportion_domain_specific_add_common.png'), width=4, height=4, units = 'in', res = 300)
print(p12)
dev.off()

p2 <- data %>% 
  mutate(disease_domain_category = factor(disease_domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain'))) %>%
  filter(disease_domain_category != 'Non-associated'& annotation != 'pLoF|missense|LC') %>%
  mutate(annotation = factor(annotation, levels= annotation_types)) %>%
  ggplot + aes(x = annotation, color = annotation, fill=annotation, alpha = disease_domain_category) +
  geom_bar(position = "fill") +
  scale_alpha_discrete(name = NULL, range=c(0.25, 1)) + 
  scale_x_discrete(labels = annotation_names) +
  annotation_color_scale + annotation_fill_scale + themes +
  labs(x = NULL, y=NULL, color = NULL, fill=NULL) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 8),
        plot.margin = unit(c(1,0,0,0), "cm")) + 
  guides(color = 'none', fill='none') 

png(paste0('~/Desktop/figure1_proportion_disease_domain_specific.png'), width=4, height=4, units = 'in', res = 300)
print(p2)
dev.off()

pp <- ggpubr::ggarrange(p12, p2, ncol=2,widths = c(1.3, 1), 
                       labels = c('(A) Phenotypic Domain', '(B) Disease Domain'),
                       font.label = list(size = 10, color = "black", face = "bold", family = NULL), hjust =0, common.legend = T)

png(paste0('~/Desktop/figure1_proportion_both_domain_specific_add_common.png'), width=6, height=3, units = 'in', res = 300)
print(pp)
dev.off()

library(gg.layers)
p <- data %>%
  filter(annotation != 'pLoF|missense|LC') %>% 
  # filter(annotation == 'pLoF') %>%
  group_by(annotation, domain_category) %>%
  dplyr::summarize(mean = mean(pLI, na.rm = TRUE),
                   sem = sd(pLI, na.rm = T)/sqrt(n()),                    
                   prop = sum(pLI > 0.9, na.rm = T) / n(), 
                   cnt = n()) %>%
  mutate(sd = 1.96* sqrt(prop* (1-prop)/cnt)) %>%
  mutate(domain_category = factor(domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain'))) %>%
  mutate(annotation = factor(annotation, levels= annotation_types)) %>%
  ggplot + geom_pointrange(aes(x = domain_category, y = prop, ymax = prop + sd, ymin = prop -sd, color = annotation)) + 
  labs(x = NULL, y = 'pLI > 0.9 (%)') +
  annotation_color_scale + annotation_fill_scale + themes +
  scale_y_continuous(label=percent) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(~annotation, labeller = labeller(annotation=annotation_names), scale='free')
png(paste0('~/Desktop/figure1_pLI_3_anno_over_0.9.png'), width=6, height=3, units = 'in', res = 300)
print(p)
dev.off()


p <- data %>%
  filter(annotation != 'pLoF|missense|LC') %>% 
  # filter(annotation == 'pLoF') %>%
  group_by(annotation, disease_domain_category) %>%
  dplyr::summarize(mean = mean(pLI, na.rm = TRUE), 
                   sem = sd(pLI, na.rm = T)/sqrt(n()),                    
                   prop = sum(pLI > 0.9, na.rm = T) / n(), 
                   cnt = n()) %>%
  mutate(sd = 1.96* sqrt(prop* (1-prop)/cnt)) %>%
  mutate(disease_domain_category = factor(disease_domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain'))) %>%
  mutate(annotation = factor(annotation, levels= annotation_types)) %>%
  ggplot + 
  geom_pointrange(aes(x = disease_domain_category, y = prop, ymax = prop + sd, ymin = prop -sd, color = annotation)) + 
  labs(x = NULL, y = 'pLI > 0.9 (%)') +
  annotation_color_scale + annotation_fill_scale + themes +
  scale_y_continuous(label=percent) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(~annotation, labeller = labeller(annotation=annotation_names), scale='free')
p
png(paste0('~/Desktop/figure1_pLI_3_anno_disease_domain_over_0.9.png'), width=6, height=3, units = 'in', res = 300)
print(p)
dev.off()

p <- data %>%
  filter(annotation != 'pLoF|missense|LC') %>% 
  group_by(annotation, domain_category) %>%
  dplyr::summarize(mean = mean(oe_lof_upper_bin, na.rm = TRUE), 
                   sem = sd(oe_lof_upper_bin, na.rm = T)/sqrt(n()),                    
                   prop = sum(oe_lof_upper_bin == 0, na.rm = T) / n(), 
                   cnt = n()) %>%
  mutate(sd = 1.96* sqrt(prop* (1-prop)/cnt)) %>%
  mutate(domain_category = factor(domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain'))) %>%
  mutate(annotation = factor(annotation, levels= annotation_types)) %>%
  ggplot + 
  # aes(x = domain_category, y = oe_lof_upper_bin, color = annotation) + 
  # geom_boxplot() + 
  geom_pointrange(aes(x = domain_category, y = prop, ymax = prop + sd, ymin = prop -sd, color = annotation)) + 
  labs(x = NULL, y = 'oe_lof_upper_bin = 0 (%)') +
  annotation_color_scale + annotation_fill_scale + themes +
  scale_y_continuous(label=percent) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(~annotation, labeller = labeller(annotation=annotation_names), scale='free')
png(paste0('~/Desktop/figure1_oe_lof_upper_bin_3_anno_bin_0.png'), width=6, height=3, units = 'in', res = 300)
print(p)
dev.off()

p <- data %>%
  filter(annotation != 'pLoF|missense|LC') %>% 
  group_by(annotation, disease_domain_category) %>%
  dplyr::summarize(mean = mean(oe_lof_upper_bin, na.rm = TRUE), 
                   sem = sd(oe_lof_upper_bin, na.rm = T)/sqrt(n()),                    
                   prop = sum(oe_lof_upper_bin == 0, na.rm = T) / n(), 
                   cnt = n()) %>%
  mutate(sd = 1.96* sqrt(prop* (1-prop)/cnt)) %>%
  mutate(disease_domain_category = factor(disease_domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain'))) %>%
  mutate(annotation = factor(annotation, levels= annotation_types)) %>%
  ggplot + 
  # aes(x = domain_category, y = oe_lof_upper_bin, color = annotation) + 
  # geom_boxplot() + 
  geom_pointrange(aes(x = disease_domain_category, y = prop, ymax = prop + sd, ymin = prop -sd, color = annotation)) + 
  labs(x = NULL, y = 'oe_lof_upper_bin = 0 (%)') +
  annotation_color_scale + annotation_fill_scale + themes +
  scale_y_continuous(label=percent) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(~annotation, labeller = labeller(annotation=annotation_names), scale='free')
png(paste0('~/Desktop/figure1_oe_lof_upper_bin_3_anno_disease_domain_bin_0.png'), width=6, height=3, units = 'in', res = 300)
print(p)
dev.off()

p <- data %>%
  filter(annotation != 'pLoF|missense|LC') %>% 
  # group_by(annotation, domain_category) %>%
  # dplyr::summarize(mean = mean(oe_lof_upper_bin, na.rm = TRUE), 
  #                  sem = sd(oe_lof_upper_bin, na.rm = T)/sqrt(n()),                    
  #                  prop = sum(oe_lof_upper_bin == 0, na.rm = T) / n(), 
  #                  cnt = n()) %>%
  # mutate(sd = 1.96* sqrt(prop* (1-prop)/cnt)) %>%
  mutate(disease_domain_category = factor(disease_domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain'))) %>%
  mutate(annotation = factor(annotation, levels= annotation_types)) %>%
  ggplot + 
  aes(x = disease_domain_category, y = oe_lof_upper, color = annotation) +
  geom_boxplot() +
  # geom_pointrange(aes(x = domain_category, y = prop, ymax = prop + sd, ymin = prop -sd, color = annotation)) + 
  labs(x = NULL, y = 'LOEUF') +
  annotation_color_scale + annotation_fill_scale + themes +
  # scale_y_continuous(label=percent) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(~annotation, labeller = labeller(annotation=annotation_names), scale='free')
png(paste0('~/Desktop/figure1_loeuf_3_anno_disease_domain.png'), width=6, height=3, units = 'in', res = 300)
print(p)
dev.off()

p <- data %>%
  filter(annotation != 'pLoF|missense|LC') %>% 
  # group_by(annotation, domain_category) %>%
  # dplyr::summarize(mean = mean(oe_lof_upper_bin, na.rm = TRUE), 
  #                  sem = sd(oe_lof_upper_bin, na.rm = T)/sqrt(n()),                    
  #                  prop = sum(oe_lof_upper_bin == 0, na.rm = T) / n(), 
  #                  cnt = n()) %>%
  # mutate(sd = 1.96* sqrt(prop* (1-prop)/cnt)) %>%
  mutate(domain_category = factor(domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain'))) %>%
  mutate(annotation = factor(annotation, levels= annotation_types)) %>%
  ggplot + 
  aes(x = domain_category, y = oe_lof_upper, color = annotation) +
  geom_boxplot() +
  # geom_pointrange(aes(x = domain_category, y = prop, ymax = prop + sd, ymin = prop -sd, color = annotation)) + 
  labs(x = NULL, y = 'LOEUF') +
  annotation_color_scale + annotation_fill_scale + themes +
  # scale_y_continuous(label=percent) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(~annotation, labeller = labeller(annotation=annotation_names), scale='free')
png(paste0('~/Desktop/figure1_loeuf_3_anno.png'), width=6, height=3, units = 'in', res = 300)
print(p)
dev.off()

p <- data %>%
  filter(annotation != 'pLoF|missense|LC') %>% 
  group_by(annotation, domain_category) %>%
  dplyr::summarize(mean = mean(cds_length, na.rm = TRUE),
                   sem = sd(cds_length, na.rm = T)/sqrt(n())) %>%
  mutate(domain_category = factor(domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain'))) %>%
  mutate(annotation = factor(annotation, levels= annotation_types)) %>%
  ggplot + 
  # aes(x = domain_category, y = cds_length, color = annotation) +
  # geom_boxplot2() +
  geom_pointrange(aes(x = domain_category, y = mean, ymax = mean + sem, ymin = mean -sem, color = annotation)) + 
  labs(x = NULL, y = 'CDS length') +
  annotation_color_scale + annotation_fill_scale + themes +
  scale_y_log10(label=comma) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(~annotation, labeller = labeller(annotation=annotation_names), scale='free')
png(paste0('~/Desktop/figure1_cds_length.png'), width=6, height=3, units = 'in', res = 300)
print(p)
dev.off()

p <- data %>%
  filter(annotation != 'pLoF|missense|LC') %>% 
  group_by(annotation, disease_domain_category) %>%
  dplyr::summarize(mean = mean(cds_length, na.rm = TRUE),
                   sem = sd(cds_length, na.rm = T)/sqrt(n())) %>%
  mutate(disease_domain_category = factor(disease_domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain'))) %>%
  mutate(annotation = factor(annotation, levels= annotation_types)) %>%
  ggplot + 
  # aes(x = domain_category, y = cds_length, color = annotation) +
  # geom_boxplot2() +
  geom_pointrange(aes(x = disease_domain_category, y = mean, ymax = mean + sem, ymin = mean -sem, color = annotation)) + 
  labs(x = NULL, y = 'CDS length') +
  annotation_color_scale + annotation_fill_scale + themes +
  scale_y_log10(label=comma) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(~annotation, labeller = labeller(annotation=annotation_names), scale='free')
png(paste0('~/Desktop/figure1_cds_length_disease_domain.png'), width=6, height=3, units = 'in', res = 300)
print(p)
dev.off()

p < - data %>%
  filter(annotation != 'pLoF|missense|LC' & cds_length < 30000) %>% 
  ggplot + aes(x = cds_length, y = total_variants, color = annotation) +
  geom_point() + labs(x = 'CDS length', y = 'Total variants') + 
  geom_smooth(method = 'lm', lty = 2, lwd =1) +
  annotation_color_scale + annotation_fill_scale + themes
png(paste0('~/Desktop/cds_length_disease_domain.png'), width=6, height=3, units = 'in', res = 300)
print(p)
dev.off()

data %>%
  filter(annotation != 'pLoF|missense|LC' & cds_length < 30000) %>% 
  ggplot + aes(x = cds_length, y = n_sig_gene, color = annotation) +
  geom_point() + labs(x = 'CDS length', y = 'N sig') + 
  geom_smooth(method = 'lm', lty = 2, lwd =1) +
  annotation_color_scale + annotation_fill_scale + themes

data %>%
  filter(annotation != 'pLoF|missense|LC' & cds_length < 30000) %>% 
  ggplot + aes(x = total_variants, y = n_sig_gene, color = annotation) +
  geom_point() + labs(x = 'Total variants', y = 'N sig') + 
  geom_smooth(method = 'lm', lty = 2, lwd =1) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, labeller = labeller(annotation = annotation_names), scale = 'free')




figure1 <- data %>%
  filter(annotation != 'pLoF|missense|LC') %>%
  filter(n_sig_gene > 1) %>%
  mutate(annotation = factor(annotation, levels = annotation_types),
         interval = get_freq_interval(CAF)) %>%
  ggplot + aes(x = n_sig_gene, color = annotation, fill=annotation) +
  labs(x = 'N associations', y = 'N genes') +
  # geom_density(stat = 'count') +
  geom_histogram(binwidth = 1, alpha=0.5, position = 'dodge', stat ='count')  +
  #scale_y_log10(label = comma) +
  # scale_x_log10(label = comma) +
  annotation_color_scale2 + 
  annotation_fill_scale2 + 
  themes + theme_classic() + 
  theme(legend.position = 'top',
        legend.direction = 'horizontal') +
  #   geom_text_repel(data = gene_name_label, aes(x = n_phewas_sig, y = 20, label = gene_symbol, color = annotation), size = 3, show.legend = FALSE)  +
  # geom_line(data = pois_data, aes(x = x, y=y, color=annotation)) +
  # geom_point(data = pois_data, aes(x = x, y=y, color=annotation)) +
  facet_grid(~annotation, labeller = label_type, scale = 'free')

figure1
figure2 <- data %>%
  filter(annotation != 'pLoF|missense|LC') %>%
  # filter(n_phewas_sig > 0) %>%
  mutate(annotation = factor(annotation, levels = annotation_types2),
         interval = get_freq_interval(CAF))  %>%
  group_by(interval, annotation) %>%
  mutate(interval = factor(interval, levels = c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'),
                           labels = c('(0.01%, 0.1%]', '(0.1%, 1%]', '(1%, 10%]', paste0('(10%, ', bquote("\U221E"), ' )') )))%>%
  summarize(n_pleiotropy = sum(n_sig_gene > 1)) %>%
  ggplot + aes(x = interval, y = n_pleiotropy, label=NULL, fill=annotation, color=annotation) + 
  geom_bar(stat='identity', alpha = .5) +
  labs(x = 'CAF interval', y = 'N pleiotropic genes') + 
  annotation_color_scale + annotation_fill_scale + themes + theme_classic() + theme(axis.text.x = element_text(size = 7)) +
  facet_grid(~annotation, scale = 'free', labeller = label_type)
figure2



figure = ggpubr::ggarrange(figure1, NULL, figure2, labels = c('(A) Number of associations per gene', '', '(B) Number of genes with association > 1 across CAF intervals'), nrow=3, common.legend=TRUE, vjust = 0, hjust = 0, font.label = list(size = 10, color = "black", face = "bold", family = NULL),heights = c(0.18, 0.02, 0.18))
figure
png(paste0(figure_path,'main_fig1.png'), height = 5, width = 7.5, units = 'in', res = 300)
print(figure)
dev.off()


figure3 <- data %>%
  filter(annotation != 'pLoF|missense|LC') %>%
  filter(n_domain_associated > 0 & n_sig_gene >1) %>%
  mutate(annotation = factor(annotation, levels = annotation_types),
         interval = get_freq_interval(CAF)) %>%
  ggplot + aes(x = n_domain_associated, color = annotation, fill=annotation) +
  labs(x = 'N domain associated', y = 'N genes') +
  geom_histogram(stat ='count')  +
  annotation_color_scale2 + 
  annotation_fill_scale2 + 
  themes + theme_classic() + 
  theme(legend.position = 'none',
        plot.margin = unit(c(1,0,0,0), "cm")) +
  facet_grid(~annotation, labeller = label_type, scale = 'free')

figure4 <- data %>%
  filter(annotation != 'pLoF|missense|LC') %>%
  filter(n_disease_domain_associated > 0 & n_sig_gene >1) %>%
  mutate(annotation = factor(annotation, levels = annotation_types),
         interval = get_freq_interval(CAF)) %>%
  ggplot + aes(x = n_disease_domain_associated, color = annotation, fill=annotation) +
  labs(x = 'N disease domain associated', y = 'N genes') +
  geom_histogram(stat ='count')  +
  annotation_color_scale2 + 
  annotation_fill_scale2 + 
  themes + theme_classic() + 
  theme(legend.position = 'none',
        plot.margin = unit(c(1,0,0,0), "cm")) +
  facet_grid(~annotation, labeller = label_type, scale = 'free')

figure = ggpubr::ggarrange(figure3, figure4, labels = c('(A) Phenotypic domains', '(B) Disease domains'), nrow=2, hjust = 0, font.label = list(size = 10, color = "black", face = "bold", family = NULL))
figure
png('~/Desktop/figure1_n_genes_n_domains.png', height = 5, width = 7.5, units = 'in', res = 300)
print(figure)
dev.off()



## logistic regression
library(pROC)
library(plotROC)
data <- data %>% filter(annotation == 'pLoF')
logistic_reg <- function(type, cut, write=FALSE){
  print(table(data$annotation, data$domain_category))
  sub_data <- data %>%
    filter(domain_category %in% c(type, 'Non-associated')) %>%
    mutate(y = if_else(domain_category == type, 1, 0)) %>%
    dplyr::select(y, oe_lof_upper, cds_length, domain_category) %>%
    filter(complete.cases(.)) %>%
    distinct(.)
  
  model <- glm(y ~ oe_lof_upper + cds_length ,family=binomial(link='logit'), 
               data=sub_data)
  
  sub_data <- sub_data %>%
    mutate(prediction = predict(model),
           residual = residuals(model))
  auc <- roc(sub_data$y, sub_data$prediction)$auc
  
  p1 <- sub_data %>%
    filter(prediction < cut) %>%
    ggplot + aes(x = prediction, y = residual,colour = domain_category) +
    geom_point( )+ 
    geom_hline(yintercept = 0, lty = 2) +
    scale_color_manual(breaks=c(type, 'Non-associated'), 
                       labels=c(type, 'Non-associated'), 
                       values = c( '#cc3311', '#004488')) +
    labs(x = 'Prediction', y = 'Residual', color = NULL) + themes
  
  p2 <- sub_data %>%
    ggplot + aes(d = y, m = prediction) + 
    labs(x = 'False positive rate', y = 'True positive rate') + 
    geom_roc(n.cuts = 0) +
    geom_abline() + 
    annotate(x = Inf, y = 0, label = paste('AUC:', round(auc,3)), geom='text', hjust = 1) + themes
  
  pp <- ggpubr::ggarrange(p1, p2, ncol=2)

  model1 <- sub_data %>%
    do(tidy(glm(y ~ oe_lof_upper + cds_length ,family=binomial(link='logit'), data=.)))
  print(model1)
  model2 <- sub_data %>%
    do(glance(glm(y ~ oe_lof_upper + cds_length ,family=binomial(link='logit'), data=.)))
  print(model2)
  
  if(write){
    png(paste0('~/Desktop/logistic_regression_', type,'.png'), width=7.5, height=2.5, units = 'in', res = 300)
    print(pp)
    dev.off()
    write_csv(model1, paste0('~/Desktop/logistic_regression_', type,'_1.csv'))
    write_csv(model2, paste0('~/Desktop/logistic_regression_', type,'_2.csv'))
  }
  
  
}


logistic_reg2 <- function(type, cut){
  sub_data <- data %>%
    filter(domain_category %in% c(type, 'Non-associated')) %>%
    mutate(y = if_else(domain_category == type, 1, 0),
           pLI_cat = as.numeric(pLI > 0.9)) %>%
    dplyr::select(y, pLI_cat, cds_length, domain_category) %>%
    filter(complete.cases(.)) %>%
    distinct(.)
  
  model <- glm(y ~ pLI_cat + cds_length ,family=binomial(link='logit'), 
               data=sub_data)
  
  sub_data <- sub_data %>%
    mutate(prediction = predict(model),
           residual = residuals(model))
  auc <- roc(sub_data$y, sub_data$prediction)$auc
  
  p1 <- sub_data %>%
    filter(prediction > cut) %>%
    ggplot + aes(x = prediction, y = residual,colour = domain_category) +
    geom_point( )+ 
    geom_hline(yintercept = 0, lty = 2) +
    scale_color_manual(breaks=c(type, 'Non-associated'), 
                       labels=c(type, 'Non-associated'), 
                       values = c( '#cc3311', '#004488')) +
    labs(x = 'Prediction', y = 'Residual', color = NULL) + themes
  
  p2 <- sub_data %>%
    ggplot + aes(d = y, m = prediction) + 
    labs(x = 'False positive rate', y = 'True positive rate') + 
    geom_roc(n.cuts = 0) +
    geom_abline() + 
    annotate(x = Inf, y = 0, label = paste('AUC:', round(auc,3)), geom='text', hjust = 1) + themes
  
  pp <- ggpubr::ggarrange(p1, p2, ncol=2)
  
  png(paste0('~/Desktop/logistic_regression_', type,'.png'), width=7.5, height=2.5, units = 'in', res = 300)
  print(pp)
  dev.off()
  
  
  model1 <- sub_data %>%
    do(tidy(glm(y ~ pLI_cat + cds_length,family=binomial(link='logit'), data=.)))
  # write_csv(model1, paste0('~/Desktop/logistic_regression_', type,'_1.csv'))
  print(model1)
  
  model2 <- sub_data %>%
    do(glance(glm(y ~ pLI_cat + cds_length,family=binomial(link='logit'), data=.)))
  # write_csv(model2, paste0('~/Desktop/logistic_regression_', type,'_2.csv'))
  print(model2)
  
}

logistic_reg('Multi-domain', cut = -4)
logistic_reg('Domain-specific', cut = -6)
logistic_reg('Trait-specific', cut = -6)

logistic_reg2('Multi-domain', cut = -4)
logistic_reg2('Domain-specific', cut = -6)
logistic_reg2('Trait-specific', cut = -6)

load_tx_summary = function(expression_cutoff=0.3) {
  tx_expression = read_delim('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/data/GTEx.v7.median_expression_per_tx_per_tissue.021018.tsv.bgz', delim = '\t') %>%
    dplyr::select(-c(v, agg_expression, Bladder, Brain_Spinalcord_cervicalc_1_, Brain_Substantianigra,
              Cervix_Ectocervix,Cervix_Endocervix, FallopianTube, Kidney_Cortex,
              MinorSalivaryGland, Uterus, Ovary,Testis, Vagina,
              Cells_EBV_transformedlymphocytes, Cells_Transformedfibroblasts, Prostate)) %>%
    rename(transcript = transcript_id)
  
  tx_melt = tx_expression %>%
    gather('tissue', 'expression', -transcript, -gene_id)
  
  tx_melt %>%
    group_by(transcript, gene_id) %>%
    summarize(mean_expression=mean(expression, na.rm=T),
              brain_mean_expression=mean(expression * ifelse(grepl('Brain', tissue), 1, NA), na.rm=T),
              n_tissues_expressed=sum(expression > expression_cutoff, na.rm=T),
              max_expression=max(expression, na.rm=T),
              max_expressed_tissue=tissue[which.max(expression)]
    ) %>% ungroup %>% 
    return
}


load_sig_data <- function(){
  data <- read_csv(paste0(data_path, 'pleiotropy_2024_gene_burden_sig_cnt_summary_annotated.csv'))
  
  data <- data %>%
    mutate(pleiotropy = if_else(n_sig_gene > 1, 'Pleiotropic', 'Non-pleiotropic'),
           domain_category = case_when(
             n_sig_gene == 0 ~ 'Non-associated',
             n_sig_gene == 1 ~ 'Trait-specific',
             n_sig_gene > 1 & n_domain_associated == 1 ~ 'Domain-specific',
             n_sig_gene > 1 & n_domain_associated > 1 ~ 'Multi-domain'
           ),
           disease_domain_category = case_when(
             n_disease_associated == 0 ~ 'Non-associated',
             n_disease_associated == 1 ~ 'Trait-specific',
             n_disease_associated > 1 & n_disease_domain_associated == 1 ~ 'Domain-specific',
             n_disease_associated > 1 & n_disease_domain_associated > 1 ~ 'Multi-domain'
           ),
           
    )
  return(data)
}

generate_expression_data_plot <-  function() {
  all_tx = transcript_data %>%
    left_join(tx_summary) %>%
    rename(transcript_oe_lof_upper_bin = oe_lof_upper_bin) %>%
    left_join(gene_data %>% 
                dplyr::select(gene_id, oe_lof_upper_bin)) 
  all_tx <- sig_data %>%
    merge(., all_tx, by.x = c('gene_symbol', 'gene_id'), by.y = c('gene', 'gene_id')) %>%
    mutate(CAF_bin = factor(get_caf_interval(CAF), levels = caf_types))
  
  return(all_tx)
}

tx_summary <- load_tx_summary()
sig_data <- load_sig_data()
transcript_data <- read_delim('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/data/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz', delim = '\t')
gene_data <- read_delim('~/gene_lists/lists/data/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', delim = '\t')
tx_colors <- c('all' = 'dodgerblue1',
              'canonical' = 'navy')

get_caf_interval = function(freq){
  interval = case_when(
    freq <= 1e-4  ~ '[0, 0.0001]',
    freq <= 1e-3  ~ '(0.0001, 0.001]',
    freq <= 1e-2  ~ '(0.001, 0.01]',
    freq <= 1e-1  ~ '(0.01, 0.1]',
    freq > 1e-1  ~ '(0.1, )',
  )
  return(interval)
}
caf_types = c('[0, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, )')
caf_names = c(paste0('CAF:[0 , 0.0001]'), 'CAF:(0.0001, 0.001]', 'CAF:(0.001, 0.01]', 'CAF:(0.01, 0.1]', paste0('CAF:(0.1, ', bquote("\U221E"), ' )'))
names(caf_names) = caf_types


plot_data <- generate_expression_data_plot()%>%
  filter(!only_canonical | canonical)%>%
  filter(annotation != 'pLoF|missense|LC') %>%
  mutate(annotation = factor(annotation, levels = annotation_types),
         domain_category = factor(domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain')))
summary_plot_data <- plot_data %>%
  group_by(domain_category, annotation, CAF_bin) %>%
  summarize(mean_metric = mean(n_tissues_expressed, na.rm=T),
            sd_metric = sd(n_tissues_expressed, na.rm=T),
            n = n(),
            sem_metric = 1.96 * sd_metric / sqrt(n))%>%
  mutate(annotation = factor(annotation, levels = annotation_types),
         domain_category = factor(domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain')),
         CAF_bin = factor(CAF_bin, levels = caf_types)
         )
p <- plot_data  %>%
  ggplot +
  labs(x = NULL, y='Number of tissues where\ncanonical transcript is expressed')+
  geom_violin(aes(x = CAF_bin, group=interaction(CAF_bin, annotation), y = n_tissues_expressed, fill = annotation),
              data = plot_data,
              alpha = 0.2, color = F) +
  geom_pointrange(aes(x = CAF_bin, y = mean_metric,
                      ymin = mean_metric - sem_metric,
                      ymax = mean_metric + sem_metric, color = annotation),
                  data = summary_plot_data)+ 
  geom_text(aes(x = CAF_bin, y = 40,
                label = paste('N:', n), color = annotation), size = 2,
                  data = summary_plot_data)+
  themes + annotation_color_scale + annotation_fill_scale + 
  facet_grid(annotation~domain_category, labeller = labeller(annotation=annotation_names, CAF_bin=caf_names), scale='free') +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  guides(color = 'none', fill = 'none')
png(paste0('~/Desktop/canonical_expression.png'), width=7.5, height=5, units = 'in', res = 300)
print(p)
dev.off()


summary_plot_data <- plot_data %>%
  group_by(domain_category, annotation) %>%
  summarize(mean_metric = mean(n_tissues_expressed, na.rm=T),
            sd_metric = sd(n_tissues_expressed, na.rm=T),
            n = n(),
            sem_metric = 1.96 * sd_metric / sqrt(n))%>%
  mutate(annotation = factor(annotation, levels = annotation_types),
         domain_category = factor(domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain'))
  )
p <- plot_data  %>%
  ggplot +
  labs(x = NULL, y='Number of tissues where\ncanonical transcript is expressed')+
  geom_violin(aes(x = domain_category, group=interaction(domain_category, annotation), y = n_tissues_expressed, fill = annotation),
              data = plot_data,
              alpha = 0.2, color = F) +
  geom_pointrange(aes(x = domain_category, y = mean_metric,
                      ymin = mean_metric - sem_metric,
                      ymax = mean_metric + sem_metric, color = annotation),
                  data = summary_plot_data)+ 
  geom_text(aes(x = domain_category, y = 40,
                label = paste('N:', n), color = annotation), size = 2,
            data = summary_plot_data)+
  themes + annotation_color_scale + annotation_fill_scale + 
  facet_grid(~annotation, labeller = labeller(annotation=annotation_names, CAF_bin=caf_names), scale='free') +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  guides(color = 'none', fill = 'none')
png(paste0('~/Desktop/canonical_expression_overall.png'), width=6, height=3, units = 'in', res = 300)
print(p)
dev.off()


plot_data <- generate_expression_data_plot()%>%
  filter(annotation != 'pLoF|missense|LC') %>%
  mutate(annotation = factor(annotation, levels = annotation_types),
         domain_category = factor(domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain')),
         )
summary_plot_data_1 <- plot_data %>%
  group_by(domain_category, annotation, CAF_bin) %>%
  summarize(mean_metric = mean(n_tissues_expressed, na.rm=T),
            sd_metric = sd(n_tissues_expressed, na.rm=T),
            n = n(),
            sem_metric = 1.96 * sd_metric / sqrt(n),
            transcripts = 'all')%>%
  mutate(annotation = factor(annotation, levels = annotation_types),
         domain_category = factor(domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain')),
         CAF_bin = factor(CAF_bin, levels = caf_types)
  )
summary_plot_data_2 <- plot_data %>%
  filter(canonical)%>%
  group_by(domain_category, annotation, CAF_bin) %>%
  summarize(mean_metric = mean(n_tissues_expressed, na.rm=T),
            sd_metric = sd(n_tissues_expressed, na.rm=T),
            n = n(),
            sem_metric = 1.96 * sd_metric / sqrt(n),
            transcripts = 'canonical')%>%
  mutate(annotation = factor(annotation, levels = annotation_types),
         domain_category = factor(domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain')),
         CAF_bin = factor(CAF_bin, levels = caf_types)
  )
summary_data <- rbind(summary_plot_data_1, summary_plot_data_2)


get_n_tissue_interval = function(n_tissues){
  interval = case_when(
    n_tissues == 0  ~ '0 tissue',
    n_tissues == 1  ~ '1 tissue',
    n_tissues <= 2  ~ '2 tissues',
    n_tissues <= 5  ~ '3-5 tissues',
    n_tissues <= 10  ~ '6-10 tissues',
    n_tissues <= 15  ~ '10-15 tissues',
    n_tissues <= 25  ~ '16-25 tissues',
    n_tissues <= 35  ~ '26-35 tissues',
    n_tissues > 35  ~ '>35 tissues',
  )
  return(interval)
}

p <- plot_data  %>%
  filter(!is.na(n_tissues_expressed)) %>%
  mutate(n_tissue_bins = factor(get_n_tissue_interval(n_tissues_expressed), 
                                levels = rev(c('0 tissue', '1 tissue', '2 tissues', '3-5 tissues', '6-10 tissues', '10-15 tissues', '16-25 tissues', '26-35 tissues', '>35 tissues')))) %>%
  group_by(annotation, n_domain_associated) %>%
  mutate(total = n()) %>%
  group_by(n_tissue_bins, annotation, n_domain_associated, total)  %>%
  dplyr::summarize(prop = n()/total) %>%
  distinct() %>%
  # filter(n_tissue_bins != '0 tissue') %>%
  ggplot + aes(x = n_domain_associated, y = prop, color = n_tissue_bins, fill=n_tissue_bins) + 
  labs(x = 'Number of associated domains', y='Proportion of genes')+
  geom_bar(stat = "identity") + 
  themes + 
  scale_color_brewer(name = NULL, palette = 'RdYlBu') +
  scale_fill_brewer(name = NULL, palette = 'RdYlBu') +
  facet_grid(~annotation, labeller = labeller(annotation=annotation_names, CAF_bin=caf_names), scale='free') +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) 
png(paste0('~/Desktop/N_tissue_expressed.png'), width=6, height=3, units = 'in', res = 300)
print(p)
dev.off()
