source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/constants.R')
library(dplyr)
library(readr)

pheno_domain <- read_delim(paste0(data_path, 'pheno', '_domain_level_gene_sig_burden_500k.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character())) %>%
  filter(!is.na(pheno_group))

pheno_info <- read_delim(paste0(data_path, 'wrap_up_results_final_599_phenotypes.txt.bgz'), delim = '\t')
table(pheno_info$pheno_group)
table(pheno_info$disease_group)



p1 <- pheno_info %>%
  group_by(pheno_group) %>%
  dplyr::summarize(cnt = n()) %>%
  arrange(desc(pheno_group)) %>%
  mutate(lab.ypos = cumsum(cnt) -0.5*cnt) %>%
  ggplot + aes(x =pheno_group, y = cnt, fill = pheno_group) +
  ylim(0, 250) +
  labs(x =NULL, y=NULL) +
    geom_bar(position=position_dodge(width=10), stat="identity", width = 1) +
    geom_text(aes(y = cnt + 15, label = cnt, color = pheno_group), size = 3, show.legend = FALSE)+
    scale_fill_manual(name = NULL, values = pheno_group_color) +
   scale_color_manual(name = NULL, values = pheno_group_color) +
    themes+ theme(legend.position = 'None',
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.text.x = element_blank(),
                  plot.margin = unit(c(1,1,0,0), "cm")) + coord_flip()

p2 <- pheno_info %>%
  filter(!is.na(disease_group)) %>%
  # mutate(disease_group = factor(disease_group, levels = names(icd_names), labels = icd_names))  %>%
  group_by(disease_group) %>%
  dplyr::summarize(cnt = n()) %>%
  arrange(desc(disease_group)) %>%
  mutate(lab.ypos = cumsum(cnt) -0.5*cnt) %>%
  ggplot + aes(x = disease_group, y = cnt, fill = disease_group) +
  labs(x =NULL, y=NULL) +
  geom_bar(position=position_dodge(width=10), stat="identity", width = 1) +
  scale_fill_manual(name = NULL, values = icd_colors) +
  scale_color_manual(name = NULL, values = icd_colors) +
  scale_x_discrete(labels = rev(icd_names)) +
  ylim(0, 30) +
  geom_text(aes(y = cnt + 3, label = cnt, color = disease_group), size = 3, show.legend = FALSE)+
  themes + theme(legend.position = 'None',
                 axis.line.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 plot.margin = unit(c(1,1,1,1), "cm")) + coord_flip()

p3 <- complex_pheno_upset(group_type = 'pheno', test_type = 'burden', max_size = 80, save = save, height = 4, width = 10) +
  ggtitle('(C) Number of genes with associations across combinations of phenotypic domains') +
  theme(plot.title = element_text(face = 'bold', size =9), )
p3

p4 <- complex_pheno_upset(group_type = 'icd', test_type = 'burden', max_size = 4, save = save, height = 4, width = 10)+
  ggtitle('(C) Number of genes with associations across combinations of disease sub-domains') +
  theme(plot.title = element_text(face = 'bold', size =9), )
p4
data_599 <- read_csv(paste0(data_path, 'pleiotropy_2024_gene_burden_sig_cnt_summary_annotated.csv'))
gene_info <- read_delim('~/Dropbox (Partners HealthCare)/github_repo/ukb_exomes/data/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', delim = '\t') %>%
  dplyr::select(gene_symbol = gene, gene_id, pLI, oe_lof_upper_bin, cds_length, oe_lof_upper)

data_599 <- data_599 %>%
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

p51 <- data_599 %>%
  mutate(domain_category = factor(domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain'))) %>%
  filter(domain_category != 'Non-associated') %>%
  # filter(annotation != 'pLoF|missense|LC') %>%
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

annotation_types3 = c('Watanabe et al.2019', annotation_types2)
annotation_names3 = c('Watanabe et al.2019', annotation_names2)
names(annotation_names3) = annotation_types3
colors3 = c(colors2, 'Watanabe et al.2019' = '#1F77B4')
annotation_fill_scale3 = scale_fill_manual(name = 'Annotation', values = colors3, breaks = annotation_types3, labels = annotation_names3)
annotation_color_scale3 = scale_color_manual(name = 'Annotation', values = colors3, breaks = annotation_types3, labels = annotation_names3)

p51_v2 <- data_599 %>%
  dplyr::filter(domain_category != 'Non-associated') %>%
  # filter(annotation != 'pLoF|missense|LC') %>%
  dplyr::group_by(annotation) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::group_by(domain_category, annotation, n) %>%
  dplyr::summarize(cnt = n()) %>%
  dplyr::mutate(prop = cnt/n) %>%
  dplyr::select(-n, -cnt) %>%
  rbind(., data.frame(domain_category = c('Trait-specific', 'Domain-specific', 'Multi-domain'),
                      annotation = rep('Watanabe et al.2019', 3),
                      prop = c(0.1880, 0.1403, 0.6717))) %>%
  mutate(domain_category = factor(domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain'))) %>%

  mutate(annotation = factor(annotation, levels= annotation_types3[c(1,2,5,3,4)])) %>%
  ggplot + aes(x = annotation, y = prop, color = annotation, fill=annotation, alpha = domain_category) +
  geom_bar(stat = "identity") +
  scale_alpha_discrete(name = NULL, range=c(0.25, 1)) +
  scale_x_discrete(labels = annotation_names3) +
  annotation_color_scale3 + annotation_fill_scale3 + themes +
  labs(x = NULL, y=NULL, color = NULL, fill=NULL) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust=1),
        # legend.position = 'right', legend.direction = 'vertical',
        plot.margin = unit(c(1,1,0,0), "cm")) +
  guides(color = 'none', fill='none')

p52 <- data_599 %>%
  mutate(disease_domain_category = factor(disease_domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain'))) %>%
  filter(disease_domain_category != 'Non-associated') %>%
  # filter(annotation != 'pLoF|missense|LC') %>%
  mutate(annotation = factor(annotation, levels= annotation_types2[c(1,4,2,3)])) %>%
  ggplot + aes(x = annotation, color = annotation, fill=annotation, alpha = disease_domain_category) +
  geom_bar(position = "fill") +
  scale_alpha_discrete(name = NULL, range=c(0.25, 1)) +
  scale_x_discrete(labels = annotation_names2) +
  annotation_color_scale2 + annotation_fill_scale2 + themes +
  labs(x = NULL, y=NULL, color = NULL, fill=NULL) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust=1),
        # legend.position = 'right', legend.direction = 'vertical',
        plot.margin = unit(c(1,0,0,0), "cm")) +
  guides(color = 'none', fill='none')


pLI_data <- read_csv(paste0(data_path, 'pleiotropy_2024_gene_burden_sig_cnt_summary_annotated.csv'))
gene_info <- read_delim('~/Dropbox (Partners HealthCare)/github_repo/ukb_exomes/data/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', delim = '\t') %>%
  dplyr::select(gene_symbol = gene, gene_id, pLI, oe_lof_upper_bin, cds_length, oe_lof_upper)

pLI_data <- pLI_data %>%
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

  ) %>%
  mutate(pLI_high = if_else(pLI > 0.9, 1, 0))

label_type = labeller(annotation = annotation_names2)
annotation_color_scale2 = scale_color_manual(name = 'Annotation', values = colors2, breaks = annotation_types2[c(1,2,4,3)], labels = annotation_names2[c(1,2,4,3)])
p6 <- pLI_data %>%
  # filter(annotation != 'pLoF|missense|LC') %>%
  group_by(annotation, domain_category) %>%
  dplyr::summarize(mean = mean(oe_lof_upper, na.rm = TRUE),
                   prop = sum(oe_lof_upper_bin == 0, na.rm = T) / n(),
                   cnt = n()) %>%
  dplyr::mutate(sd = 1.96* sqrt(prop* (1-prop)/cnt)) %>%
  dplyr::mutate(domain_category = factor(domain_category, levels = c('Non-associated', 'Trait-specific', 'Domain-specific', 'Multi-domain'))) %>%
  dplyr::mutate(annotation = factor(annotation, levels= annotation_types2[c(1,2,4,3)])) %>%
  ggplot + geom_pointrange(aes(x = domain_category, y = prop, ymax = prop + sd, ymin = prop -sd, color = annotation)) +
  labs(x = NULL, y = 'Constrained genes (%)') +
  annotation_color_scale2+ themes +
  scale_y_continuous(label=percent) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(~annotation, labeller = label_type, scale='free', nrow = 1)


p <- ggpubr::ggarrange(p1+
                         theme(plot.margin = unit(c(1.5, 0.2, 1.5,0), "cm")), p51_v2 +
                         theme(plot.margin = unit(c(0.3, 0, 0,0), "cm")), ncol=2,widths = c(1, 1.05),
                            labels = c('(A) Number of phenotypes across phenotypic domains',
                                       '(B) Genes with associations across phenotypic domains (%)'),
                            font.label = list(size = 9, color = "black", face = "bold", family = NULL), hjust =0)

p6 <- p6 +
  ggtitle('(D) Proportion of constrained genes v.s. N associations') +
  theme(plot.title = element_text(face = 'bold', size =9, hjust = -0.125), )
# pp <- p/ p3 /p6 + plot_layout(heights = unit(c(3,2,1), c("null","null","null")), widths = unit(c(2,2,0.5), c("null","null","null")))

pp <- p/p3
png(paste0(figure_path, 'figureS13/figureS13_upper.png'), width=7, height=6, units = 'in', res = 300)
print(pp)
dev.off()

p_upper <- ggplot_pdf(image_read(paste0(figure_path, 'figureS13/figureS13_upper.png')))
pp <- ggpubr::ggarrange(p_upper, p6, nrow=2,heights = c(4, 2), font.label = list(size = 9, color = "black", face = "bold", family = NULL), hjust =0)
png(paste0(figure_path, 'figureS13.png'), width=7, height=9, units = 'in', res = 300)
print(pp)
dev.off()


p <- ggpubr::ggarrange(p2+
                         theme(plot.margin = unit(c(0.5, 0.2, 0.5,0), "cm")), p52 +
                         theme(plot.margin = unit(c(0.3, 0, 0,0), "cm")), ncol=2, widths = c(1, 1.05),
                       labels = c('(A) Number of diseases across disease domains',
                                  '(B) Genes with associations across disease domains (%)'),
                       font.label = list(size = 9, color = "black", face = "bold", family = NULL), hjust =0)
p <- p/ p4

png(paste0(figure_path, 'figureS14.png'), width=7, height=6, units = 'in', res = 300)
print(p)
dev.off()
