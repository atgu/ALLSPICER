source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')
library(emojifont)

# ALB results AC 5
## Leave one out panels
raw_results_500k <- read_csv('~/Desktop/ALLSPICE/continuous_ALL_AC_5_burden_syn_var_corr_500k_2024_results_corr.csv')

results_500k <- modify_results_table(raw_results_500k, 'burden', '500k')

leave_one_out_results <- read_csv(paste0('~/Desktop/ALLSPICE/continuous_ALL_AC_5_burden_syn_var_leave_one_out_500k_2024_ALB_results.csv'))
top_hits <- results_500k %>% filter(pvalue < 1e-6)
gene_annts <- top_hits %>% merge(., results_500k, by = colnames(results_500k)[c(1, 3:5, 7:9, 11)], suffixes = c('.x', ''), all.x = T) %>% select("gene", "annotation", "pheno1", "description1", "pheno2", "description2")
gene_top_hits <- top_hits %>% filter(gene == gene_name)

gene_lab <- 'ALB'
pheno1_lab <- 'continuous_30600_both_sexes__irnt'
pheno2_lab <- 'continuous_30680_both_sexes__irnt'
tmp <- gene_top_hits[1, ]
gene_name <- tmp[,'gene']
phenocode1 <- tmp[,'pheno1']
phenocode2 <- tmp[,'pheno2']
pheno1_name <- tmp$description1
pheno2_name <- tmp$description2
threshold <- 0.05
c_hat <- results_500k %>%
  filter(gene == gene_lab & pheno1 == pheno1_lab & pheno2 == pheno2_lab) %>%
  select(annotation, c_hat)

sub_info <- alb_info  %>% filter(AC < 5) %>%
  filter(gene == gene_lab & phenoname %in% c(pheno1_lab, pheno2_lab))
wide_info <- sub_info %>%
  mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
  filter(annotation %in% c('pLoF', 'missense|LC', 'synonymous')) %>%
  mutate(
    BETA_adjusted = sqrt(2*AF*(1-AF))*BETA
  ) %>%
  pivot_wider(names_from = phenoname,  values_from = c('BETA', 'Pvalue'), id_cols = c('locus', 'alleles','AF', 'gene', 'annotation')) %>%
  mutate(significance = case_when(
    get(paste0('Pvalue_',phenocode1)) <= threshold & get(paste0('Pvalue_',phenocode2)) <= threshold ~ 'Both',
    get(paste0('Pvalue_',phenocode1)) > threshold & get(paste0('Pvalue_',phenocode2)) <= threshold ~ pheno2_name,
    get(paste0('Pvalue_',phenocode1)) <= threshold & get(paste0('Pvalue_',phenocode2)) > threshold ~ pheno1_name,
    get(paste0('Pvalue_',phenocode1)) > threshold & get(paste0('Pvalue_',phenocode2)) > threshold ~ 'None',
  )) %>%
  mutate(significance = factor(significance, levels = c('Both', pheno1_name, pheno2_name, 'None')),
         annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
  merge(., leave_one_out_results %>% filter(pheno1 == phenocode1 & pheno2 == phenocode2) %>%
          select(locus, alleles, annotation, leave_one_out_pvalue = pvalue, leave_one_out_c_hat=c_hat), by = c('locus', 'alleles', 'annotation'), all.x=T) %>%
  merge(., raw_results_500k %>% filter(gene==gene_name & pheno1==phenocode1 & pheno2==phenocode2) %>%
          mutate(annotation = factor(annotation, levels=annotation_types)) %>% select(annotation, pvalue, c_hat), by = 'annotation') %>%
  mutate(
    magnitude_change_p = if_else(is.na(leave_one_out_pvalue), 0, abs(log10(leave_one_out_pvalue/if_else(pvalue==0, 1e-320, pvalue)))),
    magnitude_change_c = if_else(is.na(leave_one_out_c_hat), 0, abs(log10(leave_one_out_c_hat/c_hat))),
    p_change_direction = if_else(leave_one_out_pvalue > pvalue, 'less significant', 'more significant')
    ) %>%
  mutate(annotation = factor(annotation, levels=annotation_types)) %>%
  filter(complete.cases(.))
p2 <- wide_info %>%
  ggplot +
  aes(x=get(paste0('BETA_',phenocode1)), y=get(paste0('BETA_',phenocode2)), color = annotation)  +
  geom_point(aes(pch = significance, size = magnitude_change_c)) +
  # geom_abline(data = c_hat %>%
  #               mutate(annotation = factor(annotation, levels=annotation_types)), aes(slope = 1/c_hat, intercept = 0, color = annotation), lwd =0.5) +
  geom_vline(xintercept = 0, lty=2, lwd = 0.25) +
  geom_hline(yintercept = 0, lty=2, lwd = 0.25) +
  annotation_color_scale + annotation_fill_scale +
  labs(x=paste0(pheno1_name), y=pheno2_name, title = NULL) +
  scale_shape_manual(name=paste0('Nominal significance (', threshold, ')'), breaks = c('Both', pheno1_name, pheno2_name, 'None'), values=c("\u25CF", "\u25D0","\u25D1", "\u25CB")) +
  scale_size(range = c(2.5, 8)) +

  geom_text_repel(data = raw_results_500k %>% filter(gene==gene_name & pheno1==phenocode1 & pheno2==phenocode2) %>%
                    mutate(annotation = factor(annotation, levels=annotation_types)), aes(x=2, y= -2, label=formatC(pvalue, format = "e", digits = 2), color=annotation), vjust = 1, size =5)+
  facet_wrap(~annotation, labeller = label_type) +
  guides(alpha = guide_legend(order = 1),
         shape = guide_legend(order = 2),
         color = "none", size = 'none') +
  theme(legend.position = 'top',
        legend.box = 'vertical',
        plot.margin = unit(c(0.5,0,0,0), "cm"))

p1 <- wide_info %>%
  ggplot +
  aes(x=get(paste0('BETA_',phenocode1)), y=get(paste0('BETA_',phenocode2)), color = annotation)  +
  geom_point(aes(pch = significance, size = magnitude_change_p, alpha=p_change_direction)) +
  # geom_abline(data = c_hat %>%
  #               mutate(annotation = factor(annotation, levels=annotation_types)), aes(slope = 1/c_hat, intercept = 0, color = annotation), lwd =0.5) +
  geom_vline(xintercept = 0, lty=2, lwd = 0.25) +
  geom_hline(yintercept = 0, lty=2, lwd = 0.25) +
  annotation_color_scale + annotation_fill_scale +
  labs(x=paste0(pheno1_name), y=pheno2_name, title = NULL) +
  scale_shape_manual(name=paste0('Nominal significance (', threshold, ')'), breaks = c('Both', pheno1_name, pheno2_name, 'None'), values=c("\u25CF", "\u25D0","\u25D1", "\u25CB")) +
  scale_size(range = c(2.5, 8)) +
  scale_alpha_discrete(name = 'Pvalue change direction', range = c(1, 0.2)) +
  geom_text_repel(data = raw_results_500k %>% filter(gene==gene_name & pheno1==phenocode1 & pheno2==phenocode2) %>%
                    mutate(annotation = factor(annotation, levels=annotation_types)), aes(x=2, y= -2, label=formatC(pvalue, format = "e", digits = 2), color=annotation), vjust = 1, size =5)+
  facet_wrap(~annotation, labeller = label_type)  +
  guides(alpha = guide_legend(order = 1),
         shape = guide_legend(order = 2),
         color = "none", size = 'none') +
  theme(legend.position = 'top',
        legend.box = 'vertical',
        plot.margin = unit(c(0.5,0,0,0), "cm"))

figure = ggpubr::ggarrange(p1, p2, nrow=2, heights = c(0.1, 0.1),
                           labels = c('(A) Magnitude of P-value change leaving each variant out when running ALLSPICE',
                             '(B) Magnitude of c_hat change leaving each variant out when running ALLSPICE'), hjust = 0,
                           font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'), common.legend = TRUE
)

png(paste0(figure_path,'figureS18.png'), height = 5, width = 8, units = 'in', res = 300)
print(figure)
dev.off()
