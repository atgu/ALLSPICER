source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')

raw_results_500k <- read_pleiotropy_results('burden', '500k') 
results_500k <- modify_results_table(raw_results_500k, 'burden', '500k')

figure_real_data_qq <- function(results, name=NULL, save=TRUE){
  results <- results %>% 
    group_by(annotation) %>%
    arrange(pvalue) %>%
    add_count() %>%
    mutate(observed = -log10(pvalue), 
           rank = order(pvalue), 
           expected = -(log10(rank / (n+1))),
           annotation = factor(annotation, levels = annotation_types))
  mx <- as.numeric(max(max(results[results$pvalue>0, 'observed']), 
                       max(results[results$pvalue>0, 'expected'])))
  
  figure <- results %>% 
    ggplot+ aes(y=observed,x=expected, color = annotation, label = gene) + 
    geom_point(alpha = 0.5) + 
    geom_abline(intercept = 0, slope = 1) +
    geom_hline(yintercept = -log10(4.23e-6), lty=2) +
    labs(x=expression(Expected -log[10](p)), y=expression(Observed -log[10](p)), color = 'Annotation') +
    annotation_color_scale + annotation_fill_scale + 
    xlim(0, mx) +
    ylim(0, mx) +
    themes +
    facet_grid(~annotation, labeller = label_type)
  # geom_text_repel(max.overlaps = 2)
  
  if(save){
    png(paste0(figure_path, name, "_qqplot.png"), width=5, height=3.5, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

p1 <- figure_real_data_qq(results_500k %>% filter(n_cases1 > 300000 & n_cases2 > 300000), save = F)
p2 <- results_500k %>%
  mutate(annotation = factor(annotation, levels = annotation_types)) %>%
  ggplot + 
  aes(x = corr, y = -log10(pvalue), color = annotation, size = n_var) +
  labs(x = 'Phenotypic correlation', y = expression(-log[10](p)), size = 'Number of variants', color = 'Annotation') +
  geom_point(alpha = 0.5) +
  annotation_color_scale + 
  # scale_y_log10() + 
  geom_hline(yintercept = -log10(0.05/11810), lty = 2) + 
  # geom_vline(xintercept = c(-0.1, 0.1, 0.8, 1), lty = 2) + 
  scale_size_continuous(range = c(0.01, 4)) + 
  facet_grid(~annotation, labeller = labeller(annotation = annotation_names)) + themes + theme(legend.position = 'top')

figure <- ggpubr::ggarrange(p1 +
                              theme(axis.title = element_text(face = 'plain', size = 11), 
                                    plot.margin = unit(c(1,0,0,0.5), "cm")), 
                            p2+
                              theme(axis.title = element_text(face = 'plain', size = 11), 
                                    plot.margin = unit(c(0.7,0,0,0.5), "cm")), 
                            labels = c('(A) QQ plots of ALLSPICE test results across high-quality phenotypes', '(B) Relationship between phenotypic correlation and ALLSPICE pvalue'),
                            nrow=2, vjust = 2, hjust = 0, font.label = list(size = 10, color = "black", face = "bold", family = NULL), heights = c(0.2, 0.18), common.legend = T)
png(paste0(figure_path,'figureS13.png'), height = 6, width = 8, units = 'in', res = 300)
print(figure)
dev.off()
