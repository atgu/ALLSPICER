setwd('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/')
source('~/Dropbox (Partners HealthCare)/ukbb_exomes_pleiotropy/R/constants.R')
if (system.file(package = "inlmisc", lib.loc = .libPaths()) == "")
  utils::install.packages("inlmisc", dependencies = TRUE)

test = 'burden'
label = if_else(test == 'burden', '', '_skato')
## Test results
result <- read_csv(paste0(result_path, paste0('continuous_af1e_4_ac_13_2_p_5e_1_500k', label,'_results.csv')))
beta <- read_csv(paste0(result_path, paste0('continuous_af1e_4_ac_13_2_p_5e_1_500k', label,'_beta.csv')))
gene <- read_delim(paste0(data_path, 'corr_testing_skato_gene_500k.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character()))

pheno_info <- distinct(gene %>% select(phenocode, description))
result <- result %>%
  merge(., pheno_info, by.x = 'pheno1', by.y = 'phenocode') %>%
  merge(., pheno_info, by.x = 'pheno2', by.y = 'phenocode')

beta <- beta %>%
  merge(., pheno_info, by.x = 'pheno1', by.y = 'phenocode') %>%
  merge(., pheno_info, by.x = 'pheno2', by.y = 'phenocode')
name = paste0('continuous_af1e_4_ac_13_2_p_5e_1_500k_', test)

# distinct phenotype pairs
pheno <- result %>% select(pheno1,pheno2)
pheno <- distinct(pheno)
a <- rbind(pheno %>% select(a = pheno1, b = pheno2), pheno %>% select(a = pheno2, b = pheno1))
a <- as.data.frame(table(a$a, a$b))
b <- a %>% filter(Freq > 1)
nrow(pheno) - nrow(b)
# distinct phenotypes
length(c(unique(unique(pheno$pheno1), unique(pheno$pheno2))))
# distinct genes
length(unique(result$gene))
# test summary
table(result$annotation, result$sig_burden)

figure_real_data_qq <- function(results, name=NULL, save=TRUE){
  results <- results %>% 
    filter(n_var>1) %>%
    filter(annotation != 'all') %>%
    group_by(annotation, sig_gene) %>%
    arrange(pvalue) %>%
    add_count() %>%
    mutate(observed = -log10(pvalue), 
           rank = order(pvalue), 
           expected = -(log10(rank / (n+1))),
           annotation = factor(annotation, levels = annotation_types))
  
  figure <- results %>% 
    ggplot+ aes(y=observed,x=expected, color = annotation, label = gene) + 
    geom_point(alpha = 0.5) + 
    geom_abline(intercept = 0, slope = 1) +
    labs(x="Expected -log10(p)", y="Observed -log10(p)", color = 'Annotation') +
    annotation_color_scale + annotation_fill_scale + 
    # xlim(0,max(results$observed)) +
    # ylim(0,max(results$observed)) +
    themes +
    facet_grid(~sig_gene, labeller = label_type) +
    geom_text_repel(max.overlaps = 5)
  
  if(save){
    png(paste0(name, "_qqplot.png"), width=8, height=6, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}
results <- result %>% filter(n_var > 2 & sig_burden == 2 & gene != 'ANK1' & gene != 'NUBP2'  ) %>% 
  mutate(sig_gene = sig_burden)
figure_real_data_qq(results ,paste0(figure_path, name, '_n_var2_filtered'))
figure_real_data_qq(all_result %>% filter(pvalue >0 & n_var > 1 & sig_burden == 2 & annotation != 'all') ,paste0(name, '_n_var1_filtered'), save=T)
figure_real_data_qq(all_result %>% filter(pvalue >0 & n_var > 1 & sig_burden == 2 & annotation != 'all' & !(gene%in% c('SLC39A8', 'STAB1', 'VCAN'))) ,paste0(name, '_n_var1_filtered_remove_slc39a8_stab1_vcan'))


## plof missense analysis
data <- read_delim(paste0(data_path, 'plof_mis_analysis_gene_wise_', test,'.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character()))
sum(data[paste0('combined_sig_', test)] > 0)
sum(data[paste0('combined_sig_', test, '_categorical')] > 0)
sum(data[paste0('combined_sig_', test, '_continuous')] > 0)
sum(data[paste0('combined_sig_', test, '_icd10')] > 0)
table(data[paste0('combined_sig_', test)])
table(data[paste0('combined_sig_', test, '_categorical')])
table(data[paste0('combined_sig_', test, '_continuous')])
table(data[paste0('combined_sig_', test, '_icd10')])

test = 'skato'
lof_mis_table <- read_delim(paste0(data_path, 'plof_mis_only_table_', test,'.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character()))
p_cut <- if_else(test == 'skato', 2.5e-7, 6.7e-7)
sum <- lof_mis_table %>%
  pivot_wider(c(1:2, 4:9), names_from = 'annotation', names_prefix='pvalue_', values_from=if_else(test=='skato', 'Pvalue', 'Pvalue_Burden')) %>%
  filter(`pvalue_pLoF|missense|LC` < p_cut & `pvalue_pLoF` > p_cut &  `pvalue_missense|LC` > p_cut) %>%
  mutate(trait_type = if_else(trait_type %in% c('icd10', 'icd_first_occurrence'), 'icd10', trait_type)) 

write_csv(sum, paste0(result_path, 'plof_mis_only_table_', test,'.csv'))





sum_data <- data %>%
  pivot_longer(6:7, names_prefix = 'combined_sig_') %>%
  mutate(test = str_split(name, '_') %>% map_chr(., 1),
         trait_type = str_split(name, '_') %>% map_chr(., 2),) %>%
  group_by(name, value, trait_type) %>%
  summarise(cnt = n()) 

data %>%
  filter(value > 0) %>%
  ggplot +
  aes(x = trait_type, y = value, color = trait_type, fill = trait_type) +
  geom_boxplot() +
  themes + theme_classic() +
  # scale_y_log10() + 
  trait_color_scale + trait_fill_scale + 
  facet_grid(~ name, scale = 'free')

## ALB figure
save_beta_figure <- function(beta, gene, annotations, pheno1, pheno2, save=TRUE, freq_filter=TRUE){
  sub <- beta %>%
    mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
    filter(
      (gene == gene) &
        (annotation %in% annotations) &
        (pheno1 == pheno1) &
        (pheno2 == pheno2)) %>%
    mutate(
      annotation = factor(annotation, levels = annotation_types, labels = annotation_names), 
      BETA_adjusted = sqrt(2*pop_AF*(1-pop_AF))*BETA
    )%>% 
    select(locus, alleles, phenocode, pop_AF, pop_AC, gene, annotation, BETA, BETA_adjusted, Pvalue)
  
  pvalue_1 <- sub %>% filter(phenocode == pheno1) %>% select(1:2, 4:7, 10)
  pvalue_2 <- sub %>% filter(phenocode == pheno2) %>% select(1:2, 4:7, 10)
  sub <- sub %>% 
    pivot_wider(names_from = phenocode, names_prefix = 'beta_adjusted_', values_from = BETA_adjusted, id_col = c('locus', 'alleles','pop_AF', 'pop_AC', 'gene', 'annotation')) 
  sub <- sub %>% 
    # pivot_wider(names_from = phenocode, names_prefix = 'pvalue_', values_from = Pvalue, id_col = c('locus', 'alleles','AF', 'AC', 'gene', 'annotation')) %>%
    merge(., pvalue_1, by = c('locus', 'alleles','pop_AF', 'pop_AC', 'gene', 'annotation')) %>%
    merge(., pvalue_2, by = c('locus', 'alleles','pop_AF', 'pop_AC', 'gene', 'annotation'))
  
  figure <- sub %>% 
    # mutate(annotation = factor(annotation, levels = annotation_types, labels = annotation_names)) %>%
    mutate(p_size = (-log(Pvalue.x) + -log(Pvalue.y))/2)%>%
    ggplot + 
    aes(x = get(paste0('beta_adjusted_', pheno1)), 
        y = get(paste0('beta_adjusted_', pheno2)), 
        color = annotation, 
        size=p_size, 
        label=locus) +
    geom_point() + 
    geom_abline(slope = 1, intercept = 0, lty =2) + 
    labs(x = 'Albumin', 
         y = 'Calcium', 
         # title = gene, 
         size = expression(bold(-log(Pvalue)[mean]))) + 
    theme_classic() +
    themes + theme(axis.title = element_text(size = 12),
                   axis.text = element_text(size = 8.5),
                   strip.text = element_text(size = 12),
                   legend.position = 'top',
                   legend.title = element_text(size = 12),
                   legend.background = element_blank()
    ) + 
    scale_size_continuous(breaks = c(2,5,10,20)) +
    annotation_color_scale + annotation_fill_scale +
    guides(color = 'none') +
    facet_wrap(~annotation) 
  if(save){
    png(paste0(figure_path, pheno1,'_',pheno2,'_', gene, if_else(freq_filter, '_freq_filtered', '_no_freq_filter'),'_beta_', tranche,'.png'), height = 3, width = 7.5, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}
tranche = '500k'
alb <- read_delim(paste0(data_path, 'alb_albumin_calcium_var_hgvsp_', tranche,'.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character()))
beta = alb
gene = 'ALB'
annotations = c('pLoF', 'missense|LC', 'synonymous')
pheno1 = '30600'
pheno2 = '30680'
save_beta_figure(alb, 'ALB',c('pLoF', 'missense|LC', 'synonymous') ,'30600', '30680', save=T, freq_filter = F)

alb %>% filter(pop_AF < 1e-4) %>%
  save_beta_figure(., 'ALB',c('pLoF', 'missense|LC', 'synonymous') ,'30600', '30680', save=T, freq_filter = T)


## number of pleiotropy variant figure
pleiotropic_cnt_bin_figure(data_type = 'gene', 'burden')
pleiotropic_cnt_bin_figure(data_type = 'gene', 'skato')
pleiotropic_cnt_bin_figure(data_type = 'var', 'burden')
pleiotropic_cnt_bin_figure(data_type = 'var', 'skato')

## number of pleiotropy gene figure
pleiotropic_gene_cnt_figure('gene', 'burden')
pleiotropic_gene_cnt_figure('gene', 'skato')

## gene list figure
gene_list_pleiotropy_figure(test = 'burden', fig_type = 'prop_pleiotropy', save_plot = T)
gene_list_pleiotropy_figure(test = 'skato', fig_type = 'prop_pleiotropy', save_plot = T)
gene_list_pleiotropy_figure(test = 'burden', fig_type = 'mean_prop_categorical', save_plot = T)
gene_list_pleiotropy_figure(test = 'skato', fig_type = 'mean_prop_categorical', save_plot = T)

test = 'skato'
data <- load_ukb_file(paste0('gene_sig_cnt_filtered_', test, '_', tranche,'.txt.bgz'), subfolder = 'analysis/')
sub <- data %>% select(1:3,16, 41:44) %>% filter(all_sig_pheno_cnt > 0 )
sub <- data %>% select(1:3,16, 41:44) %>% filter(all_sig_pheno_cnt > 14 & annotation != 'pLoF|missense|LC')
write_csv(sub, paste0(result_path, 'most_pleiotropic_gene_over_20_', test,'.csv'))

##### Domain level grouping
complex_pheno_upset(group_type = 'biomarker', test_type = 'skato', max_size = 20, save = T, height = 6, width = 15)
complex_pheno_upset(group_type = 'biomarker', test_type = 'burden', max_size = 20, save = T, height = 6, width = 15)
complex_pheno_upset(group_type = 'blood', test_type = 'skato', max_size = 100, save = T, height = 3, width = 8)
complex_pheno_upset(group_type = 'blood', test_type = 'burden', max_size = 75, save = T, height = 3, width = 8)
complex_pheno_upset(group_type = 'pheno', test_type = 'skato', max_size = 60, save = T, height = 6, width = 12)
complex_pheno_upset(group_type = 'pheno', test_type = 'burden', max_size = 35, save = T, height = 6, width = 12)

# complex_pheno_upset(group_type = 'icd', test_type = 'skato', max_size = 1000, save = F, height = 10, width = 10)
# complex_pheno_upset(group_type = 'icd', test_type = 'burden', max_size = 1000, save = T, height = 10, width = 10)

complex_pheno_upset(group_type = 'both', test_type = 'skato', max_size = 16, save = T, height = 6, width = 12)
complex_pheno_upset(group_type = 'both', test_type = 'burden', max_size = 8, save = T, height = 6, width = 12)

### ICD color palette
icd_color_schema(save_plot = T)



