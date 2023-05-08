setwd('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/')
source('~/Dropbox (Partners HealthCare)/ukbb_exomes_pleiotropy/R/constants.R')

TEST = 'skato'
# TEST = 'burden'
all_gene_sig <- read_delim(paste0(data_path, 'gene_sig_', TEST,'_all.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character()))
## Method (1): correlation test results
### high power phenotypes with > 300k cases + only AF < 0.0001 filter
result <- read_csv(paste0(result_path, paste0('continuous_pop_AF_1e_4_500k_syn_var_corr_', TEST,'_results.csv')))
beta <- read_csv(paste0(result_path, paste0('continuous_pop_AF_1e_4_500k_syn_var_corr_', TEST,'_beta.csv')))
gene <- read_delim(paste0(data_path, 'corr_testing_', TEST,'_gene_n_cases_over_300k.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character()))

gene <- gene %>%
  mutate(annotation = factor(annotation, levels = annotation_types, labels = annotation_names))
result2 <- result %>%
  select(-sig_gene.x, -sig_gene.y, -sig_gene) %>%
  merge(., gene, by.x = c('annotation', 'pheno1', 'gene'), 
        by.y = c('annotation', 'phenocode', 'gene_symbol'), all.x = T) %>%
  merge(., gene, by.x = c('annotation', 'pheno2', 'gene'), 
        by.y = c('annotation', 'phenocode', 'gene_symbol'), all.x = T) %>%
  mutate(sig_gene.x = if_else(Pvalue.x <  2.5e-7, 1, 0),
         sig_gene.y = if_else(Pvalue.y <  2.5e-7, 1, 0),) %>%
  mutate(sig_gene = sig_gene.x + sig_gene.y) %>%
  filter(!is.na(sig_gene.x) & !is.na(sig_gene.y))

### all continuous phenotypes + stringent filters on variant
result <- read_csv(paste0(result_path, paste0('continuous_af1e_4_ac_13_2_p_5e_1_500k', label,'_results.csv')))
beta <- read_csv(paste0(result_path, paste0('continuous_af1e_4_ac_13_2_p_5e_1_500k', label,'_beta.csv')))
gene <- read_delim(paste0(data_path, 'corr_testing_skato_gene_500k.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character()))

pheno_info <- distinct(gene %>% select(phenocode, description, n_cases))
result <- result %>%
  merge(., pheno_info, by.x = 'pheno1', by.y = 'phenocode') %>%
  merge(., pheno_info, by.x = 'pheno2', by.y = 'phenocode')

beta <- beta %>%
  merge(., pheno_info, by.x = 'pheno1', by.y = 'phenocode') %>%
  merge(., pheno_info, by.x = 'pheno2', by.y = 'phenocode')
name = paste0('continuous_pop_AF_1e_4_500k_high_power_pheno_', TEST)

syn_result <- result2 %>% 
  filter(annotation == 'Synonymous' & pvalue < 0.05 & n_var >1 & sig_gene ==2) %>%
  select(pheno1, pheno2, annotation, gene, pvalue, corr) %>% 
  merge(., beta %>%
          mutate(annotation = factor(annotation, levels = annotation_types, labels = annotation_names)), by = c('pheno1', 'pheno2', 'annotation', 'gene'), all.x = T) %>%
  mutate(pheno_tmp2 = if_else(pheno1 > pheno2, pheno2, pheno1),
         pheno_tmp1 = if_else(pheno1 > pheno2, pheno1, pheno2)) %>%
  mutate(pheno1 = if_else(pheno_tmp1 > pheno_tmp2, pheno_tmp2, pheno1),
         pheno2 = if_else(pheno_tmp1 > pheno_tmp2, pheno_tmp1, pheno2),)
syn_genes = unique(syn_result$gene)
syn_gene = 'NUBP2'
for(syn_gene in syn_genes){
  figure <- syn_result %>%
    filter(gene == syn_gene)%>%
    mutate(b1_adjusted = sqrt(2*AF*(1-AF))*b1, 
           b2_adjusted = sqrt(2*AF*(1-AF))*b2) %>% 
    ggplot + aes(x = b1, y = b2) +
    geom_point() +
    #geom_abline(slope = 1, intercept = 0, lty =2) + 
    labs(title = syn_gene) +
    #ylim(-0.045, 0.03) + 
    #facet_grid(pheno1~pheno2) + 
    facet_wrap(pheno1~pheno2, nrow = 4) +
    theme(strip.text.y = element_text(angle = 0)) + 
    geom_text(aes(label=paste0(' corr:',round(corr,2), '\n', 'pvalue:', format(pvalue, scientific=T, digit=2))), x = -Inf, y= Inf, hjust = -0.05, vjust = 1, size =3)
    # geom_text(aes(label=paste0(' corr:',round(corr,2), '\n', 'pvalue:', format(pvalue, scientific=T, digit=2))), x = Inf, y= -Inf, hjust = 1, vjust = -0.3, size =3)
    #geom_text(aes(label=paste0(' corr:',round(corr,2), '\n', 'pvalue:', format(pvalue, scientific=T, digit=2))), x = -Inf, y= -Inf, hjust = -0.05, vjust = -0.3, size =3)
  
  # png(paste0(figure_path, syn_gene, 'synonymous_sig_pheno_pair_beta_', tranche,'.png'),
  #     height = 2*4^length(unique(syn_result%>%
  #                                   filter(gene == syn_gene) %>%
  #                                   select(pheno1))),
  #     width = 2*4^length(unique(syn_result%>%
  #                                 filter(gene == syn_gene) %>%
  #                                 select(pheno2))), units = 'in', res = 300)
  png(paste0(figure_path, syn_gene, 'synonymous_sig_pheno_pair_beta_', tranche,'.png'), height = 8, width = 20, units = 'in', res = 300)
  print(figure)
  dev.off()
  
}


figure_real_data_qq <- function(results, name=NULL, save=TRUE){
  results <- results %>% 
    # filter(n_var>1) %>%
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
results <- result2 %>% 
  distinct() %>%
  #filter(!(gene%in% c('PIEZO1', 'ZFAT')) & pvalue >0) %>%
  # mutate(sig_gene = sig_burden) %>%
 #  filter(!(gene%in% c('NUBP2', 'ANK1')& annotation == 'synonymous')) %>%
 # filter(pvalue > 0 ) %>%
  filter(sig_gene == 2) %>%
  filter(n_var > 1 & annotation != 'all') 
# number of triplets
nrow(results)
# distinct phenotype pairs
pheno <- results %>% select(pheno1,pheno2) %>%
  mutate(pheno_tmp2 = if_else(pheno1 > pheno2, pheno2, pheno1),
         pheno_tmp1 = if_else(pheno1 > pheno2, pheno1, pheno2)) %>%
  mutate(pheno1 = if_else(pheno_tmp1 > pheno_tmp2, pheno_tmp2, pheno1),
         pheno2 = if_else(pheno_tmp1 > pheno_tmp2, pheno_tmp1, pheno2),) %>%
  select(pheno1, pheno2)
pheno <- distinct(pheno) 
nrow(pheno)
# distinct phenotypes
length(c(unique(unique(pheno$pheno1), unique(pheno$pheno2))))
# distinct genes
length(unique(results$gene))
# test summary
table(results$annotation, results$sig_gene)

figure_real_data_qq <- function(results, name=NULL, save=TRUE){
  results <- results %>% 
    # filter(annotation != 'all') %>%
    # group_by(annotation, sig_burden) %>%
    group_by(annotation) %>%
    arrange(pvalue) %>%
    add_count() %>%
    mutate(observed = -log10(pvalue), 
           rank = order(pvalue), 
           expected = -(log10(rank / (n+1))),)
           # annotation = factor(annotation, levels = annotation_types))
  mx <- as.numeric(max(max(results[results$pvalue>0, 'observed']), 
                       max(results[results$pvalue>0, 'expected'])))
  
  figure <- results %>% 
    ggplot+ aes(y=observed,x=expected, color = annotation, label = gene) + 
    geom_point(alpha = 0.5) + 
    geom_abline(intercept = 0, slope = 1) +
    labs(x="Expected -log10(p)", y="Observed -log10(p)", color = 'Annotation') +
    annotation_color_scale + annotation_fill_scale + 
    xlim(0, mx) +
    ylim(0, mx) +
    themes 
  # facet_grid(~sig_burden, labeller = label_type, scale = 'free') +
  # geom_text_repel(max.overlaps = 5)
  
  if(save){
    png(paste0(name, "_qqplot.png"), width=4, height=3, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}
name = paste0('continuous_pop_AF_1e_4_500k_syn_var_corr_', TEST)
figure_real_data_qq(results ,paste0(figure_path, name), save = F)
figure_real_data_qq(results %>% filter(n_var > 1 & pvalue > 0) ,paste0(figure_path, name, '_n_var1_p_filtered'), save = T)
figure_real_data_qq(results %>% filter(n_var > 1 & !(gene%in% c('NUBP2', 'ANK1')& annotation == 'synonymous') ) ,paste0(figure_path, name, '_n_var1_filtered_remove_ank1_nubp2_syn'), save = T)
figure_real_data_qq(results %>% filter(n_var > 1 & !(gene%in% c('NUBP2', 'ANK1')& annotation == 'synonymous') & pvalue >0) ,paste0(figure_path, name, '_n_var1_p_filtered_remove_ank1_nubp2_syn'), save = T)

figure_real_data_qq(results %>% filter(n_var > 1 & !(gene%in% c('PIEZO1', 'ZFAT')) & pvalue >0) ,paste0(figure_path, name, '_n_var1_pvalue_filtered_remove_ank1_nubp2'), save = F)


## QCed results
gene <- load_ukb_file(paste0('gene_sig_cnt_filtered_', test, '_', tranche,'.txt.bgz'), subfolder = 'analysis/') %>%
  select(1:3, 16, 41:44)
var <- load_ukb_file(paste0('var_sig_cnt_filtered_', test, '_', tranche,'.txt.bgz'), subfolder = 'analysis/')


## Gene point check
gene <- read_delim(paste0(data_path, 'gene_sig_', test,'_all.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character()))

## Variant point check
var <- read_delim(paste0(data_path, 'pleiotropic_var_point_check.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character())) %>%
  filter(Pvalue < 8e-9)


## ATM check
atm <- rbind(read_delim(paste0(data_path, 'atm_sig_phenos_skato.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character())) %>% 
               mutate(test = 'skato'),
             read_delim(paste0(data_path, 'atm_sig_phenos_burden.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character())) %>% 
               mutate(test = 'burden') %>% 
               rename(Pvalue = Pvalue_Burden))

## Domain level association check
both_data <- formatted_data('both', 'burden', write = F)
both_data <- formatted_data('both', 'skato', write = T)
icd_data <- formatted_data('icd', 'skato', write=F)
icd_data <- formatted_data('icd', 'burden', write=T)
pheno_data1 <- formatted_data('pheno', 'burden', write = F)
pheno_data2 <- formatted_data('pheno', 'skato', write = F)

pheno_data1$sum <-rowSums(pheno_data1[,6:13])

## Domain level gene examples check
test = 'burden'
gene <- read_delim(paste0(data_path, 'gene_sig_pheno_domains_', test,'.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character()))


## pLoF + missense genes
genes <- c('MPO', 'PKD1', 'STAB1', 'ABCC2', 'ANK1', 'CREB3L3', 'DGAT2', 'ABCA7', 
           'PTPRH', 'JAK2', 'MYH9', 'ATP11C', 'PIEZO1', 'GHSR', 'NEURL2', 'SAG', 
           'G6PC', 'LOXL2', 'CD36', 'GMPR', 'FBN2', 'SETD1B', 'ATM', 'PDE3B', 'NPR2', 
           'NLRP3', 'NOSTRIN', 'COQ4', 'PCSK9', 'THBS3', 'CKAP2L', 'CHD3', 'PKHD1', 
           'INSR', 'WDR6, CHEK2', 'INSC', 'COL27A1', 'MICA', 'TRIM40', 'CNPY2')

sub_gene <- gene %>%
  filter(gene_symbol %in% genes) %>%
  filter(annotation %in% c('pLoF', 'missense|LC')) %>%
  select(1,2,3,5,7,9,12) %>%
  pivot_wider(id_cols = c(1:2, 4:6), names_from = 'annotation', values_from = 'Pvalue_Burden') %>%
  filter(is.na(pLoF) | is.na(`missense|LC`))
write_csv(sub_gene, paste0(result_path, 'allelic_series_pLoF_missense_burden.csv'))         
