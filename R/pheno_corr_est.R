setwd('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/')
source('~/Dropbox (Partners HealthCare)/ukbb_exomes_pleiotropy/R/constants.R')
source('~/Dropbox (Partners HealthCare)/ukbb_exomes_pleiotropy/R/simulations.R')

pheno_corr <- read_delim(paste0(data_path, 'correlation_table_phenos_500k.txt.bgz'), delim = '\t',
                         col_types = cols(i_pheno = col_character(), j_pheno = col_character())) %>%
  select(3:5) %>%
  mutate(corr=entry)

phenotypes <- c("21001", "30000", "30010", "30080", "30510", "50" )

var_data <- read_delim(paste0(data_path, 'corr_estimate_burden_syn_var_9phenos.txt.bgz'), delim = '\t') %>%
  mutate(SD_BETA  = BETA * sqrt(2*AF*(1-AF))) %>% 
  filter(phenocode %in% phenotypes) %>%
  filter(AF < 0.0001)
sub_corr <- pheno_corr %>% filter((i_pheno %in% phenotypes) & (j_pheno %in% phenotypes) & (i_pheno != j_pheno))
pheno_info <- distinct(var_data %>% select(phenocode, n_cases, description) %>% filter(phenocode %in% phenotypes))

phenopair <-  c("21001", "30000")
get_real_pheno_pair_data <- function(var_data, phenopair){
  sub <- var_data %>% 
    filter(phenocode %in% phenopair) %>% 
    pivot_wider(id_col = c('locus', 'alleles', 'annotation'), names_from = 'phenocode', values_from = 'BETA')
  # sub <- var_data %>% 
  #   select(locus, alleles, annotation, AF, AC, gene) %>% 
  #   unique() %>% merge(., sub, by = c('locus', 'alleles', 'annotation')) %>% 
  #   # filter(complete.cases(.)) %>%
  #   mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation))
  return(sub)
}

corr_result <- data.frame()
for(i in (1:(length(phenotypes)-1))){
  for(j in ((i+1):(length(phenotypes)))){
    print(c(i, j))
    phenopair <- phenotypes[c(i,j)]
    corr <- c(sub_corr %>% filter((i_pheno == phenopair[1]) & (j_pheno == phenopair[2])) %>% select(corr))
    sub <- get_real_pheno_pair_data(var_data, phenopair)
    cor_test <- cor.test(unlist(sub[, phenopair[1]]), unlist(sub[, phenopair[2]]))
    tmp_result <- c(pheno1 = phenopair[1], pheno2 = phenopair[2], true_corr = corr, est_corr = cor_test$estimate, est_pvalue = cor_test$p.value )
    corr_result <- rbind(corr_result, tmp_result)
  }
}
corr_result <- corr_result %>%
  merge(., pheno_info, by.x = 'pheno1', by.y = 'phenocode') %>%
  merge(., pheno_info, by.x = 'pheno2', by.y = 'phenocode', suffixes = c("1", "2"))

write_csv(corr_result, paste0(result_path, 'corr_est_syn_vars_9pheno_per_sd_af_1e_4.csv'))

figure <- corr_result %>%
  ggplot + aes(x = true_corr.corr, y = est_corr.cor, label = paste0('(', pheno1, ',', pheno2, ')')) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, lty =2) + 
  labs(x = "True phenotypic correlation", y = "Estimated correlation \n (synonymous per-s.d. effect sizes)") + 
  geom_text_repel(direction = 'y', size = 2.5)

png(paste0(figure_path, 'corr_est_syn_vars_9pheno_per_sd_af_1e_4.png'), height = 4, width = 5, units = 'in', res = 300)
print(figure)
dev.off()


## Deviated phenotypes
pheno_corr <- read_delim(paste0(data_path, 'corr_estimate_syn_var_deviation02_500k.txt.bgz'), delim = '\t',
                         col_types = cols(i_phenocode = col_character(), j_phenocode = col_character())) %>%
  select(1:3, 5, 12) %>%
  mutate(corr=entry)

result_450k_old <- read_csv(paste0(result_path, paste0('continuous_pop_AF_1e_4_500k_high_power_pheno_burden_results.csv'))) %>% 
  filter(n_var > 1 & sig_gene ==2 ) %>% filter(annotation != 'all') %>%
  merge(., pheno_corr, by.x = c('pheno1', 'pheno2'), by.y = c('i_phenocode', 'j_phenocode')) %>%
  select(1:9, 16)

result_450k_new <- read_csv(paste0(result_path, paste0('continuous_pop_AF_1e_4_500k_syn_var_corr_burden_results.csv'))) %>% 
  filter(n_var > 1 & sig_gene ==2 ) %>% filter(annotation != 'all') %>%
  merge(., pheno_corr, by.x = c('pheno1', 'pheno2'), by.y = c('i_phenocode', 'j_phenocode')) %>%
  select(1:8, 'entry.x')

pheno_corr_comparison <- result_450k_new %>%
  merge(., result_450k_old, by = c('pheno1', 'pheno2', 'gene', 'annotation'), suffixes = c('_new', '_old')) %>%
  mutate(corr_diff = abs(entry.x_new - entry.x_old),
         pvalue_rate = pvalue_new/pvalue_old)

pheno_corr_comparison %>%
  ggplot + aes(x = pvalue_old, y = pvalue_new, size = corr_diff, color = pheno1) + 
  geom_point(position = position_dodge(), alpha = 0.5) +
  scale_size_continuous(range = c(1,10))
