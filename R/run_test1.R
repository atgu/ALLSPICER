source('~/Dropbox (Partners HealthCare)/ukbb_exomes_pleiotropy/R/simulations.R')
source('~/Dropbox (Partners HealthCare)/ukbb_exomes_pleiotropy/R/constants.R')
setwd('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/')

TEST = 'burden' # or 'skato'
# TEST = 'skato' # or 'skato'
CURRENT_TRANCHE = '500k'

pheno_corr <- read_delim(paste0(data_path, 'corr_estimate_syn_var_full_', CURRENT_TRANCHE, '.txt.bgz'), delim = '\t',
                         col_types = cols(i_phenocode = col_character(), j_phenocode = col_character())) %>%
  select(1:3, 5, 12) %>%
  mutate(corr=entry)
# var_file <- paste0(data_path, 'corr_testing_skato_var_n_cases_over_300k.txt.bgz')
# gene_file <- paste0('~/Downloads/corr_testing_burden_gene_n_cases_over_300k.txt.bgz')
var_file <- paste0('~/Downloads/corr_testing_', TEST, '_var_', CURRENT_TRANCHE, '_new_AF.txt.bgz')
gene_file <- paste0('~/Downloads/corr_testing_', TEST, '_gene_', CURRENT_TRANCHE, '.txt.bgz')
print(paste("reading variant data:", var_file))
var_data <- fread(cmd = paste0('gunzip -cq ', var_file))
print(paste("reading gene data:", gene_file))
gene_data <- fread(cmd = paste0('gunzip -cq ', gene_file))
name <- paste0('continuous_new_AF_1e_4_', TEST, '_syn_var_corr_', CURRENT_TRANCHE, '_oct2022_rerun')

print('Formatting gene data...')
if(TEST == 'burden'){
  gene_data <- gene_data %>%
    mutate(Pvalue = Pvalue_Burden)
}

# if(CURRENT_TRANCHE == '500k'){
#   gene_data <- gene_data %>%
#     mutate(n_cases = n_cases_defined)
# }

gene_data <- gene_data %>%
  filter(annotation != 'pLoF|missense|LC') %>%
  select(gene_symbol, annotation, phenocode, Pvalue, n_cases = n_cases_defined) %>%
  mutate(sig_gene = if_else(Pvalue< if_else(TEST == 'burden', 6.7e-7, 2.5e-7), 1, 0))

print('Formatting variant data...')
sub_data <- var_data %>%
  mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation))
genes <- unique(sub_data %>% select(gene))

test_results <- data.frame()
test_beta <- data.frame()
for(k in genes$gene){
  print(k)
  sub <- sub_data %>% filter(gene==k)
  sub$phenocode <- as.character(sub$phenocode)
  gene_sig <- gene_data %>%
    filter(gene_symbol==k & sig_gene == 1)
  phenos <- unique(c(unique(gene_sig$phenocode)))
  if(length(phenos)==0) next
  test_data <- get_real_data(sub, gene_data, pheno_corr, phenos)
  test <- get_test_result(test_data$sub,
                          r = test_data$corr,
                          pheno_corr = pheno_corr,
                          pheno_list = phenos,
                          n_ind = test_data$n_ind,
                          gene=k)
  results <- test$results
  beta <- test$beta
  if(!is.null(results)){
    gene_result <- get_gene_level_data(data = gene_data, gene = k, phenolist = phenos)
    test_full <- add_gene_level_data(gene_result, results)

    test_results <- rbind(test_results, test_full)
    test_beta <- rbind(test_beta, beta)
  }

}


write_data_fig_real(test_results, test_beta, name)




