source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/allelic_heterogeneity_test/utils.R')
source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')
setwd('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/')

TEST = 'burden' # or 'skato'
# TEST = 'skato' # or 'skato'
CURRENT_TRANCHE = '500k'

pheno_corr <- read_delim(paste0(data_path, 'corr_estimate_syn_var_full_', CURRENT_TRANCHE, '.txt.bgz'), delim = '\t',
                         col_types = cols(i_phenocode = col_character(), j_phenocode = col_character()))
pheno_corr<- pheno_corr %>%
  mutate(i_coding = if_else(is.na(i_coding), '', i_coding),
         j_coding = if_else(is.na(j_coding), '', j_coding),) %>%
  mutate(corr=entry,
         i_phenoname= paste0(i_trait_type, '_', i_phenocode, '_',  i_pheno_sex, '_',  i_coding, '_',  i_modifier),
         j_phenoname= paste0(j_trait_type, '_', j_phenocode, '_',  j_pheno_sex, '_',  j_coding, '_',  j_modifier),
         )
var_file <- paste0('~/Downloads/corr_testing_burden_var_0.01_500k.txt.bgz')
gene_file <- '~/Downloads/corr_testing_burden_gene_500k.txt.bgz'
print(paste("reading variant data:", var_file))
var_data <- fread(cmd = paste0('gunzip -cq ', var_file)) 
print(paste("reading gene data:", gene_file))
gene_data <- fread(cmd = paste0('gunzip -cq ', gene_file))
name <- 'continuous_ALL_AF_1e_4_burden_syn_var_c_hat_500k_2024'

print('Formatting gene data...')
gene_data <- gene_data %>%
  mutate(Pvalue = Pvalue_Burden) %>%
  mutate(phenocode = as.character(phenocode)) %>%
  mutate(coding = if_else(is.na(coding), '', coding)) %>%
  mutate(phenoname = paste0(trait_type, '_', phenocode, '_',  pheno_sex, '_',  coding, '_',  modifier)) %>%
  filter(annotation != 'pLoF|missense|LC') %>%
  select(gene_symbol, annotation, trait_type, phenocode, pheno_sex, coding, modifier, description, Pvalue, n_cases = n_cases_defined, phenoname) %>%
  mutate(sig_gene = if_else(Pvalue< 2.5e-6, 1, 0))

print('Formatting variant data...')
sub_data <- var_data %>%
  mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
  mutate(coding = if_else(is.na(coding), '', coding)) %>%
  mutate(phenoname = paste0(trait_type, '_', phenocode, '_',  pheno_sex, '_',  coding, '_',  modifier))
genes <- unique(sub_data %>% select(gene)) %>% unlist()

test_results <- data.frame()
test_beta <- data.frame()
genes2 <- genes
for(k in genes2){
  output_path <- paste0('~/Desktop/ALLSPICE/results/syn_c_hat/', str_replace(name, 'ALL', k), '_results.csv') 
  if(file.exists(output_path)){
    test_full <- read_csv(output_path) %>%
      filter(sig_gene == 2)
    test_results <- rbind(test_results, test_full)
    next
  }
  print(k)
  sub <- sub_data %>% filter(gene==k)
  sub$phenocode <- as.character(sub$phenocode)
  gene_sig <- gene_data %>%
    filter(gene_symbol==k & sig_gene == 1)
  phenos <- gene_sig %>%
    select(trait_type, phenocode, pheno_sex, coding, modifier, phenoname) %>%
    distinct() 
  print(paste0('Number of unique phenotypes:', nrow(phenos)))
  if(length(phenos)<2) next
  test_data <- get_real_data(sub, gene_data, pheno_corr, phenos)
  if(nrow(test_data$sub)==0) next
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
    # test_beta <- rbind(test_beta, beta)
    result_path <- '~/Desktop/ALLSPICE/results/syn_c_hat/'
    write_data_fig_real(test_full, str_replace(name, 'ALL', k), save_fig = F)
  }

}

result_path <- '~/Desktop/ALLSPICE/'
write_data_fig_real(test_results %>% distinct(), name, save_fig = F)




