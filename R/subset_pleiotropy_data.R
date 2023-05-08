source('~/Dropbox (Partners HealthCare)/ukbb_exomes_pleiotropy/R/simulations.R')
source('~/Dropbox (Partners HealthCare)/ukbb_exomes_pleiotropy/R/constants.R')
setwd('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/')

read_pleiotropy_subset <- function(test, tranche, type){
  data <- read_delim(paste0('~/Downloads/corr_testing_', test, '_', type,'_', tranche, '.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character())) %>%
  mutate(annotation = factor(annotation, levels = annotation_types, labels = annotation_names))
  if(test == 'burden') data$Pvalue = data$Pvalue_Burden
  return(data)
}

genename <- 'CD36'
annt <- 'pLoF'
phenocode1 <- '30070'
phenocode2 <- '30110'

# var_300k_sub <- read_pleiotropy_subset('burden', '300k', 'var') %>%
#   filter(gene == genename & annotation == annotation & pheno1 == phenocode1 & pheno2 == phenocode2)
#
# var_500k_sub <- read_pleiotropy_subset('burden', '500k', 'var') %>%
#   filter(gene == genename & annotation == annotation & pheno1 == phenocode1 & pheno2 == phenocode2)
#
# write_csv(var_300k_sub, paste0(result_path, 'var_subset_', genename, '_', phenocode1, '_', phenocode2, '_300k.csv'))
# write_csv(var_500k_sub, paste0(result_path, 'var_subset_', genename, '_', phenocode1, '_', phenocode2, '_500k.csv'))


gene_300k_sub <- read_pleiotropy_subset('burden', '300k', 'gene') %>%
  filter(gene_symbol == genename & annotation == annt & phenocode %in% c(phenocode1, phenocode2))

gene_500k_sub <- read_pleiotropy_subset('burden', '500k', 'gene') %>%
  filter(gene_symbol == genename & annotation == annt & phenocode %in% c(phenocode1, phenocode2))

write_csv(gene_300k_sub, paste0(result_path, 'gene_subset_', tolower(genename), '_', tolower(annt), '_', phenocode1, '_', phenocode2, '_300k.csv'))
write_csv(gene_500k_sub, paste0(result_path, 'gene_subset_', tolower(genename), '_', tolower(annt), '_', phenocode1, '_', phenocode2, '_500k.csv'))

gene_300k_sub <- read_csv(paste0(result_path, 'gene_subset_', tolower(genename), '_', tolower(annt), '_', phenocode1, '_', phenocode2, '_300k.csv'))
gene_500k_sub <- read_csv(paste0(result_path, 'gene_subset_', tolower(genename), '_', tolower(annt), '_', phenocode1, '_', phenocode2, '_500k.csv'))

mean_n_cases_300k <- mean(gene_300k_sub$n_cases)
mean_n_cases_500k <- mean(gene_500k_sub$n_cases_defined)
TRANCHE = '500k'
var_data <- read_delim(paste0(PATH, 'var_subset_cd36_', TRANCHE,'.csv'), delim='\t')
temp <- get_real_data(var_data, phenolist)
write_csv(temp, paste0('~/Downloads/var_subset_cd36_', TRANCHE,'_modified_for_test.csv'))