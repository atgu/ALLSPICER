### Simulation
simulation_results <- ALLSPICE_simulation(n_ind=10000, n_var=100, c=0.6, r=0.5, pi=0.5, sigma=1, mle = TRUE, null=TRUE)
print(simulation_results[[1]])


simulation_results <- ALLSPICE_simulation(n_ind=10000, n_var=100, c=0.6, r=0.5, pi=0.5, sigma=1, mle = TRUE, null=FALSE)
print(simulation_results[[1]])


### real data
# data <- read_delim('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/data/alb_albumin_calcium_var_hgvsp_500k.txt.bgz')
# data <- data %>%
#   pivot_wider(names_from = phenocode,  values_from = c('BETA', 'Pvalue'), id_cols = c('locus', 'alleles','pop_AF', 'gene', 'annotation'))
# write_csv(data, '~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/ALLSPICE/tests/testthat/example_data_ALB_Albumin_Calcium.csv')

data <- read_csv('example_data_ALB_Albumin_Calcium.csv')
data1 <- data %>% dplyr::filter(annotation == 'pLoF')
print(head(data1))
result1 <- ALLSPICE(data= data1, pheno_corr = 0.56, n_ind = 400000, gene = 'ALB-pLoF', pheno1 = 'Albumin', pheno2 = 'Calcium', beta1_field = 'BETA_30600', beta2_field = 'BETA_30680', af_field = 'pop_AF')
print(result1)
data2 <- data %>% dplyr::filter(annotation%in% c('missense', 'LC'))
result2 <- ALLSPICE(data= data2, pheno_corr = 0.56, n_ind = 400000, gene = 'ALB-missense', pheno1 = 'Albumin', pheno2 = 'Calcium', beta1_field = 'BETA_30600', beta2_field = 'BETA_30680', af_field = 'pop_AF')
print(result2)
