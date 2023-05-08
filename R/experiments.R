setwd('~/OneDrive/2020Work/pleiotropy_manuscript/')
source('./Rscripts/simulation.R')

##### NULL Simulation - true value of c ##############
n_var <- c(5, 20, 100)
c <- c(0, 0.2, 0.4, 0.6, 0.8, 1, 5/4, 5/3, 5/2, 5, 100)
r <- c(-1, -0.8, -0.5, 0, 0.5, 0.8, 1)
pi <- 0.5
n_ind <- 1000

sigma <- 1
null_sim1 <- pheno_corr_sims(100, n_ind, n_var, c, r, pi, sigma, mle = FALSE, write = TRUE, output = result_path, name = 'null_sim100_par_combo1_sigma1')
save_sim_fig(data = null_sim1, name =paste0(figure_path, 'null_sim100_par_combo1_sigma1_full'), mle = FALSE)
null_sim1 <- list(results = read_csv(paste0(result_path, 'simulation/null_sim100_par_combo1_sigma1', '_test_results.csv')),
                  beta = read_csv(paste0(result_path, 'simulation/null_sim100_par_combo1_sigma1', '_beta.csv')))

a = null_sim1$results %>% 
  filter(r %in% c(0, 0.5, 0.8, 1) & c %in% c(0, 0.2, 0.4, 0.6, 0.8, 1)) %>%
  group_by(r,c,n_var) %>%
  summarise(typeIerror = sum(pvalue<0.05)/n())

p = a %>% 
  ggplot + 
  aes(y = typeIerror, x = r, color = factor(r)) + 
  geom_point() +
  geom_hline(yintercept = 0.05, lty=2) +
  scale_x_continuous(breaks = c(0, 0.5, 0.8, 1)) +
  labs(x = "Phenotypic correlation", y = "Type I Error", color = 'Phenotypic correlation') +
  scale_color_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
  scale_fill_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
  themes +
  facet_grid(n_var~c) 
png(paste0(figure_path, 'null_sim100_par_combo1_sigma1', "_typeIerror.png"), width=6, height=4, units = 'in', res = 300)
print(p)
dev.off()

sigma <- 0.1
null_sim2 <- pheno_corr_sims(100, n_ind, n_var, c, r, pi, sigma, mle = FALSE, write = TRUE, output = result_path, name = 'null_sim100_par_combo2_sigma0.1')
save_sim_fig(data = null_sim2, name =paste0(figure_path, 'null_sim100_par_combo2_sigma0.1'), mle = FALSE)

sigma <- 0.01
null_sim3 <- pheno_corr_sims(100, n_ind, n_var, c, r, pi, sigma, mle = FALSE, write = TRUE, output = result_path, name = 'null_sim100_par_combo3_sigma0.01')
save_sim_fig(data = null_sim3, name =paste0(figure_path, 'null_sim100_par_combo3_sigma0.01'), mle = FALSE)
null_sim3 <- list(results = read_csv(paste0(result_path, 'simulation/null_sim100_par_combo3_sigma0.01', '_test_results.csv')),
                  beta = read_csv(paste0(result_path, 'simulation/null_sim100_par_combo3_sigma0.01', '_beta.csv')))

##### NULL Simulation - MLE of c ##############
sigma <- 1
sim1 <- pheno_corr_sims(100, n_ind, n_var, c, r, pi, sigma, mle = TRUE, write = TRUE, output = result_path, name = 'sim100_par_combo1_sigma1')
save_sim_fig(data = sim1, name =paste0(figure_path, 'sim100_par_combo1_sigma1'), mle = TRUE)

sigma <- 0.1
sim2 <- pheno_corr_sims(100, n_ind, n_var, c, r, pi, sigma, mle = TRUE, write = TRUE, output = result_path, name = 'sim100_par_combo2_sigma0.1')
save_sim_fig(data = sim2, name =paste0(figure_path, 'sim100_par_combo2_sigma0.1'), mle = TRUE)

sigma <- 0.01
sim3 <- pheno_corr_sims(100, n_ind, n_var, c, r, pi, sigma, mle = TRUE, write = TRUE, output = result_path, name = 'sim100_par_combo3_sigma0.01')
save_sim_fig(data = sim3, name =paste0(figure_path, 'sim100_par_combo3_sigma0.01'), mle = TRUE)

##### Power Analysis #####


##### REAL DATA ##############
pheno_corr <- read_delim('subset_data_pheno_correlation_bgz.txt', delim = '\t',
                         col_types = cols(i_pheno = col_character(), j_pheno = col_character())) %>% select(3:5)
# icd <- read_delim('pleiotropy_corr_testing_icd.txt.bgz', delim='\t', col_types = cols(phenocode = col_character()))

##### ATM 
atm <- read_delim('pleiotropy_corr_testing_atm.txt.bgz', delim='\t', col_types = cols(phenocode = col_character()))
# (1) 6 associated phenotypes
atm_pheno6 <- c('30040', '30050', '30090', '30120', '30260', '30270')
atm6 <- get_real_data(atm, atm_pheno6)
atm_test6 <- get_test_result(data = atm6$sub, pheno_corr = atm6$corr, pheno_list = atm_pheno6, n_ind =atm6$n_ind, gene = 'ATM')
write_data_fig_real(atm_test6$results, atm_test6$beta, 'atm_6_biomarkers', save_file = F, save_fig = T)
print_lambda_summary(atm_test6, atm6$sub)

# (2) all 63 biomarker phenotypes
atm_pheno <- unique(atm$phenocode)
atm_pheno63 <- atm_pheno[which(nchar(atm_pheno)==5)]
atm63 <- get_real_data(atm, atm_pheno63)
atm_test63 <- get_test_result(atm63$sub, atm63$corr, atm_pheno63, n_ind = atm63$n_ind)
write_data_fig(atm_test63, 'atm_var_test_63_biomarkers')
print_lambda_summary(atm_test63, atm63$sub)

# (3) 61 biomarker phenotypes removing two extreme ones
atm_pheno61 <- atm_pheno63[which(!(atm_pheno63 %in% c('30820', '30800')))]
atm61 <- get_real_data(atm, atm_pheno61)
atm_test61 <- get_test_result(atm61$sub, atm61$corr, atm_pheno61,  n_ind =atm61$n_ind )
write_data_fig(atm_test61, 'atm_var_test_61_biomarkers')
print_lambda_summary(atm_test61, atm61$sub)

# (4) non-associated biomarkers
atm_pheno55 <- atm_pheno61[which(!(atm_pheno61 %in% atm_pheno6))]
atm55 <- get_real_data(atm, atm_pheno55)
atm_test55 <- get_test_result(atm55$sub, atm55$corr, atm_pheno55, n_ind = atm55$n_ind)
write_data_fig(atm_test55, 'atm_var_test_55_biomarkers')
print_lambda_summary(atm_test55, atm55$sub)

# (5) random phenotypes
atm_rp <- read_delim('pleiotropy_corr_testing_atm_rp.txt.bgz', delim='\t', col_types = cols(phenocode = col_character()))
atm_rp <- atm_rp %>%
  filter((annotation %in% c('missense', 'synonymous') & AC>18) | annotation == 'pLoF') %>%
  filter(AF<0.0001 & phenocode == 'random_continuous') %>%
  mutate(modifier = if_else(is.na(modifier), 1, modifier),
         coding = paste0(phenocode, '_', coding))
pheno_list <- unique(atm_rp$coding)
rp_test <- data.frame()
rp_beta <- data.frame()
for(i in c(0.1, 0.2, 0.5, 1)){
  rp_data <- atm_rp %>%
    filter(modifier == i) %>%
    get_real_data_rp(.)
  temp <- get_test_result_rp(rp_data$sub, pheno_corr=0, pheno_list  = pheno_list, gene = 'ATM', n_ind = rp_data$n_ind)
  tmp_test <- temp$results %>% mutate(h = i)
  tmp_beta <- temp$beta %>% mutate(h = i)
  rp_test <- rbind(rp_test, tmp_test)
  rp_beta <- rbind(rp_beta, tmp_beta)
}

write_data_fig_rp(rp_test, rp_beta, name = 'atm_rp')


##### SLC39A8
slc39a8 <- read_delim('pleiotropy_corr_testing_slc39a8.txt.bgz', delim='\t', col_types = cols(phenocode = col_character()))
# (1)    associated phenotypes
slc39a8_pheno30 <- c('25883', '25882','25890', '25891', '25916', '25917', '25915', '25912', '25914', '25359', '25032', '25893',
                   '25913','25358', '25406', '25910', '25919', '25894','25407','25895','25881','25033','25409','25896',
                   '25892', '25899', '25911', '25166', '25880', '25909')
slc39a830 <- get_real_data(slc39a8, slc39a8_pheno30)
slc39a8_test30 <- get_test_result(slc39a830$sub, slc39a830$corr, slc39a8_pheno30, n_ind = slc39a830$n_ind)
write_data_fig(slc39a8_test30, 'slc39a8_var_test_30_mri')
print_lambda_summary(slc39a8_test30, slc39a830$sub)
# (2) 100 random MRI phenotypes
slc39a8_pheno100 <- unique(slc39a8$phenocode)[sample(1:872, 100, replace=F)]
slc39a8100 <- get_real_data(slc39a8, slc39a8_pheno100)
slc39a8_test100 <- get_test_result(slc39a8100$sub, slc39a8100$corr, slc39a8_pheno100, n_ind = slc39a8100$n_ind)
write_data_fig(slc39a8_test100, 'slc39a8_var_test_100_mri')
print_lambda_summary(slc39a8_test100, slc39a8100$sub)
# (3) 100 random non-associated phenotypes - with phenotypes driving extreme lambdas removed 
slc39a8_pheno <- unique(slc39a8$phenocode)
slc39a8_pheno_ex <- c(c('25455', '25448', '25358', '25029', '25398', '25352'), slc39a8_pheno30)
slc39a8_pheno100_v2 <- slc39a8_pheno[which(!(slc39a8_pheno %in% slc39a8_pheno_ex))][sample(1:842, 100, replace=F)]
slc39a8100_v2 <- get_real_data(slc39a8, slc39a8_pheno100_v2)
slc39a8_test100_v2 <- get_test_result(slc39a8100_v2$sub, slc39a8100_v2$corr, slc39a8_pheno_ex, n_ind = slc39a8100_v2$n_ind)
write_data_fig(slc39a8_test100_v2, 'slc39a8_var_test_100_mri_v2')
print_lambda_summary(slc39a8_test100_v2, slc39a8100_v2$sub)





########### SIMULATION I #############
n_var <- c(5, 20, 100)
c <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
r <- c(0, 0.5, 0.8)
pi <- 0.5
sigma <- 1
n_ind <- 1000
## AF FIXED FOR EACH SCENARIO
simulation1 = pheno_corr_sims(100, n_ind, n_var, c, r, pi, sigma)
result =simulation1$results
b1 = simulation$b1
b2 = simulation$b2
print_lambda_summary2(result)
write_csv(result, 'simulation_var_test_100_A_fixed_set1_sigma1.csv')
write_csv(b1, 'simulation_beta1_100_A_fixed_set1_sigma1.csv')
write_csv(b2, 'simulation_beta2_100_A_fixed_set1_sigma1.csv')

summary(result[result$n_var==5, 'lambda'])
summary(result[result$n_var==20, 'lambda'])
summary(result[result$n_var==100, 'lambda'])

summary(result[result$n_var==5&result$r==0, 'lambda'])
summary(result[result$n_var==20&result$r==0, 'lambda'])
summary(result[result$n_var==100&result$r==0, 'lambda'])
summary(result[result$n_var==5&result$r==0.5, 'lambda'])
summary(result[result$n_var==20&result$r==0.5, 'lambda'])
summary(result[result$n_var==100&result$r==0.5, 'lambda'])
summary(result[result$n_var==5&result$r==0.8, 'lambda'])
summary(result[result$n_var==20&result$r==0.8, 'lambda'])
summary(result[result$n_var==100&result$r==0.8, 'lambda'])

## AF DIFFERS AT EACH SIMULATION
# simulation_rand = pheno_corr_sims_A_rand(100, n_ind, n_var, c, r, pi, sigma)
# result =simulation_rand$results
# b1 = simulation_rand$b1
# b2 = simulation_rand$b2
# print_lambda_summary2(result)
# write_csv(result, 'simulation_var_test_100_A_rand.csv')
# write_csv(b1, 'simulation_beta1_100_A_rand.csv')
# write_csv(b2, 'simulation_beta2_100_A_rand.csv')

fig_c_hat <- result %>%
  ggplot() + 
  geom_density(aes(x = c_hat, color = factor(r)), alpha = 0.5) + themes +
  geom_vline(aes(xintercept=c), lty=2)+
  labs(x = expression(bold(hat(c))), color = 'Phenotypic correlation', y=NULL) +
  xlim(c(-2,2)) +
  scale_color_brewer(palette = 'Dark2') +
  facet_grid(n_var~c)


png('simulation100_c_hat_A_fixed_set1_sigma1.png', height = 6, width = 9, units = 'in', res = 300)
print(fig_c_hat)
dev.off()

fig_lambda <- result %>%
  ggplot() + aes(x = lambda, color = factor(r), fill = factor(r))+
  geom_density(alpha = 0.5) + themes +
  labs(x = expression(bold(Lambda)), y=NULL) +
  scale_color_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
  scale_fill_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
  facet_grid(n_var~c)

png('simulation100_lambda_A_fixed_set1_sigma1.png', height = 6, width = 9, units = 'in', res = 300)
print(fig_lambda)
dev.off()

fig_pvalue <- result %>%
  ggplot() + aes(x = pvalue, color = factor(r), fill = factor(r))+
  geom_histogram(alpha = 0.5, position = 'identity') + themes +
  geom_vline(xintercept=0.05, lty=2) +
  labs(x = 'P-value', y=NULL) +
  scale_color_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
  scale_fill_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
  facet_grid(n_var~c)

png('simulation100_pvalue_A_fixed_set1_sigma1.png', height = 6, width = 9, units = 'in', res = 300)
print(fig_pvalue)
dev.off()

########### SIMULATION 2 #############
n_var <- c(5, 20, 100)
c <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
r <- c(0, 0.5, 0.8)
pi <- 0.5
sigma <- 0.1
n_ind <- 1000
## AF FIXED FOR EACH SCENARIO
simulation2 = pheno_corr_sims(100, n_ind, n_var, c, r, pi, sigma)
result =simulation2$results
b1 = simulation2$b1
b2 = simulation2$b2
print_lambda_summary2(result)
write_csv(result, 'simulation_var_test_100_A_fixed_set2.csv')
write_csv(b1, 'simulation_beta1_100_A_fixed_set2.csv')
write_csv(b2, 'simulation_beta2_100_A_fixed_set2.csv')

summary(result[result$n_var==5, 'c_hat'])
summary(result[result$n_var==20, 'c_hat'])
summary(result[result$n_var==100, 'c_hat'])

summary(result[result$n_var==5&result$r==0, 'lambda'])
summary(result[result$n_var==20&result$r==0, 'lambda'])
summary(result[result$n_var==100&result$r==0, 'lambda'])
summary(result[result$n_var==5&result$r==0.5, 'lambda'])
summary(result[result$n_var==20&result$r==0.5, 'lambda'])
summary(result[result$n_var==100&result$r==0.5, 'lambda'])
summary(result[result$n_var==5&result$r==0.8, 'lambda'])
summary(result[result$n_var==20&result$r==0.8, 'lambda'])
summary(result[result$n_var==100&result$r==0.8, 'lambda'])

fig_c_hat <- result %>%
  ggplot() + 
  geom_density(aes(x = c_hat, color = factor(r)), alpha = 0.5) + themes +
  geom_vline(aes(xintercept=c), lty=2)+
  labs(x = expression(bold(hat(c))), color = 'Phenotypic correlation', y=NULL) +
  xlim(c(-2,2)) +
  scale_color_brewer(palette = 'Dark2') +
  facet_grid(n_var~c)


png('simulation100_c_hat_A_fixed_set2_sigma0.1.png', height = 6, width = 9, units = 'in', res = 300)
print(fig_c_hat)
dev.off()

fig_lambda <- result %>%
  ggplot() + aes(x = lambda, color = factor(r), fill = factor(r))+
  geom_density(alpha = 0.5) + themes +
  labs(x = expression(bold(Lambda)), y=NULL) +
  scale_color_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
  scale_fill_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
  facet_grid(n_var~c)

png('simulation100_lambda_A_fixed_set2_sigma0.1.png', height = 6, width = 9, units = 'in', res = 300)
print(fig_lambda)
dev.off()

fig_pvalue <- result %>%
  ggplot() + aes(x = pvalue, color = factor(r), fill = factor(r))+
  geom_histogram(alpha = 0.5, position = 'identity') + themes +
  geom_vline(xintercept=0.05, lty=2) +
  labs(x = 'P-value', y=NULL) +
  scale_color_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
  scale_fill_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
  facet_grid(n_var~c)

png('simulation100_pvalue_A_fixed_set2_sigma0.1.png', height = 6, width = 9, units = 'in', res = 300)
print(fig_pvalue)
dev.off()

########### SIMULATION 3 #############
n_var <- c(5, 20, 100)
c <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
r <- c(0, 0.5, 0.8)
pi <- 0.5
sigma <- 0.01
n_ind <- 1000

simulation = pheno_corr_sims(100, n_ind, n_var, c, r, pi, sigma)
result =simulation$results
b1 = simulation$b1
b2 = simulation$b2
print_lambda_summary2(result)
write_csv(result, 'simulation_var_test_100_A_fixed_set3_sigma0.01.csv')
write_csv(b1, 'simulation_beta1_100_A_fixed_set3_sigma0.01.csv')
write_csv(b2, 'simulation_beta2_100_A_fixed_set3_sigma0.01.csv')

summary(result[result$n_var==5, 'c_hat'])
summary(result[result$n_var==20, 'c_hat'])
summary(result[result$n_var==100, 'c_hat'])

summary(result[result$n_var==5&result$r==0, 'lambda'])
summary(result[result$n_var==20&result$r==0, 'lambda'])
summary(result[result$n_var==100&result$r==0, 'lambda'])
summary(result[result$n_var==5&result$r==0.5, 'lambda'])
summary(result[result$n_var==20&result$r==0.5, 'lambda'])
summary(result[result$n_var==100&result$r==0.5, 'lambda'])
summary(result[result$n_var==5&result$r==0.8, 'lambda'])
summary(result[result$n_var==20&result$r==0.8, 'lambda'])
summary(result[result$n_var==100&result$r==0.8, 'lambda'])

fig_c_hat <- result %>%
  ggplot() + 
  geom_density(aes(x = c_hat, color = factor(r)), alpha = 0.5) + themes +
  geom_vline(aes(xintercept=c), lty=2)+
  labs(x = expression(bold(hat(c))), color = 'Phenotypic correlation', y=NULL) +
  xlim(c(-2,2)) +
  scale_color_brewer(palette = 'Dark2') +
  facet_grid(n_var~c)


png('simulation100_c_hat_A_fixed_set3_sigma0.01.png', height = 6, width = 9, units = 'in', res = 300)
print(fig_c_hat)
dev.off()

fig_lambda <- result %>%
  ggplot() + aes(x = lambda, color = factor(r), fill = factor(r))+
  geom_density(alpha = 0.5) + themes +
  labs(x = expression(bold(Lambda)), y=NULL) +
  scale_color_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
  scale_fill_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
  facet_grid(n_var~c)

png('simulation100_lambda_A_fixed_set3_sigma0.01.png', height = 6, width = 9, units = 'in', res = 300)
print(fig_lambda)
dev.off()

fig_pvalue <- result %>%
  ggplot() + aes(x = pvalue, color = factor(r), fill = factor(r))+
  geom_histogram(alpha = 0.5, position = 'identity') + themes +
  geom_vline(xintercept=0.05, lty=2) +
  labs(x = 'P-value', y=NULL) +
  scale_color_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
  scale_fill_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
  facet_grid(n_var~c)

png('simulation100_pvalue_A_fixed_set3_sigma0.01.png', height = 6, width = 9, units = 'in', res = 300)
print(fig_pvalue)
dev.off()



#### Test derive c_hat
n_var <- 100
n_ind <- 10000
pi <- 0.5
sigma <- 1
c <- 0.5
r <- seq(-100, 100, 0.01)
AC <- get_ac_mat(n_var)
A <- get_af_mat(AC, n_ind)
X <- get_geno_mat(AC, n_ind)
b <- get_true_beta(n_ind, n_var, c, pi, sigma)
Y <- get_pheno_pair(b, X, r)
b_hat <- get_beta_hat(Y, X, A, n_ind)
b1_hat <- matrix(b_hat[1, ], nrow = 1)
b2_hat <- matrix(b_hat[2, ], nrow = 1)
c_hat <- get_c_hat(b1_hat, b2_hat,A , r)
# lambda <- get_likelihood_test_stats(n_ind, b1_hat, b2_hat, A, r)
eq <- function(r){
  Y <- get_pheno_pair(b, X, r)
  b_hat <- get_beta_hat(Y, X, A, n_ind)
  b1_hat <- matrix(b_hat[1, ], nrow = 1)
  b2_hat <- matrix(b_hat[2, ], nrow = 1)
  c_hat <- get_c_hat(b1_hat, b2_hat,A , r)
  return((n_ind/(c_hat^2-2*c_hat*r+1))*c((b1_hat - c_hat*b2_hat) %*% A %*% t(b1_hat - c_hat*b2_hat)))
  }
x = seq(-10,10, 0.01)
y = sapply(x,eq)
x[which(y==min(y))]
plot(x, y, xlab = 'r', ylab = expression(Lambda))

u <- b2_hat %*% A %*% t(b1_hat - r*b2_hat)
v <- (b2_hat-b1_hat) %*% A %*% t(b2_hat+b1_hat)
w <- b1_hat %*% A %*% t(r*b1_hat - b2_hat)
(-v+ sqrt(v^2-4*u*w))/(2*u)
(-v- sqrt(v^2-4*u*w))/(2*u)
u


eq <- function(x){(n_ind/(x^2-2*x*r+1))*c((b1_hat - x*b2_hat) %*% A %*% t(b1_hat - x*b2_hat))}
x = seq(-100,100, 0.01)
y = sapply(x,eq)
x[which(y==min(y))]
plot(x, y, xlab = 'c', ylab = expression(Lambda))

data <- slc39a8200$sub
corr <- slc39a8200$corr
A <- diag(data$AF)
b1_hat <- matrix(data[,100],nrow=1)
b2_hat <- matrix(data[,152], nrow=1)
pheno1 <- colnames(data)[100]
pheno2 <- colnames(data)[152]
r <- unlist(corr[(corr$i_pheno==pheno1 &corr$j_pheno==pheno2),'corr'])
