source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')
source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/simulations.R')
library("ggsci")

##### NULL Simulation - true value of c ##############
# null_sim11 <- read_csv(paste0(result_path, 'null_sim100_par_combo1_sigma1_test_results.csv'))
# null_sim_p11 <- save_null_sim_figures(null_sim11 %>% filter(pi == 0.5), pi_value = 0.5, name = 'figureM1', save = T)
# null_sim_p11 <- save_null_sim_figures(null_sim11 %>% filter(pi == 0.8), pi_value = 0.8, name = 'figureS3', save = T)
# 
# null_sim11 <- read_csv(paste0(result_path, 'null_sim100_par_combo2_sigma0.1_test_results.csv'))
# null_sim_p11 <- save_null_sim_figures(null_sim11 %>% filter(pi == 0.5), pi_value = 0.5, name = 'figureS1', save = T)
# 
# null_sim11 <- read_csv(paste0(result_path, 'null_sim100_par_combo3_sigma0.01_test_results.csv'))
# null_sim_p11 <- save_null_sim_figures(null_sim11 %>% filter(pi == 0.5), pi_value = 0.5, name = 'figureS2', save = T)



##### NULL Simulation - MLE of c ##############
null_mle_sim11 <- read_csv(paste0(result_path, 'null_mle_sim100_par_combo1_sigma1_test_results.csv'))
null_mle_sim_p11 <- save_null_sim_figures(null_mle_sim11 %>% filter(pi == 0.5), pi_value = 0.5, name = 'ext_figure2', save = T)
null_mle_sim_p11 <- save_null_sim_figures(null_mle_sim11 %>% filter(pi == 0.8), pi_value = 0.8, name = 'figureS8', save = T)

null_mle_sim12 <- read_csv(paste0(result_path, 'null_mle_sim100_par_combo2_sigma0.1_test_results.csv'))
null_mle_sim_p12 <- save_null_sim_figures(null_mle_sim12 %>% filter(pi == 0.5), pi_value = 0.5, name = 'figureS6', save = T)
null_mle_sim_p12 <- save_null_sim_figures(null_mle_sim12 %>% filter(pi == 0.8), pi_value = 0.8, name = 'figureS9', save = T)

null_mle_sim13 <- read_csv(paste0(result_path, 'null_mle_sim100_par_combo3_sigma0.01_test_results.csv'))
null_mle_sim_p13 <- save_null_sim_figures(null_mle_sim13 %>% filter(pi == 0.5), pi_value = 0.5, name = 'figureS7', save = T)
null_mle_sim_p13 <- save_null_sim_figures(null_mle_sim13 %>% filter(pi == 0.8), pi_value = 0.8, name = 'figureS10', save = T)

##### Power Analysis #####
alt_sim <- read_csv(paste0(result_path, 'alt_sim100_par_combo1_sigma1_test_results.csv'))
alt_sim_p <- save_alt_sim_figure(alt_sim , name = 'ext_figure3', save = T)

alt_sim <- read_csv(paste0(result_path, 'alt_sim100_par_combo2_sigma0.1_test_results.csv'))
alt_sim_p <- save_alt_sim_figure(alt_sim , name = 'figureS11', save = T)

alt_sim <- read_csv(paste0(result_path, 'alt_sim100_par_combo3_sigma0.01_test_results.csv'))
alt_sim_p <- save_alt_sim_figure(alt_sim , name = 'figureS12', save = T)


alt_beta <- read_csv(paste0(result_path, 'alt_sim100_par_combo1_sigma1_beta.csv'))
hist(alt_beta$b1_true_value)
alt_beta <- read_csv(paste0(result_path, 'alt_sim100_par_combo2_sigma0.1_beta.csv'))
hist(alt_beta$b1_true_value)
alt_beta <- read_csv(paste0(result_path, 'alt_sim100_par_combo3_sigma0.01_beta.csv'))
hist(alt_beta$b1_true_value)