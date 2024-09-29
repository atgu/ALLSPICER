source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')
source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/simulations.R')

##### NULL Simulation - true value of c ##############
sigma <- 1
null_sim1 <- pheno_corr_sims(100, all_n_ind, all_n_var, all_c, all_r, all_pi, sigma, mle = FALSE, write = TRUE, output = result_path, name = 'null_sim100_par_combo1_sigma1')

sigma <- 0.1
null_sim2 <- pheno_corr_sims(100, all_n_ind, all_n_var, all_c, all_r, all_pi, sigma, mle = FALSE, write = TRUE, output = result_path, name = 'null_sim100_par_combo2_sigma0.1')

sigma <- 0.01
null_sim3 <- pheno_corr_sims(100, all_n_ind, all_n_var, all_c, all_r, all_pi, sigma, mle = FALSE, write = TRUE, output = result_path, name = 'null_sim100_par_combo3_sigma0.01')


##### NULL Simulation - MLE of c ##############
sigma <- 1
null_mle_sim1 <- pheno_corr_sims(100, all_n_ind, all_n_var, all_c, all_r, all_pi, sigma, mle = TRUE, write = TRUE, output = result_path, name = 'null_mle_sim100_par_combo1_sigma1')

sigma <- 0.1
null_mle_sim2 <- pheno_corr_sims(100, all_n_ind, all_n_var, all_c, all_r, all_pi, sigma, mle = TRUE, write = TRUE, output = result_path, name = 'null_mle_sim100_par_combo2_sigma0.1')

sigma <- 0.01
null_mle_sim3 <- pheno_corr_sims(100, all_n_ind, all_n_var, all_c, all_r, all_pi, sigma, mle = TRUE, write = TRUE, output = result_path, name = 'null_mle_sim100_par_combo3_sigma0.01')


##### Power Analysis #####
sigma <- 1
alt_sim1 <- pheno_corr_sims(100, all_n_ind, all_n_var, c=NULL, all_r, all_pi, sigma, mle = FALSE, null=FALSE, write = TRUE, output = result_path, name = 'alt_sim100_par_combo1_sigma1')

sigma <- 0.1
alt_sim2 <- pheno_corr_sims(100, all_n_ind, all_n_var, c=NULL, all_r, all_pi, sigma, mle = FALSE, null=FALSE, write = TRUE, output = result_path, name = 'alt_sim100_par_combo2_sigma0.1')

sigma <- 0.01
alt_sim3 <- pheno_corr_sims(100, all_n_ind, all_n_var, c=NULL, all_r, all_pi, sigma, mle = FALSE, null=FALSE, write = TRUE, output = result_path, name = 'alt_sim100_par_combo3_sigma0.01')