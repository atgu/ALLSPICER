source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/simulations.R')
source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/constants.R')
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
null_mle_sim_p11 <- save_null_sim_figures(null_mle_sim11 %>% dplyr::filter(pi == 0.5), pi_value = 0.5, name = 'ext_figure2', save = T)
null_mle_sim_p11 <- save_null_sim_figures(null_mle_sim11 %>% dplyr::filter(pi == 0.8), pi_value = 0.8, name = 'figureS9', save = T)

null_mle_sim12 <- read_csv(paste0(result_path, 'null_mle_sim100_par_combo2_sigma0.1_test_results.csv'))
null_mle_sim_p12 <- save_null_sim_figures(null_mle_sim12 %>% dplyr::filter(pi == 0.5), pi_value = 0.5, name = 'figureS7', save = T)
null_mle_sim_p12 <- save_null_sim_figures(null_mle_sim12 %>% dplyr::filter(pi == 0.8), pi_value = 0.8, name = 'figureS10', save = T)

null_mle_sim13 <- read_csv(paste0(result_path, 'null_mle_sim100_par_combo3_sigma0.01_test_results.csv'))
null_mle_sim_p13 <- save_null_sim_figures(null_mle_sim13 %>% dplyr::filter(pi == 0.5), pi_value = 0.5, name = 'figureS8', save = T)
null_mle_sim_p13 <- save_null_sim_figures(null_mle_sim13 %>% dplyr::filter(pi == 0.8), pi_value = 0.8, name = 'figureS11', save = T)

##### correlated-beta Analysis
alt_sim <- read_csv(paste0(result_path, 'alt_sim100_par_combo1_sigma1_test_results.csv'))
correlated_sim <- read_csv(paste0(result_path, 'correlated_sim100_par_combo1_sigma1_test_results.csv'))
nonlinear_sim <- read_csv(paste0(result_path, 'nonlinear_sim100_par_combo1_sigma1_test_results.csv'))
figure = ggpubr::ggarrange(
  save_alt_sim_figure(alt_sim, name = 'ext_figure3A', save = F) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(plot.margin = unit(c(0.5, 0, 0.5, 0), 'cm')),
  save_alt_sim_figure(correlated_sim, name = 'ext_figure3B', save = F)+
    theme(plot.margin = unit(c(0.5, 0, 0.5, 0), 'cm')),
  save_alt_sim_figure(nonlinear_sim, name = 'ext_figure3C', save = F)+
    theme(plot.margin = unit(c(0.5, 0, 0.5, 0), 'cm')),
  ncol=1, font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'), common.legend = T,
  labels = c('(A) Independent effect sizes',
             '(B) Correlated effect sizes',
             '(C) Effect sizes in non-linear relationship'), vjust = 0, hjust = 0

)
png(paste0(figure_path,'ext_figure3.png'), height = 8, width = 7.5, units = 'in', res = 300)
print(figure)
dev.off()

alt_sim <- read_csv(paste0(result_path, 'alt_sim100_par_combo2_sigma0.1_test_results.csv'))
correlated_sim<- read_csv(paste0(result_path, 'correlated_sim100_par_combo2_sigma0.1_test_results.csv'))
nonlinear_sim <- read_csv(paste0(result_path, 'nonlinear_sim100_par_combo2_sigma0.1_test_results.csv'))
figure = ggpubr::ggarrange(
  save_alt_sim_figure(alt_sim, name = 'figureS11A', save = F) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(plot.margin = unit(c(0.5, 0, 0.5, 0), 'cm')),
  save_alt_sim_figure(correlated_sim, name = 'figureS11B', save = F)+
    theme(plot.margin = unit(c(0.5, 0, 0.5, 0), 'cm')),
  save_alt_sim_figure(nonlinear_sim, name = 'figureS11C', save = F)+
    theme(plot.margin = unit(c(0.5, 0, 0.5, 0), 'cm')),
  ncol=1, font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'), common.legend = T,
  labels = c('(A) Independent effect sizes',
             '(B) Correlated effect sizes',
             '(C) Effect sizes in non-linear relationship'), vjust = 0, hjust = 0

)
png(paste0(figure_path,'figureS12.png'), height = 8, width = 7.5, units = 'in', res = 300)
print(figure)
dev.off()

alt_sim <- read_csv(paste0(result_path, 'alt_sim100_par_combo3_sigma0.01_test_results.csv'))
correlated_sim <- read_csv(paste0(result_path, 'correlated_sim100_par_combo3_sigma0.01_test_results.csv'))
nonlinear_sim <- read_csv(paste0(result_path, 'nonlinear_sim100_par_combo3_sigma0.01_test_results.csv'))
figure = ggpubr::ggarrange(
  save_alt_sim_figure(alt_sim, name = 'figureS12A', save = F) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(plot.margin = unit(c(0.5, 0, 0.5, 0), 'cm')),
  save_alt_sim_figure(correlated_sim, name = 'figureS12B', save = F)+
    theme(plot.margin = unit(c(0.5, 0, 0.5, 0), 'cm')),
  save_alt_sim_figure(nonlinear_sim, name = 'figureS12C', save = F)+
    theme(plot.margin = unit(c(0.5, 0, 0.5, 0), 'cm')),
  ncol=1, font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'), common.legend = T,
  labels = c('(A) Independent effect sizes',
             '(B) Correlated effect sizes',
             '(C) Effect sizes in non-linear relationship'), vjust = 0, hjust = 0

)
png(paste0(figure_path,'figureS13.png'), height = 8, width = 7.5, units = 'in', res = 300)
print(figure)
dev.off()
