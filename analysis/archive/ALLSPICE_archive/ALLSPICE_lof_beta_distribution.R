source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/allelic_heterogeneity_test/utils.R')
source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')
setwd('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/')

var_file <- paste0('~/Downloads/corr_testing_burden_var_0.01_500k.txt.bgz')
var_data <- fread(cmd = paste0('gunzip -cq ', var_file)) %>% filter(gene == 'TET2' & annotation =='pLoF')

phenotypes <- unname(unlist(unique(var_data$description)))

for(i in 1:36){
  if(i < 36){
    tmp_phenotypes <- phenotypes[(10*i-9):(10*i)]
  }else{
    tmp_phenotypes <- phenotypes[(10*i-9):length(phenotypes)]
  }
  
  p <- var_data %>%
    filter(description %in% tmp_phenotypes) %>%
    ggplot + aes(x = description, y = BETA) +
    labs(x = NULL, y = 'TET2 pLoF BETA') + 
    geom_boxplot(lwd = 0.2, size = 1, outlier.size = 0.1) + themes + 
    theme(axis.text.x = element_text(angle = 90, size = 4, hjust = 1))
  
  png(paste0('~/Desktop/tet2_pLoF_beta_distribution/AF_1e_4/TET2_pLoF_effect_size_distribution_af_1e-4_', i,'.png'), height = 4, width = 6, units = 'in', res = 300)
  print(p)
  dev.off()
}


for(i in 1:36){
  if(i < 36){
    tmp_phenotypes <- phenotypes[(10*i-9):(10*i)]
  }else{
    tmp_phenotypes <- phenotypes[(10*i-9):length(phenotypes)]
  }
  
  p <- var_data %>%
    filter(AC < 5) %>%
    filter(description %in% tmp_phenotypes) %>%
    ggplot + aes(x = description, y = BETA) +
    labs(x = NULL, y = 'TET2 pLoF BETA') + 
    geom_boxplot(lwd = 0.2, size = 1, outlier.size = 0.1) + themes + 
    theme(axis.text.x = element_text(angle = 90, size = 4, hjust = 1))
  
  png(paste0('~/Desktop/tet2_pLoF_beta_distribution/AC_5/TET2_pLoF_effect_size_distribution_ac_below_5_', i,'.png'), height = 4, width = 6, units = 'in', res = 300)
  print(p)
  dev.off()
}


for(i in 1:36){
  print(paste0('i:', i))
  print((10*i-9))
  print(10*i)
  }

