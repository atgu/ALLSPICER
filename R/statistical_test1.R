setwd('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/')
source('~/Dropbox (Partners HealthCare)/ukbb_exomes_pleiotropy/R/constants.R')
source('~/Dropbox (Partners HealthCare)/ukbb_exomes_pleiotropy/R/simulations.R')

pheno_corr <- read_delim(paste0(data_path, 'correlation_table_phenos_500k.txt.bgz'), delim = '\t',
                         col_types = cols(i_pheno = col_character(), j_pheno = col_character())) %>%
  select(3:5) %>%
  mutate(corr=entry)
TEST = 'skato'
######### real data
get_real_data <- function(var_data, gene_data, phenolist){
  corr <- pheno_corr %>% filter((i_pheno %in% phenolist) & (j_pheno %in% phenolist) & (i_pheno != j_pheno))
  n_ind <- gene_data %>% filter(phenocode %in% phenolist) %>% select(phenocode, n_cases) %>% distinct()
  sub <- var_data %>% 
    filter(phenocode %in% phenolist) %>% 
    pivot_wider(id_col = c('locus', 'alleles', 'annotation'), names_from = 'phenocode', values_from = 'BETA')
  sub <- var_data %>% 
    select(locus, alleles, annotation, AF, AC, gene) %>% 
    unique() %>% merge(., sub, by = c('locus', 'alleles', 'annotation')) %>% 
    # filter(complete.cases(.)) %>%
    mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation))
  return(list(corr=corr, sub=sub, n_ind=n_ind))
}

single_test <- function(data, pheno_corr, pheno1, pheno2, n_ind, sig_level, gene){
  n_ind <- n_ind %>% filter(phenocode %in% c(pheno1, pheno2)) %>% summarise(mean=mean(n_cases))
  n_ind <- as.numeric(n_ind$mean)
  sub <- data %>% select(1:5, pheno1, pheno2) %>% filter(complete.cases(.))
  sub <- as.data.frame(sub)
  if(nrow(sub)>1){
    A <- 2*diag(sub[af_type])
  }else{
    A <- as.matrix(2*sub[af_type])
  }
  b1_hat <- t(as.matrix(sub[,pheno1]))
  b2_hat <- t(as.matrix(sub[,pheno2]))
  r <- pheno_corr
  c_hat <- get_c_hat(b1_hat, b2_hat, A, r)
  lambda <- get_likelihood_test_stats(n_ind, r, b1_hat, b2_hat, c_hat, A)
  pvalue <- 1 - pchisq(lambda, length(b1_hat)-1)
  results <- data.frame(pheno1, pheno2, c_hat, lambda, pvalue, gene, length(b1_hat))
  colnames(results) <- c('pheno1', 'pheno2', 'c_hat', 'lambda', 'pvalue', 'gene','n_var')
  beta <- data.frame(b1 = t(b1_hat), b2 = t(b2_hat), AF = c(diag(A))) %>% 
    mutate(pheno1 = pheno1, pheno2  = pheno2)
  return(list(results = results, beta = beta))
}

var_test <- function(data, pheno_corr, pheno_list, n_ind, gene, sig_level=0.05){
  results <- data.frame()
  beta <- data.frame()
  pheno_list <- colnames(data)[-(1:6)]
  n <- n_ind
  
  for(i in (1: (length(pheno_list)-1))+6){
    for(j in (i: (length(pheno_list)+6))){
      if(j == i) next
      if(length(pheno_list)<2) break
      if(nrow(data) == 0) break
      pheno1 <- colnames(data)[i]
      pheno2 <- colnames(data)[j]
      # n_ind2 <- n_ind %>% filter(phenocode %in% c(pheno1, pheno2)) %>% summarise(mean=mean(n_cases))
      n_ind <- n %>% filter(phenocode %in% c(pheno1, pheno2)) %>% summarise(mean=mean(n_cases))
      n_ind <- floor(n_ind$mean)
      
      sub <- data %>% select(1:6, pheno1, pheno2) %>% filter(complete.cases(.))
      sub <- as.data.frame(sub)
      if(nrow(sub) == 0){next}
      if(nrow(sub)>1){
        A <- 2*diag(sub$AF)
      }else{
        A <- as.matrix(2*sub$AF)
      }
      b1_hat <- t(as.matrix(sub[,pheno1]))
      b2_hat <- t(as.matrix(sub[,pheno2]))
      r <- ifelse(is.numeric(pheno_corr), pheno_corr, c(unlist(pheno_corr[((pheno_corr$i_pheno==pheno1) &(pheno_corr$j_pheno==pheno2)),'corr'])))
      c_hat <- get_c_hat(b1_hat, b2_hat, A, r)
      lambda <- get_likelihood_test_stats(n_ind, r, b1_hat, b2_hat, c_hat, A)
      pvalue <- 1 - pchisq(lambda, length(b1_hat)-1)
      temp <- data.frame(pheno1, pheno2, c_hat, lambda, pvalue, gene, length(b1_hat))
      beta_temp <- data.frame(b1 = unname(t(b1_hat)), b2 = unname(t(b2_hat)), AF = c(diag(A))) %>% 
        mutate(pheno1 = pheno1, pheno2  = pheno2, gene = gene)
      results <- rbind(results, temp)
      beta <- rbind(beta, beta_temp)
    }
  }
  if(nrow(results)>0){
    colnames(results) <- c('pheno1', 'pheno2', 'c_hat', 'lambda', 'pvalue', 'gene', 'n_var')
    return(list(results = results, beta = beta))
  }
}

get_test_result <- function(data, r, pheno_list, n_ind, gene){
  full_results <- data.frame()
  full_beta <- data.frame()
  for(i in names(table(data$annotation))){
    test <- var_test(data %>% filter(annotation == i), r, pheno_list, n_ind, gene) 
    if(is.null(test)) next
    results <- test$results %>% mutate(annotation = i)
    beta <- test$beta %>% mutate(annotation = i)
      full_results <- rbind(full_results, results)
      full_beta <- rbind(full_beta, beta)
    }
  if(nrow(full_results)>0){
    full_results <- full_results %>% 
      merge(., pheno_corr, by.x = c("pheno1", "pheno2"), by.y = c("i_pheno", "j_pheno")) %>% 
      mutate(annotation = factor(annotation, levels = annotation_types, labels = annotation_names))
    return(list(results = full_results, beta = full_beta))
    }
  }

## ALB figure
save_beta_figure <- function(beta, gene, annotations, pheno1, pheno2, name, save=TRUE){
  sub <- beta %>%
    mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
    filter(
      (gene == gene) &
        (annotation %in% annotations) &
        (pheno1 == pheno1) &
        (pheno2 == pheno2)) %>%
    mutate(
      annotation = factor(annotation, levels = annotation_types), 
      BETA_adjusted = sqrt(2*AF*(1-AF))*BETA
    )%>% 
    select(locus, alleles, phenocode, AF, AC, gene, annotation, BETA, BETA_adjusted, Pvalue)
  
  pvalue_1 <- sub %>% filter(phenocode == pheno1) %>% select(1:7, 10)
  pvalue_2 <- sub %>% filter(phenocode == pheno2) %>% select(1:7, 10)
  sub <- sub %>% 
    pivot_wider(names_from = phenocode, names_prefix = 'beta_adjusted_', values_from = BETA_adjusted, id_col = c('locus', 'alleles','AF', 'AC', 'gene', 'annotation')) %>%
    merge(., pvalue_1, by = c('locus', 'alleles','AF', 'AC', 'gene', 'annotation')) %>%
    merge(., pvalue_2, by = c('locus', 'alleles','AF', 'AC', 'gene', 'annotation'))
  
  figure <- sub %>% 
    mutate(annotation = factor(annotation, levels = annotation_types)) %>%
    mutate(p_size = (-log(Pvalue.x) + -log(Pvalue.y))/2)%>%
    ggplot + 
    aes(x = get(paste0('beta_adjusted_', pheno1)), 
        y = get(paste0('beta_adjusted_', pheno2)), 
        color = annotation, 
        size=p_size, 
        label=locus) +
    geom_point() + 
    geom_abline(slope = 1, intercept = 0, lty =2) + 
    labs(x = 'Albumin', 
         y = 'Calcium', 
         # title = gene, 
         size = expression(bold(-log(Pvalue)[mean]))) + 
    theme_classic() +
    themes + theme(axis.title = element_text(size = 12),
                   axis.text = element_text(size = 8.5),
                   strip.text = element_text(size = 12),
                   legend.position = 'right',
                   legend.title = element_text(size = 10),
                   legend.background = element_blank()
    ) + 
    scale_size_continuous(breaks = c(2,5,10,20)) +
    annotation_color_scale + annotation_fill_scale +
    guides(color = 'none') +
    facet_wrap(~annotation, labeller = label_type, scales ='free') + 
    geom_text_repel(max.overlaps = 3, size=2.5)
  if(save){
    png(paste0(figure_path, gene, '_', pheno1, '_', pheno2, '_', name,'_beta_', tranche,'.png'), height = 3, width = 9, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

TEST = 'burden'
var_file <- paste0('~/Downloads/corr_testing_', TEST, '_var_n_cases_over_300k.txt.bgz')
gene_file <- paste0('~/Downloads/corr_testing_', TEST, '_gene_n_cases_over_300k.txt.bgz')

p_filter = FALSE
gene_name = 'ALB'
pheno_pair = c('30600', '30680')
alb_result <- data.frame()
# var_data <- fread(cmd = paste0('gunzip -cq ', var_file))
var_data <- read_delim(paste0(data_path, tolower(gene_name), '_albumin_calcium_var_hgvsp_', tranche,'.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character()))
# var_data <- read_delim(paste0(data_path, tolower(gene_name), '_bmd_var_hgvsp_', tranche,'.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character()))
if(p_filter){
  var_data <- var_data %>%
    filter(Pvalue < 0.5)
}
gene_data <- fread(cmd = paste0('gunzip -cq ', gene_file)) %>% filter(gene_symbol == gene_name)
test_data <- get_real_data(var_data, gene_data, pheno_pair)


for(r in seq(0, 0.9, 0.1)){
  test_result <- get_test_result(test_data$sub, r, pheno_pair, test_data$n_ind, gene_name)
  alb_result <- rbind(alb_result, test_result$results %>%
                        mutate( r = r, note = 'All variants'))
}
# test_result <- get_test_result(test_data$sub, test_data$corr, pheno_pair, test_data$n_ind, gene_name)
# alb_result <- test_result$results %>%
#   mutate(note = 'All variants')
# save_beta_figure(var_data, gene_name, c('pLoF', 'missense|LC', 'synonymous'), pheno_pair[1], pheno_pair[2], paste0('all_variants', if_else(p_filter, '_p_5e_1', '')))

# AF < 0.0001
test_data <- get_real_data(var_data %>% filter(AF < 0.0001), gene_data, pheno_pair)
for(r in seq(0, 0.9, 0.1)){
  test_result <- get_test_result(test_data$sub, r, pheno_pair, test_data$n_ind, gene_name)
  alb_result <- rbind(alb_result, test_result$results %>%
                        mutate( r = r, note = 'AF < 0.0001'))
}
# test_result <- get_test_result(test_data$sub, test_data$corr, pheno_pair, test_data$n_ind, gene_name)
# alb_result <- rbind(alb_result, test_result$results %>% mutate(note = 'AF < 0.0001'))
# save_beta_figure(var_data %>% filter(AF < 0.0001), gene_name, c('pLoF', 'missense|LC', 'synonymous'), pheno_pair[1], pheno_pair[2], paste0('AF_below_1e_4', if_else(p_filter, '_p_5e_1', '')))


# AF < 0.0001 & AC >5
test_data <- get_real_data(var_data %>% filter(AF < 0.0001 &  AC >5), gene_data, pheno_pair)
for(r in seq(0, 0.9, 0.1)){
  test_result <- get_test_result(test_data$sub, r, pheno_pair, test_data$n_ind, gene_name)
  alb_result <- rbind(alb_result, test_result$results %>%
                        mutate( r = r, note = 'AF < 0.0001 & AC > 5'))
}
# test_result <- get_test_result(test_data$sub, test_data$corr, pheno_pair, test_data$n_ind, gene_name)
# alb_result <- rbind(alb_result, test_result$results %>% mutate(note = 'AF < 0.0001 & AC > 5'))
# save_beta_figure(var_data %>% filter(AF < 0.0001  &  AC >5), gene_name, c('pLoF', 'missense|LC', 'synonymous'), pheno_pair[1], pheno_pair[2], paste0('AF_below_1e_4_AC_over_5', if_else(p_filter, '_p_5e_1', '')))

# AF < 0.0001 & AC >13
test_data <- get_real_data(var_data %>% filter(AF < 0.0001 &  AC >13), gene_data, pheno_pair)
for(r in seq(0, 0.9, 0.1)){
  test_result <- get_test_result(test_data$sub, r, pheno_pair, test_data$n_ind, gene_name)
  alb_result <- rbind(alb_result, test_result$results %>%
                        mutate( r = r, note = 'AF < 0.0001 & AC > 13'))
}
# test_result <- get_test_result(test_data$sub, test_data$corr, pheno_pair, test_data$n_ind, gene_name)
# alb_result <- rbind(alb_result, test_result$results %>% mutate(note = 'AF < 0.0001 & AC > 13'))
# save_beta_figure(var_data %>% filter(AF < 0.0001  &  AC >13), gene_name, c('pLoF', 'missense|LC', 'synonymous'), pheno_pair[1], pheno_pair[2], paste0('AF_below_1e_4_AC_over_13', if_else(p_filter, '_p_5e_1', '')))

# AC < 13
test_data <- get_real_data(var_data %>% filter(AC < 13), gene_data, pheno_pair)
for(r in seq(0, 0.9, 0.1)){
  test_result <- get_test_result(test_data$sub, r, pheno_pair, test_data$n_ind, gene_name)
  alb_result <- rbind(alb_result, test_result$results %>%
                        mutate( r = r, note = 'AC < 13'))
}
# test_result <- get_test_result(test_data$sub, test_data$corr, pheno_pair, test_data$n_ind, gene_name)
# alb_result <- rbind(alb_result, test_result$results %>% mutate(note = 'AC < 13'))
# save_beta_figure(var_data %>% filter(AC < 13), gene_name, c('pLoF', 'missense|LC', 'synonymous'), pheno_pair[1], pheno_pair[2], paste0('AC_below_13', if_else(p_filter, '_p_5e_1', '')))

# AC < 5
test_data <- get_real_data(var_data %>% filter(AC < 5), gene_data, pheno_pair)
for(r in seq(0, 0.9, 0.1)){
  test_result <- get_test_result(test_data$sub, r, pheno_pair, test_data$n_ind, gene_name)
  alb_result <- rbind(alb_result, test_result$results %>%
                        mutate( r = r, note = 'AC < 5'))
}
# test_result <- get_test_result(test_data$sub, test_data$corr, pheno_pair, test_data$n_ind, gene_name)
# alb_result <- rbind(alb_result, test_result$results %>% mutate(note = 'AC < 5'))
# save_beta_figure(var_data %>% filter(AC < 5), gene_name, c('pLoF', 'missense|LC', 'synonymous'), pheno_pair[1], pheno_pair[2], paste0('AC_below_5', if_else(p_filter, '_p_5e_1', '')))

# singletons
test_data <- get_real_data(var_data %>% filter(AC == 1), gene_data, pheno_pair)
for(r in seq(0, 0.9, 0.1)){
  test_result <- get_test_result(test_data$sub, r, pheno_pair, test_data$n_ind, gene_name)
  alb_result <- rbind(alb_result, test_result$results %>%
                        mutate( r = r, note = 'Singletons'))
}
# test_result <- get_test_result(test_data$sub, test_data$corr, pheno_pair, test_data$n_ind, gene_name)
# alb_result <- rbind(alb_result, test_result$results %>% mutate(note = 'singletons'))
# save_beta_figure(var_data %>% filter(AC == 1), gene_name, c('pLoF', 'missense|LC', 'synonymous'), pheno_pair[1], pheno_pair[2], paste0('singletons', if_else(p_filter, '_p_5e_1', '')))
write_csv(alb_results, paste0(result_path, tolower(gene_name), '_stat_test1_all_scenarios_switch_r_', if_else(p_filter, '_p_5e_1', ''),'.csv'))
write_csv(alb_result, paste0(result_path, tolower(gene_name), '_stat_test1_all_scenario', if_else(p_filter, '_p_5e_1', ''),'.csv'))
t <- read_csv(paste0(result_path, 'alb_stat_test1_all_scenario_no_p_filter.csv'))
t <- t %>%
  mutate(r = corr,
         note = if_else(note == 'singletons', 'Singletons', note)) %>%
  select(1:10, r, note)
alb_results <- rbind(t, alb_result)
gene_name = 'ALB'
p_filter = FALSE
figure <- alb_results %>%
  mutate(note = factor(note, levels = c('Singletons', 'AC < 5', 'AC < 13', 'AF < 0.0001', 'AF < 0.0001 & AC > 13', 'AF < 0.0001 & AC > 5', 'All variants'))) %>%
  ggplot +  aes(x =r, y = -log10(pvalue), color = annotation, size = log10(n_var)) +
  geom_point(position = 'dodge', alpha = 0.8) +
  geom_vline(aes(xintercept = corr), lty = 2) + 
  geom_hline(yintercept = -log10(0.05), lty = 2, color = 'red')+
  labs(size = 'Number of variants') +
  scale_size_continuous(breaks = c(1, log10(50), 2, log10(200)), labels = c(10, 50, 100, 200),range = c(2,3)) + 
  scale_x_continuous(breaks = seq(0, 0.9, 0.1), labels =seq(0, 0.9, 0.1) ) +
  annotation_color_scale +
  facet_grid(~note, scale = 'free_y') +
  theme(axis.text.x = element_text(angle=45))

png(paste0(figure_path, tolower(gene_name), '_stat_test1_all_scenarios_switch_r_', if_else(p_filter, '_p_5e_1', ''), tranche,'.png'), height = 4, width = 15, units = 'in', res = 300)
print(figure)
dev.off()

## 300k vs. 450k test results
gene_300k <- read_delim('~/OneDrive/2020Work/Pleiotropy/pleiotropy_subdata_corr_testing_all_continuous_pheno_var_p5e_1_af1e_4_ac9_gene.txt.bgz', delim='\t', col_types = cols(phenocode = col_character()))
result_300k <- read_csv('~/OneDrive/2020Work/Pleiotropy/all_af1e_4_ac_9_p_5e_1_full_results.csv') %>% 
  mutate(sig_gene = sig_burden) %>%
  filter(n_var > 1 & sig_burden ==2 ) %>% filter(annotation != 'all')
gene_450k <- read_delim(paste0(data_path, 'corr_testing_skato_gene_500k.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character()))
result_450k <- read_csv(paste0(result_path, paste0('continuous_af1e_4_ac_13_2_p_5e_1_500k_results.csv'))) %>% 
   mutate(sig_gene = sig_burden) %>%
  filter(n_var > 1 & sig_burden ==2 ) %>% filter(annotation != 'all')

pheno_info <- distinct(gene_450k %>% select(phenocode, description, n_cases_defined))
result_450k <- result_450k %>%
  merge(., pheno_info, by.x = 'pheno1', by.y = 'phenocode') %>%
  merge(., pheno_info, by.x = 'pheno2', by.y = 'phenocode')
name = 'continuous_af1e_4_ac_9_p_5e_1_high_qual_burden_300k'
figure_real_data_qq(result_300k %>% filter(pvalue >0 & n_var > 1 & sig_burden == 2 & annotation != 'all'  & !(gene%in% c('SLC39A8', 'STAB1', 'VCAN'))) ,paste0(name, '_n_var1'), save=T)

result_450k_2 <- result_450k %>%
  merge(., result_300k, by.x = c('pheno1', 'pheno2', 'gene', 'annotation'), 
        by.y = c('pheno1', 'pheno2', 'gene', 'annotation')) %>%
  mutate(pvalue =pvalue.x)

result_300k_2 <- result_450k %>%
  merge(., result_300k, by.x = c('pheno1', 'pheno2', 'gene', 'annotation'), 
        by.y = c('pheno1', 'pheno2', 'gene', 'annotation')) %>%
  mutate(pvalue =pvalue.y)
name = 'continuous_af1e_4_ac_9_p_5e_1_high_qual_burden'
figure_real_data_qq(result_450k_2 %>% filter(n_var.x > 1 & sig_burden.x == 2 & annotation != 'all' & pvalue >0 ) ,paste0(name, '_450k_subset_n_var1_p_filtered'), save=T)
figure_real_data_qq(result_300k_2 %>% filter(n_var.x > 1 & sig_burden.x == 2 & annotation != 'all' & pvalue >0 ) ,paste0(name, '_300k_subset_n_var1_p_filtered'), save=T)
