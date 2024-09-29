source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')

## Table S1: Poisson dispersion test 
data <- read.csv(paste0(data_path, 'gene_phewas_burden_sig_count_239.csv'), sep = '\t')
sum_data <- data %>%
  dplyr::group_by(annotation) %>%
  dplyr::summarize(p_sig = sum(n_phewas_sig > 0)/n(),
                   p_pleiotropy = sum(n_phewas_sig > 1)/n(),
                   p_pleiotropy_sig = sum(n_phewas_sig > 1)/sum(n_phewas_sig > 0),
                   mean_n_sig = mean(n_phewas_sig),
                   var_n_sig = var(n_phewas_sig),
                   stat = sum((n_phewas_sig-mean_n_sig)^2/mean_n_sig),
                   cnt = n()) 
sum_data <- sum_data %>%
  mutate(
    chisq = (cnt-1)*var_n_sig/mean_n_sig,
  ) %>%
  mutate(
    p_chisq = pchisq(chisq, cnt-1, lower.tail = F),
  )
write_csv(sum_data, paste0(result_path, 'tableS1_poisson_dispersion_test.csv'))

## Table S2: logistic regression of pleiotropic or not ~ CAF + CDS

# Gene power
data_239 <- read.csv(paste0(data_path, 'gene_phewas_burden_sig_count_239.csv'), sep = '\t') %>%
  mutate(pleiotropy = n_phewas_sig > 1)
gene_info <- read_delim('~/gene_lists/lists/data/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', delim = '\t') %>%
  dplyr::select(gene_symbol = gene, gene_id, pLI, oe_lof_upper_bin, cds_length, oe_lof_upper)

data_239 <- data_239 %>%
  merge(., gene_info, by = c('gene_symbol', 'gene_id'), all.x = T)

data_239 %>% 
  group_by(annotation) %>%
  do(tidy(glm(pleiotropy ~ CAF + cds_length ,family=binomial(link='logit'), data=.)))

# All genes
glm <- glm(pleiotropy ~ CAF + cds_length ,family=binomial(link='logit'), 
           data=data_239)
summary(glm)

# Genes with at least one association
data_239 <- data_239 %>%
  filter(n_phewas_sig > 0) 
glm <- glm(pleiotropy ~ CAF + cds_length ,family=binomial(link='logit'), 
           data=data_239)
summary(glm)

# Table S2: Group by annotation + genes with at least one association
gene_info <- read_delim('~/gene_lists/lists/data/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', delim = '\t') %>%
  dplyr::select(gene_symbol = gene, gene_id, pLI, oe_lof_upper_bin, cds_length, oe_lof_upper)

data_239 <- read.csv(paste0(data_path, 'gene_phewas_burden_sig_count_239.csv'), sep = '\t')%>%
  mutate(pleiotropy = n_phewas_sig > 1) %>%
  filter(annotation != 'pLoF|missense|LC') %>%
  filter(n_phewas_sig > 0) %>%
  merge(., gene_info, by = c('gene_symbol', 'gene_id'), all.x = T)
data_239%>% 
  group_by(annotation) %>%
  do(tidy(glm(pleiotropy ~ CAF + cds_length ,family=binomial(link='logit'), data=.)))
data_239 %>%
  mutate(interval = get_freq_interval(CAF)) %>%
  group_by(annotation, interval) %>%
  dplyr::summarize(prop = sum(pleiotropy)/n())

data_239%>% 
  group_by(annotation) %>%
  do(tidy(glm(pleiotropy ~ CAF,family=binomial(link='logit'), data=.)))
data_239 %>%
  mutate(interval = get_freq_interval(CAF)) %>%
  group_by(annotation, interval) %>%
  dplyr::summarize(prop = sum(pleiotropy)/n(),
                   total = n())

data_239 %>%
  group_by(annotation) %>%
  do(tidy(wilcox.test(CAF ~ pleiotropy, data = ., alternative = "two.sided")))

data_239 %>%
  group_by(annotation) %>%
  dplyr::summarize(tidy(cor.test(as.numeric(pleiotropy), CAF, method = "spearman")))


## Table S3: comparison to Watanabe et al. 2019
# `%+%` <- function(x, y)  mapply(sum, x, y, MoreArgs = list(na.rm = TRUE))
# data <- read_delim(paste0(data_path, 'pleiotropy_2024_gene_burden_sig_cnt_summary.txt.bgz'), delim = '\t') %>%
#   mutate(n_domain_associated = ((n_pheno_group_sig_Biomarkers > 0) %+% 
#                                   (n_pheno_group_sig_Brain > 0) %+%
#                                   (n_pheno_group_sig_Diet > 0) %+% 
#                                   (n_pheno_group_sig_Diseases > 0) %+%
#                                   (n_pheno_group_sig_Mental > 0) %+%
#                                   (n_pheno_group_sig_Physical > 0))) %>%
#   mutate(n_disease_domain_associated = 
#            ((n_disease_group_sig_A > 0) %+% (n_disease_group_sig_H1 > 0) %+% (n_disease_group_sig_M > 0) %+%
#               (n_disease_group_sig_C > 0) %+% (n_disease_group_sig_H2 > 0) %+% (n_disease_group_sig_N > 0) %+%
#               (n_disease_group_sig_D > 0) %+% (n_disease_group_sig_I > 0) %+% (n_disease_group_sig_O > 0) %+%
#               (n_disease_group_sig_E > 0) %+% (n_disease_group_sig_J > 0) %+% (n_disease_group_sig_Q > 0)%+%
#               (n_disease_group_sig_F > 0) %+% (n_disease_group_sig_K > 0) %+% (n_disease_group_sig_R > 0) %+%
#               (n_disease_group_sig_G > 0) %+% (n_disease_group_sig_L > 0)),
#          n_disease_associated = 
#            ((n_disease_group_sig_A ) %+% (n_disease_group_sig_H1 ) %+% (n_disease_group_sig_M ) +
#               (n_disease_group_sig_C ) %+% (n_disease_group_sig_H2 ) %+% (n_disease_group_sig_N ) +
#               (n_disease_group_sig_D ) %+% (n_disease_group_sig_I ) %+% (n_disease_group_sig_O ) +
#               (n_disease_group_sig_E ) %+% (n_disease_group_sig_J ) %+% (n_disease_group_sig_Q )+
#               (n_disease_group_sig_F ) %+% (n_disease_group_sig_K ) %+% (n_disease_group_sig_R ) +
#               (n_disease_group_sig_G ) %+% (n_disease_group_sig_L )))
# write_csv(data, paste0(data_path, 'pleiotropy_2024_gene_burden_sig_cnt_summary_annotated.csv'))
data <- read_csv(paste0(data_path, 'pleiotropy_2024_gene_burden_sig_cnt_summary_annotated.csv'))
N_table <- data  %>%
  group_by(annotation) %>%
  dplyr::summarize(
    total_gene = n(),
    associated = sum(n_sig_gene > 0, na.rm=T),
    pleiotropic = sum(n_sig_gene > 1, na.rm=T),
    multi_domain = sum(n_domain_associated > 1, na.rm=T),
    domain_specific = sum(n_domain_associated == 1 & n_sig_gene > 1, na.rm=T),
    trait_specific = sum(n_sig_gene == 1, na.rm=T),
    non_associated = sum(n_sig_gene == 0, na.rm=T))
print(N_table)
p_table <- data   %>%
  group_by(annotation) %>%
  dplyr::summarize(
    p_total_gene = n()/n(),
    p_associated = percent(sum(n_sig_gene > 0, na.rm=T)/n(), accuracy = 0.01),
    p_pleiotropic = percent(sum(n_sig_gene > 1, na.rm=T)/sum(n_sig_gene > 0, na.rm=T), accuracy = 0.01),
    p_multi_domain = percent(sum(n_domain_associated > 1, na.rm=T)/sum(n_sig_gene > 0, na.rm=T), accuracy = 0.01),
    p_domain_specific = percent(sum(n_domain_associated == 1 & n_sig_gene > 1, na.rm=T)/sum(n_sig_gene > 0, na.rm=T), accuracy = 0.01),
    p_trait_specific = percent(sum(n_sig_gene == 1, na.rm=T)/sum(n_sig_gene > 0, na.rm=T), accuracy = 0.01),
    p_non_associated = percent(sum(n_sig_gene == 0, na.rm=T)/n(), accuracy = 0.01))
print(p_table)
table <- N_table %>%
  merge(., p_table, by = 'annotation') %>%
  mutate(annotation = factor(annotation, levels=annotation_types2)) %>%
  arrange(annotation)
print(table)
write_csv(table, paste0(result_path, 'tableS3_comparison_watanabe_et_al.csv'))

## Table S4: proportion of pleiotropic genes at disease domain level 
data <- read_csv(paste0(data_path, 'pleiotropy_2024_gene_burden_sig_cnt_summary_annotated.csv'))
N_table <- data  %>%
  group_by(annotation) %>%
  dplyr::summarize(
    total_gene = n(),
    disease_associated = sum(n_disease_associated > 0, na.rm=T),
    pleiotropic = sum(n_disease_associated > 1, na.rm=T),
    multi_domain = sum(n_disease_domain_associated > 1, na.rm=T),
    domain_specific = sum(n_disease_domain_associated == 1 & n_disease_associated > 1, na.rm=T),
    disease_specific = sum(n_disease_associated == 1, na.rm=T))
print(N_table)
p_table <- data   %>%
  group_by(annotation) %>%
  dplyr::summarize(
    p_total_gene = n()/n(),
    p_associated = percent(sum(n_disease_associated > 0, na.rm=T)/n(), accuracy = 0.01),
    p_pleiotropic = percent(sum(n_disease_associated > 1, na.rm=T)/sum(n_disease_associated > 0, na.rm=T), accuracy = 0.01),
    p_multi_domain = percent(sum(n_disease_domain_associated > 1, na.rm=T)/sum(n_disease_associated > 0, na.rm=T), accuracy = 0.01),
    p_domain_specific = percent(sum(n_disease_domain_associated == 1 & n_disease_associated > 1, na.rm=T)/sum(n_disease_associated > 0, na.rm=T), accuracy = 0.01),
    p_trait_specific = percent(sum(n_disease_associated == 1, na.rm=T)/sum(n_disease_associated > 0, na.rm=T), accuracy = 0.01))
print(p_table)
table <- N_table %>%
  merge(., p_table, by = 'annotation') %>%
  mutate(annotation = factor(annotation, levels=annotation_types2))%>%
  arrange(annotation)
print(table)
write_csv(table, paste0(result_path, 'tableS4_pleiotropic_genes_disease_domain.csv'))

## Table S5: correlation of number of associations among phenotypic domains
pheno_domain <- read_delim(paste0(data_path, 'pheno', '_domain_level_gene_sig_burden_', tranche,'.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character())) %>%
  filter(!is.na(pheno_group))

pheno_by_domain <- pheno_domain %>%
  pivot_wider(., id_cols = colnames(pheno_domain)[1:6], values_from = pheno_group_sig_cnt, names_from = pheno_group)

pheno_by_domain_bi <- pheno_domain %>%
  pivot_wider(., id_cols = colnames(pheno_domain)[1:6], values_from = pheno_group_sig, names_from = pheno_group)

cor.test(pheno_by_domain$Biomarkers, pheno_by_domain$`Physical measures`)
cor(pheno_by_domain[, 7:12])
cor(pheno_by_domain_bi[, 7:12])


## Table S6: logistic regression on domain-level association status ~ loeuf + CDS
data <- read_csv(paste0(data_path, 'pleiotropy_2024_gene_burden_sig_cnt_summary_annotated.csv'))
gene_info <- read_delim('~/gene_lists/lists/data/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', delim = '\t') %>%
  dplyr::select(gene_symbol = gene, gene_id, pLI, oe_lof_upper_bin, cds_length, oe_lof_upper)

data <- data %>%
  merge(., gene_info, by = c('gene_symbol', 'gene_id'), all.x = T) %>%
  mutate(pleiotropy = if_else(n_sig_gene > 1, 'Pleiotropic', 'Non-pleiotropic'),
         domain_category = case_when(
           n_sig_gene == 0 ~ 'Non-associated',
           n_sig_gene == 1 ~ 'Trait-specific',
           n_sig_gene > 1 & n_domain_associated == 1 ~ 'Domain-specific',
           n_sig_gene > 1 & n_domain_associated > 1 ~ 'Multi-domain'
         ),
         disease_domain_category = case_when(
           n_disease_associated == 0 ~ 'Non-associated',
           n_disease_associated == 1 ~ 'Trait-specific',
           n_disease_associated > 1 & n_disease_domain_associated == 1 ~ 'Domain-specific',
           n_disease_associated > 1 & n_disease_domain_associated > 1 ~ 'Multi-domain'
         ),
         
  )

library(pROC)
library(plotROC)
data <- data %>% filter(annotation == 'pLoF')
logistic_reg <- function(type, cut, write=FALSE){
  print(table(data$annotation, data$domain_category))
  sub_data <- data %>%
    filter(domain_category %in% c(type, 'Non-associated')) %>%
    mutate(y = if_else(domain_category == type, 1, 0)) %>%
    dplyr::select(y, oe_lof_upper, cds_length, domain_category) %>%
    filter(complete.cases(.)) %>%
    distinct(.)
  
  model <- glm(y ~ oe_lof_upper + cds_length ,family=binomial(link='logit'), 
               data=sub_data)
  
  sub_data <- sub_data %>%
    mutate(prediction = predict(model),
           residual = residuals(model))
  auc <- roc(sub_data$y, sub_data$prediction)$auc
  
  p1 <- sub_data %>%
    filter(prediction < cut) %>%
    ggplot + aes(x = prediction, y = residual,colour = domain_category) +
    geom_point( )+ 
    geom_hline(yintercept = 0, lty = 2) +
    scale_color_manual(breaks=c(type, 'Non-associated'), 
                       labels=c(type, 'Non-associated'), 
                       values = c( '#cc3311', '#004488')) +
    labs(x = 'Prediction', y = 'Residual', color = NULL) + themes
  
  p2 <- sub_data %>%
    ggplot + aes(d = y, m = prediction) + 
    labs(x = 'False positive rate', y = 'True positive rate') + 
    geom_roc(n.cuts = 0) +
    geom_abline() + 
    annotate(x = Inf, y = 0, label = paste('AUC:', round(auc,3)), geom='text', hjust = 1) + themes
  
  pp <- ggpubr::ggarrange(p1, p2, ncol=2)
  
  model1 <- sub_data %>%
    do(tidy(glm(y ~ oe_lof_upper + cds_length ,family=binomial(link='logit'), data=.)))
  print(model1)
  model2 <- sub_data %>%
    do(glance(glm(y ~ oe_lof_upper + cds_length ,family=binomial(link='logit'), data=.)))
  print(model2)
  
  if(write){
    png(paste0('~/Desktop/logistic_regression_', type,'.png'), width=7.5, height=2.5, units = 'in', res = 300)
    print(pp)
    dev.off()
    write_csv(model1, paste0('~/Desktop/logistic_regression_', type,'_1.csv'))
    write_csv(model2, paste0('~/Desktop/logistic_regression_', type,'_2.csv'))
  }
  
  
}
logistic_reg('Multi-domain', cut = -4)
logistic_reg('Domain-specific', cut = -6)
logistic_reg('Trait-specific', cut = -6)