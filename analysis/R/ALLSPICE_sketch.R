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
genes2 <- c('ALB')
cutoff <- 5
for(k in genes2){
  output_path <- paste0('~/Desktop/ALLSPICE/results/syn_c_hat/', str_replace(name, 'ALL', k), '_results.csv') 
  if(file.exists(output_path)){
    test_full <- read_csv(output_path) %>%
      filter(sig_gene == 2)
    test_results <- rbind(test_results, test_full)
    next
  }
  print(k)
  sub <- sub_data %>% filter(gene==k & AC <cutoff)
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
  
  get_c_hat <- function(b1_hat, b2_hat, A, r){
    u <- c(b2_hat %*% A %*% t(b1_hat - r* b2_hat))
    v <- c((b2_hat - b1_hat) %*% A %*% t(b2_hat + b1_hat))
    w <- c(b1_hat %*% A %*% t(r* b1_hat- b2_hat))
    c1 <- c((-v + sqrt(v^2-4*u*w))/(2*u))
    c2 <- c((-v - sqrt(v^2-4*u*w))/(2*u))
    c_hat <- if_else(u>0, max(c1,c2), min(c1,c2))
    results <- c(u, v, w, c1, c2, c_hat)
    names(results) <- c('u', 'v', 'w', 'c1', 'c2', 'c_hat')
    return(results)
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
        n_ind <- n %>%
          filter(phenoname %in% c(pheno1, pheno2)) %>%
          summarise(mean=mean(n_cases))
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
        r <- ifelse(is.numeric(pheno_corr), pheno_corr, c(unlist(pheno_corr[((pheno_corr$i_phenoname==pheno1) &(pheno_corr$j_phenoname==pheno2)),'corr'])))
        phenoname1 = pheno1
        phenoname2 = pheno2
        gene_name = gene
        # c_hat <- get_c_hat(b1_hat, b2_hat, A, r)
        c_hat <- 1.3103
        lambda <- get_likelihood_test_stats(n_ind, r, b1_hat, b2_hat, c_hat, A)
        # lambda <- get_likelihood_test_stats(n_ind, r, b1_hat, b2_hat, c_hat['c_hat'], A)
        # lambda2 <- get_likelihood_test_stats(n_ind, r, b1_hat, b2_hat, c_hat['c2'], A)
        # lambda1 <- get_likelihood_test_stats(n_ind, r, b1_hat, b2_hat, c_hat['c1'], A)
        print(paste0('lambda (c_hat): ', lambda))
        # print(paste0('c_1:', c_hat[2]))
        # print(paste0('lambda (c_1): ', get_likelihood_test_stats(n_ind, r, b1_hat, b2_hat, c_hat[2], A)))
        # print(paste0('c_1:', c_hat[3]))
        # print(paste0('lambda (c_2): ', get_likelihood_test_stats(n_ind, r, b1_hat, b2_hat, c_hat[3], A)))
        pvalue <- 1 - pchisq(as.numeric(lambda), length(b1_hat)-1)
        temp <- data.frame(pheno1, pheno2, c_hat, lambda, pvalue, gene, length(b1_hat))
        # temp <- data.frame(pheno1, pheno2, c_hat['u'], c_hat['v'], c_hat['w'], c_hat['c1'], c_hat['c2'], c_hat['c_hat'], lambda, lambda1, lambda2, pvalue, gene, length(b1_hat))
        beta_temp <- data.frame(locus=sub$locus, alleles= sub$alleles,b1 = unname(t(b1_hat)), b2 = unname(t(b2_hat)), AF = c(diag(A))) %>%
          mutate(pheno1 = pheno1, pheno2  = pheno2, gene = gene)
        results <- rbind(results, temp)
        beta <- rbind(beta, beta_temp)
      }
    }
    if(nrow(results)>0){
      # colnames(results) <- c('pheno1', 'pheno2', 'u', 'v', 'w', 'c1', 'c2', 'c_hat', 'lambda', 'lambda1', 'lambda2', 'pvalue', 'gene', 'n_var')
      colnames(results) <- c('pheno1', 'pheno2','c_hat', 'lambda','pvalue', 'gene', 'n_var')
      return(list(results = results, beta = beta))
    }
  }
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
    result_path <- '~/Desktop/ALLSPICE/results/'
    name <- paste0('continuous_ALL_AC_', cutoff,'_burden_syn_var_syn_c_hat_500k_2024')
    write_data_fig_real(test_full, str_replace(name, 'ALL', k), save_fig = F)
  }

}

result_path <- '~/Desktop/ALLSPICE/'
write_data_fig_real(test_results %>% distinct(), name, save_fig = F)


## Check lambda differences
### AF < 1e-4 original
data1 <- read_csv('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/results/continuous_final_AF_1e_4_burden_syn_var_corr_500k_2023_results.csv')%>%
  filter(pheno1 == 'continuous_30600_both_sexes__irnt' & pheno2 == 'continuous_30680_both_sexes__irnt' & gene == 'ALB') %>%
  select(pheno1, pheno2, annotation, gene, n_var, c_hat, lambda, pvalue)
### AF < 1e-4 synonymous c_hat
data2 <- read_csv('~/Desktop/ALLSPICE/continuous_ALL_AF_1e_4_burden_syn_var_c_hat_500k_2024_results_corr.csv')%>%
  filter(pheno1 == 'continuous_30600_both_sexes__irnt' & pheno2 == 'continuous_30680_both_sexes__irnt'& gene == 'ALB')%>%
  select(pheno1, pheno2, annotation, gene, n_var, c_hat, lambda, pvalue)
### AF < 1e-4 leave-one-out
data3 <- read_csv('~/Desktop/ALLSPICE/ALLSPICE_leave_one_out_1e-4_example_ALB_results.csv') %>%
  filter(pheno1 == 'continuous_30600_both_sexes__irnt' & pheno2 == 'continuous_30680_both_sexes__irnt'& gene == 'ALB')%>%
  select(pheno1, pheno2, annotation, gene, n_var, c_hat, lambda, pvalue, locus, alleles)
### AF < 1e-4 leave-one-out synounymous c_hat
data4 <- read_csv('~/Desktop/ALLSPICE/continuous_ALL_AC_80_burden_syn_var_leave_one_out_syn_c_hat_500k_2024_ALB_results.csv') %>%
  filter(pheno1 == 'continuous_30600_both_sexes__irnt' & pheno2 == 'continuous_30680_both_sexes__irnt'& gene == 'ALB')%>%
  select(pheno1, pheno2, annotation, gene, n_var, c_hat, lambda, pvalue, locus, alleles)

### AC < 5 original
data5 <- read_csv('~/Desktop/ALLSPICE/continuous_ALL_AC_5_burden_syn_var_corr_500k_2024_results_corr.csv')%>%
  filter(pheno1 == 'continuous_30600_both_sexes__irnt' & pheno2 == 'continuous_30680_both_sexes__irnt'& gene == 'ALB')%>%
  select(pheno1, pheno2, annotation, gene, n_var, c_hat, lambda, pvalue)
### AC <5 synonymous c_hat
data6<- read_csv('~/Desktop/ALLSPICE/results/continuous_ALB_AC_5_burden_syn_var_syn_c_hat_500k_2024_results.csv')%>%
  filter(pheno1 == 'continuous_30600_both_sexes__irnt' & pheno2 == 'continuous_30680_both_sexes__irnt'& gene == 'ALB')%>%
  select(pheno1, pheno2, annotation, gene, n_var, c_hat, lambda, pvalue)
### AC < 5 leave-one-out
data7 <- read_csv('~/Desktop/ALLSPICE/continuous_ALL_AC_5_burden_syn_var_leave_one_out_500k_2024_ALB_results.csv') %>%
  filter(pheno1 == 'continuous_30600_both_sexes__irnt' & pheno2 == 'continuous_30680_both_sexes__irnt'& gene == 'ALB')%>%
  select(pheno1, pheno2, annotation, gene, n_var, c_hat, lambda, pvalue, locus, alleles)
### AC < 5 leave-one-out synounymous c_hat
data8 <- read_csv('~/Desktop/ALLSPICE/continuous_ALL_AC_5_burden_syn_var_leave_one_out_syn_c_hat_500k_2024_ALB_results.csv') %>%
  filter(pheno1 == 'continuous_30600_both_sexes__irnt' & pheno2 == 'continuous_30680_both_sexes__irnt'& gene == 'ALB')%>%
  select(pheno1, pheno2, annotation, gene, n_var, c_hat, lambda, pvalue, locus, alleles)

data <- merge(merge(
  merge(data1, data2, by = c('pheno1', 'pheno2', 'annotation', 'gene'), suffix = c('_af1e_4_original', '_af1e_4_syn_c_hat')),
  merge(data3, data4, by = c('pheno1', 'pheno2', 'annotation', 'gene', 'locus', 'alleles'), suffix = c('_af1e_4_leaveoneout_original', '_af1e_4_leaveoneout_syn_c_hat')),
  by = c('pheno1', 'pheno2', 'annotation', 'gene'), all.y = T),
  merge(
    merge(data5, data6, by = c('pheno1', 'pheno2', 'annotation', 'gene'), suffix = c('_ac5_original', '_ac5_syn_c_hat')),
    merge(data7, data8, by = c('pheno1', 'pheno2', 'annotation', 'gene', 'locus', 'alleles'), suffix = c('_ac5_leaveoneout_original', '_ac5_leaveoneout_syn_c_hat')),
    by = c('pheno1', 'pheno2', 'annotation', 'gene'), all.y = T),
  by = c('pheno1', 'pheno2', 'annotation', 'gene', 'locus', 'alleles'))


data <- data %>%
  mutate(af1e_4_lambda_diff = lambda_af1e_4_original - lambda_af1e_4_leaveoneout_original,
         af1e_4_lambda_diff_c_hat = lambda_af1e_4_syn_c_hat - lambda_af1e_4_leaveoneout_syn_c_hat,
         ac5_lambda_diff = lambda_ac5_original - lambda_ac5_leaveoneout_original,
         ac5_lambda_diff_c_hat = lambda_ac5_syn_c_hat - lambda_ac5_leaveoneout_syn_c_hat,
) %>% filter(annotation == 'missense|LC')


#### results plotting etc.
raw_results_500k <- read_pleiotropy_results('burden', '500k') 
# pheno_corr <- read_delim(paste0(data_path, 'corr_estimate_syn_var_full_500k.txt.bgz'), delim = '\t',
#                          col_types = cols(i_phenocode = col_character(), j_phenocode = col_character()))
# raw_results_500k <- read_csv( '~/Desktop/ALLSPICE/continuous_ALL_AF_1e_4_burden_syn_var_c_hat_500k_2024_results.csv')
# raw_results_500k <- raw_results_500k %>%
#   merge(., pheno_corr %>% 
#           mutate(i_phenoname = paste0(i_trait_type, '_', i_phenocode, '_',  i_pheno_sex, '_',  i_coding, '_',  i_modifier), 
#                  j_phenoname = paste0(j_trait_type, '_', j_phenocode, '_',  j_pheno_sex, '_',  j_coding, '_',  j_modifier),
#                  corr = entry) %>%
#           select(i_phenoname, j_phenoname, corr), by.x = c('pheno1', 'pheno2'), by.y = c('i_phenoname', 'j_phenoname'), all.x = T)
results_500k <- modify_results_table(raw_results_500k, 'burden', '500k')

top_triplets <- read_delim(paste0(data_path, "top_significant_triplets.txt.bgz"), delim='\t', col_types = cols(phenocode = col_character())) %>% 
  mutate(coding = if_else(is.na(coding), '', coding)) %>%
  mutate(phenoname = paste0(trait_type, '_', phenocode, '_',  pheno_sex, '_',  coding, '_',  modifier)) %>%
  mutate(BETA_adjusted = sqrt(2*AF*(1-AF))*BETA)


# p <- results_500k %>%
#   ggplot + aes(x = c_hat, y = -log10(pvalue), color = annotation) + 
#   geom_vline(xintercept = 1, lty=2) +
#   geom_point() + themes + 
#   xlim(-10, 10) + 
#   annotation_color_scale + 
#   annotation_fill_scale 
# output_figure(p, 'syn_c_hat_pvalue', height = 3, width = 8)
# 
p <- results_500k %>%
  # filter(pvalue>0.05) %>%
  mutate(annotation = factor(annotation, levels = annotation_types)) %>%
  ggplot + aes(x = corr, y = c_hat, color = annotation, size = n_var) +
  geom_point() + themes +
  geom_abline(slope = 1, intercept = 0) +
  ylim(-10, 10) +
  annotation_color_scale +
  annotation_fill_scale +
  scale_alpha(range = c(0.5, 1)) +
  facet_grid(~annotation)
# output_figure(p, 'syn_c_hat_pheno_corr', height = 3, width = 8)
# 
# p <- results_500k %>%
#   mutate(annotation = factor(annotation, levels = annotation_types)) %>%
#   ggplot + aes(x =  c_hat, color = annotation) + 
#   geom_density() + themes + 
#   xlim(-5, 5) + 
#   annotation_color_scale + 
#   annotation_fill_scale +
#   scale_alpha(range = c(0.5, 1)) + 
#   facet_grid(~annotation, scale = 'free')
# output_figure(p, 'c_hat_dist_zoomed', height = 3, width = 8)
# 
# 
# p <- results_500k %>%
#   ggplot + aes(x = n_var, y = -log10(pvalue), color = annotation) + 
#   geom_point() + themes + 
#   annotation_color_scale + 
#   annotation_fill_scale 
# output_figure(p, 'nvar_pvalue_syn_s_hat.png', height = 3, width = 8)
# 
# p <- results_500k %>%
#   mutate(annotation = factor(annotation, levels = annotation_types)) %>%
#   ggplot + aes(x = corr, y = -log10(pvalue), color = annotation, size = n_var) + 
#   geom_point() + themes + 
#   annotation_color_scale + 
#   annotation_fill_scale +
#   facet_grid(~annotation)
# output_figure(p, 'pvalue_pheno_corr_syn_c_hat', height = 3, width = 8)
# 
# p <- results_500k %>%
#   mutate(annotation = factor(annotation, levels = annotation_types)) %>%
#   ggplot + aes(x =  -log10(pvalue), color = annotation) + 
#   geom_density() + themes + 
#   geom_vline(xintercept = -log10(0.05/11810), lty=2) + 
#   annotation_color_scale + 
#   annotation_fill_scale +
#   scale_alpha(range = c(0.5, 1)) + 
#   facet_grid(~annotation, scale = 'free')
# output_figure(p, 'pvalue_dist_zoomed', height = 3, width = 8)
# 
# top_hits <- results_500k %>% filter(pvalue < 1e-6 & sig_gene==2) 
# gene_annts <- top_hits %>% merge(., results_500k, by = colnames(results_500k)[c(1, 3:5, 7:9, 11)], suffixes = c('.x', ''), all.x = T) %>% select("gene", "annotation", "pheno1", "description1", "pheno2", "description2")
# for(i in 1:nrow(top_hits)){
#   tmp <- top_hits[i, ]
#   c_hat<- results_500k %>%
#     filter(gene == tmp[,'gene'] & pheno1 == tmp[,'pheno1'] & pheno2 == tmp[,'pheno2']) %>%
#     select(annotation, c_hat)
#   # if(!dir.exists(paste0('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/for_siwei/files_for_siwei/', tmp[,'gene'], '/'))) next
#   annotations <- unlist(merge(tmp, gene_annts, by = c("gene", "pheno1", "description1", "pheno2", "description2"), suffixes = c('.x', '')) %>% select(annotation))
#   figure_path <- '~/Desktop/'
#   figure_beta_triplets_all_annt(data=top_triplets  %>% filter(AF < 1e-4), 
#                                 gene=tmp[,'gene'], 
#                                 phenocode1=tmp[,'pheno1'], 
#                                 phenocode2=tmp[,'pheno2'], 
#                                 c_hat = c_hat,
#                                 tmp[,'description1'], 
#                                 tmp[,'description2'], 
#                                 figure_name ='beta_figure', 
#                                 save=T, 
#                                 threshold = 1e-2)
# }
# 
# 
# gene_name <- 'APOB'
# leave_one_out_results <- read_csv(paste0('~/Desktop/ALLSPICE/ALLSPICE_leave_one_out_1e-4_example_', gene_name,'_results.csv'))
# pext <- read_delim('~/Downloads/tx_annotation_coding_variants_APOB.txt.bgz', delim = '\t')
# top_hits <- results_500k %>% filter(pvalue < 1e-6 & sig_gene==2) 
# gene_annts <- top_hits %>% merge(., results_500k, by = colnames(results_500k)[c(1, 3:5, 7:9, 11)], suffixes = c('.x', ''), all.x = T) %>% select("gene", "annotation", "pheno1", "description1", "pheno2", "description2")
# gene_top_hits <- top_hits %>% filter(gene == gene_name)
# for(i in 1:nrow(gene_top_hits)){
#   tmp <- gene_top_hits[i, ]
#   gene_name <- tmp[,'gene']
#   phenocode1 <- tmp[,'pheno1']
#   phenocode2 <- tmp[,'pheno2']
#   pheno1_name <- tmp$description1
#   pheno2_name <- tmp$description2
#   c_hat<- results_500k %>%
#     filter(gene == gene_name & pheno1 ==phenocode1  & pheno2 == phenocode2) %>%
#     select(annotation, c_hat)
#   annotations <- unlist(merge(tmp, gene_annts, by = c("gene", "pheno1", "description1", "pheno2", "description2"), suffixes = c('.x', '')) %>% select(annotation))
#   figure_path <- '~/Desktop/'
#   sub_info <- top_triplets  %>% filter(AF < 1e-4) %>%
#     filter(gene == gene_name & phenoname %in% c(phenocode1, phenocode2))
#   wide_info <- sub_info %>%
#     mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
#     filter(annotation %in% c('pLoF', 'missense|LC', 'synonymous')) %>%
#     mutate(
#       BETA_adjusted = sqrt(2*AF*(1-AF))*BETA
#     ) %>%
#     pivot_wider(names_from = phenoname,  values_from = c('BETA', 'Pvalue'), id_cols = c('locus', 'alleles','AF', 'gene', 'annotation')) %>%
#     mutate(significance = case_when(
#       get(paste0('Pvalue_',phenocode1)) <= threshold & get(paste0('Pvalue_',phenocode2)) <= threshold ~ 'Both',
#       get(paste0('Pvalue_',phenocode1)) > threshold & get(paste0('Pvalue_',phenocode2)) <= threshold ~ pheno2_name,
#       get(paste0('Pvalue_',phenocode1)) <= threshold & get(paste0('Pvalue_',phenocode2)) > threshold ~ pheno1_name,
#       get(paste0('Pvalue_',phenocode1)) > threshold & get(paste0('Pvalue_',phenocode2)) > threshold ~ 'None',
#     )) %>%
#     mutate(significance = factor(significance, levels = c('Both', pheno1_name, pheno2_name, 'None')),
#            annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
#     merge(., leave_one_out_results %>% filter(pheno1 == phenocode1 & pheno2 == phenocode2) %>% 
#             select(locus, alleles, annotation, leave_one_out_pvalue = pvalue, leave_one_out_c_hat=c_hat), by = c('locus', 'alleles', 'annotation'), all.x=T) %>%
#     merge(., raw_results_500k %>% filter(gene==gene_name & pheno1==phenocode1 & pheno2==phenocode2) %>%
#             mutate(annotation = factor(annotation, levels=annotation_types)) %>% select(annotation, pvalue, c_hat), by = 'annotation') %>%
#     mutate(
#       # magnitude_change_p = if_else(is.na(leave_one_out_pvalue), 0, abs(log10(leave_one_out_pvalue/if_else(pvalue==0, 1e-320, pvalue)))),
#       magnitude_change_p = if_else(is.na(leave_one_out_pvalue) | annotation == 'synonymous', 0, abs(leave_one_out_pvalue-pvalue)),
#       magnitude_change_c = if_else(is.na(leave_one_out_c_hat), 0, abs(log10(leave_one_out_c_hat/c_hat)))) %>%
#     mutate(annotation = factor(annotation, levels=annotation_types)) %>% 
#     merge(., pext %>% select(-annotation), by = c('locus', 'alleles'))
#   figure <- wide_info %>%
#     ggplot + 
#     aes(x=get(paste0('BETA_',phenocode1)), y=get(paste0('BETA_',phenocode2)), color = annotation)  +
#     # aes(x=sqrt(2*AF*(1-AF))*pvalue_info[,paste0('beta_',phenocode1)], y=sqrt(2*AF*(1-AF))*pvalue_info[,paste0('beta_',phenocode2)], color = annotation, size = -log(mean_p)) +
#     # labs(x=pheno1_name, y=pheno2_name, title =paste0(gene_name, '-', annotation) ) + 
#     geom_point(aes(pch = significance, size = mean_proportion)) +
#     geom_abline(data = c_hat %>%
#                   mutate(annotation = factor(annotation, levels=annotation_types)), aes(slope = 1/c_hat, intercept = 0, color = annotation), lwd =0.5) +
#     geom_vline(xintercept = 0, lty=2, lwd = 0.25) + 
#     geom_hline(yintercept = 0, lty=2, lwd = 0.25) + 
#     annotation_color_scale + annotation_fill_scale + 
#     labs(x=paste0(pheno1_name), y=pheno2_name, title = NULL) + 
#     scale_shape_manual(name=paste0('Nominal significance (', threshold, ')'), breaks = c('Both', pheno1_name, pheno2_name, 'None'), values=c("\u25CF", "\u25D0","\u25D1", "\u25CB")) + 
#     scale_size(range = c(2, 8)) + 
#     scale_alpha(range = c(0.2, 1)) + 
#     geom_text_repel(data = raw_results_500k %>% filter(gene==gene_name & pheno1==phenocode1 & pheno2==phenocode2) %>%
#                       mutate(annotation = factor(annotation, levels=annotation_types)), aes(x=2, y= -2, label=formatC(pvalue, format = "e", digits = 2), color=annotation), vjust = 1, size =5)+
#     facet_wrap(~annotation, labeller = label_type) +
#     guides(color = "none", size = 'none') +
#     theme(legend.position = 'top') 
#     # geom_point(data = wide_info %>% filter(locus == 'chr2:21006019'), aes(x = get(paste0('BETA_',phenocode1)), y=get(paste0('BETA_',phenocode2))), size=5, pch=1, color='#6CA6CB') +
#     # geom_text_repel(data = wide_info %>% filter(locus == 'chr2:21006019'), aes(x = get(paste0('BETA_',phenocode1)), y=get(paste0('BETA_',phenocode2)), label = 'chr2:21006019:CA:C'), size=3, color='#6CA6CB', face = 'bold', vjust = 1, hjust = 0.5)
#   
#   png(paste0(figure_path, '/ALLSPICE/',gene_name, '_', str_split(phenocode1, '_')[[1]][2], '_', str_split(phenocode2, '_')[[1]][2],'_1e_4_pext.png'), height = 3, width = 8, units = 'in', res = 300)
#   print(figure)
#   dev.off()
# }
