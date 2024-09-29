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
var_data <- var_data %>% filter(gene %in% c('ALB', 'ALPL', 'APOB'))
print(paste("reading gene data:", gene_file))
gene_data <- fread(cmd = paste0('gunzip -cq ', gene_file)) %>% filter(gene_symbol %in% c('ALB', 'ALPL', 'APOB'))
name <- 'ALLSPICE_leave_one_out_1e-4_example'

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
genes2 <- genes
genes2 <- c('ALB')
cutoff <- 5
for(k in genes2){
  output_path <- paste0('~/Desktop/ALLSPICE_results/', str_replace(name, 'ALL', k), '_results.csv') 
  if(file.exists(output_path)){
    test_full <- read_csv(output_path) %>%
      filter(sig_gene == 2)
    test_results <- rbind(test_results, test_full)
    next
  }
  print(k)
  sub <- sub_data %>% filter(gene==k & AC < cutoff)
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
  leave_one_out_results <- data.frame()

  for(anno in c('pLoF', 'missense|LC', 'synonymous')){
    if(anno == 'missense|LC' & k == 'APOB') next
      sub_anno_data <- test_data$sub %>% filter(annotation == anno)
      for(i in 1:nrow(sub_anno_data)){
        print(paste0(i, '/', nrow(sub_anno_data)))
        
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
              c_hat <- get_c_hat(b1_hat, b2_hat, A, r)
              lambda <- get_likelihood_test_stats(n_ind, r, b1_hat, b2_hat, c_hat, A)
              pvalue <- 1 - pchisq(as.numeric(lambda), length(b1_hat)-1)
              temp <- data.frame(pheno1, pheno2, c_hat, lambda, pvalue, gene, length(b1_hat))
              beta_temp <- data.frame(locus=sub$locus, alleles= sub$alleles,b1 = unname(t(b1_hat)), b2 = unname(t(b2_hat)), AF = c(diag(A))) %>%
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
        test <- get_test_result(sub_anno_data[-i,],
                                r = test_data$corr,
                                pheno_corr = pheno_corr,
                                pheno_list = phenos,
                                n_ind = test_data$n_ind,
                                gene=k)
        tmp_results <- test$results %>% mutate(locus = sub_anno_data[i,] %$% locus,
                                               alleles = sub_anno_data[i,] %$% alleles,
                                               )
        leave_one_out_results <- rbind(leave_one_out_results, tmp_results)
    }


  }
  


  if(!is.null(leave_one_out_results)){
    gene_result <- get_gene_level_data(data = gene_data, gene = k, phenolist = phenos)
    leave_one_out_results <- add_gene_level_data(gene_result, leave_one_out_results)
    # test_beta <- rbind(test_beta, beta)
    result_path <- '~/Desktop/ALLSPICE/'
    name <- paste0('continuous_ALL_AC_', cutoff,'_burden_syn_var_leave_one_out_500k_2024')
    write_data_fig_real(results =leave_one_out_results, name =paste0(name, '_', k), save_fig = F)
  }
}


