setwd('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/')
source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')

all_n_var <- c(5, 10, 20, 100)
all_c <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
all_r <- c(-1, -0.8, -0.5, -0.2, -0.1, 0, 0.1, 0.2,  0.5, 0.8, 1)
# all_pi <- c(0, 0.2, 0.5, 0.8, 1)
all_pi <- c(0.5, 0.8)
all_n_ind <- 1000

# Basic simulation formulas
get_single_geno <- function(cnt, n_ind){
  if(n_ind < cnt) stop('n_ind should be greater than cnt')
  geno <- sample(c(rep(1, cnt),rep(0, n_ind-cnt)), n_ind, replace=FALSE)
  return(geno)
}

get_geno_mat <- function(AC, n_ind){
  X <- t(as.matrix(sapply(diag(AC), get_single_geno, n_ind = n_ind)))
  return(X)
}

get_ac_mat <- function(n_var, max_cnt = 100){
  AC <- diag(floor(runif(n = n_var, min = 1, max =max_cnt)),ncol = n_var, nrow = n_var)
  return(AC)
}

get_af_mat <-function(AC, n_ind){
  A <- AC/n_ind
  return(A)
}

get_true_beta <- function(n_ind, n_var, c, pi, sigma, null=TRUE){
  b2 <- rbinom(n_var, 1, pi) * rnorm(n_var, 0, sigma)
  b1 <- rbinom(n_var, 1, pi) * rnorm(n_var, 0, sigma)
  if(null){
    b1 <- c * b2
  }
  b <- matrix(c(b1, b2), nrow = 2, byrow = T)
  return(b)
}

get_pheno_pair <- function(b, X, r){
  R <- matrix(c(1, r, r, 1), nrow = 2, byrow = T)
  MU <- b %*% X
  Y <- apply(MU, 2, mvtnorm::rmvnorm, n = 1, sigma = R)
  return(Y)
}

get_beta_hat <- function(Y, X, A, n_ind){
  b_hat <- Y %*% t(X) %*% solve(A)/n_ind
  return(b_hat)
}

get_mle_beta <- function(b1_hat, b2_hat, c, r, null=TRUE){
  if(null){
    b2_mle <- c((c*b1_hat - c*r*b2_hat - r*b1_hat + b2_hat)/(c^2 - 2*c*r + 1))
    b1_mle <- c* b2_mle
  }else{
    b2_mle <- b2_hat
    b1_mle <- b1_hat
  }
  b_mle <- matrix(c(b1_mle, b2_mle), nrow = 2, byrow = T)
  return(b_mle)
}

get_c_hat <- function(b1_hat, b2_hat, A, r){
  u <- c(b2_hat %*% A %*% t(b1_hat - r* b2_hat))
  v <- c((b2_hat - b1_hat) %*% A %*% t(b2_hat + b1_hat))
  w <- c(b1_hat %*% A %*% t(r* b1_hat- b2_hat))
  c1 <- c((-v + sqrt(v^2-4*u*w))/(2*u))
  c2 <- c((-v - sqrt(v^2-4*u*w))/(2*u))
  c_hat <- if_else(u>0, max(c1,c2), min(c1,c2))
  return(c_hat)
}

get_likelihood_test_stats <- function(n_ind, r, b1_hat, b2_hat, c, A){
  lambda <- c((n_ind/(c^2-2*c*r+1))*((b1_hat - c*b2_hat) %*% A %*% t(b1_hat - c*b2_hat)))
  return(lambda)
}

######## Simulations 
pheno_corr_test <- function(n_ind, n_var, c, r, pi, sigma, mle = TRUE, null=TRUE){
  AC <- get_ac_mat(n_var)
  A <- get_af_mat(AC, n_ind)
  X <- get_geno_mat(AC, n_ind)
  b <- get_true_beta(n_ind, n_var, c, pi, sigma, null=null)
  Y <- get_pheno_pair(b, X, r)
  b_hat <- get_beta_hat(Y, X, A, n_ind)
  b1_hat <- matrix(b_hat[1, ], nrow = 1)
  b2_hat <- matrix(b_hat[2, ], nrow = 1)
  c_hat <- get_c_hat(b1_hat, b2_hat, A, r)
  b_mle <- get_mle_beta(b1_hat, b2_hat, c, r, null=null)
  df <- n_var-1
  if(null){
    c_hat <- if_else(mle, c_hat, c)
    df <- if_else(mle, df, n_var)
  }
  lambda <- get_likelihood_test_stats(n_ind, r, b1_hat, b2_hat, c_hat, A)
  pvalue <- pchisq(lambda, df, lower.tail = FALSE)
  results <- c(c_hat = c_hat, lambda = lambda, pvalue = pvalue)
  b1 <- data.frame(cbind(b[1,], b_hat[1,], b_mle[1,]) %>% set_colnames(c('b1_true_value', 'b1_sum_stats', 'b1_mle')))
  b2 <- data.frame(cbind(b[2,], b_hat[2,], b_mle[2,]) %>% set_colnames(c('b2_true_value', 'b2_sum_stats', 'b2_mle')))
  beta <- cbind(b1, b2, AF = c(diag(A)))
  return(list(results = results, beta = beta))
}

pheno_corr_sims <- function(times, n_ind, n_var, c, r, pi, sigma, mle = TRUE, null=TRUE, write = TRUE, output=NULL, name=NULL){
  if(null){
    par_set <- expand.grid(n_ind, n_var, c, r, pi, sigma)
    colnames(par_set) <- c('n_ind', 'n_var', 'c', 'r', 'pi', 'sigma')
  }else{
    par_set <- expand.grid(n_ind, n_var, r, pi, sigma)
    colnames(par_set) <- c('n_ind', 'n_var', 'r', 'pi', 'sigma')
  }
  par_set <- par_set %>% mutate(par_set = rownames(par_set))
  results <- data.frame()
  beta <- data.frame()
  AF <- data.frame()
  for(i in 1:nrow(par_set)){
    print(i)
    print(par_set[i,])
    n_var <- par_set[i, 'n_var']
    n_ind <- par_set[i, 'n_ind']
    if(null){
      c <- par_set[i, 'c']
    }
    r <- par_set[i, 'r']
    pi <- par_set[i, 'pi']
    sigma <- par_set[i, 'sigma']
    temp <- replicate(n = times, pheno_corr_test(n_ind, n_var, c, r, pi, sigma, mle, null), simplify = FALSE)
    results <- rbind(results, map_dfr(temp, function(x) x$results, .id = 'simulation')  %>% mutate(par_set = i))
    beta <- rbind(beta, map_df(temp, function(x) x$beta, .id = 'simulation') %>% mutate(par_set = i))
  }
  results <- results %>% merge(., par_set, by = 'par_set')
  beta <- beta %>% merge(., par_set, by = 'par_set')
  if(write){
    write_csv(results, file=paste0(output, name, '_test_results.csv'))
    write_csv(beta, file=paste0(output, name, '_beta.csv'))
  }
  return(list(beta = beta, results=results))
}

######### Simulation summary
figure_sim_test_qq <- function(test_data, name=NULL, save=TRUE){
  test_data <- test_data %>% 
    group_by(par_set) %>%
    arrange(pvalue) %>%
    add_count() %>%
    mutate(observed = -log10(pvalue), 
           rank = order(pvalue), 
           expected = -(log10(rank / (n+1))))
  
  figure <- test_data %>% 
    mutate(n_var = factor(n_var, levels = c(5, 20, 100), labels= paste('N(var) = ', c(5, 20, 100))),
           c = factor(c, levels = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels= paste('c = ', c(0, 0.2, 0.4, 0.6, 0.8, 1)))) %>%
    ggplot+ aes(y = observed, x = expected, color = factor(r)) + 
    geom_point() + 
    geom_abline(intercept = 0, slope = 1) +
    labs(x = "Expected -log10(p)", y = "Observed -log10(p)", color = 'Phenotypic correlation') +
    scale_color_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
    scale_fill_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
    xlim(0,max(test_data$observed, test_data$expected)) +
    ylim(0,max(test_data$observed, test_data$expected)) +
    themes +
    facet_grid(n_var~c) 
  
  if(save){
    png(paste0(name, "_qqplot.png"), width=6, height=4.5, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figure_sim_test_corr <- function(test_data, name=NULL, save=TRUE){
  figure <- test_data %>%
    ggplot() + 
    geom_density(aes(x = c_hat, color = factor(r)), alpha = 0.5) + themes +
    geom_vline(aes(xintercept=c), lty=2)+
    labs(x = expression(bold(hat(c))), color = 'Phenotypic correlation', y=NULL) +
    xlim(c(-2,2)) +
    scale_color_brewer(palette = 'Dark2') +
    facet_grid(n_var~c)
  if(save){
    png(paste0(name, '_corr.png'), height = 6, width = 9, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figure_sim_test_c_hat_pvalue <- function(test_data, name=NULL, save=TRUE){
  figure <- test_data %>%
    mutate(max_c =  pmax(abs(c_hat), abs(1/c_hat))) %>% 
    ggplot + 
    aes(y = pvalue, x = max_c, color = factor(r)) + 
    labs(y = 'pvalue', x = 'Maximum of the absolute values of c and 1/c') +
    scale_x_log10(label = comma) +
    geom_point() + 
    theme_classic() + themes + 
    theme(axis.text.x = element_text(size = 6)) +
    scale_color_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
    facet_grid(n_var~c)
  if(save){
    png(paste0(name, '_c_hat_vs_pvalue.png'), height = 6, width = 9, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figure_sim_test_beta <- function(test_data, name=NULL, af_adjust = TRUE, save=TRUE){
  if(af_adjust){
    test_data <- test_data %>% mutate(mutiplier = sqrt(2*AF*(1-AF)))
  }else{
    test_data$mutiplier <- 1
  }
  lab <- if_else(af_adjust, ' - AF adjusted', '')
  figure <- test_data %>%
    ggplot + 
    aes(y = b1_sum_stats*mutiplier, x = b2_sum_stats*mutiplier, color = factor(r)) + 
    labs(y = expression(beta[1]), x = expression(beta[2]), title = paste('Summary statistics', lab)) +
    geom_point() + 
    theme_classic() + themes + 
    geom_abline(aes(slope = c, intercept = 0), lty=2)+
    scale_color_brewer(name = 'Phenotypic correlation', palette = 'Dark2') +
    facet_grid(n_var~c)
  if(save){
    png(paste0(name, '_beta_sum_stats.png'), height = 6, width = 9, units = 'in', res = 300)
    print(figure)
    dev.off()
  } 
  return(figure)
}


save_sim_fig <- function(data, name, mle = TRUE, save_file=TRUE, save_fig=TRUE){
  result <- data$results %>%
    filter(r %in% c(0, 0.5, 0.8, 1) & c %in% c(0, 0.2, 0.4, 0.6, 0.8, 1))
  beta <- data$beta %>%
    filter(r %in% c(0, 0.5, 0.8, 1) & c %in% c(0, 0.2, 0.4, 0.6, 0.8, 1))
  figure_sim_test_qq(result, name, save = save_fig)
  figure_sim_test_beta(beta, name, af_adjust = FALSE, save = save_fig)
  figure_sim_test_beta(beta, paste0(name, '_af_adjusted'), af_adjust = TRUE, save = save_fig)
  if(mle){
    figure_sim_test_corr(test_data=result, name, save = save_fig)
    figure_sim_test_c_hat_pvalue(test_data=result, name, save=save_fig)
  }
} 


######### random phenotype
# Test
get_real_data_rp <- function(data){
  n_ind <- data %>% summarise(mean = mean(n_cases))
  sub <- data %>% 
    pivot_wider(id_col = c('locus', 'alleles', 'annotation'), names_from = 'phenocode', values_from = 'BETA')
  sub <- data %>% 
    select(locus, alleles, annotation, AF, AC, gene) %>% 
    unique() %>% merge(., sub, by = c('locus', 'alleles', 'annotation'))
  return(list(sub=sub, n_ind=floor(n_ind$mean)))
}

get_test_result_rp <- function(data, pheno_corr, pheno_list, gene, n_ind){
  all_test <- var_test(data, pheno_corr = pheno_corr, pheno_list = pheno_list, gene = gene,  n_ind = n_ind)
  lof_test <- var_test(data %>% filter(annotation == 'pLoF'), pheno_corr = pheno_corr, pheno_list = pheno_list, gene = gene,  n_ind = n_ind)
  mis_test <- var_test(data %>% filter(annotation == 'missense|LC'), pheno_corr = pheno_corr, pheno_list = pheno_list, gene = gene,  n_ind = n_ind)
  syn_test <- var_test(data %>% filter(annotation == 'synonymous'), pheno_corr = pheno_corr, pheno_list = pheno_list, gene = gene,  n_ind = n_ind)

  all_results <- all_test$results %>% mutate(annotation = 'all')
  lof_results <- lof_test$results %>% mutate(annotation = 'pLoF')
  mis_results <- mis_test$results %>% mutate(annotation = 'missense|LC')
  syn_results <- syn_test$results %>% mutate(annotation = 'synonymous')
  results <- rbind(all_results, lof_results, mis_results, syn_results) %>% mutate(annotation = factor(annotation, levels = annotation_types))
  
  all_beta <- all_test$beta %>% mutate(annotation = 'all')
  lof_beta <- lof_test$beta %>% mutate(annotation = 'pLoF')
  mis_beta <- mis_test$beta %>% mutate(annotation = 'missense|LC')
  syn_beta <- mis_test$beta %>% mutate(annotation = 'synonymous')
  beta <- rbind(all_beta, lof_beta, mis_beta, syn_beta) %>% mutate(annotation = factor(annotation, levels = annotation_types))
  return(list(results = results, beta = beta))
}

# Summary
figure_rp_corr <- function(results, name=NULL, save=TRUE){
  figure = results %>%
    ggplot() + 
    geom_boxplot(aes(y  = c_hat, x = annotation, color = annotation, fill=annotation), alpha=0.8) +
    labs(x = 'Annotation', y = expression(bold(hat(c)))) +
    ylim(-2,2) +
    theme_classic() + themes + annotation_color_scale + annotation_fill_scale + 
    scale_x_discrete(labels = annotation_names) +
    facet_wrap(~h, labeller = label_type)
  if(save){
    png(paste0(name, '_corr.png'), height = 6, width = 9, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figure_rp_qq <- function(results, name=NULL, save=TRUE){
  results <- results %>% 
    filter(pvalue>0) %>%
    group_by(annotation, h) %>%
    arrange(pvalue) %>%
    add_count() %>%
    mutate(observed = -log10(pvalue), 
           rank = order(pvalue), 
           expected = -(log10(rank / (n+1))),
           annotation = factor(annotation, levels = annotation_types))
  
  figure <- results %>% 
    ggplot+ aes(y=observed,x=expected, color = annotation) + 
    geom_point() + 
    geom_abline(intercept = 0, slope = 1) +
    labs(x="Expected -log10(p)", y="Observed -log10(p)", color = 'Annotation') +
    annotation_color_scale + annotation_fill_scale + 
    xlim(0,max(results$observed)) +
    ylim(0,max(results$observed)) +
    themes +
    facet_grid(h~annotation, labeller = label_type) 
  
  if(save){
    png(paste0(name, "_qqplot.png"), width=8, height=6, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figure_rp_c_hat_pvalue <- function(results, name=NULL, save=TRUE){
  figure <- results %>%
    mutate(max_c =  pmax(abs(c_hat), abs(1/c_hat))) %>% 
    ggplot + 
    aes(y = pvalue, x = max_c, color = annotation) + 
    labs(y = 'pvalue', x = 'Maximum of the absolute values of c and 1/c') +
    scale_x_log10(label = comma) +
    geom_point() + 
    theme_classic() + themes + annotation_color_scale + annotation_fill_scale + 
    theme(axis.text.x = element_text(size = 6)) +
    facet_grid(~h)
  if(save){
    png(paste0(name, '_c_hat_vs_pvalue.png'), height = 6, width = 9, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figure_rp_beta <- function(beta, name=NULL, af_adjust = TRUE, save=TRUE){
  if(af_adjust){
    beta <- beta %>% mutate(mutiplier = sqrt(2*AF*(1-AF)))
  }else{
    beta$mutiplier <- 1
  }
  lab <- if_else(af_adjust, ' - AF adjusted', '')
  figure <- beta %>%
    mutate(annotation = factor(annotation, levels = annotation_types)) %>%
    ggplot + 
    aes(y = b1*mutiplier, x = b2*mutiplier, color = annotation) + 
    labs(y = expression(beta[1]), x = expression(beta[2]), title = paste('Summary statistics', lab)) +
    geom_point() + 
    theme_classic() + themes + annotation_color_scale + annotation_fill_scale + 
    scale_x_discrete(labels = annotation_names) +
    theme(axis.text.x = element_text(size = 6)) +
    facet_grid(annotation~h)
  if(save){
    png(paste0(name, '_beta_sum_stats.png'), height = 6, width = 9, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

write_data_fig_rp <- function(results, beta, name, save_file = TRUE, save_fig = TRUE){
  if(save_file){
    write_csv(results, file=paste0(name, '_results.csv'))
    write_csv(beta, file=paste0(name, '_beta.csv'))
  }
  if(save_fig){
    figure_rp_qq(results, name, save = save_fig)
    figure_rp_beta(beta, name, af_adjust = F, save = save_fig)
    figure_rp_corr(results, name, save = save_fig)
    figure_rp_c_hat_pvalue(results, name, save=save_fig)
  }
}

######### real data
get_real_data <- function(var_data, gene_data, pheno_corr, phenolist){
  corr <- pheno_corr %>%
    filter(i != j) %>%
    filter((i_phenoname %in% phenolist$phenoname) & (j_phenoname %in% phenolist$phenoname))
  n_ind <- gene_data %>%
    filter(phenoname %in% phenolist$phenoname) %>%
    select(trait_type, phenocode, pheno_sex, coding, modifier, description, n_cases, phenoname) %>%
    distinct()
  sub <- var_data %>%
    filter(phenoname %in% phenolist$phenoname) %>%
    pivot_wider(., id_cols = c("locus", "alleles", 'annotation', 'AC', 'AF', 'gene'),
                names_from = c(phenoname),
                values_from = BETA
                ) %>%
    mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation))
  return(list(corr=corr, sub=sub, n_ind=n_ind))
}

single_test <- function(data, pheno_corr, pheno1, pheno2, n_ind, sig_level, gene){
  # n_ind <- n_ind %>% filter(phenocode %in% c(pheno1, pheno2)) %>% summarise(mean=mean(n_cases))
  # n_ind <- as.numeric(n_ind$mean)
  sub <- data %>% select(1:5, pheno1, pheno2) %>% filter(complete.cases(.))
  sub <- as.data.frame(sub)
  if(nrow(sub)>1){
    A <- 2*diag(sub$AF)
  }else{
    A <- as.matrix(2*sub$AF)
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
      r <- ifelse(is.numeric(pheno_corr), pheno_corr, c(unlist(pheno_corr[((pheno_corr$i_phenocode==pheno1) &(pheno_corr$j_phenocode==pheno2)),'corr'])))
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

get_test_result <- function(data, r, pheno_corr, pheno_list, n_ind, gene){
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
      merge(., pheno_corr, by.x = c("pheno1", "pheno2"), by.y = c("i_phenocode", "j_phenocode")) %>% 
      mutate(annotation = factor(annotation, levels = annotation_types, labels = annotation_names))
    return(list(results = full_results, beta = full_beta))
  }
}

## Real data summary
figure_real_data_qq <- function(results, name=NULL, save=TRUE){
  results <- results %>% 
    filter(n_var>1) %>%
    filter(annotation != 'all') %>%
    group_by(annotation, sig_gene) %>%
    arrange(pvalue) %>%
    add_count() %>%
    mutate(observed = -log10(pvalue), 
           rank = order(pvalue), 
           expected = -(log10(rank / (n+1))),
           annotation = factor(annotation, levels = annotation_types))
  
  figure <- results %>% 
    ggplot+ aes(y=observed,x=expected, color = annotation, label = gene) + 
    geom_point(alpha = 0.5) + 
    geom_abline(intercept = 0, slope = 1) +
    labs(x="Expected -log10(p)", y="Observed -log10(p)", color = 'Annotation') +
    annotation_color_scale + annotation_fill_scale + 
    # xlim(0,max(results$observed)) +
    # ylim(0,max(results$observed)) +
    themes +
    facet_grid(~sig_gene, labeller = label_type) +
    geom_text_repel(max.overlaps = 5)
  
  if(save){
    png(paste0(name, "_qqplot.png"), width=8, height=6, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figure_real_data_corr <- function(results, name=NULL, save=TRUE){
  figure = results %>%
    mutate(max_c =  pmax(abs(c_hat), abs(1/c_hat))) %>% 
    ggplot() + aes(x = corr, y = max_c, color = annotation) +
    geom_point(alpha=0.8) +
    labs(x = 'Phenotypic correlation', y = 'Maximum of the absolute values of c and 1/c') +
    scale_y_log10(label = comma) + 
    # xlim(c(-2,2)) +
    # ylim(c(-2,2)) +
    theme_classic() + themes + annotation_color_scale +  
    facet_wrap(~annotation, nrow=2, labeller = label_type)
  if(save){
    png(paste0(name, '_corr.png'), height = 6, width = 9, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}



figure_real_data_c_hat_pvalue <- function(results, name=NULL, save=TRUE){
  figure <- results %>%
    mutate(max_c =  pmax(abs(c_hat), abs(1/c_hat))) %>% 
    ggplot + 
    aes(y = pvalue, x = max_c, color = annotation) + 
    labs(y = 'pvalue', x = 'Maximum of the absolute values of c and 1/c') +
    scale_x_log10(label = comma) +
    geom_point() + 
    theme_classic() + themes + annotation_color_scale + annotation_fill_scale + 
    theme(axis.text.x = element_text(size = 6)) +
    facet_grid(~annotation, labeller = label_type)
  if(save){
    png(paste0(name, '_c_hat_vs_pvalue.png'), height = 6, width = 9, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figure_real_data_beta <- function(beta, name=NULL, af_adjust = TRUE, save=TRUE){
  if(af_adjust){
    beta <- beta %>% mutate(mutiplier = sqrt(2*AF*(1-AF)))
  }else{
    beta$mutiplier <- 1
  }
  lab <- if_else(af_adjust, ' - AF adjusted', '')
  figure <- beta %>%
    mutate(annotation = factor(annotation, levels = annotation_types)) %>%
    ggplot + 
    aes(y = b1*mutiplier, x = b2*mutiplier, color = annotation) + 
    labs(y = expression(beta[1]), x = expression(beta[2]), title = paste('Summary statistics', lab)) +
    geom_point() + 
    theme_classic() + themes + annotation_color_scale + annotation_fill_scale + 
    scale_x_discrete(labels = annotation_names) +
    theme(axis.text.x = element_text(size = 6)) +
    facet_wrap(~annotation, nrow = 2)
  if(save){
    png(paste0(name, '_beta_sum_stats.png'), height = 6, width = 9, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

write_data_fig_real <- function(results, beta, name, save_file = TRUE, save_fig = TRUE){
  if(save_file){
    write_csv(results, file=paste0(result_path, name, '_results.csv'))
    write_csv(beta, file=paste0(result_path, name, '_beta.csv'))
  }
  if(save_fig){
    figure_real_data_qq(results, name, save = save_fig)
    # figure_real_data_beta(beta, name, af_adjust = F, save = save_fig)
    figure_real_data_corr(results, name, save = save_fig)
    figure_real_data_c_hat_pvalue(results, name, save=save_fig)
  }
}

### irnt functions
get_gene_level_data <- function(data, gene, phenolist){
  gene_result <- data %>% 
    filter(phenoname %in% phenolist$phenoname) %>%
    filter(gene_symbol==gene) %>%
    select(annotation, phenocode, sig_gene)
    # select(annotation, phenocode, description, sig_gene)
  gene_all <- gene_result %>% 
    # group_by(phenocode, description) %>% 
    group_by(phenocode) %>% 
    summarise(sig_gene = if_else(sum(sig_gene, na.rm = T)>0, 1, 0)) %>%
    mutate(annotation= 'all') %>%
    select(colnames(gene_result))
  gene_result <- rbind(gene_result, gene_all)
  return(gene_result)
}

add_gene_level_data <- function(gene_level_data, variant_test_data){
  test_full <- variant_test_data %>% 
    merge(., gene_level_data, by.x = c('annotation', 'pheno1'), 
          by.y = c('annotation', 'phenocode'), all = T) %>%
    merge(., gene_level_data, by.x = c('annotation', 'pheno2'), 
          by.y = c('annotation', 'phenocode'), all = T) %>%
    mutate(sig_gene.x = if_else(is.na(sig_gene.x), 0, sig_gene.x),
           sig_gene.y = if_else(is.na(sig_gene.y), 0, sig_gene.y),) %>%
    mutate(sig_gene = sig_gene.x + sig_gene.y) %>%
    filter(!is.na(pheno1) & !is.na(pheno2))
  return(test_full)
}


# Pipeline
run_test <- function(var_data, gene_data, p_filter, AF_upper, AC_lower, random_pheno, name){
  ## Formatting
  gene_data <- gene_data %>% 
    filter(trait_type == 'continuous') %>% 
    select(gene_symbol, annotation, phenocode, Pvalue_Burden) %>%
    mutate(sig_gene = if_else(Pvalue_Burden<6.7e-7, 1, 0))
  sub_data <- var_data %>% 
    filter(((annotation %in% c('missense', 'synonymous') & AC > AC_lower) | annotation == 'pLoF')) %>%
    filter(trait_type == 'continuous') %>% 
    filter(AF < AF_upper & Pvalue < p_filter) %>%
    mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation))
  genes <- unique(sub_data %>% select(gene))
  
  if(random_pheno){
    sub_data <- sub_data %>%
      filter( phenocode == 'random_continuous') %>%
      mutate(modifier = if_else(is.na(modifier), 1, modifier),
             phenocode = paste0(phenocode, '_', coding))
  }
  
  test_results <- data.frame()
  test_beta <- data.frame()
  for(k in genes$gene){
    print(k)
    sub <- sub_data %>% filter(gene==k)
    sub$phenocode <- as.character(sub$phenocode)
    gene_sig <- gene_data %>% 
      filter(gene_symbol==k & sig_gene == 1) 
    if(random_pheno){
      for(i in c(0.1, 0.2, 0.5, 1)){
        tmp_data <- sub %>%
          filter(modifier == i) %>%
          get_real_data_rp(.)
        phenos <- unique(sub$phenocode[sub$modifier == i])
        temp <- get_test_result_rp(tmp_data$sub, pheno_corr=0, pheno_list  = phenos, gene = k, n_ind = tmp_data$n_ind)
        tmp_test <- temp$results %>% mutate(h = i)
        tmp_beta <- temp$beta %>% mutate(h = i)
        test_results <- rbind(test_results, tmp_test)
        test_beta <- rbind(test_beta, tmp_beta)
      }
      gene_result <- get_gene_level_data(data = gene_data, gene = k, phenolist = phenos)
      test_full <- add_gene_level_data(gene_result, test_results)
    }else{
      phenos <- unique(c(unique(gene_sig$phenocode)))
      if(length(phenos)==0) next
      test_data <- get_real_data(sub, phenos)
      test <- get_test_result(test_data$sub, 
                              r = test_data$corr,
                              pheno_corr = test_data$corr, 
                              pheno_list = phenos, 
                              n_ind = test_data$n_ind, 
                              gene=k)
      results <- test$results
      beta <- test$beta
      if(!is.null(results)){
        gene_result <- get_gene_level_data(data = gene_data, gene = k, phenolist = phenos)
        test_full <- add_gene_level_data(gene_result, results)
        
        test_results <- rbind(test_results, test_full)
        test_beta <- rbind(test_beta, beta)
      }
      
    }
  }
  if(random_pheno){
    write_data_fig_rp(test_results, test_beta, name)
  }else{
    write_data_fig_real(test_results, test_beta, name)
  }
  return(test_results)
}
