source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')


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
    return(list(results = full_results, beta = full_beta))
  }
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
    select(annotation, phenoname, sig_gene)
  return(gene_result)
}

add_gene_level_data <- function(gene_level_data, variant_test_data){
  test_full <- variant_test_data %>% 
    merge(., gene_level_data, by.x = c('annotation', 'pheno1'), 
          by.y = c('annotation', 'phenoname'), all = T) %>%
    merge(., gene_level_data, by.x = c('annotation', 'pheno2'), 
          by.y = c('annotation', 'phenoname'), all = T) %>%
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
