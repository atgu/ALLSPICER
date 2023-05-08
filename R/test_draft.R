setwd('~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/')
source('~/Dropbox (Partners HealthCare)/ukbb_exomes_pleiotropy/R/simulations.R')


ll_function <- function(A, b_hat, R, b, n_ind){
  term1 <- sum(diag(A%*%t(b_hat)%*%solve(R)%*%b))
  term2 <- sum(diag(A%*%t(b)%*%solve(R)%*%b))
  ll <- n_ind*term1 - n_ind * term2 /2
  return(ll)
}

get_ind_beta <- function(n_var, c, sigma, null=TRUE){
  ind1 <- sample(c(0,1), n_var, replace = T)
  ind2 <- abs(1-ind1)
  b2 <- ind2 * rnorm(n_var, 0, sigma)
  b1 <- ind1 * rnorm(n_var, 0, sigma)
  if(!null){
    b2 <- rnorm(n_var, 0, sigma)
    b1 <- c * b2
  }
  b <- matrix(c(b1, b2), nrow = 2, byrow = T)
  return(b)
}

get_cos_similarity <- function(b){
  b1 <- b[1, ]
  b2 <- b[2, ]
  prod <- c(b1 %*% b2)
  cos <- prod/(sqrt(sum(b1^2))*(sqrt(sum(b2^2))))
  return(cos)
}

get_cos_derivation <- function(b, index){
  index0 <- if_else(index==1, 2, 1)
  b1 <- b[1, ]
  b2 <- b[2, ]
  prod <- c(b1 %*% b2)
  b1_dis <- sqrt(sum(b1^2))
  b2_dis <- sqrt(sum(b2^2))
  deriv <- b[index0, ]/(b1_dis*b2_dis) - b[index, ]*prod/(sqrt(sum(b[index, ]^2))^3 * sqrt(sum(b[index0, ]^2)))
  return(deriv)
}

get_mle_beta_gd <- function(b0, b_hat, A, n_ind, r, learn_rate, ll_conv_threshold, cos_conv_threshold, max_iter, tuning_par, plot = T) {
  converged <- F
  iterations <- 0

  b <- b0
  cos <- get_cos_similarity(b)
  sig <- if_else((abs(cos) < cos_conv_threshold), 0, sign(cos))
  R <- matrix(c(1, r, r, 1), nrow = 2, byrow = T)
  max_ll <- ll_function(A = A, b_hat = b_hat, R = R, b = b, n_ind = n_ind) - tuning_par * abs(cos)

  logliks <- c()
  cosines <- c()
  b_lst <- list()
  while((converged == F) & (iterations <= max_iter)) {
    ## Implement the gradient descent algorithm
    deriv_1 <- (n_ind*diag(A)/(1-r^2))* (b_hat[1,] - r* b_hat[2,] - b[1,] + r*b[2,]) - tuning_par * sig * get_cos_derivation(b, 1) # constraint term
    deriv_2 <- (n_ind*diag(A)/(1-r^2))* (b_hat[2,] - r* b_hat[1,] - b[2,] + r*b[1,]) - tuning_par * sig * get_cos_derivation(b, 2) # constraint term
    deriv <- matrix(c(deriv_1, deriv_2), nrow = 2, byrow = T)
    b <- b + learn_rate * deriv
    cos <- get_cos_similarity(b)
    sig <- if_else((abs(cos) < cos_conv_threshold), 0, sign(cos))
    ll <- ll_function(A = A, b_hat = b_hat, R = R, b = b, n_ind = n_ind) - tuning_par * abs(cos)
    if(is.nan(ll)) {
      converged = T
      print(paste("Warning:not converged in", iterations, "iterations"))
      # print(paste("Maximum likelihood:", max_ll))
      if(plot)plot(1:length(logliks), logliks)
      if(plot) plot(1:length(cosines), cosines)
      next
    }
    logliks <- c(logliks, ll)
    cosines <- c(cosines, cos)
    b_lst <- c(b_lst, list(b))
    iterations = iterations + 1
    b <- b_lst[[which(logliks == max(logliks))[1]]]
    if(abs(ll - max_ll) < ll_conv_threshold & abs(cos) < cos_conv_threshold) {
      converged = T
      print(paste("Iteration:", iterations))
    }
    max_ll  <- max(max_ll, ll)
  }
  print('cosine similarity of b MLE: ')
  print(get_cos_similarity(b))
  if(plot)plot(1:length(logliks), logliks + tuning_par*abs(cosines))
  if(plot) plot(1:length(cosines), cosines)
  return(list(b = b, convergence = converged, iterations = iterations, cosines=cosines, ll = logliks))
}

# beta point check
beta_point_check <- function(n_var = 10,
                             n_ind =1000,
                             r = 0.6,
                             sigma = 1,
                             learn_rate=0.001,
                             ll_conv_threshold=0.0001,
                             cos_conv_threshold=0.001,
                             max_iter = 1000,
                             null = TRUE,
                             c = 2,
                             tuning_par){
  AC <- get_ac_mat(n_var)
  A <- get_af_mat(AC, n_ind)
  X <- get_geno_mat(AC, n_ind)
  if(null){c <- NULL}else{c <- c}
  b0 <- get_ind_beta(n_var=n_var,c = c,sigma=sigma, null=null)
  Y <- get_pheno_pair(b0, X, r)
  R <- matrix(c(1, r, r, 1), nrow = 2, byrow = T)
  b_hat <- get_beta_hat(Y, X, A, n_ind)
  ll_alt <- ll_function(A = A, b_hat = b_hat, R = R, b = b_hat, n_ind = n_ind)
  gd_results <- get_mle_beta_gd(b0 = b0,
                                b_hat = b_hat,
                                A = A,
                                n_ind =n_ind,
                                r = r,
                                learn_rate=learn_rate,
                                ll_conv_threshold=ll_conv_threshold,
                                cos_conv_threshold=cos_conv_threshold,
                                max_iter = max_iter,
                                tuning_par = tuning_par,
                                plot = T)
  b_gd <- gd_results$b
  plot(b_hat[1,], b_hat[2,], col='red', xlab = 'b1', ylab = 'b2')
  points(b_gd[1,], b_gd[2,], col='blue')
  ll_null <- ll_function(A = A, b_hat = b_hat, R = R, b = b_gd, n_ind = n_ind)
  lambda <- 2 * (ll_alt - ll_null)
  pvalue <- pchisq(lambda, 1, lower.tail = FALSE)
  print(paste('Pvalue:', pvalue))
  return(gd_results)
}

gd <- beta_point_check(n_var = 100, n_ind =1000,
                             r = 0.6,
                             sigma = 1,
                             learn_rate=0.001,
                             ll_conv_threshold=0.0001,
                             cos_conv_threshold=0.001,
                             max_iter = 1000,
                             null = TRUE,
                             c = 2,
                             tuning_par = 0)


pheno_ind_test <- function(n_ind, n_var, r, sigma, tuning_par, learn_rate=1e-3, ll_conv_threshold=1e-3, cos_conv_threshold=1e-4, c = NULL, max_iter = 1000, null = TRUE){
  AC <- get_ac_mat(n_var)
  A <- get_af_mat(AC, n_ind)
  X <- get_geno_mat(AC, n_ind)
  not_ind <- TRUE

  b0 <- get_ind_beta(n_var=n_var,c = c,sigma=sigma, null=null)
  Y <- get_pheno_pair(b0, X, r)
  R <- matrix(c(1, r, r, 1), nrow = 2, byrow = T)
  b_hat <- get_beta_hat(Y, X, A, n_ind)
  ll_alt <- ll_function(A = A, b_hat = b_hat, R = R, b = b_hat, n_ind = n_ind)

  gd_results <- get_mle_beta_gd(b0 = b0, b_hat = b_hat, A = A, n_ind =n_ind, r = r, learn_rate=learn_rate, ll_conv_threshold=ll_conv_threshold, cos_conv_threshold=cos_conv_threshold, max_iter = max_iter, tuning_par = tuning_par, plot = F)
  b_mle <- gd_results$b
  convergence <- gd_results$convergence
  iterations <- gd_results$iterations
  ll_null <- ll_function(A = A, b_hat = b_hat, R = R, b = b_mle, n_ind = n_ind)

  lambda <- 2 * (ll_alt - ll_null)
  cosine_mle <- get_cos_similarity(b_mle)
  pvalue <- pchisq(lambda, 1, lower.tail = FALSE)
  results <- c(lambda = lambda, pvalue = pvalue, ll_null = ll_null, ll_alt = ll_alt, convergence = convergence, iterations=iterations, cosine_mle = cosine_mle)
  b1 <- data.frame(cbind(b0[1,], b_hat[1,], b_mle[1,]) %>% set_colnames(c('b1_true_value', 'b1_sum_stats', 'b1_mle')))
  b2 <- data.frame(cbind(b0[2,], b_hat[2,], b_mle[2,]) %>% set_colnames(c('b2_true_value', 'b2_sum_stats', 'b2_mle')))
  beta <- cbind(b1, b2, AF = c(diag(A)))
  return(list(results = results, beta = beta))
}

pheno_ind_sims <- function(times, n_ind, n_var, r, sigma, c, tuning_par, learn_rate, ll_conv_threshold, cos_conv_threshold, max_iter, null=TRUE, write = TRUE, output=NULL, name=NULL){
  par_set <- expand.grid(n_ind, n_var, r, sigma, tuning_par, learn_rate, ll_conv_threshold, cos_conv_threshold)
  colnames(par_set) <- c('n_ind', 'n_var', 'r', 'sigma', 'tuning_par', 'learn_rate', 'll_conv_threshold', 'cos_conv_threshold')

  par_set <- par_set %>% mutate(par_set = rownames(par_set))
  results <- data.frame()
  beta <- data.frame()
  AF <- data.frame()
  for(i in 1:nrow(par_set)){
    print(i)
    print(par_set[i,])
    n_var <- par_set[i, 'n_var']
    n_ind <- par_set[i, 'n_ind']
    r <- par_set[i, 'r']
    sigma <- par_set[i, 'sigma']
    tuning_par <- par_set[i, 'tuning_par']
    learn_rate <- par_set[i, 'learn_rate']
    ll_conv_threshold <- par_set[i, 'll_conv_threshold']
    cos_conv_threshold <- par_set[i, 'cos_conv_threshold']
    temp <- replicate(n = times, pheno_ind_test(n_ind=n_ind, n_var=n_var, r=r, sigma=sigma, tuning_par = tuning_par,
                                                learn_rate=learn_rate, ll_conv_threshold=ll_conv_threshold, cos_conv_threshold=cos_conv_threshold, max_iter = max_iter, c=c, null = null), simplify = FALSE)
    results <- rbind(results, map_dfr(temp, function(x) x$results, .id = 'simulation')  %>% mutate(par_set = i))
    beta <- rbind(beta, map_df(temp, function(x) x$beta, .id = 'simulation') %>% mutate(par_set = i))
    while(sum(results$iterations > 0 & results$par_set == i) < times){
      temp <- replicate(n = times, pheno_ind_test(n_ind=n_ind, n_var=n_var, r=r, sigma=sigma, tuning_par = tuning_par,
                                                  learn_rate=learn_rate, ll_conv_threshold=ll_conv_threshold, cos_conv_threshold=cos_conv_threshold, max_iter = max_iter, c=c, null = null), simplify = FALSE)
      results <- rbind(results, map_dfr(temp, function(x) x$results, .id = 'simulation')  %>% mutate(par_set = i))
      beta <- rbind(beta, map_df(temp, function(x) x$beta, .id = 'simulation') %>% mutate(par_set = i))
    }

  }
  results <- results %>% merge(., par_set, by = 'par_set') %>% filter(iterations > 0 )
  beta <- beta %>% merge(., par_set, by = 'par_set') %>% merge(., results %>% select('par_set', 'simulation'), by = c('par_set', 'simulation'))
  if(write){
    write_csv(results, file=paste0(output, name, '_test_results.csv'))
    write_csv(beta, file=paste0(output, name, '_beta.csv'))
  }
  return(list(beta = beta, results=results))
}

figure_sim_test_qq <- function(test_data, name=NULL, save=TRUE){
  test_data <- test_data %>%
    group_by(par_set) %>%
    arrange(pvalue) %>%
    add_count() %>%
    mutate(observed = -log10(pvalue),
           rank = order(pvalue),
           expected = -(log10(rank / (n+1))))

  figure <- test_data %>%
    ggplot+ aes(y = observed, x = expected, color = factor(r)) +
    geom_point(size = 1) +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = "Expected -log10(p)", y = "Observed -log10(p)", color = 'Phenotypic correlation') +
    scale_color_brewer(name = 'Phenotypic correlation', palette = 'Spectral') +
    scale_fill_brewer(name = 'Phenotypic correlation', palette = 'Spectral') +
    xlim(0, 5) +
    ylim(0, 5) +
    themes  +
    facet_grid(tuning_par~n_var)

  if(save){
    png(paste0(name, "_qqplot.png"), width=8, height=8, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

# null simulation results
n_var = c(5, 10, 20, 50, 100)
n_ind =1000
r = seq(0, 0.9, 0.1)
sigma = 1
learn_rate = c(1e-3)
ll_conv_threshold=c(1e-4)
cos_conv_threshold=c(1e-3)
# tuning_par = c(10, 12, 14, 16, 20, 24, 28)
tuning_par = c(10, 14, 28)
results <- pheno_ind_sims(times = 20,
                          n_ind=n_ind,
                          n_var=n_var,
                          r=r,
                          sigma=sigma,
                          tuning_par = tuning_par,
                          learn_rate=learn_rate,
                          ll_conv_threshold=ll_conv_threshold,
                          cos_conv_threshold=cos_conv_threshold,
                          max_iter = 1000,
                          c = NULL,
                          null = TRUE,
                          write = TRUE,
                          output = result_path, name = 'ind_test_null_sim20_max_pen_ll_final_selected_ver0516')
null_result <- results$results
null_b <- results$beta

null_result %>%
  figure_sim_test_qq(., name = paste0(figure_path,'ind_test_null_sim20_cosine_max_ll_min_cos_final_selected_LR_1e_3'), save = T)

n_var = c(5, 10, 20, 50, 100)
n_ind =1000
r = seq(0, 0.9, 0.1)
sigma = 1
learn_rate = c(1e-3)
ll_conv_threshold=c(1e-2)
cos_conv_threshold=c(1e-2)
tuning_par = c(500, 1000, 2000, 5000)
results <- pheno_ind_sims(times = 20,
                          n_ind=n_ind,
                          n_var=n_var,
                          r=r,
                          sigma=sigma,
                          tuning_par = tuning_par,
                          learn_rate=learn_rate,
                          ll_conv_threshold=ll_conv_threshold,
                          cos_conv_threshold=cos_conv_threshold,
                          max_iter = 1000,
                          c = 2,
                          null = FALSE,
                          write = TRUE,
                          output = result_path, name = 'ind_test_alt_sim20_max_pen_ll_final_selected_ver0516')
alt_result <- results$results
alt_b <- results$beta

alt_result %>%
  figure_sim_test_qq(., name = paste0(figure_path,'ind_test_alt_sim20_max_pen_ll_final_selected'), save = T)

