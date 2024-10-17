#' ALLSPICE
#'
#' ALLSPICE (ALLelic Spectrum of Pleiotropy Informed Correlated Effects)
#'
#' @param data Input data with number of rows indicating number of variants, three columns are required:
#' 1) effect sizes of variants for phenotype 1, 2) effect sizes of variants for phenotype 2, 3) allele frequency of variants
#' Note: this should include variants from ONE gene that is associated with the two phenotypes,
#' preferably of the SAME functional category after being filtered to variants with allele frequency below a certain threshold (e.g. 1e-4)
#' @param pheno_corr phenotypic correlation between the two phenotypes being tested
#' @param n_ind total number of individuals
#' @param gene name of the gene being tested, default `GENENAME`
#' @param pheno1 descriptive name of phenotype 1, default `PHENO1`
#' @param pheno2 descriptive name of phenotype 2, default `PHENO2`
#' @param beta1_field field name for effect sizes of variants on phenotype 1, default `BETA1`
#' @param beta2_field field name for effect sizes of variants on phenotype 2, default `BETA2`
#' @param af_field field name for allele frequencies of variants, default `AF`
#'
#' @return A list of summary statistics from ALLSPICE test
#' including phenotype names, gene names, MLE of slope c, ALLSPICE test statistic - lambda, pvalue from a chi-square distribution, total number of variants being tested
#' @examples
#' data <- data.frame(x = rnorm(10), y = rnorm(10), z = runif(10, 0,1))
#' ALLSPICE(data,pheno_corr=0.5,n_ind=10000,beta1_field='x',beta2_field='y',af_field='z')
#' @import readr
#' @import mvtnorm
#' @importFrom dplyr filter select if_else
#' @importFrom magrittr %>% set_colnames
#' @importFrom stats complete.cases pchisq runif rbinom rnorm
#' @export

ALLSPICE <- function(data, pheno_corr, n_ind, gene='GENENAME', pheno1='PHENO1', pheno2='PHENO2', beta1_field = 'BETA1', beta2_field = 'BETA2', af_field = 'AF'){
  fields <- colnames(data)
  if(!(beta1_field %in% fields) | !(beta2_field %in% fields) | !(af_field %in% fields)) {
    stop("Effect sizes of variants on the two phenotypes or Allele frequency is missing from the data, please check")
  }
  data <- format_ALLSPICE_data(data=data, beta1_field = beta1_field, beta2_field = beta2_field, af_field = af_field)
  data <- data %>% dplyr::select('beta1', 'beta2', 'AF')
  data <- data[complete.cases(data), ]
  if(nrow(data)>1){
    A <- 2*diag(data$AF)
  }else{
    A <- as.matrix(2*data$AF)
  }
  b1_hat <- t(as.matrix(data$beta1))
  b2_hat <- t(as.matrix(data$beta2))
  r <- pheno_corr
  c_hat <- get_c_hat(b1_hat, b2_hat, A, r)
  lambda <- get_likelihood_test_stats(n_ind, r, b1_hat, b2_hat, c_hat, A)
  pvalue <- 1 - pchisq(lambda, length(b1_hat)-1)
  results <- data.frame(pheno1, pheno2, gene, c_hat, lambda, pvalue, length(b1_hat))
  colnames(results) <- c('pheno1', 'pheno2', 'gene', 'c_hat', 'lambda', 'pvalue','n_var')
  return(results)
}


#' ALLSPICE_simulation
#'
#' Simulate data and run ALLSPICE
#'
#' @param n_ind total number of individuals
#' @param n_var total number of variants
#' @param c slope between the two sets of variant effect sizes, only applicable when `null` == TRUE
#' @param r phenotypic correlation between the two phenotypes
#' @param pi probability of variant of having no effect on the phenotype
#' @param sigma variance of the two sets of effect sizes
#' @param mle whether to use MLE of c to compute the test statistic, use true c value if FALSE
#' @param null whether to simulate data under the null hypothesis (no linear relationship) or the alternative hypothesis
#'
#' @return A list of two pieces of results:
#' 1) ALLSPICE test results
#' 2) effect size table: true effect size simulated, effect size estimate from linear model, effect size estimated from MLE
#' @examples
#' ALLSPICE_simulation(n_ind=10000, n_var=100, c=0.6, r=0.5, pi=0.5, sigma=1, mle = TRUE, null=TRUE)
#'
#' @export

ALLSPICE_simulation <- function(n_ind, n_var, c, r, pi, sigma, mle = TRUE, null=TRUE){
  AC <- get_ac_mat(n_var)
  A <- get_af_mat(AC, n_ind)
  X <- get_geno_mat(AC, n_ind)
  b <- get_true_beta(n_var, c, pi, sigma, null=null)
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
  results <- c(n_ind = n_ind, n_var = n_var, c = c, r=r, pi=pi, sigma =sigma, c_hat = c_hat, lambda = lambda, pvalue = pvalue)
  b1 <- data.frame(cbind(b[1,], b_hat[1,], b_mle[1,]) %>% set_colnames(c('b1_true_value', 'b1_sum_stats', 'b1_mle')))
  b2 <- data.frame(cbind(b[2,], b_hat[2,], b_mle[2,]) %>% set_colnames(c('b2_true_value', 'b2_sum_stats', 'b2_mle')))
  beta <- cbind(b1, b2, AF = c(diag(A)))
  return(list(results = results, beta = beta))
}
