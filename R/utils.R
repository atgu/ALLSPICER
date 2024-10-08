#' get_single_geno
#'
#' simulation function: simulate genotype information for one locus, where `cnt` samples out of `n_ind` has the mutation
#'
#' @param cnt number of individuals with the mutation
#' @param n_ind total number of individuals
#'
#' @return A binary vector representing the genotype information of `n_ind` individuals for a particular locus, where `cnt` entries has value 1.
#' @examples
#' geno <- get_single_geno(cnt = 100, n_ind = 10000)
#' @import readr
#' @import mvtnorm
#' @importFrom dplyr filter select if_else
#' @importFrom magrittr %>% set_colnames
#' @importFrom stats complete.cases pchisq runif rbinom rnorm
#' @export
get_single_geno <- function(cnt, n_ind){
  if(n_ind < cnt) stop('n_ind should be greater than cnt')
  geno <- sample(c(rep(1, cnt),rep(0, n_ind-cnt)), n_ind, replace=FALSE)
  return(geno)
}

#' get_geno_mat
#'
#' simulation function: simulate genotype information for a set of loci with allele counts `AC`
#'
#' @param AC allele counts of loci (length `m`)
#' @param n_ind total number of indicitions
#'
#' @return An `n_ind`x`m` matrix of genotype information of `n_ind` individuals and `m` variants
#' @examples
#' geno_mat <- get_geno_mat(AC = c(20, 50, 10, 1, 5), n_ind = 10000)
#' @export

get_geno_mat <- function(AC, n_ind){
  X <- t(as.matrix(sapply(diag(AC), get_single_geno, n_ind = n_ind)))
  return(X)
}


#' get_ac_mat
#'
#' simulation function: simulate allele count information for `n_var` variants, with a maximum allele count `max_cnt`
#'
#' @param n_var total number of variants
#' @param max_cnt maximum allele count, default 100
#'
#' @return A `n_var`x`n_var` diagnal matrix of allele count information for `n_var` variants
#' @examples
#' ac_mat <- get_ac_mat(n_var=100, max_cnt = 100)
#' @export

get_ac_mat <- function(n_var, max_cnt = 100){
  AC <- diag(floor(runif(n = n_var, min = 1, max =max_cnt)),ncol = n_var, nrow = n_var)
  return(AC)
}


#' get_af_mat
#'
#' simulation function: compute allele frequency information variants with allele counts stored in diagonal matrix `AC` from a population of sample size `n_ind`
#'
#' @param AC a diagonal matrix of allele count information for all variants
#' @param n_ind total number of individuals in the population
#'
#' @return A `n_var`x`n_var` diagnal matrix of allele frequency information for `n_var` (dimension of `AC`) variants
#' @examples
#' af_mat <- get_af_mat(AC = c(20, 50, 10, 1, 5), n_ind = 10000)
#' @export

get_af_mat <-function(AC, n_ind){
  A <- AC/n_ind
  return(A)
}


#' get_true_beta
#'
#' simulation function: simulate true effect size information of `n_var` variants for two phenotypes
#'
#' @param n_var total number of variants
#' @param c slope between the two sets of variant effect sizes, only applicable when `null` == TRUE
#' @param pi probability of variant of having no effect on the phenotype
#' @param sigma variance of the two sets of effect sizes
#' @param null whether to simulate data under the null hypothesis (no linear relationship) or the alternative hypothesis
#'
#' @return A 2x`n_var` matrix of effect size information for `n_var` variants (first row corresponds to the first phenotype, second row corresponds to the second phenotype)
#' @examples
#' true_beta <- get_true_beta(n_var=100, c=0.6, pi=0.5, sigma=1, null=TRUE)
#' @export

get_true_beta <- function(n_var, c, pi, sigma, null=TRUE){
  b2 <- rbinom(n_var, 1, pi) * rnorm(n_var, 0, sigma)
  b1 <- rbinom(n_var, 1, pi) * rnorm(n_var, 0, sigma)
  if(null){
    b1 <- c * b2
  }
  b <- matrix(c(b1, b2), nrow = 2, byrow = T)
  return(b)
}


#' get_pheno_pair
#'
#' simulation function: simulate true phenotype values of a pair of phenotypes
#'
#' @param b true effect size matrix of variants on the two phenotypes
#' @param X genotype matrix
#' @param r phenotypic correlation between the two phenotypes
#'
#' @return A 2x`n_ind` matrix of phenotype information (first row corresponds to the first phenotype, second row corresponds to the second phenotype)
#' @examples
#' AC <- get_ac_mat(n_var=100)
#' X <- get_geno_mat(AC, n_ind=10000)
#' b <- get_true_beta(n_var=100, c=0.6, pi=0.5, sigma=1, null=TRUE)
#' Y <- get_pheno_pair(b=b, X=X, r=0.5)
#' @export

get_pheno_pair <- function(b, X, r){
  R <- matrix(c(1, r, r, 1), nrow = 2, byrow = T)
  MU <- b %*% X
  Y <- apply(MU, 2, mvtnorm::rmvnorm, n = 1, sigma = R)
  return(Y)
}


#' get_beta_hat
#'
#' simulation function: compute effect sizes estimated form linear regression model
#'
#' @param Y phenotype information
#' @param X genotype information
#' @param A Allele frequency information
#' @param n_ind total number of individuals
#' @return A 2x`n_var` matrix of estimated effect size information (first row corresponds to the first phenotype, second row corresponds to the second phenotype)
#' @examples
#' AC <- get_ac_mat(n_var=100)
#' A <- get_af_mat(AC=AC, n_ind=10000)
#' X <- get_geno_mat(AC, n_ind=10000)
#' b <- get_true_beta(n_var=100, c=0.6, pi=0.5, sigma=1, null=TRUE)
#' Y <- get_pheno_pair(b=b, X=X, r=0.5)
#' b_hat <- get_beta_hat(Y=Y, X=X, A=A, n_ind=10000)
#' @export

get_beta_hat <- function(Y, X, A, n_ind){
  b_hat <- Y %*% t(X) %*% solve(A)/n_ind
  return(b_hat)
}


#' get_mle_beta
#'
#' ALLSPICE function: compute the effect size estimates that maximize the likelihood (maximum likelihood estimate - MLE) conditioning on c
#'
#' @param b1_hat estimated effect size of the first phenotype across all variants
#' @param b2_hat estimated effect size of the second phenotype across all variants
#' @param c slope between the two sets of variant effect sizes, only applicable when `null` == TRUE
#' @param r phenotypic correlation between the two phenotypes
#' @param null whether to simulate data under the null hypothesis (no linear relationship) or the alternative hypothesis
#' @return A 2x`n_var` matrix of MLE estimated effect size information (first row corresponds to the first phenotype, second row corresponds to the second phenotype)
#' @examples
#' AC <- get_ac_mat(n_var=100)
#' A <- get_af_mat(AC=AC, n_ind=10000)
#' X <- get_geno_mat(AC, n_ind=10000)
#' b <- get_true_beta(n_var=100, c=0.6, pi=0.5, sigma=1, null=TRUE)
#' Y <- get_pheno_pair(b=b, X=X, r=0.5)
#' b_hat <- get_beta_hat(Y=Y, X=X, A=A, n_ind=10000)
#' b1_hat <- matrix(b_hat[1, ], nrow = 1)
#' b2_hat <- matrix(b_hat[2, ], nrow = 1)
#' b_mle <- get_mle_beta(b1_hat=b1_hat, b2_hat=b2_hat, c=0.6, r=0.5, null=TRUE)
#' @export

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

#' get_c_hat
#'
#' ALLSPICE function: compute the slope `c` that maximize the likelihood (maximum likelihood estimate - MLE)
#'
#' @param b1_hat estimated effect size of the first phenotype across all variants
#' @param b2_hat estimated effect size of the second phenotype across all variants
#' @param A Allele frequency information
#' @param r phenotypic correlation between the two phenotypes
#'
#' @return the MLE of slope between two sets of effect sizes
#' @examples
#' AC <- get_ac_mat(n_var=100)
#' A <- get_af_mat(AC=AC, n_ind=10000)
#' X <- get_geno_mat(AC, n_ind=10000)
#' b <- get_true_beta(n_var=100, c=0.6, pi=0.5, sigma=1, null=TRUE)
#' Y <- get_pheno_pair(b=b, X=X, r=0.5)
#' b_hat <- get_beta_hat(Y=Y, X=X, A=A, n_ind=10000)
#' b1_hat <- matrix(b_hat[1, ], nrow = 1)
#' b2_hat <- matrix(b_hat[2, ], nrow = 1)
#' c_hat <- get_c_hat(b1_hat=b1_hat, b2_hat=b2_hat, A=A, r=0.5)
#' @export

get_c_hat <- function(b1_hat, b2_hat, A, r){
  u <- c(b2_hat %*% A %*% t(b1_hat - r* b2_hat))
  v <- c((b2_hat - b1_hat) %*% A %*% t(b2_hat + b1_hat))
  w <- c(b1_hat %*% A %*% t(r* b1_hat- b2_hat))
  c1 <- c((-v + sqrt(v^2-4*u*w))/(2*u))
  c2 <- c((-v - sqrt(v^2-4*u*w))/(2*u))
  c_hat <- if_else(u>0, max(c1,c2), min(c1,c2))
  return(c_hat)
}


#' get_likelihood_test_stats
#'
#' ALLSPICE function: compute the maximum likelihood ratio of the ALLSPICE test statistic
#'
#' @param n_ind total number of individuals
#' @param r phenotypic correlation between the two phenotypes
#' @param b1_hat estimated effect size of the first phenotype across all variants
#' @param b2_hat estimated effect size of the second phenotype across all variants
#' @param c MLE of the slope between the two sets of variant effect sizes
#' @param A Allele frequency information
#'
#' @return A single numeric value representing the test statistic of ALLSPICE (maximum likelihood ratio)
#' @examples
#' AC <- get_ac_mat(n_var=100)
#' A <- get_af_mat(AC=AC, n_ind=10000)
#' X <- get_geno_mat(AC, n_ind=10000)
#' b <- get_true_beta(n_var=100, c=0.6, pi=0.5, sigma=1, null=TRUE)
#' Y <- get_pheno_pair(b=b, X=X, r=0.5)
#' b_hat <- get_beta_hat(Y=Y, X=X, A=A, n_ind=10000)
#' b1_hat <- matrix(b_hat[1, ], nrow = 1)
#' b2_hat <- matrix(b_hat[2, ], nrow = 1)
#' c_hat <- get_c_hat(b1_hat=b1_hat, b2_hat=b2_hat, A=A, r=0.5)
#' lambda <- get_likelihood_test_stats(n_ind=10000, r=0.5, b1_hat=b1_hat, b2_hat=b2_hat, c=c_hat, A=A)
#' @export

get_likelihood_test_stats <- function(n_ind, r, b1_hat, b2_hat, c, A){
  lambda <- c((n_ind/(c^2-2*c*r+1))*((b1_hat - c*b2_hat) %*% A %*% t(b1_hat - c*b2_hat)))
  return(lambda)
}

#' format_ALLSPICE_data
#'
#' data formatting function: format raw data to be loaded into ALLSPICE
#'
#' @param data raw input data
#' @param beta1_field field name of effect size for the first phenotype
#' @param beta2_field field name of effect size for the second phenotype
#' @param af_field field name of allele frequency information
#'
#' @return a data frame containing effect sizes of variants on two phenotypes and their allele frequency information
#' @examples
#' data <- data.frame(x = rnorm(10), y = rnorm(10), z = runif(10, 0,1))
#' data <- format_ALLSPICE_data(data=data, beta1_field = 'x', beta2_field = 'y', af_field = 'z')
#' @export

format_ALLSPICE_data <- function(data, beta1_field, beta2_field, af_field){
  data <- as.data.frame(data)
  data$beta1 <- data[[beta1_field]]
  data$beta2 <- data[[beta2_field]]
  data$AF <- data[[af_field]]
  return(data)
}


