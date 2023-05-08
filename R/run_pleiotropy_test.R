###### Functions #####
# Basic simulation equations
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

single_test <- function(data, pheno_corr, pheno1, pheno2, n_ind, gene){
  sub <- data %>% select(1:6, pheno1, pheno2) %>% filter(complete.cases(.))
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
  return(results)
}


###### Preparation ########
# Gene: CD36
# Annotation: pLoF
# phenotype1: 30070 - Red blood cell (erythrocyte) distribution width
# phenotype2: 30110 - Platelet distribution width
# phenotypic correlation:
genename <- 'CD36'
phenocode1 <- '30070'
phenocode2 <- '30110'
phenolist <- c(phenocode1, phenocode2)
phenocorr_300k <- -0.01419300
phenocorr_500k <- -0.01566100
mean_n_cases_300k <- 273539.5
mean_n_cases_500k <- 473357

# User specify:
TRANCHE <- '300k' # or 300k
PATH <- '~/Downloads/'

###### Run Test #############
n_ind <- if_else(TRANCHE == '300k', mean_n_cases_300k, mean_n_cases_500k)
pheno_corr <- if_else(TRANCHE == '300k', phenocorr_300k, phenocorr_500k)
var_data <- read_csv(paste0('~/Downloads/var_subset_cd36_', TRANCHE,'_modified_for_test.csv'))
single_test(var_data, pheno_corr, phenocode1, phenocode2, n_ind, genename)