source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/ALLSPICE/ALLSPICE_utils.R')

############## Setup Arguments ##############
option_list <- list(
  make_option(c("-a","--alpha"), type="double", default=5e-8,
              help="defult significant level [default= %default]", metavar="double"),
  make_option(c("-b","--bp_col"), type="character",
              help="bp column [default= %default]", metavar="character"),
  make_option(c("-c","--chr_col"), type="character",
              help="chromosome column [default= %default]", metavar="character"),
  make_option(c("-d","--dn_cut"), type="double", default=0,
              help="Downsample cutoff [default= %default], -log10(p-value) below this cutoff value are down-sampled", metavar="double"),
  make_option(c("-e","--exp_pcol"), type="character", default="Pvalue_expected_log10",
              help="expected pvalue column [default= %default]", metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-g", "--gclambda"), type="double", default=NULL,
              help="LambdaGC", metavar="character"),
  make_option(c("-i","--locus"), type="character", default = 'locus',
              help="if given then chr:bp:ref:alt OR chr:bp identifier assumed and chr and bp are read from there [default= %default]", metavar="character"),
  make_option(c("-l", "--log10p"), type="logical", default=FALSE,
              help="whether the p-values are -log10 or not [default= %default]", metavar="logical"),
  make_option(c("-n", "--name"), type="character", default='gene_symbol',
              help="label field to print on manhattan plot [default= %default]", metavar="character"),
  make_option(c("-m", "--multi_pheno"), type="character",
              help="phenotype column [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-p","--pval_col"), type="character", default="Pvalue_log10",
              help="pvalue column [default= %default]", metavar="character"),
  make_option(c("-q","--no_qqplot"), type="logical", default = FALSE,
              help="Making Manhattan plots only, ignore QQ plots [default= %default]", metavar="logical"),
  make_option(c("-s","--siglim"), type="integer", default=20,
              help="-upper limit for the number of significant bps to highlight [default= %default]", metavar="integer"),
  make_option(c("-t","--top"), type="integer", default=3,
              help="-number of top significant bps to highlight [default= %default]", metavar="integer"),
  make_option(c("-w","--af_col"), type="character", default='AF_Allele2',
              help="allele frequency column [default= %default]", metavar="character")
);


############## Load Arguments ##############
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser, positional_arguments=0);

print(str(opt))
alpha <- opt$options$alpha
dn_cut <- opt$options$dn_cut
top <- opt$options$top
siglim <- opt$options$siglim
pcol <- opt$options$pval_col
no_qqplot <- opt$options$no_qqplot
af_col <- opt$options$af_col
log_p <- opt$options$log10p
locus_field <- opt$options$locus
name <- opt$options$name
exp_pcol <- opt$options$exp_pcol
lambda_gc <- opt$options$gclambda

file <- opt$options$file
file <- '~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/data/corr_testing_burden_var_0.01_500k.txt.bgz'
print(paste("Reading file:", file))

data <- fread(paste0('gunzip -cq ~/Downloads/corr_testing_burden_var_0.01_500k.txt.bgz'))
data <- read_delim(file, delim = '\t')

test_data <- data %>%
  filter(gene == 'GIGYF1' & annotation == 'pLoF')


if(endswith(file, 'gz')){
  data <- fread(paste0('gunzip -cq ', file))
}else{
  data <- read_delim(file, delim = '\t')
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