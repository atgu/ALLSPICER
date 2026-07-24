source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/constants.R')
library(dplyr)
library(stringr)

pheno_summary <- read_gcs_fast('gs://ukb-diverse-pops/wlu/genebass_notebooks/allpice_phenotype_summary.tsv.bgz')
print(table(pheno_summary$exclusion_reason))
write_tsv(pheno_summary, '~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/data/allspice_phenotype_summary.tsv')
