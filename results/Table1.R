source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')

full_pheno_AS_table <- read_delim(paste0(data_path, 'test_allelic_series_genes.txt.bgz'), delim='\t')
full_pheno_AS_table_wide <- full_pheno_AS_table %>%
  dplyr::select(-total_variants) %>%
  pivot_wider(names_from = 'annotation', values_from = c('Pvalue_Burden', 'BETA_Burden')) %>%
  arrange(gene_symbol) 
write_csv(full_pheno_AS_table_wide, paste0(result_path, 'allelic_series_full_table.csv'))

mis_lof_combined <- read_delim(paste0(data_path, 'gene_lof_mis_sig_burden.txt.bgz'), delim = '\t',col_types = cols(phenocode = col_character()))
allelic_series <- mis_lof_combined %>% filter(mis_lof_combined$lof_not_mis_cnt > 0 & mis_lof_combined$mis_not_lof_cnt > 0) 

gene_point_check <- read.csv(paste0(data_path, 'gene_point_check.csv'), sep = ',')
gene_view <- gene_point_check %>% dplyr::select(2:4, 14:15, 17)


allelic_series_genes_1 <- allelic_series %>%
  filter(lof_mis_diff_1) %>%
  dplyr::select(gene_symbol) %>%
  merge(., gene_view %>% filter(annotation %in% c('pLoF', 'missense|LC')), by = 'gene_symbol') %>%
  group_by(gene_symbol, annotation) %>%
  mutate(phenos = paste0(description, collapse = ", ")) %>%
  dplyr::select(-description, - gene_id, -description_more, -pheno_group) %>%
  distinct()  %>% 
  group_by(annotation, gene_symbol) %>%
  pivot_wider(., id_cols = gene_symbol, names_from = annotation, values_from =phenos) %>%
  dplyr::select(gene_symbol, pLoF, Missense = `missense|LC`)
allelic_series_genes_1
# write_csv(allelic_series_genes_1, paste0(result_path, 'allelic_series_genes_1.csv'))


allelic_series_genes_2 <- allelic_series %>%
  filter(lof_mis_diff_2) %>%
  dplyr::select(gene_symbol) %>%
  merge(., gene_view %>% filter(annotation %in% c('pLoF', 'missense|LC') & pheno_group %in% c('Biomarkers', 'Diseases')), by = 'gene_symbol') %>%
  group_by(gene_symbol, annotation) %>%
  mutate(phenos = paste0(description, collapse = ", ")) %>%
  dplyr::select(-description, - gene_id, -description_more, -pheno_group) %>%
  distinct() %>%
  group_by(annotation, gene_symbol) %>%
  summarise(ColB = paste0(phenos, collapse = "")) %>%
  pivot_wider(., id_cols = gene_symbol, names_from = annotation, values_from =ColB) %>%
  dplyr::select(gene_symbol, pLoF, Missense = `missense|LC`)
allelic_series_genes_2
# write_csv(allelic_series_genes_2, paste0(result_path, 'allelic_series_genes_2.csv'))
