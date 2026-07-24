source('~/Dropbox (Partners HealthCare)/github_repo/ALLSPICER/analysis/R/constants.R')

### Section 2
data_599 <- read_csv(paste0(data_path, 'pleiotropy_2024_gene_burden_sig_cnt_summary_annotated.csv')) %>%
  filter(n_sig_gene > 1 & annotation %in% c('pLoF', 'missense|LC'))
data_239 <- read.csv(paste0(data_path, 'gene_phewas_burden_sig_count_239.csv'), sep = '\t') %>%
  filter(n_phewas_sig > 1 & annotation %in% c('pLoF', 'missense|LC'))

table(data_599$annotation)
length(unique(data_599$gene_id))
table(data_239$annotation)
length(unique(data_239$gene_symbol))

### Table 1
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

## phenotype correlation
full_pheno_AS_table <- read_delim(paste0(data_path, 'test_allelic_series_genes.txt.bgz'), delim='\t')
pheno_corr <- fread('gunzip -c ~/Downloads/corr_estimate_syn_var_full_500k.txt.bgz')
# copy the data from ~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/data/'

genes <- c('ATM', 'PKD1', 'SLC34A3', 'TMPRSS6', 'FAM234A', 'JAK2', 'MICA')
phenotypes <- list(
  ATM = c('2453', 'C25', 'C50', 'C77', 'C78', 'C79', '130660', '30260', '30040', '30050'),
  PKD1 = c('30880', '30500', '30020', '30670', '30720', '30610', '30130', '30700'),
  SLC34A3 = c('132036', '30670', '30720', '30700', '30810'),
  TMPRSS6 = c('130622', '30060', '30270', '30750', '30020', '30070', '30030', '30040', '30050'),
  FAM234A = c('30740', '2443'),
  JAK2 = c('30150', '30210', '130670', '30100', '30070', '30770', '30080', '30090'),
  MICA = c('30120', '130696')
)


truncate_sentence <- function(x, n = 3) {
  vapply(x, function(s) {
    if (is.na(s) || !nzchar(s)) return(s)
    # normalize any Unicode spaces to regular spaces, then squish
    s <- str_replace_all(s, "\\p{Z}+", " ")
    s <- str_squish(s)
    words <- str_split(s, " ", simplify = FALSE)[[1]]
    if (length(words) > n) {
      paste(paste(words[seq_len(n)], collapse = " "), "...")
    } else {
      s
    }
  }, character(1))
}

plot_pheno_corr <- function(pheno_corr, phenotypes, gene_lab){
  sub_anno_info <- full_pheno_AS_table %>%
    filter(gene_symbol == gene_lab & Pvalue_Burden < 2.5e-6) %>%
    select(phenocode, gene_symbol, annotation) %>%
    mutate(
      annotation = factor(annotation, levels = annotation_types, labels = annotation_names)
    )%>%
    group_by(phenocode) %>%
    mutate(n = n()) %>%
    mutate(annotation = if_else(n == 2, 'Both', annotation)) %>%
    distinct()
  if(gene_lab == 'MICA'){
    sub_pheno_corr <- pheno_corr %>%
      filter((i_phenocode %in% phenotypes[gene_lab][[1]] | (i_phenocode == '20002' & i_coding == '1226')) &
             (j_phenocode %in% phenotypes[gene_lab][[1]] | (j_phenocode == '20002' & j_coding == '1226'))) %>%
      merge(., sub_anno_info, by.x = 'i_phenocode', by.y = 'phenocode', all.x=T) %>%
      merge(., sub_anno_info, by.x = 'j_phenocode', by.y = 'phenocode', all.x=T) %>%
      mutate(
        i_description = if_else(i_description =='Non-cancer illness code, self-reported', 'Hypothyroidism/myxoedema', i_description),
        j_description = if_else(j_description =='Non-cancer illness code, self-reported', 'Hypothyroidism/myxoedema', j_description)
      )
  }else if(gene_lab == 'PKD1'){
    sub_pheno_corr <- pheno_corr %>%
      filter((i_phenocode %in% phenotypes[gene_lab][[1]] | (i_phenocode == '20002' & i_coding == '1065')) &
               (j_phenocode %in% phenotypes[gene_lab][[1]] | (j_phenocode == '20002' & j_coding == '1065'))) %>%
      merge(., sub_anno_info, by.x = 'i_phenocode', by.y = 'phenocode', all.x=T) %>%
      merge(., sub_anno_info, by.x = 'j_phenocode', by.y = 'phenocode', all.x=T) %>%
      mutate(
        i_description = if_else(i_description =='Non-cancer illness code, self-reported', 'Hypertension', i_description),
        j_description = if_else(j_description =='Non-cancer illness code, self-reported', 'Hypertension', j_description)
      )
  }else{
    sub_pheno_corr <- pheno_corr %>%
      filter(i_phenocode %in% phenotypes[gene_lab][[1]] & j_phenocode %in% phenotypes[gene_lab][[1]]) %>%
      merge(., sub_anno_info, by.x = 'i_phenocode', by.y = 'phenocode', all.x=T) %>%
      merge(., sub_anno_info, by.x = 'j_phenocode', by.y = 'phenocode', all.x=T)
  }


  p <- sub_pheno_corr  %>%
    mutate(
      phenotype1_grp = interaction(truncate_sentence(str_to_sentence(i_description)), annotation.x, sep = " | "),
      phenotype2_grp = interaction(truncate_sentence(str_to_sentence(j_description)),annotation.y, sep = " | ")
    ) %>%
    ggplot() +
    geom_tile(aes(x = phenotype1_grp, y = phenotype2_grp, fill = entry)) +
    # geom_text(aes(x = phenotype1_grp, y = phenotype2_grp, label = round(entry,2)), size=4) +
    coord_equal() +
    scale_fill_gradient2(name = 'Correlation', mid="#FBFEF9",low="#0C6291",high="#A63446") +
    labs(x =NULL, y = NULL) + theme_classic() +
    themes +
    theme_classic()  +
    theme(legend.title = element_text(face = 'plain', size = 11),
          axis.title = element_text(face = 'plain', size = 11),
          plot.margin = unit(c(0.5,0,0,0), "cm"),
          axis.text = element_text(family = "mono", size = 7),
          axis.text.x = element_text(angle = 30, hjust = 1))
  return(p)
}

p1 <- plot_pheno_corr(pheno_corr, phenotypes, 'ATM')
p2 <- plot_pheno_corr(pheno_corr, phenotypes, 'PKD1')
p3 <- plot_pheno_corr(pheno_corr, phenotypes, 'SLC34A3')
p4 <- plot_pheno_corr(pheno_corr, phenotypes, 'TMPRSS6')
p5 <- plot_pheno_corr(pheno_corr, phenotypes, 'FAM234A')
p6 <- plot_pheno_corr(pheno_corr, phenotypes, 'JAK2')
p7 <- plot_pheno_corr(pheno_corr, phenotypes, 'MICA')


figure = ggpubr::ggarrange(ggpubr::ggarrange(p1, p2,
                                             labels = c('(A) ATM', '(B) PKD1'),
                                             nrow=1, hjust=0, widths = c(1, 1), align = 'h',
                                             font.label = list(size = 10, color = "black", face = "bold.italic", family = 'Arial')) +
                             theme(plot.margin = unit(c(0.5,0,0,0), "cm")),
                           ggpubr::ggarrange(p4, p6,
                                             labels = c('(C) TMPRSS6', '(D) JAK2'),
                                             nrow=1, hjust=0, widths = c(1, 1), align = 'h',
                                             font.label = list(size = 10, color = "black", face = "bold.italic", family = 'Arial'))+
                             theme(plot.margin = unit(c(0.5,0,0,0), "cm")),
                           ggpubr::ggarrange(p3, p5, p7,
                                             labels = c('(E) SLC34A3', '(F) FAM234A', '(G) MICA'),
                                             nrow=1, hjust=0, widths = c(0.95,1.05 , 1), align = 'h',
                                             font.label = list(size = 10, color = "black", face = "bold.italic", family = 'Arial'))+
                             theme(plot.margin = unit(c(0.7,0,0,0), "cm")),
                           ncol=1, hjust = 0, heights = c(1,1, 0.8),
                           font.label = list(size = 10, color = "black", face = "bold.italic", family = 'Arial')
)
png(paste0(figure_path,'figureS12.png'), height = 10, width = 12, units = 'in', res = 300)
print(figure)
dev.off()

