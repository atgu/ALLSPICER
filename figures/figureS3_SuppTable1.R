source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')
##### Process haploinsufficient genes
# haploinsufficient_genes <- data.frame() 
# for(name in c('haploinsufficiency_moderate_curated_2016.tsv', 'haploinsufficiency_severe_curated_2016.tsv', 'haploinsufficiency_mild_curated_2016.tsv')){
#   gene_list =  as.data.frame(fread(paste0('~/gene_lists/lists/', name), quote="", header = F))
#   haploinsufficient_genes <- rbind(haploinsufficient_genes, gene_list)
# }
# haploinsufficient_genes <- unique(haploinsufficient_genes)
# write_tsv(haploinsufficient_genes, '~/gene_lists/lists/haploinsufficient.tsv', quote = NULL, col_names = FALSE)

##### Process ClinGen haploinsufficient genes 
# clinGen1 <- fread(paste0('~/gene_lists/lists/clingen_level3_genes_2015_02_27.tsv'), quote="", header = F)
# clinGen2 <- fread(paste0('~/gene_lists/lists/clingen_level3_genes_2018_09_13.tsv'), quote="", header = F)
# sum(clinGen1$V1 %in% clinGen2$V1)

data_599 <- read.csv(paste0(data_path, 'gene_phewas_burden_sig_count_599.csv'), sep = '\t') %>%
  mutate(interval = get_freq_interval(freq=CAF))  %>%
  mutate(interval = factor(interval, levels = c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'),
                           labels = c('(0.01%, 0.1%]', '(0.1%, 1%]', '(1%, 10%]', paste0('(10%, ', bquote("\U221E"), ' )') )))%>% 
  filter(n_phewas_sig >= 1)

pLoF_mis_599 <- data_599 %>% filter(annotation %in% c('pLoF', 'missense|LC'))
no_combine_599 <- data_599 %>% filter(annotation != 'pLoF|missense|LC')

data_239 <- read.csv(paste0(data_path, 'gene_phewas_burden_sig_count_239.csv'), sep = '\t') %>%
  mutate(interval = get_freq_interval(freq=CAF))  %>%
  mutate(interval = factor(interval, levels = c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'),
                           labels = c('(0.01%, 0.1%]', '(0.1%, 1%]', '(1%, 10%]', paste0('(10%, ', bquote("\U221E"), ' )') )))%>% 
  filter(n_phewas_sig >= 1)

pLoF_mis_239 <- data_239 %>% filter(annotation %in% c('pLoF', 'missense|LC'))
no_combine_239 <- data_239 %>% filter(annotation != 'pLoF|missense|LC')

# test <- 'burden'
# figure <- gene_list_pleiotropy_figure(pLoF_mis, test='burden', panel = T, filter_cat=F, overwrite = T)
# figure 
# 
# png(paste0(figure_path, 'final_2024/figure3_gene_list_', test,'_2_anno_v1.png'), height = 5, width = 5, units = 'in', res = 300)
# print(figure)
# dev.off()
# 
# figure <- gene_list_pleiotropy_figure(pLoF_mis, test='burden', panel = T, filter_cat=T, overwrite = T)
# figure 
# 
# png(paste0(figure_path, 'final_2024/figure3_gene_list_', test,'_2_anno_v2.png'), height = 5, width = 5, units = 'in', res = 300)
# print(figure)
# dev.off()

# figure <- gene_list_pleiotropy_figure(no_combine, test='burden', panel = F, filter_cat=T, overwrite = T)
# figure 
# 
# png(paste0(figure_path, 'final_2024/figure3_gene_list_', test,'_3_anno_v1.png'), height = 4, width = 10, units = 'in', res = 300)
# print(figure)
# dev.off()

figure_239_full <- gene_list_pleiotropy_figure(no_combine_239, test='burden_239', panel = F, filter_cat=F, overwrite = T) + 
  theme(plot.margin = unit(c(1,1,0.5,0.2), "cm"))
figure_239_full
figure_599_full <- gene_list_pleiotropy_figure(no_combine_599, test='burden_599', panel = F, filter_cat=F, overwrite = T)+ 
  theme(plot.margin = unit(c(1,1,0.5,0.2), "cm"))
figure_599_full

figure <- ggpubr::ggarrange(figure_599_full, figure_239_full, 
                            labels = c('(A) Proportion of pleiotropic genes among 599 high-quality phenotypes', '(B) Proportion of pleiotropic genes among 239 independent phenotypes'),
                            nrow=2, common.legend=TRUE, vjust = 1.5, hjust = -0.1, font.label = list(size = 10, color = "black", face = "bold", family = NULL), legend = 'bottom', heights = c(0.2, 0.18))
png(paste0(figure_path,'figureS3.png'), height = 8, width = 8, units = 'in', res = 300)
print(figure)
dev.off()


# SuppTable 1
supptable_239 <- read_csv(paste0(result_path, 'gene_list_pleiotropy_burden_239.csv'))
colnames(supptable_239)[3] <- 'n_with_associations'
supptable_599 <- read_csv(paste0(result_path, 'gene_list_pleiotropy_burden_599.csv'))
colnames(supptable_599)[3] <- 'n_with_associations'
gene_list_names['universe'] <- 'Universe'
supptable <- merge(supptable_239, supptable_599, by = c('gene_list', 'annotation', 'n_total'), suffixes = c('.239', '.599'), all=T) %>%
  mutate(gene_list = gene_list_names[gene_list],
         annotation = factor(annotation, levels=annotation_types)) %>%
  arrange(gene_list, annotation)
write_csv(supptable, paste0(result_path, 'gene_list_pleiotropy_burden.csv'))



