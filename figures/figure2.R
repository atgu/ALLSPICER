source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')
library(emojifont)

raw_results_500k <- read_pleiotropy_results('burden', '500k') 
results_500k <- modify_results_table(raw_results_500k, 'burden', '500k')

top_triplets <- read_delim(paste0(data_path, "top_significant_triplets.txt.bgz"), delim='\t', col_types = cols(phenocode = col_character())) %>% 
  mutate(coding = if_else(is.na(coding), '', coding)) %>%
  mutate(phenoname = paste0(trait_type, '_', phenocode, '_',  pheno_sex, '_',  coding, '_',  modifier)) %>%
  mutate(BETA_adjusted = sqrt(2*AF*(1-AF))*BETA)

## ALB panel
alb_info <- read_delim(paste0(data_path, 'alb_albumin_calcium_var_hgvsp_500k.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character())) %>% 
  mutate(coding = if_else(is.na(coding), '', coding)) %>%
  mutate(phenoname = paste0(trait_type, '_', phenocode, '_',  pheno_sex, '_',  coding, '_',  modifier),
         AF = pop_AF) 
gene_lab <- 'ALB'
pheno1_lab <- 'continuous_30600_both_sexes__irnt'
pheno2_lab <- 'continuous_30680_both_sexes__irnt'
c_hat <- results_500k %>%
  filter(gene == gene_lab & pheno1 == pheno1_lab & pheno2 == pheno2_lab) %>%
  select(annotation, c_hat)

p1<- figure_beta_triplets_all_annt(alb_info  %>% filter(AF < 1e-4), 
                                   gene = gene_lab, 
                                   phenocode1 = pheno1_lab, 
                                   phenocode2 = pheno2_lab, 
                                   c_hat = c_hat, 
                                   pheno1_name = 'Albumin', 
                                   pheno2_name = 'Calcium', 'figure4', 
                                   save=F, 
                                   threshold = 0.05) + 
  theme(legend.title = element_text(face = 'plain', size = 11), 
        # legend.position = 'none',
        axis.title = element_text(face = 'plain', size = 11), 
        plot.margin = unit(c(0.5,0,0,0), "cm")) +
  geom_point(data = alb_info %>% filter(locus == 'chr4:73412072') %>%
               mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
               pivot_wider(names_from = phenocode,  values_from = c('BETA', 'Pvalue'), id_cols = c('locus', 'alleles','AF', 'gene', 'annotation')) %>%
               mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation))%>%
               mutate(annotation = factor(annotation, levels=annotation_types)), aes(x = BETA_30600, y=BETA_30680), size=5, pch=1, color='#6CA6CB') +
  geom_text_repel(data = alb_info %>% filter(locus == 'chr4:73412072') %>%
                    mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
                    pivot_wider(names_from = phenocode,  values_from = c('BETA', 'Pvalue'), id_cols = c('locus', 'alleles','AF', 'gene', 'annotation')) %>%
                    mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation))%>%
                    mutate(annotation = factor(annotation, levels=annotation_types)), aes(x = BETA_30600, y=BETA_30680, label = 'chr4:73412072:A:G'), size=3, color='#6CA6CB', face = 'bold', vjust = 1, hjust = 0.5) +
  geom_point(data = alb_info %>% filter(locus == 'chr4:73405164') %>%
               mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
               pivot_wider(names_from = phenocode,  values_from = c('BETA', 'Pvalue'), id_cols = c('locus', 'alleles','AF', 'gene', 'annotation')) %>%
               mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation))%>%
               mutate(annotation = factor(annotation, levels=annotation_types)), aes(x = BETA_30600, y=BETA_30680), size=5, pch=1, color='#6CA6CB') +
  geom_text_repel(data = alb_info %>% filter(locus == 'chr4:73405164') %>%
                    mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
                    pivot_wider(names_from = phenocode,  values_from = c('BETA', 'Pvalue'), id_cols = c('locus', 'alleles','AF', 'gene', 'annotation')) %>%
                    mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation))%>%
                    mutate(annotation = factor(annotation, levels=annotation_types)), aes(x = BETA_30600, y=BETA_30680, label = 'chr4:73405164:T:C'), size=3, color='#6CA6CB', face = 'bold', vjust = 1, hjust = 0.5)


p1

## ALPL panel
alpl_info <- read_delim(paste0(data_path, 'alpl_var_hgvsp_500k.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character())) %>% 
  mutate(coding = if_else(is.na(coding), '', coding)) %>%
  mutate(phenoname = paste0(trait_type, '_', phenocode, '_',  pheno_sex, '_',  coding, '_',  modifier))

gene_lab <- 'ALPL'
pheno1_lab <- 'continuous_30610_both_sexes__irnt'
pheno2_lab <- 'continuous_30810_both_sexes__irnt'
c_hat <- results_500k %>%
  filter(gene == gene_lab & pheno1 == pheno1_lab & pheno2 == pheno2_lab) %>%
  select(annotation, c_hat)

p2 <- figure_beta_triplets_all_annt(alpl_info  %>% filter(AF < 1e-4), 
                                    gene_name = 'ALPL', 
                                    phenocode1 = 'continuous_30610_both_sexes__irnt', 
                                    phenocode2 = 'continuous_30810_both_sexes__irnt',
                                    c_hat = c_hat,
                                    pheno1_name = 'Alkaline phosphatase', 
                                    pheno2_name = 'Phosphate', 
                                    'figureS14', save=F, 
                                    threshold = 0.05) + 
  theme(legend.title = element_text(face = 'plain', size = 11), 
        axis.title = element_text(face = 'plain', size = 11), 
        plot.margin = unit(c(0.5,0,0,0), "cm"))
p2


figure = ggpubr::ggarrange(p1, p2, nrow=2, heights = c(0.1, 0.1), labels = c('(A) Comparison of effect sizes of ALB variants on Albumin and Calcium', '(B) Comparison of effect sizes of ALPL variants on Alkaline phosphatase and Phosphate', ''), hjust = 0, 
                           font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')
)
png(paste0(figure_path,'figure2.png'), height = 5, width = 7.5, units = 'in', res = 300)
print(figure)
dev.off()
