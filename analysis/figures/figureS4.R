source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("STRINGdb", force = T)
# Requirements
# packages = c('gplots', 'hash', 'sqldf', 'plotrix', 'tidygraph', "DBI", "STRINGdb", "magrittr")
# 
# for(p in packages){
#   if(!require(p, character.only = T)){
#     install.packages(p)
#   }
# }

data_599 <- read.csv(paste0(data_path, 'gene_phewas_burden_sig_count_599.csv'), sep = '\t') %>%
  mutate(interval = get_freq_interval(freq=CAF))  %>%
  mutate(interval = factor(interval, levels = c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'),
                           labels = c('(0.01%, 0.1%]', '(0.1%, 1%]', '(1%, 10%]', paste0('(10%, ', bquote("\U221E"), ' )') )))

ppi_599 <- ppi_figure(data_599, name = 'figureS4/figureS4_ppi_599', save = T)


data_239 <- read.csv(paste0(data_path, 'gene_phewas_burden_sig_count_239.csv'), sep = '\t') %>%
  mutate(interval = get_freq_interval(freq=CAF))  %>%
  mutate(interval = factor(interval, levels = c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'),
                           labels = c('(0.01%, 0.1%]', '(0.1%, 1%]', '(1%, 10%]', paste0('(10%, ', bquote("\U221E"), ' )') )))

ppi_239 <- ppi_figure(data_239, name = 'figureS4/figure2_ppi_239', save = T)

figure <- ggpubr::ggarrange(ppi_599 +
                              theme(legend.title = element_text(face = 'plain', size = 11), 
                                    axis.title = element_text(face = 'plain', size = 11), 
                                    plot.margin = unit(c(1,0,0,0.5), "cm")), 
                            ppi_239+
                              theme(legend.title = element_text(face = 'plain', size = 11), 
                                    axis.title = element_text(face = 'plain', size = 11), 
                                    plot.margin = unit(c(0.7,0,0,0.5), "cm")), 
                            labels = c('(A) Number of protein-protein interactions across genes among 599 independent phenotypes', 
                                       '(B) Number of protein-protein interactions across genes among 239 independent phenotypes'),
                            nrow=2, vjust = 2, hjust = 0, font.label = list(size = 10, color = "black", face = "bold", family = NULL), 
                            common.legend=TRUE, heights = c(0.2, 0.18))
png(paste0(figure_path,'figureS4.png'), height = 6, width = 8, units = 'in', res = 300)
print(figure)
dev.off()