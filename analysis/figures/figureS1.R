source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')
require("reticulate")
source_python("~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/python/pickle.py")

plof_1 <- read_pickle_file(paste0(data_path,"pvalue_permutation/gene_plof_(0.01%, 0.1%]_proportion"))
plof_2 <- read_pickle_file(paste0(data_path,"pvalue_permutation/gene_plof_(0.1%, 1%]_proportion"))
plof_3 <- read_pickle_file(paste0(data_path,"pvalue_permutation/gene_plof_(1%, 10%]_proportion"))
plof_4 <- read_pickle_file(paste0(data_path,"pvalue_permutation/gene_plof_(10%, ∞)_proportion"))
mis_1 <- read_pickle_file(paste0(data_path,"pvalue_permutation/gene_missense_lc_(0.01%, 0.1%]_proportion"))
mis_2 <- read_pickle_file(paste0(data_path,"pvalue_permutation/gene_missense_lc_(0.1%, 1%]_proportion"))
mis_3 <- read_pickle_file(paste0(data_path,"pvalue_permutation/gene_missense_lc_(1%, 10%]_proportion"))
mis_4 <- read_pickle_file(paste0(data_path,"pvalue_permutation/gene_missense_lc_(10%, ∞)_proportion"))
syn_1 <- read_pickle_file(paste0(data_path,"pvalue_permutation/gene_synonymous_(0.01%, 0.1%]_proportion"))
syn_2 <- read_pickle_file(paste0(data_path,"pvalue_permutation/gene_synonymous_(0.1%, 1%]_proportion"))
syn_3 <- read_pickle_file(paste0(data_path,"pvalue_permutation/gene_synonymous_(1%, 10%]_proportion"))
syn_4 <- read_pickle_file(paste0(data_path,"pvalue_permutation/gene_synonymous_(10%, ∞)_proportion"))
# combined_1 <- read_pickle_file(paste0(data_path,"gene_plof_missense_lc_(0.01%, 0.1%]_proportion"))
# combined_2 <- read_pickle_file(paste0(data_path,"gene_plof_missense_lc_(0.1%, 1%]_proportion"))
# combined_3 <- read_pickle_file(paste0(data_path,"gene_plof_missense_lc_(1%, 10%]_proportion"))
# combined_4 <- read_pickle_file(paste0(data_path,"gene_plof_missense_lc_(10%, ∞)_proportion"))

true <-  read.csv(paste0(data_path, 'gene_phewas_burden_sig_count_239.csv'), sep = '\t')
true_sum <- true %>%
  filter(annotation != 'pLoF|missense|LC') %>%
  mutate(interval = get_freq_interval(CAF)) %>%
  group_by(annotation, interval) %>%
  summarize(p_sig = sum(n_phewas_sig > 0)/n(),
            p_pleiotropy = sum(n_phewas_sig > 1)/n(),
            p_pleiotropy_sig = sum(n_phewas_sig > 1)/sum(n_phewas_sig > 0))%>%
  mutate(annotation = factor(annotation, levels = annotation_types2),
         interval = factor(interval, levels = c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'),
                           labels = c('(0.01%, 0.1%]', '(0.1%, 1%]', '(1%, 10%]', paste0('(10%, ', bquote("\U221E"), ' )'))))

permute_data <- data.frame(prop_pleiotropy = c(unlist(plof_1), unlist(plof_2), unlist(plof_3), unlist(plof_4), unlist(mis_1), unlist(mis_2), unlist(mis_3), unlist(mis_4), unlist(syn_1), unlist(syn_2), unlist(syn_3), unlist(syn_4)), 
                           annotation = rep(c('pLoF', 'missense|LC', 'synonymous'), each =400), interval =rep(rep(c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'), each=100), times=3))
mean_permute <- permute_data %>%
  group_by(annotation, interval) %>%
  summarize(mean_permute = mean(prop_pleiotropy))%>%
  mutate(annotation = factor(annotation, levels = annotation_types2),
         interval = factor(interval, levels = c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'),
                           labels = c('(0.01%, 0.1%]', '(0.1%, 1%]', '(1%, 10%]', paste0('(10%, ', bquote("\U221E"), ' )'))))

label_type = labeller(annotation = annotation_names2)
figureB <- permute_data %>%
  mutate(annotation = factor(annotation, levels = annotation_types2),
         interval = factor(interval, levels = c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'),
                           labels = c('(0.01%, 0.1%]', '(0.1%, 1%]', '(1%, 10%]', paste0('(10%, ', bquote("\U221E"), ' )'))))%>%
  ggplot + aes(x = prop_pleiotropy, color = annotation, fill = annotation) + 
  labs(y = 'Count', x = 'Proportion of pleiotropic genes') +
  geom_histogram(alpha = 0.5, bin=100) + 
  geom_vline(data = true_sum, aes(xintercept = p_pleiotropy, color = annotation), lty=2, show.legend = F) + 
  scale_x_continuous(labels = scales::percent) +
  annotation_color_scale2 + annotation_fill_scale2 + themes + theme_classic() +
  facet_grid(interval~annotation, scale = 'free', labeller = label_type) + 
  geom_text(data = true_sum, aes(x = Inf, y = Inf, label = paste0('True: ',round(p_pleiotropy*100, 2), '%'), color =annotation), vjust = 1, hjust = 1, show.legend = F) + 
  geom_text(data = mean_permute, aes(x = Inf, y = Inf, label = paste0('Mean permuted:',round(mean_permute*100, 2), '%'), color =annotation), vjust = 3, hjust = 1, show.legend = F) + 
  theme(legend.position = 'top') +
  theme(plot.margin = unit(c(0.5,1,0.5,0.2), "cm"))
figureB

plof <- read_pickle_file(paste0(data_path, "pvalue_permutation/gene_plof_proportion"))
mis <- read_pickle_file(paste0(data_path,"pvalue_permutation/gene_missense_lc_proportion"))
syn <- read_pickle_file(paste0(data_path,"pvalue_permutation/gene_synonymous_proportion"))

true_sum <- gene_data %>%
  group_by(annotation) %>%
  dplyr::summarize(p_sig = sum(n_phewas_sig > 0)/n(),
                   p_pleiotropy = sum(n_phewas_sig > 1)/n(),
                   p_pleiotropy_sig = sum(n_phewas_sig > 1)/sum(n_phewas_sig > 0))%>%
  mutate(annotation = factor(annotation, levels = annotation_types))

permute_data <- data.frame(prop_pleiotropy = c(unlist(plof), unlist(mis), unlist(syn)), annotation = rep(c('pLoF', 'missense|LC', 'synonymous'), each =100))
mean_permute <- permute_data %>%
  group_by(annotation) %>%
  summarize(mean_permute = mean(prop_pleiotropy))%>%
  mutate(annotation = factor(annotation, levels = annotation_types))


label_type = labeller(annotation = annotation_names)
figureA <- permute_data %>%
  mutate(annotation = factor(annotation, levels = annotation_types)) %>%
  ggplot + aes(x = prop_pleiotropy, color = annotation, fill = annotation) +
  labs(y = 'Count', x = 'Proportion of pleiotropic genes') +
  geom_histogram(alpha = 0.5, bin=100) +
  geom_vline(data = true_sum, aes(xintercept = p_pleiotropy, color = annotation), lty=2) +
  scale_x_continuous(labels = scales::percent) +
  annotation_color_scale + annotation_fill_scale + themes + theme_classic() +
  facet_grid(~annotation, scale = 'free', labeller = label_type) +
  geom_text(data = true_sum, aes(x = Inf, y = Inf, label = paste0('True: ',round(p_pleiotropy*100, 2), '%'), color=annotation), vjust = 1, hjust = 1, show.legend = F) +
  geom_text(data = mean_permute, aes(x = Inf, y = Inf, label = paste0('Mean permuted:',round(mean_permute*100, 2), '%'), color = annotation), vjust = 3, hjust = 1, show.legend = F) +
  guides(text = 'None') +
  theme(plot.margin = unit(c(0.5,1,0.5,0.2), "cm"))

figureA

figure = ggpubr::ggarrange(figureA, figureB, nrow=2, heights = c(0.1, 0.25), labels = c('(A) Permutation across all CAF', '(B) Permutation within CAF bins'),
                           hjust = 0, common.legend = T,
                           font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')
)


png(paste0(figure_path,'figureS1.png'), height = 8, width = 8, units = 'in', res = 300)
print(figure)
dev.off()
