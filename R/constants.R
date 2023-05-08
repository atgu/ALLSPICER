source('~/ukb_exomes/R/constants.R')
detach('package:plyr')
library(UpSetR)
library(ComplexUpset)
library(ggsci)

root <- '~/Dropbox (Partners HealthCare)/analysis/ukb_exomes_pleiotropy/'
data_path <- paste0(root, 'data/')
figure_path <- paste0(root, 'figures/')
result_path <- paste0(root, 'results/')
r_path <- paste0(root, 'R/')

annotation_types2 = c(annotation_types, 'pLoF|missense|LC')
annotation_names2 = c(annotation_names, 'pLoF|Missense')
names(annotation_names2) = annotation_types2
colors2 = c(colors, 'pLoF|missense|LC' = '#FFA600')
annotation_fill_scale2 = scale_fill_manual(name = 'Annotation', values = colors2, breaks = annotation_types2, labels = annotation_names2)
annotation_color_scale2 = scale_color_manual(name = 'Annotation', values = colors2, breaks = annotation_types2, labels = annotation_names2)

## pleiotropic genes/variant count in bins figure
pleiotropic_cnt_bin_figure <- function(data_type='gene', test, save_plot=T){
  data <- load_ukb_file(paste0(data_type, '_sig_cnt_filtered_', test, '_', tranche,'.txt.bgz'), subfolder = 'analysis/')
  annotation_types <- c(annotation_types, 'pLoF|missense|LC')
  annotation_names <- c(annotation_names, 'pLoF|Missense')
  names(annotation_names) <- annotation_types
  label_type = labeller(annotation = annotation_names)
  sum_data <- data %>%
    mutate(interval = get_freq_interval(get(if_else(data_type == 'gene', 'CAF', 'AF')))) %>%
    group_by(interval, annotation) %>% add_count() %>%
    filter(all_sig_pheno_cnt>1) %>%
    mutate(cnt_interval = case_when(
      all_sig_pheno_cnt <= 5  ~ '[2, 5]',
      all_sig_pheno_cnt <= 10  ~ '(5, 10]',
      all_sig_pheno_cnt <= 15  ~ '(10, 15]',
      all_sig_pheno_cnt <= 20  ~ '(15, 20]',
      all_sig_pheno_cnt <= 25  ~ '(20, 25]',
      all_sig_pheno_cnt <= 30  ~ '(25, 30]',
      all_sig_pheno_cnt <= 40  ~ '(30, 40]',
      all_sig_pheno_cnt <= 60  ~ '(40, 60]',
      all_sig_pheno_cnt > 60  ~  paste0('(60, ', bquote("\U221E"), ' )'),
    )) %>%
    group_by(interval, annotation, cnt_interval, n) %>% summarise(cnt = n()) %>%
    mutate(annotation = factor(annotation, levels = annotation_types2),
           interval = factor(interval, levels = rev(c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]')),
                             labels = rev(c('(0.01%, 0.1%]', '(0.1%, 1%]', '(1%, 10%]', paste0('(10%, ', bquote("\U221E"), ' )')))),
           cnt_interval = factor(cnt_interval, levels = c('[2, 5]', '(5, 10]', '(10, 15]', '(15, 20]', '(20, 25]', '(25, 30]', '(30, 40]', '(40, 60]', paste0('(60, ', bquote("\U221E"), ' )')),
                                 labels = c('[2, 5]', '(5, 10]', '(10, 15]', '(15, 20]', '(20, 25]', '(25, 30]', '(30, 40]', '(40, 60]', paste0('(60, ', bquote("\U221E"), ' )')))) 
  label_data <- sum_data %>%
    group_by(annotation, interval) %>%
    summarise(n_sig = sum(cnt))
  
  sum_data <- sum_data %>%
    merge(., label_data, by = c("annotation", 'interval')) %>%
    mutate(prop_sig = cnt/n_sig)
  
  cnts <- c('[2, 5]', '(5, 10]', '(10, 15]', '(15, 20]', '(20, 25]', '(25, 30]', '(30, 40]', '(40, 60]', '(60, ∞ )')
  cnt_colors <- inlmisc::GetColors(length(cnts), scheme = "smooth rainbow")
  names(cnt_colors) <- cnts
  test_name <- if_else(data_type=='gene', if_else(test=='skato', '(SKAT-O)', '(Burden test)'), '')
  if(data_type=='gene'){
    y_labels = c(200, 150, 100, 50, 0, 0.25, 0.5, 0.75, 1)
  }else{
    y_labels = c(4000, 3000, 2000, 1000, 0, 0.25, 0.5, 0.75, 1)
  }
  figure <- sum_data %>%
    ggplot + aes(x = interval,
                 y = prop_sig*2,
                 color = reorder(factor(cnt_interval), desc(cnt_interval)),
                 fill = reorder(factor(cnt_interval), desc(cnt_interval)),
                 label = if_else(prop_sig < 0.004, '', as.character(cnt))) +
    geom_bar(position='stack', stat='identity')+
    labs(y = paste('Number/Proportion of pleiotropy', if_else(data_type=='gene', 'genes', 'variants'), test_name), x = NULL, 
         color = paste('Number of Associations'), fill = paste('Number of Associations'))  +
    scale_color_manual(values = cnt_colors) + 
    scale_fill_manual(values = cnt_colors) + 
    scale_y_continuous(limits = c(if_else(data_type=='gene', if_else(test == 'burden', -0.35, -0.52), -0.4), 2.001), 
                       expand = c(0, 0), 
                       breaks = c(-0.4, -0.3, -0.2, -0.1, 0, 0.5, 1, 1.5, 2), 
                       labels =y_labels) +
    facet_grid(annotation~., switch = "y", scales='free_x', labeller = label_type) + 
    coord_flip()+    theme(panel.spacing = unit(0, "lines"),
                           panel.grid = element_blank(),
                           panel.border = element_blank(),
                           strip.placement = "outside",
                           strip.text.y = element_text(angle = 0),
                           # strip.text.y.left = element_text(angle = 0)
                           strip.text = element_text(face = 'bold', size = 22),
                           axis.text= element_text(size = 20, face = 'bold'),
                           legend.title = element_text(face = 'bold', size = 25),
                           legend.text = element_text(size = 18),
                           axis.text.x = element_text(vjust = 1, hjust = 0.95),
                           axis.text.y = element_blank(),
                           axis.title = element_text(size = 25, face = 'bold'),
                           axis.ticks.y = element_blank(), 
                           axis.line.y= element_blank(),
                           legend.position="top", legend.box="horizontal", legend.margin=margin()
    ) +
    geom_text(size = 7, position = position_stack(vjust = 0.5), color = 'white') + 
    geom_text(aes(y = -0.25, label=paste0(n_sig, ' /', n), color=NULL), size = 8) +
    geom_label(data = label_data, aes(y = 0, label=interval, fill = NULL,color = NULL), size = 8, fontface='bold', alpha = 0) +
    geom_bar(data = label_data, aes(x = interval, y = -n_sig/if_else(data_type=='gene', 500,10000), label=NULL, fill=NULL, color=NULL), stat='identity', alpha = .5) +
    guides(fill=guide_legend(nrow=1,byrow=TRUE),
           color=guide_legend(nrow=1,byrow=TRUE)) 
  if(save_plot){
    png(paste0(figure_path, data_type, '_cnt_all_pheno_', test,'_filtered_by_interval_bar_binned.png'), height =12, width = 20, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

## pleiotropic gene count figure
pleiotropic_gene_cnt_figure <- function(data_type, test, save_plot=T){
  data <- load_ukb_file(paste0(data_type,'_sig_cnt_filtered_', test, '_', tranche,'.txt.bgz'), subfolder = 'analysis/')
  annotation_types2 <- c(annotation_types, 'pLoF|missense|LC')
  annotation_names2 <- c(annotation_names, 'pLoF|Missense')
  names(annotation_names2) <- annotation_types2
  label_type = labeller(annotation = annotation_names2)
  sum_data <- data %>%
    # filter(annotation != 'pLoF|missense|LC') %>%
    mutate(interval = get_freq_interval(get(if_else(data_type == 'gene', 'CAF', 'AF')))) %>%
    group_by(interval, annotation) %>% add_count() %>%
    group_by(interval, annotation, all_sig_pheno_cnt, n) %>% summarise(cnt = n()) %>%
    mutate(prop = cnt/n) %>%
    mutate(sd = sqrt(prop*(1-prop))/n,
           annotation = factor(annotation, levels = annotation_types2),
           interval = factor(interval, levels = rev(c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]')),
                             labels = rev(c('(0.01%, 0.1%]', '(0.1%, 1%]', '(1%, 10%]', if_else(data_type == 'gene', paste0('(10%, ', bquote("\U221E"), ' )'), '(0.1, 1]' ))))) %>%
    filter(all_sig_pheno_cnt>1)
  label_data <- sum_data %>%
    group_by(annotation, interval) %>%
    summarise(n_sig = sum(cnt)) 
  
  sum_data <- sum_data %>%
    merge(., label_data, by = c("annotation", 'interval')) %>%
    mutate(prop_sig = cnt/n_sig)
  
  cnts <- sort(unique(sum_data$all_sig_pheno_cnt))
  cnt_colors <- inlmisc::GetColors(if_else(data_type == 'gene', 33, 87), scheme = "smooth rainbow")
  names(cnt_colors) <- as.character(2:(if_else(data_type == 'gene', 33, 87)+1))
  cnt_colors <- cnt_colors[as.character(cnts)]
  test_name <- if_else(test=='skato', 'SKAT-O', 'Burden test')
  figure <- sum_data %>%
    ggplot + 
    aes(x = interval,
        y = cnt/n_sig,
        color = reorder(factor(all_sig_pheno_cnt), -all_sig_pheno_cnt),
        fill = reorder(factor(all_sig_pheno_cnt), -all_sig_pheno_cnt),
        label = if_else(cnt == 1, '', as.character(cnt))) +
    geom_bar(position='stack', stat='identity') +
    # geom_hline(yintercept = c(0)) + 
    labs(y = paste('Number/Proportion of pleiotropy genes (', test_name,')'), x = NULL, 
         color = paste('Number of Associations'), fill = paste('Number of Associations'))  +
    scale_color_manual(values = cnt_colors) + 
    scale_fill_manual(values = cnt_colors) + 
    scale_y_continuous(limits = c(if_else(test == 'burden', -0.35, -0.52), 1.001), expand = c(0, 0), breaks = c(-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.25, 0.5, 0.75, 1), labels = c(250, 200, 150, 100, 50, 0, 0.25, 0.5, 0.75, 1)) +
    facet_grid(annotation~., switch = "y", scales='free_x', labeller = label_type) + 
    # theme_classic() +
    coord_flip()+
    theme(panel.spacing = unit(0, "lines"),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          strip.placement = "outside",
          strip.text.y = element_text(angle = 0),
          # strip.text.y.left = element_text(angle = 0)
          strip.text = element_text(face = 'bold', size = 22),
          axis.text= element_text(size = 20, face = 'bold'),
          legend.title = element_text(face = 'bold', size = 25),
          legend.text = element_text(size = 18),
          axis.text.x = element_text(vjust = 1, hjust = 0.95),
          axis.text.y = element_blank(),
          axis.title = element_text(size = 25, face = 'bold'),
          axis.ticks.y = element_blank(), 
          axis.line.y= element_blank(),
          legend.position="right", legend.box="horizontal", legend.margin=margin()
    ) +
    geom_text(size = 8, position = position_stack(vjust = 0.5), color = 'white') + 
    geom_text(aes(y = if_else(test =='burden', -0.25, -0.45), label=paste0(n_sig, ' /', n)), size = 8, fontface='bold') +
    geom_label(data = label_data, aes(y = 0, label=interval, color = NULL), size = 8, fontface='bold', fill = 'white', alpha = 0) +
    geom_bar(data = label_data, aes(x = interval, y = -n_sig/500, label=NULL, fill=NULL, color=NULL), stat='identity', alpha = .5) +
    guides(fill=guide_legend(nrow=5,byrow=TRUE),
           color=guide_legend(nrow=5,byrow=TRUE)) 
  if(save_plot){
    png(paste0(figure_path, data_type, '_cnt_all_pheno_', test,'_filtered_by_interval_bar_normalized.png'), height =12, width = 22, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

## gene list figure
get_gene_list_sumstats <- function(gene_list, pleiotropic_genes, list_name){
  gene_data = gene_list %>%
    mutate(n_total = n()) %>%
    merge(., pleiotropic_genes, by = 'gene_symbol', all.x = T)  %>%
    add_count(annotation)
  sum_table = gene_data %>% 
    filter(!is.na(annotation)) %>%
    group_by(annotation, n_total, n) %>%
    summarise(n_pleiotropy = sum(sig_cnt > 1),
              # sum_prop_con_sig = sum(prop_con_sig * (sig_cnt > 1), na.rm = T),
              # sum_prop_icd_sig = sum(prop_icd_sig * (sig_cnt > 1), na.rm = T),
              # sum_prop_cat_sig = sum(prop_cat_sig * (sig_cnt > 1), na.rm = T),
              ) %>% 
    mutate(prop_pleiotropy = n_pleiotropy/n,
           # mean_prop_con_sig = sum_prop_con_sig/n_pleiotropy, 
           # mean_prop_icd_sig = sum_prop_icd_sig/n_pleiotropy, 
           # mean_prop_cat_sig = sum_prop_cat_sig/n_pleiotropy, 
           gene_list = list_name,
           sd = sqrt(prop_pleiotropy*(1-prop_pleiotropy)/n)) 
  # %>%
  #   mutate(mean_prop_con_sig = if_else(is.nan(mean_prop_con_sig), NA, mean_prop_con_sig),
  #          mean_prop_icd_sig = if_else(is.nan(mean_prop_icd_sig), NA, mean_prop_icd_sig),
  #          mean_prop_cat_sig = if_else(is.nan(mean_prop_cat_sig), NA, mean_prop_cat_sig))
  return(as.data.frame(sum_table))
}

gene_list_pleiotropy_figure <- function(test, fig_type, save_plot=T){
  gene_sig_after <- load_ukb_file(paste0('gene_sig_cnt_filtered_', test, '_', tranche,'.txt.bgz'), subfolder = 'analysis/')
  gene_info <- load_ukb_file('gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', subfolder = 'analysis/')
  gene_dd <- load_ukb_file('forKonrad_sig31kDNM_consensus_genes_2021_01_12.txt', subfolder = 'analysis/') %>% select(gene_symbol = symbol) 
  
  pleiotropic_genes <- gene_sig_after %>%
    distinct(gene_id, 
             gene_symbol, 
             annotation, 
             sig_cnt = all_sig_pheno_cnt, 
             con_sig_cnt = get(paste0('continuous_sig_pheno_cnt_', test)),
             icd_sig_cnt = get(paste0('icd10_sig_pheno_cnt_', test)),
             cat_sig_cnt = get(paste0('categorical_sig_pheno_cnt_', test)),) %>%
    mutate(prop_con_sig = con_sig_cnt/sig_cnt,
           prop_icd_sig = icd_sig_cnt/sig_cnt,
           prop_cat_sig = cat_sig_cnt/sig_cnt)
  constrained = gene_sig_after %>%
    merge(., gene_info[, c('gene_id', 'gene', 'oe_lof_upper_bin')], by.y = c('gene_id', 'gene'), by.x =c('gene_id', 'gene_symbol')) %>%
    filter(oe_lof_upper_bin == 0) %>% distinct(gene_symbol)
  
  files <- list.files('~/gene_lists/lists/')
  files <- files[grepl(".tsv$",files)]
  strs <- unlist(strsplit(files,".tsv"))
  results <- get_gene_list_sumstats(constrained, pleiotropic_genes, 'constrained')
  results <- rbind(results, get_gene_list_sumstats(gene_dd, pleiotropic_genes, 'developmental delay'))
  for(i in 1: length(strs)){
    gene_list =  as.data.frame(fread(paste0('~/gene_lists/lists/', files[i]), quote="", header = T))
    colnames(gene_list)[1] = 'gene_symbol'
    sum_table <- get_gene_list_sumstats(gene_list, pleiotropic_genes, strs[i])
    results <- rbind(results, data.frame(sum_table))
  }
  annotation_types2 = c(annotation_types, 'pLoF|missense|LC')
  annotation_names2 = c(annotation_names, 'pLoF|Missense')
  names(annotation_names2) = annotation_types2
  results <- results %>% 
    mutate(annotation  = factor(annotation, levels = annotation_types2))
  write_csv(results, paste0(result_path, 'gene_list_pleiotropy_', test, '.csv'))
  test_name <- if_else(test == 'skato', 'SKAT-O', 'Burden test')

  colors2 = c(colors, 'pLoF|missense|LC' = '#FFA600')
  annotation_fill_scale2 = scale_fill_manual(name = 'Annotation', values = colors2, breaks = annotation_types2, labels = annotation_names2)
  annotation_color_scale2 = scale_color_manual(name = 'Annotation', values = colors2, breaks = annotation_types2, labels = annotation_names2)
  gene_list_names = sort(c('ACMG_2_0' = 'ACMG V2.0', 'all_ad' = 'Autosomal dominant', 'all_ar' = 'Autosomal recessive', 'berg_ad' = 'Autosomal dominant-Berg (OMIM)', 'berg_ar' = 'Autosomal recessive-Berg (OMIM)', 'berg_xd' = 'X-linked dominant genes-Berg (OMIM)', 'berg_xr' = 'X-linked recessive genes-Berg (OMIM)',
                      'blekhman_ad' = 'Autosomal dominant-Blekhman (OMIM)', 'blekhman_ar' = 'Autosomal recessive-Blekhman (OMIM)', 'blekhman_x' = 'X-linked genes-Blekhman (OMIM)', 'BROCA_Cancer_Risk_Panel' = 'BROCA-Cancer risk panel', 'clingen_level3_genes_2015_02_27' = 'ClinGen haploinsufficient V1', 
                      'clingen_level3_genes_2018_09_13' = 'ClinGen haploinsufficient V2', 'clinvar_path_likelypath' = 'Clinvar (likely) pathogenic', 'constrained' = 'Constrained', 'core_essentials_hart' = 'Essential in culture', 'developmental delay' = 'Developmental Delay',  
                      'DRG_KangJ' = 'DNA repair genes-KangJ', 'DRG_WoodRD' = 'DNA repair genes-WoodRD', 'drug_targets_nelson' = 'Drug targets', 'fda_approved_drug_targets' = 'FDA approved drug targets', 'fmrp_list_gencode' = 'FMRP interactors', 'gpcr_guide' = 'GPCRs from guidetopharmacology', 'gpcr_union' = 'GPCRs all', 
                      'gpcr_uniprot' = 'GPCRs from uniprot', 'gpi_anchored' = 'GPI-anchored proteins', 'gwascatalog' = 'GWAS catalog', 'CEGv2_subset_universe' = 'Cell essential', 'NEGv1_subset_universe' = 'Cell non-essential', 'haploinsufficiency_mild_curated_2016' = 'Haploinsufficient mild',
                      'haploinsufficiency_moderate_curated_2016' = 'Haploinsufficient moderate', 'haploinsufficiency_severe_curated_2016' = 'Haploinsufficient severe', 'homozygous_lof_tolerant_twohit' = 'Homozygous lof tolerant', 'kinases' = 'Kinases', 'mgi_essential' = 'MGI essential', 
                      'natural_product_targets' = 'Natural product targets', 'olfactory_receptors' = 'Olfactory receptors', 'x-linked_clinvar' = 'X-linked genes-ClinVar'), decreasing = T)
  if(fig_type == 'prop_pleiotropy'){
    figure <- results %>%
      filter(gene_list != 'universe') %>%
      filter(n > 10 & n_pleiotropy >0) %>%
      ggplot +
      aes(x = reorder(gene_list, desc(gene_list)), y = prop_pleiotropy, ymin =  prop_pleiotropy-sd, ymax =  prop_pleiotropy+sd, color = annotation) +
      geom_pointrange(stat = "identity", position = position_dodge(width = 2), size=0.3) +
      geom_hline(data = results[results$gene_list == 'universe', ], aes(yintercept = prop_pleiotropy, color = annotation), lty = 2) +
      labs(x = NULL, y = paste0('Proportion of pleiotropic genes (', test_name,')')) + 
      scale_x_discrete(labels = gene_list_names) + 
      scale_y_continuous(labels = scales::percent_format(accuracy=1)) +
      annotation_color_scale2 + annotation_fill_scale2 + 
      coord_flip() +
      theme_classic() + themes +
      theme(axis.text.x = element_text(size = 8, angle = 0, vjust = 0.8),
            legend.position = 'none') + 
      facet_grid(~annotation, scales = 'fixed',labeller = labeller(annotation = annotation_names2)) 
  }else{
    figure <- results %>%
      filter(gene_list != 'universe') %>%
      add_count(gene_list, name = 'n_annt_group') %>%
      filter(n > 10 & n_pleiotropy >1 & n_annt_group>2) %>%
      ggplot +
      aes(x = reorder(gene_list, desc(gene_list)), y = 1 - mean_prop_con_sig, color = annotation) +
      # aes(x = gene_list, y = 1 - mean_prop_con_sig, color = annotation) +
      # aes(x =annotation, y = mean_prop_con_sig, color = annotation) +
      geom_point(stat = "identity", position = position_dodge(width = 2)) +
      geom_hline(data = results[results$gene_list == 'universe', ], aes(yintercept = 1 - mean_prop_con_sig, color = annotation), lty = 2) +
      geom_hline(aes(yintercept = 0, color = annotation)) +
      labs(x = NULL, y = paste0('Mean proportion of categorical associations among pleiotropic genes (', test_name,')')) + 
      scale_x_discrete(labels = gene_list_names) + 
      scale_y_continuous(labels = scales::percent_format(accuracy=1)) +
      annotation_color_scale2 + annotation_fill_scale2 + 
      coord_flip() +
      theme_classic() + themes +
      theme(axis.text.x = element_text(size = 7, angle = 0, vjust = 0.8),
            axis.title.x = element_text(size = 7, angle = 0, vjust = 0.8),
            legend.position = 'none') +
      # geom_text(data = results[results$gene_list != 'universe', ], aes(x = reorder(gene_list, desc(gene_list)), label = n_pleiotropy, color = annotation, y = 1-mean_prop_con_sig +0.02), size = 5) + 
      facet_grid(~annotation, scales = 'free_x',labeller = labeller(annotation = annotation_names2, gene_list = gene_list_names)) 

  }

  if(save_plot){
    png(paste0(figure_path, paste0('gene_list_', fig_type,'_', tranche, '_', test), '.png'), height = 5, width = 6.5, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

## pheno group figure
formatted_data <- function(group_type, test_type, write=T){
  data <- read_delim(paste0(data_path, group_type, '_domain_level_gene_sig_', test_type, '_', tranche,'.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character())) %>%
    filter(!is.na(pheno_group))
  data <- format_upset_data(data, group_type) %>%
    mutate(annotation = factor(annotation, levels = annotation_types2, labels = annotation_names2))
  if(write){
    write_csv(data, paste0(result_path, group_type, '_domain_level_gene_sig_', test_type, '_', tranche, '.csv'))
  }
  return(data)
}

format_upset_data <- function(data, group){
  if(group %in% c('icd', 'both')){
    data <- data %>%
      mutate(pheno_group = str_to_sentence(pheno_group)) %>%
      mutate(pheno_group_sig = if_else(pheno_group_sig, 1, 0), 
             group_name = if_else(pheno_group %in% names(icd_names), icd_names[pheno_group], pheno_group)) %>%
      as.data.frame() %>%
      pivot_wider(id_cols = colnames(.)[1:(ncol(.)-4)], names_from = 'group_name', values_from = 'pheno_group_sig', values_fn = function(x) sum(x)) 
    if(group == 'both'){
      data <- data %>%
        mutate(sum_disease = rowSums(data[,6:(ncol(data)-7)]),
               sum_biomarker = rowSums(data[,(ncol(data)-6):ncol(data)])) %>%
        filter(sum_disease != 0 & sum_biomarker!=0 )
    }

    
  }else{
    data <- data %>%
      filter(!is.na(pheno_group)) %>%
      mutate(pheno_group = str_to_sentence(pheno_group)) %>%
      mutate(pheno_group_sig = if_else(pheno_group_sig, 1, 0)
             ) %>%
      pivot_wider(id_cols = colnames(.)[1:(ncol(.)-3)], names_from = 'pheno_group', values_from = 'pheno_group_sig') %>%
      as.data.frame(.) 
  }
  data <- data %>%
    mutate(group_cnt = rowSums(data[,-c(1:5)])) %>%
    filter(group_cnt > 1) 
  print(data[,c('gene_symbol', 'annotation', 'group_cnt')])
  return(data %>% select(-group_cnt))
}

group_upset_plot <- function(data, group_type, test_type,groups, query.list, max_size, group_metadata = NULL){
  test_name = if_else(test_type == 'skato', 'SKAT-O', 'Burden test')
  if(group_type == 'icd'){
    p <- format_upset_data(data, group =group_type) %>%
      ComplexUpset::upset(., 
                          intersect = groups, name = NULL, 
                          mode = 'exclusive_intersection',
                          set_sizes = (upset_set_size() +
                                         geom_text(aes(label=..count..), hjust=-0.8, stat='count', size = 3) +
                                         theme(axis.text=element_text(size=8))
                          ),
                          guides='over',
                          matrix=(
                            intersection_matrix(
                              geom=geom_point(
                                shape='square',
                                size=2
                              ),
                              segment=geom_segment(
                                linetype=2,lwd = 0 
                              ))),
                          base_annotations=list(
                            ' '=intersection_size(
                              counts=FALSE,
                              mapping=aes(fill=annotation)
                            ) + annotation_fill_scale2+ 
                              scale_y_continuous(limits = c(0, max_size), expand = c(0, 0))+
                              theme_classic() + theme(axis.ticks.x = element_blank(),
                                                      axis.title.x = element_blank(),
                                                      axis.text.x = element_blank())
                          ),
                          width_ratio = 0.2, 
                          max_size = max_size, 
                          wrap = TRUE,
                          group_by = 'sets',
                          themes = upset_modify_themes(list(
                            'intersections_matrix'=theme(text=element_text(size=15)), 
                            'overall_sizes'=theme(axis.text.x=element_text(size=12, angle = 45))
                          )), 
                          queries=query.list
      )
  }
  if(group_type %in% c('biomarker', 'blood')){
    p <- format_upset_data(data, group =group_type) %>%
      ComplexUpset::upset(., 
                          intersect = groups, name = test_name,  
                          mode = 'exclusive_intersection',
                          set_sizes = (upset_set_size() +
                                         geom_text(aes(label=..count..), hjust=-0.8, stat='count', size = 3) +
                                         theme(axis.text=element_text(size=8))+
                                         scale_y_continuous(trans=reverse_log_trans()) 
                          ),
                          guides='over',
                          matrix=(
                            intersection_matrix(
                              geom=geom_point(
                                shape='square',
                                size=1.5
                              ),
                              segment=geom_segment(
                                linetype=3,lwd = 0.2 
                              ))),
                          base_annotations=list(
                            ' '=intersection_size(
                              counts=FALSE,
                              mapping=aes(fill=annotation)
                            ) + annotation_fill_scale2+ 
                              scale_y_continuous(limits = c(0, max_size), expand = c(0, 0))+
                              theme_classic() + theme(axis.ticks.x = element_blank(),
                                                      axis.title.x = element_blank(),
                                                      axis.text.x = element_blank())
                          ),
                          width_ratio = 0.2, 
                          max_size = max_size, 
                          wrap = TRUE,
                          group_by = 'sets',
                          themes = upset_modify_themes(list(
                            'intersections_matrix'=theme(text=element_text(size=15)), 
                            'overall_sizes'=theme(axis.text.x=element_text(size=10))
                          )), 
                          queries=query.list
      )
  }else if(group_type == 'pheno'){
    p <- format_upset_data(data, group =group_type) %>%
      ComplexUpset::upset(., 
                          intersect = groups, name = NULL, 
                          mode = 'exclusive_intersection',
                          set_sizes = FALSE,
                          guides='over',
                          matrix=(
                            intersection_matrix(
                              geom=geom_point(
                                shape='square',
                                size=1.5
                              ),
                              segment=geom_segment(
                                linetype=3,lwd = 0.2 
                              ))),
                          base_annotations=list(
                            ' '=intersection_size(
                              counts=FALSE,
                              mapping=aes(fill=annotation)
                            ) + annotation_fill_scale2+ 
                              scale_y_continuous(limits = c(0, max_size), expand = c(0, 0))+
                              theme_classic() + theme(axis.ticks.x = element_blank(),
                                                      axis.title.x = element_blank(),
                                                      axis.text.x = element_blank())
                          ),
                          width_ratio = 0.2,
                          max_size = max_size,
                          wrap = TRUE,
                          group_by = 'sets',
                          themes = upset_modify_themes(list(
                            'intersections_matrix'=theme(text=element_text(size=15)),
                            'overall_sizes'=theme(axis.text.x=element_text(size=12, angle = 45))
                          )),
                          queries=query.list
      )
  }else if(group_type == 'both'){
    p <- format_upset_data(data, group =group_type) %>%
      ComplexUpset::upset(., 
                          intersect = groups, name = test_name, sort_sets=FALSE, 
                          mode = 'exclusive_intersection',
                          set_sizes = FALSE,
                          guides='over',
                          matrix=(
                            intersection_matrix(
                              geom=geom_point(
                                shape='square',
                                size=1.5
                              ),
                              segment=geom_segment(
                                linetype=2, lwd = 0
                              ))),
                          stripes=upset_stripes(
                            mapping=aes(color=label),
                            colors=c(
                              'Biomarker'='#334195',
                              'Disease'='#b594b6'
                            ),
                            data=group_metadata
                          ), 
                          base_annotations=list(
                            ' '=intersection_size(
                              counts=FALSE,
                              mapping=aes(fill=annotation)
                            ) + annotation_fill_scale2 + 
                              scale_y_continuous(limits = c(0, max_size), expand = c(0, 0))+
                              theme_classic() + theme(axis.ticks.x = element_blank(),
                                                                                 axis.title.x = element_blank(),
                                                                                 axis.text.x = element_blank())
                          ),
                          width_ratio = 0.2,
                          max_size = max_size,
                          wrap = TRUE,
                          group_by = 'sets',
                          themes = upset_modify_themes(list(
                            'intersections_matrix'=theme(text=element_text(size=11.5)),
                            'overall_sizes'=theme(axis.text.x=element_text(size=12, angle = 45))
                          )),
                          queries=query.list
      )
    
  }
  return(p)
}

complex_pheno_upset <- function(group_type, test_type, max_size =100, save = T, height = 10, width = 20){
  data <- read_delim(paste0(data_path, group_type, '_domain_level_gene_sig_', test_type, '_', tranche,'.txt.bgz'), delim='\t', col_types = cols(phenocode = col_character())) %>%
    filter(!is.na(pheno_group))
  # data <- data %>% filter(annotation != 'pLoF|missense|LC')
  query.list <- list()
  group_metadata = NULL
  if(group_type == 'biomarker'){
    groups <- c("Bone and joint", "Cardiovascular", "Diabetes", "Hormone", "Liver", "Renal", "Blood")
    col.pal <- brewer.pal(n = 12, name = 'Set3')[c(1:6,8)]
  }else if(group_type == 'pheno'){
    keep <- as.data.frame(data %>% group_by(pheno_group) %>% summarise(sum(pheno_group_sig_cnt)))[,2] != 0
    groups <- unique(data$pheno_group)[keep]
    groups <- str_to_sentence(tolower(groups))[groups != 'Diet']
    groups <- sort(groups)
    col.pal <- pal_d3("category20")(length(groups))
  }else if(group_type == 'blood'){
    groups <- c("Red blood cells", "White blood cells", "Platelet", "Reticulocyte")
    col.pal <- brewer.pal(n = 12, name = 'Set3')[c(10, 9, 11, 12)]
  }else if(group_type == 'icd'){
    keep <- as.data.frame(data %>% group_by(pheno_group) %>% summarise(sum(pheno_group_sig_cnt)))[,2] != 0
    groups <- unique(unname(icd_names[unique(data$pheno_group)]))[keep]
    col.pal <- icd_colors[unique(data$pheno_group)][keep]
  }else if(group_type == 'both'){
    groups <- unique(data$pheno_group)
    icd_groups <- unique(unname(icd_names[groups[nchar(groups)==1]]))
    bio_groups <- c("Bone and joint", "Cardiovascular", "Diabetes", "Hormone", "Liver", "Renal", "Blood")
    groups <- c(bio_groups, icd_groups)
    col.pal <- rep(c('black'), length(groups))
    group_metadata <- data.frame(
      set=groups,
      label=c(rep('Biomarker', length(bio_groups)), rep('Disease', length(icd_groups)))
    )
  }
  names(col.pal) <- groups
  for(i in 1:length(groups)){
    query.list <- append(query.list, list(upset_query(group=groups[i], color=col.pal[groups[i]])))
    query.list <- append(query.list, list(upset_query(set=groups[i], fill=col.pal[groups[i]])))
  }
  p <- group_upset_plot(data, group_type, test_type, groups, query.list, max_size, group_metadata = group_metadata)
  if(save){
    png(paste0(figure_path, group_type, '_domain_level_gene_sig_', test_type, '_', tranche, '.png'), 
        height = height, width = width, units = 'in', res = 300)
    print(p)
    dev.off()
  }
  return(p)
}


## amino acid position 
get_amino_acid_data <- function(data){
  p_data = data %>%
    # filter(annotation == 'missense') %>%
    filter(!is.na(hgvsp)) %>%
    mutate(protein = str_split(hgvsp, ':') %>% map_chr(., 2)) %>%
    mutate(mutation = str_split(protein, '\\.') %>% map_chr(., 2)) %>%
    filter(nchar(mutation)<10) %>%
    mutate(amino_acid1 = mutation %>% str_sub(., 1, 3),
           position = mutation %>% str_extract(., "[[:digit:]]+"),
           amino_acid2 = mutation %>% str_sub(., -3, -1),)
  return(p_data)
}

save_amino_acid_figure <- function(data, gene, save = T){
  figure <- data %>%
    # filter(AF < 0.0001) %>%
    # filter(as.numeric(position)>200 & as.numeric(position)<400) %>%
    mutate(annotation = factor(if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation), levels = annotation_types))%>%
    filter(annotation != 'synonymous') %>%
    ggplot + aes(x = as.numeric(position), y = BETA, color = phenocode) + 
    geom_point(size = 0.5) + 
    geom_hline(yintercept = 0, lty=2) + 
    # ylim(-1.5,1.5) +
    # xlim(200, 400) +
    labs(x = paste0(gene, '- Amino Acid Position'))+
    geom_smooth(lwd = 0.5) +
    theme_classic() +
    scale_color_brewer(palette = 'Set1') +
    themes + theme(axis.title = element_text(size = 12),
                   axis.text = element_text(size = 8.5),
                   strip.text = element_text(size = 12)) + 
    facet_wrap(~annotation, ncol = 1, labeller = label_type)
  if(save){
    png(paste0(figure_path, gene, "_amino_acid_position.png"), width=6, height=6, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}
## icd color palette
icd_color_schema <- function(save_plot = T){
  x = rep(1:5, 5)[-25]
  y = rep(5:1, each=5)[-25]
  x_offset = rep(0.25, 24)
  y_offset = rep(0.25, 24)
  col_data = data.frame(x, y, x_offset, y_offset, icd_group = names(icd_colors), icd_colors)
  
  figure = col_data %>%
    ggplot + aes(x = x , 
                 y = y, label = paste(icd_group, '\n', icd_names[icd_group], '\n', icd_colors), fill = icd_colors) +
    geom_tile(color = icd_colors, fill=icd_colors) + theme_void() +
    scale_fill_manual(values = icd_colors)+
    geom_text(fontface='bold')
  
  if(save_plot){
    png(paste0(figure_path, 'icd_color_schema.png'), height = 7.5, width = 10, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

## heterogeneity test result

