library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(viridis)
library(matrixStats)
args <- c('~/NIH/occular_transcriptomes_paper/', '~/NIH/eyeintegration_splicing/dl_data/all_tissues_psi.tsv.gz', 
          '~/NIH/eyeintegration_splicing/sampleTableV6.tsv', 
          '~/NIH/eyeintegration_splicing/dl_data/gfc_TCONS_to_st_MSTRG.tsv.gz', 
          '~/NIH/occular_transcriptomes_paper/clean_data/V1_exon_classification_data.Rdata',
          '~/NIH/occular_transcriptomes_paper/all_tissues.combined.gtf',
          '~/NIH/occular_transcriptomes_paper/clean_data/tissue_to_colors.Rdata',
          '~/NIH/occular_transcriptomes_paper/clean_data/splicing_analysis.Rdata')
args <- commandArgs(trailingOnly = T)
wd <- args[1]
psi_file <- args[2]
sample_file <- args[3]
tc2mstrg_file <- args[4]
exon_classification_file <- args[5]
gtf_file <- args[6]
color_df <- args[7]
outfile <- args[8]
#save(args, file='/tmp/sphm_args.rdata')
setwd(wd)

tcons2mstrg <- read_tsv(tc2mstrg_file)
load(exon_classification_file)
psi_tab <- read_tsv(psi_file)
colnames(psi_tab) <- str_remove(colnames(psi_tab), '_psi')
full_gtf <- rtracklayer::readGFF(gtf_file)
load(color_df)
sample_table <- read_tsv(sample_file) %>% left_join(tissue_color_mapping_df) %>%  filter(sample %in% colnames(psi_tab))


psi_tab[is.na(psi_tab)] <- 0
set.seed(34543)
mat <- psi_tab[,-(1:4)] %>% filter(rowSums(.)< (1*(ncol(.)-3)), rowSums(.)>(.05*ncol(.))) %>% sample_n(10000) %>% .[,sample_table$sample] %>% t()

t_mat <- mat[, sample(1:ncol(mat), 1000)]
ht_opt(fast_hclust = TRUE)

tcolors <- tissue_color_mapping_df$color
names(tcolors) <- tissue_color_mapping_df$body_location
color_list <- list( body_location=tcolors)

col_list <- tissue_color_mapping_df$color
names(col_list) <- tissue_color_mapping_df$body_location
body_location <- sample_table$body_location
ha <- HeatmapAnnotation(body_location=sample_table$body_location, col = color_list,which = 'row')
hm <- Heatmap(t_mat,col = viridis(100)  ,name = 'PSI', right_annotation = ha, show_row_dend = F, show_column_dend =F,
              heatmap_height = unit(20,'cm'), heatmap_width = unit(40,'cm'), 
              show_column_names = F, show_row_names = F)
splicing_heatmap <- draw(hm)
#splicing_heatmap
save(splicing_heatmap,file = outfile)
MIN_PSI <- .1
summarise_splicing <- function(s_tissue){
    ctab <- tcons2mstrg %>% select(transcript_id, !!s_tissue) %>% filter(!is.na(.[,s_tissue]))
    exons_in_tissue <- full_gtf %>% filter(transcript_id %in% ctab$transcript_id) %>% 
        select(seqid, strand, start, end) %>% distinct
    novel_exons_in_tissue_splicing <- novel_exons_TSES %>% inner_join(exons_in_tissue) %>% 
        filter(!nv_type_rc %in% c( 'novel_TES', 'novel_TSS'))
    psi_tissue <- filter(sample_table, subtissue == s_tissue) %>% pull(sample) %>% 
        {select(psi_tab, seqid, strand, start, end, .)} %>% mutate( start=start+1) %>% 
        inner_join(novel_exons_in_tissue_splicing,.)
    meta_cols <- colnames(novel_exons_in_tissue_splicing)
    psi_only <- psi_tissue %>% select(-meta_cols)
    not_det <- psi_only %>% apply(2, is.na) %>% {rowSums(.) == ncol(.) }
    psi_tissue_det <- psi_tissue %>% filter(!not_det) %>% mutate(avg_psi= select(., -meta_cols) %>% rowMeans())
    tibble(num_const_exons=nrow(novel_exons_in_tissue_splicing), num_det_rmats=nrow(psi_tissue_det),
           num_exp_psi=psi_tissue_det %>% filter(avg_psi >= MIN_PSI ) %>% nrow)
}
subtissues <- filter(sample_table, !subtissue %in%c('Cornea_Fetal.Tissue', 'synth')) %>% pull(subtissue) %>% unique
splicing_sum <- lapply(subtissues, summarise_splicing) %>% bind_rows %>% mutate(subtissue=subtissues) %>% 
    select(subtissue,everything())



summarise_splicing_not_found <- function(s_tissue){
    ctab <- tcons2mstrg %>% select(transcript_id, !!s_tissue) %>% filter(!is.na(.[,s_tissue]))
    exons_in_tissue <- full_gtf %>% filter(transcript_id %in% ctab$transcript_id) %>% 
        select(seqid, strand, start, end) %>% distinct
    novel_exons_in_tissue_splicing <- novel_exons_TSES %>% inner_join(exons_in_tissue) %>% 
        filter(!nv_type_rc %in% c( 'novel_TES', 'novel_TSS'))
    psi_tissue <- filter(sample_table, subtissue == s_tissue) %>% pull(sample) %>% 
        {select(psi_tab, seqid, strand, start, end, .)} %>% mutate( start=start+1) %>% 
        inner_join(novel_exons_in_tissue_splicing,.)
    meta_cols <- colnames(novel_exons_in_tissue_splicing)
    psi_only <- psi_tissue %>% select(-meta_cols)
    not_det <- psi_only %>% apply(2, is.na) %>% {rowSums(.) == ncol(.) }
    psi_tissue_det <- psi_tissue %>% filter(!not_det) %>% mutate(avg_psi= select(., -meta_cols) %>% rowMeans())
    psi_tissue_det %>% filter(avg_psi < MIN_PSI ) %>% pull(nv_type_rc) %>% table 
}
subtissues <- filter(sample_table, !subtissue %in%c('Cornea_Fetal.Tissue', 'synth')) %>% pull(subtissue) %>% unique
not_det_splicing_types <- lapply(subtissues, summarise_splicing_not_found)
undetected_exons_by_event <-  not_det_splicing_types %>% do.call(rbind, .) %>% as_tibble() %>% gather(event_type, misclassed_events)

save(splicing_heatmap, splicing_sum,undetected_exons_by_event, file = outfile)
