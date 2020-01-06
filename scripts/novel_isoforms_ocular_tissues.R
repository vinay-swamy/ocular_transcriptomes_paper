library(tidyverse)
library(enrichR)
library(UpSetR)
library(matrixStats)
library(ggpubr)
# args <- c('/Volumes/data/ocular_transcriptomes_paper/',
#           '/Volumes/data/ocular_transcriptomes_pipeline/',
#           '~/NIH/ocular_transcriptomes_paper/new_data_122619/all_tissues.combined_NovelAno.gtf',
#           '/Volumes/data/ocular_transcriptomes_pipeline/data/rdata/novel_exon_classification.Rdata',
#           '/Volumes/data/ocular_transcriptomes_pipeline/sampleTableFull.tsv',
#           '~/NIH/ocular_transcriptomes_paper/new_data_122619/TCONS2MSTRG.tsv', 
#           '~/NIH/ocular_transcriptomes_paper/new_data_122619/all_tissue_quant.Rdata',
#           '/Volumes/data/ocular_transcriptomes_paper/clean_data/rdata/novel_isoforms.Rdata')



args <- commandArgs(trailingOnly = T)


save(args, file = 'testing/niot_args.rdata')
working_dir <- args[1] 
data_dir <- args[2]
gtf_file <- args[3]
novel_exon_class <- args[4]
sample_table_file <- args[5]
tcons_2_mstrg_file <- args[6]
quant_file <- args[7]
clean_outdata <- args[8]


setwd(working_dir)
load(novel_exon_class)

all_novel_tx <- c(novel_transcripts$transcript_id, novel_loci_distinct$transcript_id)
gtf <- rtracklayer::readGFF(gtf_file)
anno_tab <- gtf %>% filter(type == "transcript") %>% select(transcript_id, gene_name, oId)
tc2ms <- read_tsv(tcons_2_mstrg_file) %>% filter(!transcript_id %in% novel_loci_distinct$transcript_id)
eye_tissues <- c('Retina_Fetal.Tissue', 'Retina_Adult.Tissue', 'RPE_Fetal.Tissue', 'RPE_Adult.Tissue',
                 'Cornea_Adult.Tissue', 'Cornea_Fetal.Tissue', "ESC_Stem.Cell.Line", "Lens_Stem.Cell.Line")
brain_tissues <- tc2ms %>% select(contains('Brain')) %>% colnames
body_tissues <- tc2ms %>% select(-transcript_id, -class_code, -eye_tissues, -brain_tissues) %>% colnames 


sample_table <- read_tsv(sample_table_file)

## prep data 

det_mat <- apply(tc2ms[,-1],2, function(t) !is.na(t)) %>% cbind(tc2ms[,1], .)
det_eye <- det_mat[,eye_tissues] %>% {rowSums(.) > 0}
det_brain <- det_mat[,brain_tissues] %>% {rowSums(.) > 0}
det_body <- det_mat[,body_tissues] %>% {rowSums(.) == 0}

novel_eye_tx_det <- filter(det_mat, det_eye, transcript_id %in% all_novel_tx ) %>% select(c('transcript_id',eye_tissues))
novel_brain_tx_det <- filter(det_mat, det_brain, transcript_id %in% all_novel_tx ) %>% select(c('transcript_id',brain_tissues))
novel_body_det <- filter(det_mat, det_body, transcript_id %in% all_novel_tx ) %>% select(c('transcript_id',body_tissues))



##upset plot 
determine_specifcity <- function(x,df){
    df %>% filter(.[,x]) %>% pull(transcript_id)
}



novel_tx_by_tissue <- lapply(eye_tissues, determine_specifcity, novel_eye_tx_det)
names(novel_tx_by_tissue) <- eye_tissues
novel_eye_tx_by_tissue <- novel_tx_by_tissue


novel_tx_by_tissue[['Brain(all)']] <- novel_brain_tx_det$transcript_id
novel_tx_by_tissue[['Body(all)']] <- novel_body_det$transcript_id
#upset(fromList(novel_eye_tx_by_tissue),nintersects = 20, nsets = length(novel_eye_tx_by_tissue), order.by = 'freq')

#upset(fromList(novel_tx_by_tissue),nintersects = 20, nsets = length(novel_eye_tx_by_tissue), order.by = 'freq')
novel_eye_fetal <- novel_eye_tx_by_tissue[grepl('Fetal', names(novel_eye_tx_by_tissue))]
upset(fromList(novel_eye_fetal))


# piu plots 


sample_table_eye <- sample_table %>% filter(subtissue %in% eye_tissues)

t2g <- gtf %>% filter(type == 'transcript') %>% select(transcript_id, gene_name) %>% distinct
load(quant_file)
all_quant[is.na(all_quant)] <- 0
counts_eye <- all_quant[,c('transcript_id', sample_table_eye$sample)]
counts_eye_by_tissue <- lapply(eye_tissues,
                               function(tis) filter(sample_table_eye, subtissue == tis) %>% pull(sample) %>%
                                   {counts_eye[,c('transcript_id', .)]} %>%
                                   mutate(!!tis := rowMedians(.[,-1] %>% as.matrix)) %>%
                                   select(transcript_id, !!tis)
                                   ) %>%
    reduce(left_join) %>% left_join(t2g, .)


keep <- rowSums(counts_eye_by_tissue[,-(1:2)]) > 0
med_0_counts_eye <- counts_eye_by_tissue[!keep,]
counts_eye_by_tissue <- counts_eye_by_tissue[keep,] %>% distinct
eye_spec_tx <- novel_eye_tx_det





calc_isoform_percentage <- function(t_tissue){
    df <- counts_eye_by_tissue %>% select(transcript_id, gene_name, !!t_tissue)
    tt_sym <- as.symbol(t_tissue)
    df_gene_sums <- df %>% 
        select(-transcript_id) %>% 
        group_by(gene_name) %>%  
        summarise(gene_sums:= sum(!!tt_sym)) %>% 
        left_join(df, .) %>% 
        mutate(piu = .[[t_tissue]] / .[['gene_sums']], !!t_tissue :=NULL ) %>% 
        select(transcript_id, gene_name, !!t_tissue:=piu)
    return(df_gene_sums)
    
}

replace_nan <- function(df) {
    df_fixed <- lapply(colnames(df),function(col) pull(df, col) %>%  
                           {replace(., is.nan(.), 0)}) %>% bind_cols %>% as_tibble
    colnames(df_fixed) <- colnames(df)
    return(df_fixed)
    
}



piu_raw <- lapply(colnames(counts_eye_by_tissue)[-(1:2)], calc_isoform_percentage) %>% reduce(left_join)

piu <-replace_nan(piu_raw)



summarise_novel_transcripts_in_tissue <- function(t_tissue) {
    eye_spec_novel <- filter(eye_spec_tx, transcript_id %in% novel_transcripts$transcript_id)
    fetal <- paste0(t_tissue, '_Fetal.Tissue')
    adult <- paste0(t_tissue, '_Adult.Tissue')
    exp_in_tissue <- (eye_spec_novel[,-1] %>% select(contains(t_tissue)) %>% {rowSums(.) >0})
    not_exp_othertissues <-  eye_spec_novel[,-1] %>% select(-contains(t_tissue)) %>% {rowSums(.) == 0}
    tissue_spec_det <- eye_spec_novel[exp_in_tissue & not_exp_othertissues, ] 
    
    tissue_spec_det <- tissue_spec_det[,c('transcript_id', fetal, adult)]
    tissue_spec_det %>% pull(transcript_id) %>% unique() %>% length 
    tissue_spec_genes <-  tissue_spec_det %>% pull(transcript_id) %>% unique() %>% 
        {filter(gtf, transcript_id %in% .)} %>% pull(gene_name) %>% unique()
    fetal_exp <- tissue_spec_det %>% filter(.[,fetal]) 
    adult_exp <- tissue_spec_det %>% filter(.[,adult]) 
    both <- tissue_spec_det %>% filter(.[,adult], .[,fetal])
    # 
    # adult_piu_ns <- det_mat %>% 
    #     filter(.[,adult]) %>% 
    #     pull(transcript_id) %>% 
    #     {filter(piu, transcript_id %in% ., 
    #                  !transcript_id %in% tissue_spec_det$transcript_id,
    #                  !transcript_id %in% novel_transcripts$transcript_id ) } %>% 
    #     select(transcript_id, piu := !!adult) %>% 
    #     mutate(tissue=t_tissue, stage='adult_ns')
    # 
    # fetal_piu_ns <- det_mat %>% 
    #     filter(.[,fetal]) %>% 
    #     pull(transcript_id) %>% 
    #     {filter(piu, transcript_id %in% ., 
    #                  !transcript_id %in% tissue_spec_det$transcript_id,
    #                  !transcript_id %in% novel_transcripts$transcript_id ) } %>% 
    #     select(transcript_id, piu := !!fetal) %>% 
    #     mutate(tissue=t_tissue, stage='fetal_ns')
    
    
    
    adult_piu_spec <- tissue_spec_det %>% 
        filter(.[,adult]) %>% 
        pull(transcript_id) %>% 
        {filter(piu, transcript_id %in% .)} %>% 
        select(transcript_id, piu := !!adult) %>% 
        mutate(tissue=t_tissue, stage='adult')
    
    fetal_piu_spec <- tissue_spec_det %>% 
        filter(.[,fetal]) %>% 
        pull(transcript_id) %>% 
        {filter(piu, transcript_id %in% .)} %>% 
        select(transcript_id, piu := !!fetal) %>% 
        mutate(tissue=t_tissue, stage='fetal')
    
    
    
    fetal_locations <- gtf %>% 
        filter(transcript_id %in% fetal_exp$transcript_id) %>% 
        select(id, nv_type_rc, exon_location) %>% 
        filter(!is.na(id)) %>% 
        distinct %>% 
        pull(exon_location) %>% 
        table %>% 
        {./sum(.)}
    
    
    adult_locations <- gtf %>% 
        filter(transcript_id %in% adult_exp$transcript_id) %>% 
        select(id, nv_type_rc, exon_location) %>% 
        filter(!is.na(id)) %>% 
        distinct %>% 
        pull(exon_location) %>% 
        table %>% 
        {./sum(.)}
    list(loc_df=bind_rows(adult_locations, fetal_locations) %>% mutate(age=c('adult', 'fetal'), tissue=t_tissue) ,
         piu_df=bind_rows(adult_piu_spec, fetal_piu_spec)
    )
}

tissues <- c('Retina', 'Cornea', 'RPE')
res <- lapply(tissues, summarise_novel_transcripts_in_tissue)


location_df <- lapply(res, function(x) x[['loc_df']]) %>% bind_rows %>% gather(location, count, -age, -tissue)
piu_df <- lapply(res, function(x) x[['piu_df']]) %>% bind_rows

p <- ggboxplot(piu_df, x='stage', y='piu', color = 'tissue',
               title ='Comparison of percent isoform usage(piu) of novel \ntranscripts in fetal and adult eye tissues')+
    stat_compare_means(label.y = 1.1) +
    scale_color_manual(values = c('green', 'blue', 'red'))
facet(p,facet.by = 'tissue')+
    theme_minimal()

ggplot(location_df) + 
    geom_bar(aes(x=age, fill=location, y=count), position = 'fill', stat = 'identity' ) + 
    facet_wrap(~ tissue) + 
    ylab('percentage of novel exons') + 
    ggtitle('location of novel exons in occular tissues')
save(novel_tx_by_tissue, novel_eye_tx_by_tissue, piu_df, location_df, file = clean_outdata)





    
    
    
    
    
    
    
    