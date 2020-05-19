library(tidyverse)
library(data.table)
#library(enrichR)
library(yaml)
library(matrixStats)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--workingDir', action='store', dest='working_dir')
parser$add_argument('--dataDir', action = 'store', dest = 'data_dir')
parser$add_argument('--filesYaml', action = 'store', dest = 'files_yaml')

list2env(parser$parse_args(), .GlobalEnv)

####
# working_dir <- '/data/swamyvs/ocular_transcriptomes_paper/'
# data_dir <- '/data/swamyvs/ocular_transcriptomes_pipeline/'
# files_yaml <- '/data/swamyvs/ocular_transcriptomes_paper/files.yaml'
####
files <- read_yaml(files_yaml)
setwd(working_dir)
load(files$exon_class_rdata)

all_novel_tx <- c(novel_transcripts$transcript_id) %>% unique
gtf <- rtracklayer::readGFF(files$anno_gtf)
anno_tab <- gtf %>% filter(type == "transcript") %>% select(transcript_id, gene_name, oId)
tc2ms_all <- fread(files$tcons2mstrg) %>% as_tibble
sample_table <- read_tsv(files$sample_table)
tc2ms <- tc2ms_all %>% filter(!transcript_id %in% novel_loci_distinct$transcript_id) 
eye_tissues <- c('Retina_Fetal.Tissue', 'Retina_Adult.Tissue', 'RPE_Fetal.Tissue', 'RPE_Adult.Tissue',
                 'Cornea_Adult.Tissue', 'Cornea_Fetal.Tissue')
brain_tissues <- sample_table %>% filter(tissue == 'Brain') %>% pull(subtissue) %>% unique
body_tissues <- sample_table %>% 
    filter(!tissue%in% c(eye_tissues, brain_tissues)) %>% pull(subtissue) %>% unique




## prep data 

det_mat <- apply(tc2ms[,-(1:4)],2, function(t) t!='') %>% cbind(tc2ms[,1:4], .)
det_eye <- det_mat[,eye_tissues] %>% {rowSums(.) > 0}
det_brain <- det_mat[,brain_tissues] %>% {rowSums(.) > 0}
det_body <- det_mat[,body_tissues] %>% {rowSums(.) > 0}

novel_eye_tx_det <- filter(det_mat, det_eye, transcript_id %in% all_novel_tx ) %>% select(all_of(c('transcript_id',eye_tissues)))
novel_brain_tx_det <- filter(det_mat, det_brain, transcript_id %in% all_novel_tx ) %>% select(all_of(c('transcript_id',brain_tissues)))
novel_body_det <- filter(det_mat, det_body, transcript_id %in% all_novel_tx ) %>% select(all_of(c('transcript_id',body_tissues)))
novel_eye_exp_tx <- novel_eye_tx_det %>% inner_join(tc2ms %>% select(transcript_id, class_code),.)


##upset plot 
determine_specifcity <- function(x,df){
    df %>% filter(.[,x]) %>% pull(transcript_id)
}



novel_tx_by_tissue <- lapply(eye_tissues,function(x) determine_specifcity(x, novel_eye_tx_det))
names(novel_tx_by_tissue) <- eye_tissues
novel_eye_tx_by_tissue <- novel_tx_by_tissue


novel_tx_by_tissue[['Brain(all)']] <- novel_brain_tx_det$transcript_id
novel_tx_by_tissue[['Body(all)']] <- novel_body_det$transcript_id
#upset(fromList(novel_eye_tx_by_tissue),nintersects = 20, nsets = length(novel_eye_tx_by_tissue), order.by = 'freq')




# piu plots 


sample_table_eye <- sample_table %>% filter(subtissue %in% eye_tissues)

t2g <- gtf %>% 
    filter(type == 'transcript') %>% 
    select(transcript_id, gene_name) %>% 
    distinct
load(files$all_tissue_quant)
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



long_counts_with_codes <- counts_eye_by_tissue %>% 
    inner_join(tc2ms_all %>% select(transcript_id, class_code), .) %>% 
    gather(key = 'subtissue', value = 'counts', -transcript_id, -class_code, -gene_name) %>% 
    mutate(log2counts = log2(counts + 1), 
           transcript_novelty = case_when(transcript_id %in% all_novel_tx ~ 'novel_isoforms', 
                                         class_code == 'u' ~ 'novel Loci',
                                         class_code == '='~ 'ref_isoform'
                                         ))
# ggplot(long_counts_with_codes) +
#     geom_boxplot(aes(x = transcript_novelty, y= log2counts, fill = subtissue) )+
#     facet_wrap(~subtissue)



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


save.image('testing/novel_iso_ws.rdata')
summarise_novel_transcripts_in_tissue <- function(t_tissue) {

    
    fetal <- paste0(t_tissue, '_Fetal.Tissue')
    adult <- paste0(t_tissue, '_Adult.Tissue')
    
    exp_in_tissue <- (novel_eye_exp_tx[,-(1:2)] %>% select(contains(t_tissue)) %>% {rowSums(.) >0})
    tissue_spec_det <- novel_eye_exp_tx %>% filter(exp_in_tissue, transcript_id %in% all_novel_tx) %>% select(transcript_id, !!fetal, !!adult)
        
    
    tissue_spec_genes <-  tissue_spec_det %>% pull(transcript_id) %>% unique() %>% 
        {filter(gtf, transcript_id %in% .)} %>% pull(gene_name) %>% unique()
    fetal_exp <- tissue_spec_det %>% filter(.[,fetal]) 
    adult_exp <- tissue_spec_det %>% filter(.[,adult]) 
    both <- tissue_spec_det %>% filter(.[,adult], .[,fetal])
    
    

    
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
    
    list(loc_df <- tibble(),
         #loc_df=bind_rows(adult_locations, fetal_locations) %>% mutate(age=c('adult', 'fetal'), tissue=t_tissue) ,
         piu_df=bind_rows(adult_piu_spec, fetal_piu_spec)
    )
}

tissues <- c('Retina', 'Cornea', 'RPE')
res <- lapply(tissues, summarise_novel_transcripts_in_tissue)


location_df <- lapply(res, function(x) x[['loc_df']]) %>% bind_rows #%>% gather(location, count, -age, -tissue)
piu_df <- lapply(res, function(x) x[['piu_df']]) %>% bind_rows
save(novel_tx_by_tissue, novel_eye_tx_by_tissue, piu_df, location_df, file = files$novel_isoform_analysis_rdata)





    
    
    
    
    
    
    
    