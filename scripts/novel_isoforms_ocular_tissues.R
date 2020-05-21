library(tidyverse)
library(data.table)
library(enrichR)
library(yaml)
library(matrixStats)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--workingDir', action='store', dest='working_dir')
parser$add_argument('--dataDir', action = 'store', dest = 'data_dir')
parser$add_argument('--filesYaml', action = 'store', dest = 'files_yaml')

list2env(parser$parse_args(), .GlobalEnv)

####
working_dir <- '/data/swamyvs/ocular_transcriptomes_paper/'
data_dir <- '/data/swamyvs/ocular_transcriptomes_pipeline/'
files_yaml <- '/data/swamyvs/ocular_transcriptomes_paper/files.yaml'
####
files <- read_yaml(files_yaml)
setwd(working_dir)
load(files$exon_class_rdata)

all_novel_tx <- c(novel_transcripts$transcript_id) %>% unique
gtf <- rtracklayer::readGFF(files$anno_gtf)
anno_tab <- gtf %>% filter(type == "transcript") %>% select(transcript_id, gene_name, oId)
tc2ms_all <- fread(files$tcons2mstrg) %>% as_tibble
sample_table <- read_tsv(files$sample_table)
subtissues <- unique(sample_table$subtissue)
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

sample_table_eye <- sample_table %>% filter(subtissue %in% eye_tissues)

t2g <- gtf %>% 
    filter(type == 'transcript') %>% 
    select(transcript_id, gene_name) %>% 
    distinct
load(files$all_tissue_quant)
all_quant[is.na(all_quant)] <- 0
counts_eye <- all_quant[,c('transcript_id', sample_table_eye$sample)]
counts_by_tissue <- lapply(subtissues,
                               function(tis) filter(sample_table, subtissue == tis) %>% pull(sample) %>%
                                   {all_quant[,c('transcript_id', .)]} %>%
                                   mutate(!!tis := rowMedians(.[,-1] %>% as.matrix)) %>%
                                   select(transcript_id, !!tis)
                                   ) %>%
    reduce(full_join) %>% left_join(t2g, .)

keep <- rowSums(counts_by_tissue[,-(1:2)]) > 0
med_0_counts <- counts_by_tissue[!keep,]
counts_by_tissue <- counts_by_tissue[keep,] %>% distinct



# long_counts_with_codes <- counts_eye_by_tissue %>% 
#     inner_join(tc2ms_all %>% select(transcript_id, class_code), .) %>% 
#     gather(key = 'subtissue', value = 'counts', -transcript_id, -class_code, -gene_name) %>% 
#     mutate(log2counts = log2(counts + 1), 
#            transcript_novelty = case_when(transcript_id %in% all_novel_tx ~ 'novel_isoforms', 
#                                          class_code == 'u' ~ 'novel Loci',
#                                          class_code == '='~ 'ref_isoform'
#                                          ))
# ggplot(long_counts_with_codes) +
#     geom_boxplot(aes(x = transcript_novelty, y= log2counts, fill = subtissue) )+
#     facet_wrap(~subtissue)



calc_isoform_percentage <- function(t_tissue){
    df <- counts_by_tissue %>% select(transcript_id, gene_name, !!t_tissue)
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



piu_raw <- lapply(colnames(counts_by_tissue)[-(1:2)], calc_isoform_percentage) %>% reduce(full_join)
piu_all <-replace_nan(piu_raw)
piu <- piu_all %>% select(transcript_id, gene_name, all_of(eye_tissues) )
#### now identify primary isoforms for each tissue
appris_tags <- fread(files$ref_gtf, header= F ) %>% 
    as_tibble %>% 
    mutate(appris_tag = str_extract(V9, "appris_\\w+_\\d+"), 
           transcript_id = str_extract(V9,"transcript_id \"\\w+\\d+\\.\\d+\";") %>% str_split('"') %>% sapply(function(x)x[2])) %>% 
    filter(V3 == 'transcript') %>% select(transcript_id, appris_tag)

dntx2enst <- filter(tc2ms_all, class_code == '=') %>% select(dntx_id = transcript_id, transcript_id = refid)


appris_tags_clean <- filter(appris_tags, !is.na(appris_tag)) %>% 
    mutate(tag_level = str_split(appris_tag, '_') %>% sapply(function(x) x[3]) %>% as.numeric, 
           tag_level = case_when(grepl('alternative', appris_tag) & tag_level == 1 ~ 6, 
                                 grepl('alternative', appris_tag) & tag_level == 2 ~ 7,
                                 TRUE ~ tag_level), 
    )
table(appris_tags_clean$tag_level)
appris_in_dntx <-  inner_join(dntx2enst, appris_tags_clean)
table(appris_in_dntx$tag_level)
nrow(appris_in_dntx %>% filter(tag_level <6)) / nrow(appris_tags_clean %>% filter(tag_level <6))
tx2code <- tc2ms_all %>% select(transcript_id, class_code)
primary_isoforms_per_tissue <- piu_all %>% 
    gather(key = 'subtissue', value = 'piu', -transcript_id, -gene_name) %>% 
    group_by(subtissue, gene_name) %>% 
    summarise(primary_iso = max(piu), transcript_id  = transcript_id[which.max(piu)]) %>% 
    ungroup %>% 
    filter(primary_iso >0) %>% 
    mutate(is_appris = transcript_id %in% appris_in_dntx$dntx_id) %>% 
    inner_join(tx2code) %>% 
    mutate(primary_isoform_origin = case_when(is_appris ~ 'appris',  
                                              !is_appris & class_code == '=' ~ 'gencode', 
                                              !is_appris & class_code != '=' ~ 'dntx')) %>% 
    group_by(subtissue, primary_isoform_origin) %>% 
    summarise(count = n())


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
        mutate(tissue=t_tissue, stage='adult',subtissue =adult  )
    
    fetal_piu_spec <- tissue_spec_det %>% 
        filter(.[,fetal]) %>% 
        pull(transcript_id) %>% 
        {filter(piu, transcript_id %in% .)} %>% 
        select(transcript_id, piu := !!fetal) %>% 
        mutate(tissue=t_tissue, stage='fetal', subtissue = fetal)
    
    list(loc_df <- tibble(),
         #loc_df=bind_rows(adult_locations, fetal_locations) %>% mutate(age=c('adult', 'fetal'), tissue=t_tissue) ,
         piu_df=bind_rows(adult_piu_spec, fetal_piu_spec)
    )
}

tissues <- c('Retina', 'Cornea', 'RPE')
res <- lapply(tissues, summarise_novel_transcripts_in_tissue)


location_df <- lapply(res, function(x) x[['loc_df']]) %>% bind_rows #%>% gather(location, count, -age, -tissue)
long_eye_det <- det_mat %>% select(transcript_id, all_of(eye_tissues)) %>% 
    filter(rowSums(.[,-1]) >0) %>% 
    gather(key = 'subtissue', value = 'det', -transcript_id)

piu_df <- lapply(res, function(x) x[['piu_df']]) %>% bind_rows %>% inner_join(long_eye_det) %>% filter(det)


save(novel_tx_by_tissue, novel_eye_tx_by_tissue, piu_df, location_df, primary_isoforms_per_tissue, file = files$novel_isoform_analysis_rdata)





    
    
    
    
    
    
    
    