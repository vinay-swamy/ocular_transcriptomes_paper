library(ggpubr)
library(tidyverse)
count_num_det_ref_novel <- function(s_tissue){
    nm_col_clean <- function(col){
        raw_name <- col %>% .[.!='-'] %>% .[!grepl('ENST', .)] %>% head
        name=str_split(raw_name, '\\.\\d|:|\\|')[[1]][2] %>% str_split('_MSTRG') %>% .[[1]] %>% .[1]
        return(name)
    }
    
    gtf_file <- paste0(path_to_raw_tissues_gtfs, s_tissue, '.combined.gtf')
    tracking_file <- paste0(path_to_raw_tissues_gtfs, s_tissue, '.tracking')
    gtf <- rtracklayer::readGFF(gtf_file) %>% mutate(strand=ifelse(strand == '-', '-', '+'))
    track_tab <- read_tsv(tracking_file, col_names=F)
    names <- c('transcript_id', 'gene_id','refid','code', apply(track_tab[,-(1:4)], 2, nm_col_clean))
    
    if(any(is.na(names)) ){
        stop()
    }
    colnames(track_tab) <- names
    
    #first, identify de novo transcripts that absolutely match the refernece annotation
    num_samp <- ncol(track_tab)-4
    gtf_tx <- gtf %>% filter(type == 'transcript') 
    gffc_ref_absmatch <- gtf_tx %>% filter(class_code == '=') %>% inner_join(all_transcripts) %>%
        pull(transcript_id)
    # next remove transcripts that dont meet detection req
    det_df <- apply(track_tab[,-(1:4)],2, function(x) x!='-') %>% as.data.frame %>%  bind_cols(track_tab[,1:4],.)
    num_det_df <-det_df %>% mutate(num_det= rowSums(det_df[,-(1:4)]),
                                   frac_det= num_det / num_samp,
                                   ano_type=ifelse(transcript_id %in% gffc_ref_absmatch, 'ref', 'novel'),
                                   subtissue=s_tissue) %>%   
        select(transcript_id, gene_id, refid, code,ano_type,subtissue, num_det, frac_det) 
    return(num_det_df)
}

eye_tissues <- c('Retina_Adult.Tissue', 'Retina_Fetal.Tissue', 'RPE_Adult.Tissue', 'RPE_Fetal.Tissue', 'Cornea_Adult.Tissue', 'Cornea_Fetal.Tissue')
raw_eye_det <- lapply(eye_tissues, count_num_det_ref_novel) %>% bind_rows()

base<- ggplot(data = all_eye_tissue_det_counts)+
    geom_boxplot(aes(x=ano_type, y=frac_det)) + 
    facet_wrap(~subtissue)  + 
    ggtitle('Fraction of total samples transcripts are constructed in, no filtering') +
    theme_minimal()

read_detfile <- function(s_tissue){
    detfile <- paste0('/Volumes/data/ocular_transcriptomes_pipeline/data/misc/raw_dd/',s_tissue,'.dd.tsv.gz')
    det <- read_tsv(detfile)
    num_samp <- ncol(det)-4
    num_det <- det %>% mutate(num_det= rowSums(det[,-(1:4)]),
                              frac_det= num_det / num_samp,
                              subtissue= s_tissue,
                              ano_type = ifelse(code == '=', 'ref', 'novel')) %>% 
        select(transcript_id, gene_id, refid, code, subtissue, ano_type, frac_det, num_det)
    return(num_det)
}

gfc_eye_det <- lapply(eye_tissues, read_detfile) %>% bind_rows()





gfc_filt <- ggplot(data = all_eye_det_filt)+
    geom_boxplot(aes(x=ano_type, y=frac_det)) + 
    facet_wrap(~subtissue)  + 
    ggtitle('Fraction of total samples transcripts are constructed in, filtering by study') +
    theme_minimal()



read_final_detfile <- function(s_tissue){
    detfile <- paste0('/Volumes/data/ocular_transcriptomes_pipeline/data/misc/final_dd/',s_tissue,'.dd.tsv')
    det <- read_tsv(detfile) %>% select(transcript_id, everything())
    num_samp <- ncol(det)-4
    num_det <- det %>% mutate(num_det= rowSums(det[,-(1:4)]),
                              frac_det= num_det / num_samp,
                              subtissue= s_tissue,
                              ano_type = ifelse(code == '=', 'ref', 'novel')) %>% 
        select(transcript_id, gene_id, refid, code, subtissue, ano_type, frac_det, num_det)
    return(num_det)
}

final_eye_det <- lapply(eye_tissues, read_final_detfile) %>% bind_rows()

final <- ggplot(data = final_eye_det)+
    geom_boxplot(aes(x=ano_type, y=frac_det)) + 
    facet_wrap(~subtissue)  + 
    ggtitle('Fraction of total samples transcripts are constructed in, final') +
    theme_minimal()

ggarrange(base, gfc_filt, final, ncol=1, labels = "AUTO") 
save(raw_eye_det, gfc_eye_det, final_eye_det, file = '/Volumes/data/ocular')




