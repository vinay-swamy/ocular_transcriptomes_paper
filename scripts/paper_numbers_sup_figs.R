'
NUMBERS WE NEED TO AUTO GENERATE:
-number of samples 
-number of tissues
-final pan body transcriptome ref counts 
-final pan body transcriptoem eye counts 
-final pan eye transcriptome """""""
-final amount of unannotated sequence 
-net gain of reads 
SUP FIGS
- boxplot of reference and novel frac sample det
'
library(tidyverse)
library(RBedtools)
library(ggpubr)
args <- c('/Volumes/data/ocular_transcriptomes_paper/', 
          '/Volumes/data/ocular_transcriptomes_pipeline/',
          '/Volumes/data/ocular_transcriptomes_pipeline/sampleTableFull.tsv', 
          '/Volumes/data/ocular_transcriptomes_pipeline/data/gtfs/all_tissues.combined.gtf',
          '/Volumes/data/ocular_transcriptomes_paper/clean_data/pan_eye_txome.combined.gtf',
          '/Volumes/data/ocular_transcriptomes_pipeline/data/gtfs/raw_tissue_gtfs/',
          '/Volumes/data/ocular_transcriptomes_pipeline/rdata/all_ref_tx_exons.rdata',
          '/Volumes/data/ocular_transcriptomes_pipeline/ref/core_tight.Rdata',
          '/Volumes/data/ocular_transcriptomes_paper/clean_data/gencode_salmon_mapping_rates.tab',
          '/Volumes/data/ocular_transcriptomes_paper/clean_data/DNTX_salmon_mapping_rates.tab',
          '/Volumes/data/ocular_transcriptomes_paper/clean_data/rdata/paper_numbers.Rdata',
          '/Volumes/data/ocular_transcriptomes_paper/clean_data/rdata/sup_fig_data.Rdata')

args <- commandArgs(trailingOnly = T)

working_dir <- args[1]
data_dir <- args[2]
sample_table_file <- args[3]
pan_body_gtf_file <- args[4]
pan_eye_gtf_file <- args[5]
path_to_raw_tissues_gtfs <- args[6]
all_ref_ano <- args[7]
core_tight_file <- args[8]
gencode_mapping_rate_file <- args[9]
DNTX_mapping_rate_file <- args[10]
out_num_file <- args[11]
out_supdata_file <- args[12]

setwd(working_dir)


sample_table <- read_tsv(sample_table_file)
NUM_TOTAL_SAMPLES <- nrow(sample_table)
NUM_EYE_SAMPLES <- sample_table %>% filter(!body_location %in% c('Brain', 'Body')) %>% nrow 
NUM_BODY_SAMPLES <- sample_table %>% filter(body_location %in% c('Brain', 'Body')) %>% nrow 
NUM_BODY_TISSUES <- sample_table %>% filter(body_location %in% c('Brain', 'Body')) %>% 
    pull(subtissue) %>% unique %>% length 

pan_body_gtf <- rtracklayer::readGFF(pan_body_gtf_file)
pan_eye_gtf <- rtracklayer::readGFF(pan_eye_gtf_file)
NUM_REF_TX_BODY <- pan_body_gtf %>% filter( type == 'transcript', class_code=='=')%>% nrow 
NUM_REF_TX_EYE  <- pan_eye_gtf  %>% filter( type == 'transcript', class_code=='=')%>% nrow 

NUM_NOVEL_TX_BODY <- pan_body_gtf %>% filter( type == 'transcript', class_code!='=')%>% nrow 
NUM_NOVEL_TX_EYE  <- pan_eye_gtf  %>% filter( type == 'transcript', class_code!='=')%>% nrow 
load(all_ref_ano)
all_exon_bed <- all_exons %>% mutate(score= 888) %>% select(seqid, start, end, origin, score, strand) %>% 
    from_data_frame %>% 
    RBedtools('sort', i=.)
dntx_exon_bed <- pan_body_gtf %>% filter(type == 'exon') %>% mutate(score = 999) %>% 
    select(seqid, start, end, transcript_id, score, strand) %>% from_data_frame %>% 
    RBedtools('sort', i=.)
novel_sequence <- RBedtools('subtract', options = '-s', output = 'stdout', a=dntx_exon_bed, b=all_exon_bed) %>% 
    RBedtools('sort', output = 'stdout', i=.) %>% 
    RBedtools('merge', options = '-s -c 6 -o distinct',i=. )%>% 
    to_data_frame
NUM_TOTAL_NOVEL_SEQ <- novel_sequence %>% mutate(length=X3-X2) %>% pull(length) %>% sum 

DNTX_mapping_rates <- process_lib_size_tabs(DNTX_mapping_rate_file , 'DNTX')
gencode_mapping_rates <- process_lib_size_tabs(gencode_mapping_rate_file, 'gencode')




## SUPPLEMENTARY FIGURES 
#---- 
load(core_tight_file)
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
save(list= ls()[grepl('NUM_', ls())], file = out_num_file)
save(raw_eye_det, gfc_eye_det, final_eye_det, file = out_supdata_file)
