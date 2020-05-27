library(tidyverse)
library(RBedtools)
library(argparse)
library(yaml)

parser <- ArgumentParser()
parser$add_argument('--workingDir', action = 'store', dest = 'working_dir')
parser$add_argument('--dataDir', action = 'store', dest = 'data_dir')
parser$add_argument('--filesYaml', action = 'store', dest = 'files_yaml')
list2env(parser$parse_args(), .GlobalEnv)

###
# working_dir <- '/data/swamyvs/ocular_transcriptomes_paper/'
# files_yaml <- '/data/swamyvs/ocular_transcriptomes_paper/files.yaml'
###


setwd(working_dir)
files <- read_yaml(files_yaml)
######### base numbers 
sample_table <- read_tsv(files$sample_table)
NUM_TOTAL_SAMPLES <- nrow(sample_table)
NUM_EYE_SAMPLES <- sample_table %>% filter(!body_location %in% c('Brain', 'Body')) %>% nrow 
NUM_BODY_SAMPLES <- sample_table %>% filter(body_location %in% c('Brain', 'Body')) %>% nrow 
NUM_BODY_TISSUES <- sample_table %>% filter(body_location %in% c('Brain', 'Body') | subtissue%in%c('Lens_Stem.Cell.Line', 'ESC_Stem.Cell.Line') ) %>% 
    pull(subtissue) %>% unique %>% length 
NUM_EYE_TISSUES <- sample_table %>% filter(!(body_location %in% c('Brain', 'Body') | subtissue%in%c('Lens_Stem.Cell.Line', 'ESC_Stem.Cell.Line')) ) %>% 
    pull(subtissue) %>% unique %>% length 
load(files$core_tight_rdata)
core_tight <- core_tight %>% as_tibble %>%  filter(sample_accession %in% sample_table$sample)
NUM_STUDIES = length(unique(core_tight$study_accession))

pan_body_gtf <- rtracklayer::readGFF(files$anno_gtf)
pan_eye_gtf <- rtracklayer::readGFF(files$pan_eye_gtf)
NUM_REF_TX_BODY <- pan_body_gtf %>% filter( type == 'transcript', class_code=='=')%>% nrow 
NUM_REF_TX_EYE  <- pan_eye_gtf  %>% filter( type == 'transcript', class_code=='=')%>% nrow 

NUM_NOVEL_TX_BODY <- pan_body_gtf %>% filter( type == 'transcript', class_code!='=')%>% nrow 
NUM_NOVEL_TX_EYE  <- pan_eye_gtf  %>% filter( type == 'transcript', class_code!='=')%>% nrow 
NUM_TOTAL_TX_BODY <- NUM_NOVEL_TX_BODY + NUM_NOVEL_TX_BODY
NUM_TOTAL_TX_BODY <- NUM_NOVEL_TX_EYE + NUM_NOVEL_TX_EYE

load(files$ref_tx_exon_rdata)
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
#######ORFs
cds_origin_df <- fread(files$CDS_origin_df)
NUM_REF_ORF <-  cds_origin_df %>% select(-transcript_id) %>% distinct %>% pull(is_annotaed_cds) %>% {sum(.)}
NUM_NOVEL_ORF <- cds_origin_df %>% select(-transcript_id) %>% distinct %>% pull(is_annotaed_cds) %>% {sum(!.)}





####long read data numbers 
load(files$long_read_results_rdata)
NUM_BASE_TXOME_CONST_ACC <- acc_df_alltx %>% filter(build == 'stringtie') %>% pull(accuracy) %>% round(digits = 3)
NUM_TXOME_CONST_ACC_2000_ABOVE <- acc_df_by_length %>% filter(as.numeric(qbin) >2000) %>% pull(stringtie_acc) %>% mean %>% round(digits = 3)
NUM_TXOME_CONST_ACC_2000_BELOW <- acc_df_by_length %>% filter(as.numeric(qbin) <=2000) %>% pull(stringtie_acc) %>% mean %>% round(digits = 3)
NUM_FILT_TXOME_ACC_2000_ABOVE <- tpm_df_plotting %>% filter(filtering_type == 'mean',co == 1 ) %>% 
    filter(as.numeric(qbin) >2000) %>% pull(acc) %>% mean %>% round(digits = 3)
NUM_FILT_TXOME_ACC_2000_BELOW <- tpm_df_plotting %>% filter(filtering_type == 'mean',co == 1 ) %>% 
    filter(as.numeric(qbin) <=2000) %>% pull(acc) %>% mean %>% round(digits = 3)
#######
####### novel isoform numbers 
load(files$novel_isoform_analysis_rdata)
appris_totals <- primary_isoforms_per_tissue %>% 
    spread(key = primary_isoform_origin, value = count) %>% 
    ungroup %>% 
    mutate(total =  rowSums(.[,-1]),
           percent_appris = appris / total,
           percent_gencode = gencode / total, 
           percent_novel = total )
NUM_AVG_APPRIS <- mean(appris_totals$percent_appris)
NUM_AVG_NONAPPRIS_GENCODE = mean(appris_totals$percent_gencode)
NUM_AVG_NOVEL = mean(appris_totals$percent_novel)

####### Correlation between expression and transcript length 
load(files$gencode_quant)
subtissues <- unique(sample_table$subtissue)
ref_gtf <- rtracklayer::readGFF(files$ref_gtf)
ref_tx_lengths <- filter(ref_gtf, type == 'exon') %>% select(seqid, strand, start, end, transcript_id) %>% 
    mutate(length = end-start) %>% 
    group_by(transcript_id) %>% 
    summarise(length = sum(length) + n())

calc_txome_size <- function(t_tissue){
    samples <- filter(sample_table, subtissue == t_tissue) %>% pull(sample) 
    keep <- rowSums(gencode_quant[,samples]) != 0
    tissue_quant <- gencode_quant %>% filter(keep) %>% 
        mutate(avg_quant  = rowMeans(.[,samples])) %>% 
        inner_join(ref_tx_lengths) %>% 
        select(transcript_id, avg_quant, length)
    return(cor(tissue_quant$avg_quant, tissue_quant$length, method = 'spearman') )
    
}
res <- sapply(subtissues, calc_txome_size)
NUM_AVG_EXP_TX_LENGTH_COR <- mean(res)

save(list= ls()[grepl('NUM_', ls())], file = files$paper_numbers_rdata)

