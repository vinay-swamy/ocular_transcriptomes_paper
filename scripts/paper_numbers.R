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
library(argparse)
args <- list(
    working_dir= '/Volumes/data/ocular_transcriptomes_paper/', 
    data_dir= '/Volumes/data/ocular_transcriptomes_pipeline/',
    sample_table_file= '/Volumes/data/ocular_transcriptomes_pipeline/sampleTableFull.tsv', 
    pan_body_gtf_file= '/Volumes/data/ocular_transcriptomes_pipeline/data/gtfs/all_tissues.combined.gtf',
    pan_eye_gtf_file= '/Volumes/data/ocular_transcriptomes_paper/clean_data/pan_eye_txome.combined.gtf',
    path_to_raw_tissues_gtfs= '/Volumes/data/ocular_transcriptomes_pipeline/data/gtfs/raw_tissue_gtfs/',
    all_ref_ano= '/Volumes/data/ocular_transcriptomes_pipeline/rdata/all_ref_tx_exons.rdata',
    core_tight_file= '/Volumes/data/ocular_transcriptomes_pipeline/ref/core_tight.Rdata',
    out_num_file=  '/Volumes/data/ocular_transcriptomes_paper/clean_data/rdata/paper_numbers.Rdata')
list2env(args, .GlobalEnv)

parser <- ArgumentParser()
parser$add_argument('--workingDir', action = 'store', dest = 'working_dir')
parser$add_argument('--dataDir', action = 'store', dest = 'data_dir')
parser$add_argument('--sampleTableFile', action = 'store', dest = 'sample_table_file')
parser$add_argument('--allTissueGtf', action = 'store', dest = 'pan_body_gtf_file')
parser$add_argument('--EyeOnlyGtf', action = 'store', dest = 'pan_eye_gtf_file')
parser$add_argument('--pathToRawGtfs', action = 'store', dest = 'path_to_raw_tissues_gtf')
parser$add_argument('--allRefAno', action = 'store', dest = 'all_ref_ano')
parser$add_argument('--coreTight', action = 'store', dest = 'core_tight_file')
parser$add_argument('--outNumFile', action = 'store', dest = 'out_num_file')
list2env(parser$parse_args(), .GlobalEnv)


setwd(working_dir)


sample_table <- read_tsv(sample_table_file)
NUM_TOTAL_SAMPLES <- nrow(sample_table)
NUM_EYE_SAMPLES <- sample_table %>% filter(!body_location %in% c('Brain', 'Body')) %>% nrow 
NUM_BODY_SAMPLES <- sample_table %>% filter(body_location %in% c('Brain', 'Body')) %>% nrow 
NUM_BODY_TISSUES <- sample_table %>% filter(body_location %in% c('Brain', 'Body')) %>% 
    pull(subtissue) %>% unique %>% length 

pan_body_gtf <- rtracklayer::readGFF(pan_body_gtf_file)
pan_eye_gtf <- rtracklayer::readGFF(pan_eye_gtf_file)
NUM_TOTAL_TX_BODY <- pan_body_gtf %>% filter(type == 'transcript') %>% nrow 
NUM_TOTAL_TX_EYE <- pan_eye_gtf %>% filter(type == 'transcript') %>% nrow 

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

process_lib_size_tabs <- function(file, type){
    df_messy <- read.delim(file, sep = ' ', header = F) 
    sample <- df_messy$V1 %>% as.character %>%  str_split('/') %>% 
        sapply( function (x){idx=(which(grepl('aux_info', x)) -1);return(x[idx])} )
    total_reads <-  df_messy$V3 %>% str_split(',') %>% sapply(function(x)x[1]) %>% as.numeric 
    percent_mapped <- df_messy$V5 %>% str_split(',') %>% sapply(function(x)x[1]) %>% as.numeric %>% {./100}
    df <- tibble(sample, total_reads, percent_mapped)
    colnames(df) <- c('sample', paste0(type, c('_total_reads', '_percent_mapped')))
    return(df)
}

DNTX_mapping_rates <- process_lib_size_tabs(DNTX_mapping_rate_file , 'DNTX')
gencode_mapping_rates <- process_lib_size_tabs(gencode_mapping_rate_file, 'gencode')

inner_join(DNTX_mapping_rates, gencode_mapping_rates) %>%
    mutate(diff = (DNTX_total_reads*DNTX_percent_mapped) - (gencode_total_reads*gencode_percent_mapped) ) %>% 
    View 
    {sum(. < 0)}
    pull(diff) %>% sum
    
    median_sample_mapping_rates <- inner_join(DNTX_mapping_rates, gencode_mapping_rates) %>% 
        inner_join(select(sample_table, sample, subtissue)) %>% 
        mutate(map_diff=DNTX_percent_mapped - gencode_percent_mapped) %>% 
        group_by(subtissue) %>% 
        summarise(med_diff=median(map_diff), med_libsize=median(gencode_total_reads) ) %>% arrange(desc(med_diff)) %>% 
        mutate(med_num_reads_gained = med_diff * med_libsize) #%>% 
        left_join(tx_counts) %>% left_join(nsamp_by_tissue)
load(core_tight_file)
NUM_STUDIES <- core_tight %>% filter(sample_accession %in% sample_table$sample) %>% pull(study_accession) %>% unique %>% length 

save(list= ls()[grepl('NUM_', ls())], file = out_num_file)

