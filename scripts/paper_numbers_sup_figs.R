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

sample_table <- read_tsv(files$sample_table)
NUM_TOTAL_SAMPLES <- nrow(sample_table)
NUM_EYE_SAMPLES <- sample_table %>% filter(!body_location %in% c('Brain', 'Body')) %>% nrow 
NUM_BODY_SAMPLES <- sample_table %>% filter(body_location %in% c('Brain', 'Body')) %>% nrow 
NUM_BODY_TISSUES <- sample_table %>% filter(body_location %in% c('Brain', 'Body')) %>% 
    pull(subtissue) %>% unique %>% length 

pan_body_gtf <- rtracklayer::readGFF(files$anno_gtf)
pan_eye_gtf <- rtracklayer::readGFF(files$pan_eye_gtf)
NUM_REF_TX_BODY <- pan_body_gtf %>% filter( type == 'transcript', class_code=='=')%>% nrow 
NUM_REF_TX_EYE  <- pan_eye_gtf  %>% filter( type == 'transcript', class_code=='=')%>% nrow 

NUM_NOVEL_TX_BODY <- pan_body_gtf %>% filter( type == 'transcript', class_code!='=')%>% nrow 
NUM_NOVEL_TX_EYE  <- pan_eye_gtf  %>% filter( type == 'transcript', class_code!='=')%>% nrow 
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

DNTX_mapping_rates <- process_lib_size_tabs(files$DNTX_mr , 'DNTX')
gencode_mapping_rates <- process_lib_size_tabs(files$gencode_mr, 'gencode')

    
median_sample_mapping_rates <- inner_join(DNTX_mapping_rates, gencode_mapping_rates) %>% 
        inner_join(select(sample_table, sample, subtissue)) %>% 
        mutate(map_diff=DNTX_percent_mapped - gencode_percent_mapped) %>% 
        group_by(subtissue) %>% 
        summarise(med_diff=median(map_diff), med_libsize=median(gencode_total_reads) ) %>% arrange(desc(med_diff)) %>% 
        mutate(med_num_reads_gained = med_diff * med_libsize) #%>% 
        #left_join(tx_counts) %>% left_join(nsamp_by_tissue)

load(files$gencode_quant)
subtissues <- unique(sample_table$subtissue)
calc_txome_size <- function(t_tissue){
    samples <- filter(sample_table, subtissue == t_tissue) %>% pull(sample) 
    all_exp <- sum(rowSums(gencode_quant[,samples]) != 0)
    avg_1tpm <- sum(rowMeans(gencode_quant[,samples]) >=1)
    return(tibble(subtissue = t_tissue, all_exp = all_exp, avg_1tpm = avg_1tpm))
}

gencode_tx_size <- lapply(subtissues, calc_txome_size) %>% bind_rows()

## SUPPLEMENTARY FIGURES 
#---- 
load(files$core_tight_rdata)

save(list= ls()[grepl('NUM_', ls())], file = files$paper_numbers_rdata)
save(gencode_tx_size, file = files$sup_fig_data)
