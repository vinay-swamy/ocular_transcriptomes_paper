library(tidyverse)
library(matrixStats)
library(yaml)
library(argparse)
source('~/scripts/read_salmon.R')

parser <- ArgumentParser()
parser$add_argument('--workingDir', action='store', dest='working_dir')
parser$add_argument('--dataDir', action = 'store', dest = 'data_dir')
parser$add_argument('--filesYaml', action = 'store', dest = 'files_yaml')
#######
# working_dir <- '/data/swamyvs/ocular_transcriptomes_paper/'
# data_dir <- '/data/swamyvs/ocular_transcriptomes_pipeline/'
# files_yaml <- '/data/swamyvs/ocular_transcriptomes_paper/files.yaml'
#######
list2env(parser$parse_args(), .GlobalEnv)
files <- read_yaml(files_yaml)

setwd(data_dir)
process_lib_size_tabs <- function(file, type){
    df_messy <- read.delim(file, sep = ' ', header = F) 
    sample <- df_messy$V1 %>% 
        as.character %>%  
        str_split('/') %>% 
        sapply( function (x){idx=(which(grepl('aux_info', x)) -1);return(x[idx])} )
    total_reads <-  df_messy$V3 %>% 
        str_split(',') %>% 
        sapply(function(x)x[1]) %>%
        as.numeric 
    percent_mapped <- df_messy$V5 %>% 
        str_split(',') %>% 
        sapply(function(x)x[1]) %>% 
        as.numeric %>% {./100}
    df <- tibble(sample, total_reads, percent_mapped)
    colnames(df) <- c('sample', paste0(type, c('_total_reads', '_percent_mapped')))
    return(df)
}
sample_table <- read_tsv(files$sample_table) %>% filter(subtissue != 'synth')
t2m <- read_tsv(files$tcons2mstrg)
head(t2m)
#COUNT GTF TX COUNTS
#----




tx_counts <- read_delim(files$gtf_tx_counts, ' ' , col_names=c('raw', 'tx_count')) %>% 
    filter(!grepl('filtered', raw)) %>% 
    mutate(build=ifelse(grepl('final',raw), 'filtered','base'),
           subtissue= raw %>% 
                        str_remove("\\.combined") %>% 
                        str_split('/') %>% 
                        sapply(function(x) x[grep('\\.gtf', x) ] %>% str_remove('\\.gtf|_st\\.gtf')),
           raw=NULL)%>% 
        spread(build, tx_count) #%>%
    #select(subtissue, raw, gfc_filt, compfilt, final)
base_enst_tx_counts <- read_delim(files$base_gtf_enst_tx_counts, ' ' , col_names=c('raw', 'tx_count')) %>% 
    filter(!grepl('filtered', raw)) %>% 
    mutate(build=ifelse(grepl('final',raw), 'filtered','base'),
           subtissue= raw %>% 
               str_remove("\\.combined") %>% 
               str_split('/') %>% 
               sapply(function(x) x[grep('\\.gtf', x) ] %>% str_remove('\\.gtf|_st\\.gtf')),
           raw=NULL)%>% 
    spread(build, tx_count)
#----


# COMPARE SALMON MAPPING RATES
#----

nsamp_by_tissue <- sample_table %>% group_by(subtissue) %>% summarise(nsamp=n())

DNTX_mapping_rates <- process_lib_size_tabs(files$DNTX_mr , 'DNTX')
gencode_mapping_rates <- process_lib_size_tabs(files$gencode_mr, 'gencode')
all_sample_mapping_rate_difference <- inner_join(DNTX_mapping_rates, gencode_mapping_rates) %>% 
    mutate(mapping_rate_diff = DNTX_percent_mapped - gencode_percent_mapped) %>% 
    inner_join(select(sample_table, sample, body_location))

all_sample_mapping_rates <- sample_mapping_rates <- inner_join(DNTX_mapping_rates, gencode_mapping_rates) %>% 
    inner_join(select(sample_table, sample, body_location)) %>% 
    select(sample, body_location, contains('percent')) %>% 
    gather(key='build',value = 'mapping_rate', -body_location, -sample) %>% 
    mutate(build=ifelse(build == 'DNTX_percent_mapped', 'DNTX', 'gencode'))

load(files$gencode_quant)
subtissues <- unique(sample_table$subtissue)
calc_txome_size <- function(t_tissue){
    samples <- filter(sample_table, subtissue == t_tissue) %>% pull(sample) 
    all_exp <- sum(rowSums(gencode_quant[,samples]) != 0)
    avg_1tpm <- sum(rowMeans(gencode_quant[,samples]) >=1)
    return(tibble(subtissue = t_tissue, all_exp = all_exp, avg_1tpm = avg_1tpm))
}

gencode_tx_size <- lapply(subtissues, calc_txome_size) %>% bind_rows()



save(gencode_tx_size, all_sample_mapping_rates,all_sample_mapping_rate_difference, base_enst_tx_counts, tx_counts,  file = files$txome_stats_rdata)







