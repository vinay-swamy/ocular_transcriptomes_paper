library(tidyverse)
library(matrixStats)
library(argparse)
source('~/scripts/read_salmon.R')
# args <- c('/Volumes/data/ocular_transcriptomes_paper/',
#           '/Volumes/data/ocular_transcriptomes_pipeline/',
#           'sampleTableFull.tsv',
#           'data/gtfs/all_tissues.combined_NovelAno.gtf',
#           'data/salmon_quant/RPE_Fetal.Tissue/',
#           'data/misc/TCONS2MSTRG.tsv',
#           'RPE_Fetal.Tissue',
#           '/Volumes/data/ocular_transcriptomes_paper/clean_data/DNTX_salmon_mapping_rates.tab',
#           '/Volumes/data/ocular_transcriptomes_paper/clean_data/gencode_salmon_mapping_rates.tab',
#           '/Volumes/data/ocular_transcriptomes_paper/clean_data/all_gtfs_tx_counts.tab',
#           '/Volumes/data/ocular_transcriptomes_paper/clean_data/transcriptome_pipeline_stats.Rdata'
#           )

parser <- ArgumentParser()
parser$add_argument('--workingDir', action='store', dest='working_dir')
parser$add_argument('--dataDir', action = 'store', dest = 'data_dir')
parser$add_argument('--sampleTableFile', action = 'store', dest = 'sample_table_file')
parser$add_argument('--anoGtfFile', action = 'store', dest='ano_gtf_file')
parser$add_argument('--pathToBs', action = 'store', dest = 'path_to_bs')
parser$add_argument('--tcons2mstrgFile', action = 'store', dest = 'tcons2mstrg_file')
parser$add_argument('--tissue', action='store', dest = 'tissue')
parser$add_argument('--gencodeMapRateFile', action = 'store', dest = 'gencode_mapping_rate_file')
parser$add_argument('--DNTXMapRateFile',action = 'store', dest = 'DNTX_mapping_rate_file')
parser$add_argument('--gtfTxCountFile', action = 'store', dest = 'gtf_trasncript_count_file')
parser$add_argument('--cleanData', action = 'store', dest = 'clean_data')
list2env(parser$parse_args(), .GlobalEnv)
save.image('testing/tps_args.Rdata')
# working_dir <- args[1]
# data_dir <- args[2]
# sample_table_file <- args[3]
# ano_gtf_file <-  args[4]
# path_to_bs <- args[5]
# tcons2mstrg_file <- args[6]
# tissue <- args[7]
# gencode_mapping_rate_file <- args[8]
# DNTX_mapping_rate_file <- args[9]
# gtf_trasncript_count_file <- args[10]
# clean_data <- args[11]

setwd(data_dir)
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
sample_table <- read_tsv(sample_table_file) %>% filter(subtissue != 'synth')
t2m <- read_tsv(tcons2mstrg_file)
#PROCESS BOOTSTRAPS
#----


#---- 


#COUNT GTF TX COUNTS
#----




tx_counts <- read_delim(gtf_trasncript_count_file, ' ' , col_names=c('raw', 'tx_count')) %>% 
    filter(!grepl('filtered', raw)) %>% 
    mutate(build=ifelse(grepl('final',raw), 'filtered','base'),
           subtissue= raw %>% 
                        str_remove("\\.combined") %>% 
                        str_split('/') %>% 
                        sapply(function(x) x[grep('\\.gtf', x) ] %>% str_remove('\\.gtf|_st\\.gtf')),
           raw=NULL)%>% 
        spread(build, tx_count)# %>%
    select(subtissue, raw, gfc_filt, compfilt, final)
#----


# COMPARE SALMON MAPPING RATES
#----

nsamp_by_tissue <- sample_table %>% group_by(subtissue) %>% summarise(nsamp=n())

DNTX_mapping_rates <- process_lib_size_tabs(DNTX_mapping_rate_file , 'DNTX')
gencode_mapping_rates <- process_lib_size_tabs(gencode_mapping_rate_file, 'gencode')
median_sample_mapping_rates <- inner_join(DNTX_mapping_rates, gencode_mapping_rates) %>% 
    inner_join(select(sample_table, sample, subtissue)) %>% 
    mutate(map_diff=DNTX_percent_mapped - gencode_percent_mapped) %>% 
    group_by(subtissue) %>% 
    summarise(med_diff=median(map_diff), med_libsize=median(gencode_total_reads) ) %>% arrange(desc(med_diff)) %>% 
    mutate(med_num_reads_gained = med_diff * med_libsize) %>% 
    left_join(tx_counts) %>% left_join(nsamp_by_tissue)

all_sample_mapping_rates <- sample_mapping_rates <- inner_join(DNTX_mapping_rates, gencode_mapping_rates) %>% 
    inner_join(select(sample_table, sample, body_location)) %>% 
    select(sample, body_location, contains('percent')) %>% 
    gather(key='build',value = 'mapping_rate', -body_location, -sample) %>% 
    mutate(build=ifelse(build == 'DNTX_percent_mapped', 'DNTX', 'gencode'))

#save(all_sample_mapping_rates, file = '/Volumes/data/ocular_transcriptomes_paper/clean_data/rdata/all_tx_mappingrates.Rdata')

#----

# library(viridis)
# library(ComplexHeatmap)
# col <- viridis(100)
# gtf <- rtracklayer::readGFF(gtf_file)
# tctab <- t2m
# bool_tab <- apply(tctab[,-1], 2, function(x) as.numeric(!is.na(x))) %>% as.data.frame
# cor_tab <- cor(bool_tab,method = 'spearman')
# hm <- Heatmap(cor_tab, col = col, name = 'spearman correlation')
# construction_correlation_heatmap <- hm


#----
save(median_sample_mapping_rates,all_sample_mapping_rates , tx_counts,  file = clean_data)







