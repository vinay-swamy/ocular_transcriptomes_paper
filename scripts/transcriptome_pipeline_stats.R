library(tidyverse)
library(matrixStats)
source('~/scripts/read_salmon.R')
args <- c('/Volumes/data/occular_transcriptomes_paper/',
          '/Volumes/data/eyeintegration_splicing/',
          'sampleTableFull.tsv',
          'data/gtfs/all_tissues.combined_NovelAno.gtf',
          'data/salmon_quant/RPE_Fetal.Tissue/',
          'data/misc/TCONS2MSTRG.tsv',
          'RPE_Fetal.Tissue',
          '/Volumes/data/occular_transcriptomes_paper/clean_data/DNTX_salmon_mapping_rates.tab',
          '/Volumes/data/occular_transcriptomes_paper/clean_data/gencode_salmon_mapping_rates.tab',
          '/Volumes/data/occular_transcriptomes_paper/clean_data/all_gtfs_tx_counts.tab',
          '/Volumes/data/occular_transcriptomes_paper/clean_data/transcriptome_pipeline_stats.Rdata'
          )
args <- commandArgs(trailingOnly = T)
save(args, file='tmp/tps.args')
working_dir <- args[1]
data_dir <- args[2]
sample_table_file <- args[3]
ano_gtf_file <-  args[4]
path_to_bs <- args[5]
tcons2mstrg_file <- args[6]
tissue <- args[7]
gencode_mapping_rate_file <- args[8]
DNTX_mapping_rate_file <- args[9]
gtf_trasncript_count_file <- args[10]
clean_data <- args[11]

setwd(data_dir)
proc_boostrap <- function(file){
    BS_raw <- read_tsv(file, col_names = F)
    name <- str_split(file,'/') %>%.[[1]] %>% {.[(which(grepl('quant_bootstrap',.)) -1)]} 
    res <- tibble(transcript_id=BS_raw$X1, !!name := BS_raw[,-1] %>% as.matrix %>% {matrixStats::rowVars(.)} )
    return(res)
}
var_sum <- function(all_bs, t_tissue){
    novel_bs <- filter(all_bs, grepl(t_tissue, transcript_id))
    ref_bs <- filter(all_bs, !grepl(t_tissue, transcript_id) )
    novel_medvar <- novel_bs[,-1] %>% as.matrix() %>% rowMedians()
    ref_medvar <- ref_bs[,-1] %>% as.matrix() %>% rowMedians()
    return(list(novel_median_variance=novel_medvar, ref_median_variance=ref_medvar))
}
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
bs_files <- list.files(path = path_to_bs, pattern = 'quant_bootstraps.tsv.gz', recursive = T, full.names = T)
all_bs <- lapply(bs_files, proc_boostrap) %>% reduce(left_join)
median_bootstrap_variance <- var_sum(all_bs, t_tissue = tissue)

#---- 


#COUNT GTF TX COUNTS
#----




tx_counts <- read_delim(gtf_trasncript_count_file, ' ' , col_names=c('raw', 'tx_count')) %>% 
    mutate(build= case_when( grepl('final_gtf',raw) ~ 'final',
                             grepl('compfilt',raw) ~ 'compfilt',
                             grepl('combined',raw) ~ 'raw',
                             grepl('gfcfilt',raw) ~ 'gfc_filt'),
           subtissue= raw %>% 
                        str_remove("final_gtf|\\.compfilt|\\.combined|\\.gfcfilt") %>% 
                        str_split('/') %>% 
                        sapply(function(x) x[grep('\\.gtf', x) ] %>% str_remove('\\.gtf|_st\\.gtf')),
           raw=NULL) %>% 
        spread(build, tx_count) %>%
    select(subtissue, raw, gfc_filt, compfilt, final)
#----


# COMPARE SALMON MAPPING RATES
#----

nsamp_by_tissue <- sample_table %>% group_by(subtissue) %>% summarise(nsamp=n())

DNTX_mapping_rates <- process_lib_size_tabs(DNTX_mapping_rate_file , 'DNTX')
gencode_mapping_rates <- process_lib_size_tabs(gencode_mapping_rate_file, 'gencode')
sample_mapping_rates <- inner_join(DNTX_mapping_rates, gencode_mapping_rates) %>% 
    inner_join(select(sample_table, sample, subtissue)) %>% 
    mutate(map_diff=DNTX_percent_mapped - gencode_percent_mapped) %>% 
    group_by(subtissue) %>% 
    summarise(med_diff=median(map_diff), med_libsize=median(gencode_total_reads) ) %>% arrange(desc(med_diff)) %>% 
    mutate(med_num_reads_gained = med_diff * med_libsize) %>% 
    left_join(tx_counts) %>% left_join(nsamp_by_tissue)


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
save(median_bootstrap_variance, sample_mapping_rates, tx_counts,  file = clean_data)







