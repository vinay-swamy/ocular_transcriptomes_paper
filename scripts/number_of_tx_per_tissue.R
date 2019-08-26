library(tidyverse)
setwd('~/NIH/dev_eyeintegration_splicing/')
gtf_file <- 'data/gtfs/all_tissues.combined.gtf'
sample_table_file <- 'sampleTableDev.tsv'
pep_info_file <- 'data/seqs/pep_fasta_meta_info.tsv' 
pep_info <- read_tsv(pep_info_file, col_names = c('trdID', 'transcript_id','orf', 'type', 'length', 'score', 'misc')) %>% 
  mutate(transcript_id= str_split(transcript_id, '.p') %>% sapply(function(x) x[[1]]))
tcons2mstrg <- 'data/misc/gfc_TCONS_to_st_MSTRG.tsv' %>% read_tsv 

gtf <- rtracklayer::readGFF(gtf_file) %>% mutate(is.pc= transcript_id %in% pep_info$transcript_id) 

count_and_class_tx <- function(tx_list){
  tissues <- colnames(tcons2mstrg)[-1]
  pc <- filter(gtf, is.pc, type == 'transcript', transcript_id %in% tx_list) %>% pull(transcript_id)
  nc <- filter(gtf, !is.pc, type == 'transcript', transcript_id %in% tx_list) %>% pull(transcript_id)
  nv_to_PC <- filter(tcons2mstrg, transcript_id %in% pc ) %>% 
    {apply(.[,-1], 2 ,function(tissue) sum(!is.na(tissue))) } %>% {tibble(tissue=names(.), type='protein coding', count=.)} 
  nv_to_NC <- filter(tcons2mstrg, transcript_id %in% nc ) %>% 
    {apply(.[,-1], 2 ,function(tissue) sum(!is.na(tissue))) } %>% {tibble(tissue=names(.),type='non-coding', count=.)} 
  return(rbind(nv_to_PC, nv_to_NC))
}

novel_transcripts <- gtf %>% filter(!grepl('TCONS', gene_name), type == 'transcript', grepl('MSTRG', oId)) %>% 
  pull(transcript_id)
novel_loci <- gtf %>% filter(grepl('TCONS', gene_name), type == 'transcript') %>% pull(transcript_id)
novel_transcripts_by_tissue <- count_and_class_tx(novel_transcripts) 
novel_loci_by_tissue <- count_and_class_tx(novel_loci)
save(novel_loci_by_tissue, novel_transcripts_by_tissue, file='~/NIH/occular_transcriptomes_paper/data/novel_tx_by_tissue.Rdata')
