library(tidyverse)
# setwd('~/NIH/dev_eyeintegration_splicing/')
# gtf_file <- 'data/gtfs/all_tissues.combined.gtf'
# sample_table_file <- 'sampleTableDev.tsv'
# pep_info_file <- 'data/seqs/pep_fasta_meta_info.tsv'
# exon_classification_file <- 'rdata/novel_exon_classification.rdata'
# transdecoder_gff_file <- 'testing/transcripts.fa.transdecoder.genome.gff3'

args <- commandArgs(trailingOnly = T)
wd <- args[1]
gtf_file <- args[2]
sample_table_file <- args[3]
pep_info_file <- args[4]
exon_classification_file <- args[5]
transdecoder_gff_file <- args[6]
tcons2mstrg_file <- args[7]
outfile <- args[8]

setwd(wd)
load(exon_classification_file)
cds_gff <- rtracklayer::readGFF(transdecoder_gff_file) %>% as.data.frame %>% mutate(ID=str_extract(ID,'TCONS_[0-9]+|ENSG[0-9]+'))
pep_info <- read_tsv(pep_info_file, col_names = c('trdID', 'transcript_id','orf', 'type', 'length', 'score', 'misc')) %>%
  mutate(transcript_id= str_split(transcript_id, '.p') %>% sapply(function(x) x[[1]]))
tcons2mstrg <-  read_tsv(tcons2mstrg_file)

# now we need to classify these exons as either in the coding region of a transcript,
# a non coding region of a coding transcript
# or a noncoding transcript
gtf <- rtracklayer::readGFF(gtf_file) %>% mutate(is.pc= transcript_id %in% cds_gff$ID)
gtf_exons <- filter(gtf, type == 'exon') %>% select(seqid, strand, start, end , transcript_id, is.pc)
novel_exons_txid <- inner_join(novel_exons_TSES, gtf_exons)
cds_df <- cds_gff %>% filter(type == 'CDS') %>% select(seqid, strand, start, end) %>% distinct %>% mutate(is.CDS=T)
novel_exons_PC <- novel_exons_txid %>% left_join(cds_df) %>% mutate(is.CDS=replace_na(is.CDS, F))
# a hacky way to count cds exons


count_and_class_tx <- function(tissue){
    txs <- tcons2mstrg[,c('transcript_id', tissue)] %>% filter(!is.na(.[,2])) %>% pull(transcript_id)
    novel_exons_tissue <- filter(novel_exons_PC, transcript_id %in% txs) %>% select(transcript_id, is.pc, is.CDS) %>% distinct
    num_noncoding <- sum(!novel_exons_tissue$is.pc)
    num_nocodingchange <- filter(novel_exons_tissue, is.pc) %>% pull(is.CDS) %>%{sum(!.)}
    num_codingchange <- filter(novel_exons_tissue, is.pc) %>% pull(is.CDS) %>%{sum(.)}
    return(tibble(tissue=rep(tissue, 3),type=c('Non Coding', 'No Change', 'Coding Change'),
                  count=c(num_noncoding, num_nocodingchange, num_codingchange) ))

}
count_and_class_loci <- function(tx_list){
  tissues <- colnames(tcons2mstrg)[-1]
  pc <- filter(gtf, is.pc, type == 'transcript', transcript_id %in% tx_list) %>% pull(transcript_id)
  nc <- filter(gtf, !is.pc, type == 'transcript', transcript_id %in% tx_list) %>% pull(transcript_id)
  nv_to_PC <- filter(tcons2mstrg, transcript_id %in% pc ) %>%
    {apply(.[,-1], 2 ,function(tissue) sum(!is.na(tissue))) } %>% {tibble(tissue=names(.), type='protein coding', count=.)}
  nv_to_NC <- filter(tcons2mstrg, transcript_id %in% nc ) %>%
    {apply(.[,-1], 2 ,function(tissue) sum(!is.na(tissue))) } %>% {tibble(tissue=names(.),type='non-coding', count=.)}
  return(rbind(nv_to_PC, nv_to_NC))
}


novel_transcripts_by_tissue <- lapply(colnames(tcons2mstrg)[-1], count_and_class_tx) %>% bind_rows()
novel_loci <- novel_loci_distinct$transcript_id
novel_loci_by_tissue <- count_and_class_loci(novel_loci)
save(novel_loci_by_tissue, novel_transcripts_by_tissue, file=outfile)
