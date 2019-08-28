library(tidyverse)
# setwd('~/NIH/dev_eyeintegration_splicing/')
# gtf_file <- 'data/gtfs/all_tissues.combined.gtf'
# sample_table_file <- 'sampleTableDev.tsv'
args <- commandArgs(trailingOnly = T)
wd <- args[1]
gtf_file <- args[2]
sample_table_file <- args[3]
outfile <- args[4]
setwd(wd)
gtf <- rtracklayer::readGFF(gtf_file)
sample_table <- read_tsv(sample_table_file, col_names = c('sample', 'run', 'paired', 'tissue', 'subtissue', 'origin'))
overall_stats <- tibble(description=c('number of tissues', 'number of samples', 'number of novel transcripts', 'number of novel loci'),
       value=c(sample_table %>% pull(subtissue) %>% unique %>% length ,
               sample_table %>% pull(sample) %>% unique %>% length , 
               gtf %>% filter(!grepl('TCONS', gene_name), type == 'transcript', grepl('MSTRG', oId)) %>% nrow , 
               gtf %>% filter(grepl('TCONS', gene_name), type == 'transcript') %>% nrow
              )
      )
save(overall_stats, file=outfile)
















