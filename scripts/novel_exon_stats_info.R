library(tidyverse)

args <- commandArgs(trailingOnly = T)
wd <- args[1]
exon_classifcation_file <- args[2]
gtf_file <- args[3]
outfile <- args[4]
load(exon_classifcation_file)
setwd(wd)
gtf <- rtracklayer::readGFF(gtf_file)
novel_exon_bed <- novel_exons_TSES %>% mutate(seqid=as.character(seqid), score=999) %>% 
  select(seqid, start, end, id, score, strand)
ref_exon_bed <- all_exons %>% mutate(score= 999, seqid=as.character(seqid)) %>% 
  select(seqid, start, end, origin, score, strand)
k <- bedr(input = list(a=novel_exon_bed, b=ref_exon_bed), method = 'intersect', params = '-loj -s', check.chr = F, check.sort = F)
write_tsv(novel_exon_bed, '/tmp/novel_exons.bed', col_names = F)
write_tsv(ref_exon_bed, '/tmp/ref_exons.bed', col_names = F)
system2(command = 'bedtools', args = 'intersect -loj -s -a /tmp/novel_exons.bed -b /tmp/ref_exons.bed > /tmp/tout.bed')
res <- read_tsv('/tmp/tout.bed', col_names = F) %>% filter(X8 == -1) %>% .[,1:6]
total_exons_desc <- tibble(desc=c('Number of Alternative Exons', 'Number of Fully Novel Exons'),
                      value=c((nrow(novel_exons_TSES) - nrow(res)), nrow(res)))
rna_proc_type <-  novel_exons_TSES %>%
  mutate(broad_type=case_when(nv_type_rc == 'novel_TSS' ~ 'Novel Transcriptional Start Site',
                              nv_type_rc == 'novel_TES' ~ 'Novel Transcriptional Termination Site',
                                                      T ~ 'Alternative Splicing') ) %>% 
  pull(broad_type) %>% table %>% {tibble(desc=names(.), count=.)}
save(total_exons_desc, rna_proc_type, file = outfile)

                      