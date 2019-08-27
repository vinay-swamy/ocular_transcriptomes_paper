library(tidyverse)
library(bedr)
load('rdata/novel_exon_classification.rdata')
gtf <- rtracklayer::readGFF('data/gtfs/all_tissues.combined.gtf')
novel_exon_bed <- novel_exons_TSES %>% mutate(seqid=as.character(seqid), score=999) %>% 
  select(seqid, start, end, id, score, strand)
ref_exon_bed <- all_exons %>% mutate(score= 999, seqid=as.character(seqid)) %>% 
  select(seqid, start, end, origin, score, strand)
k <- bedr(input = list(a=novel_exon_bed, b=ref_exon_bed), method = 'intersect', params = '-loj -s', check.chr = F, check.sort = F)
write_tsv(novel_exon_bed, 'testing/novel_exons.bed', col_names = F)
write_tsv(ref_exon_bed, 'testing/ref_exons.bed', col_names = F)
system2(command = 'bedtools', args = 'intersect -loj -s -a testing/novel_exons.bed -b testing/ref_exons.bed > testing/tout.bed')
res <- read_tsv('testing/tout.bed', col_names = F) %>% filter(X8 == -1) %>% .[,1:6]
total_exons_desc <- tibble(desc=c('Number of Alternative Exons', 'Number of Fully Novel Exons'),
                      value=c((nrow(novel_exons_TSES) - nrow(res)), nrow(res)))
rna_proc_type <-  novel_exons_TSES %>%
  mutate(broad_type=case_when(nv_type_rc == 'novel_TSS' ~ 'Novel Transcriptional Start Site',
                              nv_type_rc == 'novel_TES' ~ 'Novel Transcriptional Termination Site',
                                                      T ~ 'Alternative Splicing') ) %>% 
  pull(broad_type) %>% table %>% {tibble(desc=names(.), count=.)}
save(total_exons_desc, rna_proc_type, file = '~/NIH/occular_transcriptomes_paper/data/exon_classification.Rdata')

                      