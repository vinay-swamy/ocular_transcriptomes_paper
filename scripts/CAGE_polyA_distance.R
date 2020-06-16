library(tidyverse)
library(RBedtools)
library(yaml)


files_yaml <- '/data/swamyvs/ocular_transcriptomes_paper/files.yaml'
data_dir <- '/data/swamyvs/ocular_transcriptomes_pipeline/'

files <- read_yaml(files_yaml)
setwd(data_dir)
gencode_ref <- rtracklayer::readGFF(files$anno_gtf)
gencode_starts_bed <- gencode_ref %>% 
    filter(type == 'transcript') %>% 
    select(seqid, strand, start) %>% 
    distinct %>% 
    mutate(label = 'gencode',
           score = 333, 
           end = start+1) %>% select(seqid,start, end ,label, score, strand) %>%
    from_data_frame %>% RBedtools('sort', i=.)

dntx_gtf <- rtracklayer::readGFF(files$ref_gtf)
all_dntx_starts_bed <-dntx_gtf %>% 
    filter(type == "transcript") %>%
    select(seqid, strand, start) %>% 
    distinct %>% 
    mutate(label = 'dntx',
           score = 333, 
           end = start+1) %>% select(seqid,start, end ,label, score, strand) %>%
    from_data_frame %>% RBedtools('sort', i=.)

cage_normal= RBedtools('sort', i=files$CAGE_annotation)
gencode_CAGE_phase12 <-  RBedtools('closest', options= '-s -D b -t first', a=gencode_starts_bed, b=cage_normal ) %>% to_data_frame()
dntx_CAGE_phase12 <- RBedtools('closest',  options= '-s -D b -t first', a=all_dntx_starts_bed, b=cage_normal) %>% to_data_frame

all_CAGE_phase12 <- bind_rows( gencode_CAGE_phase12%>% select(build = X4, dist = X16),
                                dntx_CAGE_phase12 %>% select(build = X4, dist = X16)) %>% 
    mutate(dist=abs(dist))

gencode_ends_bed <- gencode_ref %>% filter(type == "transcript") %>% 
    select(seqid, strand, start=end) %>%
    distinct %>% 
    mutate(label = 'gencode', score = 34, end= start+1) %>% 
    select(seqid, start, end, label, score, strand) %>% 
    from_data_frame %>% 
    RBedtools('sort', i=.)
dntx_ends_bed <-  dntx_gtf %>% filter(type == "transcript") %>% 
    select(seqid, strand, start=end) %>%
    distinct %>% 
    mutate(label = 'dntx', score = 34, end= start+1) %>% 
    select(seqid, start, end, label, score, strand) %>% 
    from_data_frame %>% 
    RBedtools('sort', i=.)

gencode_polya_closest <- RBedtools('closest', '-s -D b -t first', 
                                   a=gencode_ends_bed, b=files$polyA_annotation) %>%
    to_data_frame %>% 
    select(label = X4, dist = X18) %>%
    mutate(label = 'gencode', dist=abs(dist))
dntx_polya_closest <- RBedtools('closest', '-s -D b -t first', 
                                a=dntx_ends_bed, b=files$polyA_annotation) %>% 
    to_data_frame %>% 
    select(label = X4, dist = X18) %>% 
    mutate(label = 'dntx', dist = abs(dist))

all_polya_closest <- bind_rows(gencode_polya_closest, 
                               dntx_polya_closest)

save(all_CAGE_phase12, all_polya_closest, file =files$CAGE_polyA_rdata)

#files$CAGE_polyA_rdata <- '/data/swamyvs/ocular_transcriptomes_paper/clean_data/rdata/CAGE_poly_dist.Rdata'
