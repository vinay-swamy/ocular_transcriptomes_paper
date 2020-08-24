library(tidyverse)
library(RBedtools)
library(yaml)
library(glue)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--dataDir', action = 'store', dest = 'data_dir')
parser$add_argument('--filesYaml', action = 'store', dest = 'files_yaml')
# ####
# files_yaml <- '/data/swamyvs/ocular_transcriptomes_paper/files.yaml'
# data_dir <- '/data/swamyvs/ocular_transcriptomes_pipeline/'
# ####
list2env(parser$parse_args(), .GlobalEnv)

files <- read_yaml(files_yaml)
setwd(data_dir)
gencode_ref <- rtracklayer::readGFF(files$ref_gtf)
gencode_starts_bed <- gencode_ref %>% 
    filter(type == 'transcript') %>% 
    select(seqid, strand, start) %>% 
    distinct %>% 
    mutate(label = 'gencode',
           score = 333, 
           end = start+1) %>% select(seqid,start, end ,label, score, strand) %>%
    from_data_frame %>% RBedtools('sort', i=.)

dntx_gtf <- rtracklayer::readGFF(files$anno_gtf)
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
    mutate(abs_dist=abs(dist))


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
    mutate(label = 'gencode', abs_dist=abs(dist))
dntx_polya_closest <- RBedtools('closest', '-s -D b -t first', 
                                a=dntx_ends_bed, b=files$polyA_annotation) %>% 
    to_data_frame %>% 
    select(label = X4, dist = X18) %>% 
    mutate(label = 'dntx', abs_dist = abs(dist))

all_polya_closest <- bind_rows(gencode_polya_closest, 
                               dntx_polya_closest)

######################
######compare phylop scores 
###base comparison: gencode vs DNTX vs all
# remove scaffold chroms since they are not in phylop
gencode_exons <- gencode_ref %>% 
    filter(type == 'exon', grepl('chr', as.character(seqid))) %>% 
    mutate(score = 333) %>% 
    select(seqid, start, end, exon_id, score, strand) %>% 
    from_data_frame %>% 
    RBedtools('sort', i=.)

dntx_exons <- dntx_gtf %>% 
    filter(type == 'exon',grepl('chr', as.character(seqid))) %>%
    mutate(score = 444, 
           exon_id = paste0('DNEX_', 1:nrow(.))) %>% 
    select(seqid, start, end, exon_id, score,strand) %>% 
    from_data_frame %>% 
    RBedtools('sort', i=.)
system2('cat', glue('{unclass(gencode_exons)}  {unclass(dntx_exons)} > /tmp/gc_dntx_c.bed'))
dntx_gc_distinct <- RBedtools('sort', output = 'stdout', i='/tmp/gc_dntx_c.bed') %>% 
                    RBedtools('merge', options = '-s', i=.)

gencode_phylop_score <- RBedtools('intersect', options =  '-wa -wb -sorted ', output = 'stdout',
                                  a=gencode_exons, b=files$phylop_bed) %>% 
    RBedtools('groupby', options ='-g 1,2,3,4,5,6 -c 11 -o mean',i= .) %>%
    to_data_frame %>%
    select(exon_id=X4, mean_phylop_score=X7)

dntx_phylop_score <- RBedtools('intersect', options =  '-wa -wb -sorted ', output = 'stdout',
                               a=dntx_exons, b=files$phylop_bed) %>% 
    RBedtools('groupby', options ='-g 1,2,3,4,5,6 -c 11 -o mean',i= .) %>%
    to_data_frame %>%
    select(exon_id=X4, mean_phylop_score=X7)

all_phylop <- bind_rows(gencode_phylop_score %>% mutate(build = 'gencode'),
                        dntx_phylop_score %>% mutate(build = 'dntx'))
save(all_CAGE_phase12, all_polya_closest,all_phylop, file =files$CAGE_polyA_rdata)

#files$CAGE_polyA_rdata <- '/data/swamyvs/ocular_transcriptomes_paper/clean_data/rdata/CAGE_poly_dist.Rdata'
