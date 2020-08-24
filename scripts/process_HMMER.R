library(tidyverse)
library(data.table)
library(argparse)
library(glue)
library(yaml)

parser <- ArgumentParser()
parser$add_argument('--dataDir', action  = 'store', dest = 'data_dir')
parser$add_argument('--filesYaml', action = 'store', dest = 'file_yaml')
list2env(parser$parse_args(), .GlobalEnv)
files = read_yaml(file_yaml)
hmmer_dir = files$hmmer_dir
setwd(data_dir)
eye_tissues <- c("RPE_Fetal.Tissue",  "Retina_Adult.Tissue",  "RPE_Adult.Tissue", "Cornea_Adult.Tissue", "Cornea_Fetal.Tissue", "Retina_Fetal.Tissue")
cn=c('target_name', 'target_accession','query_name', 'query_accession','evalue', 'total_score', 'total_bias','dom_val', 'dom_score', 'dom_bias', 'exp','reg', 'clu','ov','env', 'dom','rep','inc',  'description')
seq_hits <- read_delim(glue('{hmmer_dir}/seq_hits.tsv'),skip = 3, delim=' ' , col_names = cn) %>%
    mutate(query_name=query_name %>% strsplit('\\.') %>% sapply(function(x) x[1]), 
           query_name = str_remove_all(query_name,' ')) %>% 
    filter(target_name != '#') %>% 
    mutate_if(!colnames(.)%in%c('target_name', 'target_accession','query_name', 'query_accession','description' ), as.numeric ) %>% 
    filter(evalue< .01)

conv_tab_all <- fread(files$tcons2mstrg) %>% as_tibble %>% 
    {bind_cols(.[,1:4], apply(.[,-(1:4)], 2, function(x) x!='') %>% as_tibble)} %>% 
    rename(query_name = transcript_id) %>% as_tibble 
conv_tab_eye <- conv_tab_all %>% filter(rowSums(.[,eye_tissues]) >0) %>% select(query_name, eye_tissues)
seq_hits_eye <- inner_join(seq_hits,conv_tab_eye)

domain_cn=c('target_name', 'target_accession' ,'tlen' ,'query_name', 'query_accession' ,'qlen' ,'evalue', 'total_score', 'total_bias', 
            'dom_num' ,  'total_doms' , 'c_evalue' , 'i_evalue' , 'domain_score' , 'domain_bias' , 'start_hmm_coord' ,'end_hmmcoord',
            'start_seq', 'end_seq' , 'from_env' ,'to_env',  'post_prob', 'description')
domain_hits <- read_delim(glue('{hmmer_dir}/domain_hits.tsv'),skip = 3, delim=' ' , col_names = domain_cn) %>%
    mutate(query_name=query_name %>% strsplit('\\.') %>% sapply(function(x) x[1]), 
           query_name = str_remove_all(query_name,' ')) %>% 
    filter(target_name != '#') %>% 
    mutate_if(!colnames(.)%in%c('target_name', 'target_accession','query_name', 'query_accession','description' ), as.numeric ) %>% 
    filter(!query_name %in% seq_hits_eye$query_name, evalue <.01 )
domain_hits_eye <- inner_join(domain_hits, conv_tab_eye)

save(domain_hits, domain_hits_eye, seq_hits, seq_hits_eye, file = files$hmmer_results)
