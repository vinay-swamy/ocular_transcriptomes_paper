library(tidyverse)
library(data.table)
library(matrixStats)
library(argparse)
library(yaml)
read_VEP <- function(file, tissue){
    fread(file, sep = '\t', skip = '##') %>% 
        as_tibble %>% rename(allele_id = `#Uploaded_variation`) %>% 
        filter(allele_id %in% VUS_ids) %>% 
        select(allele_id,Feature,Consequence, Extra) %>% 
        mutate(IMPACT = str_split(Extra, '=|;') %>% sapply(function(x) x[2]),
               Consequence= str_split(Consequence, ',')) %>% 
        unnest(Consequence) %>% 
        inner_join(clinsig_rank) %>% 
        select(allele_id,impact_code) %>% 
        distinct %>% 
        group_by(allele_id) %>% #select the max variant for each allele:transcript pair 
        summarise(!!tissue := max(impact_code))
}
parser <- ArgumentParser()
parser$add_argument('--dataDir', action = 'store', dest = 'data_dir')
parser$add_argument('--fileYaml', action = 'store', dest = 'file_yaml')
####
data_dir <- '/data/swamyvs/ocular_transcriptomes_pipeline/'
files_yaml <- '/data/swamyvs/ocular_transcriptomes_paper/files.yaml'
####
list2env(parser$parse_args(), .GlobalEnv)
files <- read_yaml(files_yaml)




setwd(data_dir)
library(tidyverse)
library(data.table)
library(RBedtools)
library(matrixStats)
library(ComplexHeatmap)
library(RColorBrewer)

clinvar_anno <- fread('ref/clinvar_variant_summary.txt.gz', sep = '\t') %>% 
    as_tibble %>% rename(allele_id = `#AlleleID`)

VUS_ids <- clinvar_anno %>% filter(ClinicalSignificance == "Uncertain significance") %>% pull(allele_id) %>% unique
clinsig_rank <- tibble(IMPACT = c('HIGH', 'MODERATE', 'LOW', "MODIFIER")) %>%  
    mutate(impact_code  = case_when(IMPACT == 'HIGH' ~ 3, 
                                    IMPACT == 'MODERATE' ~2, 
                                    TRUE ~ 1) )

sample_table <- read_tsv('sampleTableFullV3.tsv')
subtissues <- unique(sample_table$subtissue)

load('old_data/data/all_tissue_quant.Rdata')
sample_table <- sample_table %>% filter(sample %in% colnames(all_quant))
avg_tissue_tpm <- lapply(subtissues, function(x) filter(sample_table, subtissue == x) %>% 
                             pull(sample) %>% 
                             {select(all_quant, transcript_id, .)} %>% 
                             filter(!is.na(.[,2])) %>% 
                             mutate(avg_exp := rowMeans(.[,-1])) %>% 
                             select(transcript_id, avg_exp)
)
median_tissue_tpm <- lapply(subtissues, function(x) filter(sample_table, subtissue == x) %>% 
                                pull(sample) %>% 
                                {select(all_quant, transcript_id, .)} %>% 
                                filter(!is.na(.[,2])) %>% 
                                mutate(med_exp := rowMedians(.[,-1] %>% as.matrix)) %>% 
                                select(transcript_id, med_exp)
)

names(avg_tissue_tpm) <-names(median_tissue_tpm) <-  subtissues

read_VEP <- function(file, tissue){
    fread(file, sep = '\t', skip = '##') %>% 
        as_tibble %>% rename(allele_id = `#Uploaded_variation`) %>% 
        filter(allele_id %in% VUS_ids) %>% 
        select(allele_id,transcript_id = Feature,Consequence, Extra) %>% 
        mutate(IMPACT = str_split(Extra, '=|;') %>% sapply(function(x) x[2]),
               Consequence= str_split(Consequence, ',')) %>% 
        unnest(Consequence) %>% 
        inner_join(clinsig_rank) %>% 
        left_join(median_tissue_tpm[[tissue]]) %>% 
        left_join(avg_tissue_tpm[[tissue]]) %>%
        select(allele_id,impact_code, avg_exp, med_exp) %>% 
        distinct %>% 
        group_by(allele_id) %>% #select the max variant for each allele:transcript pair 
        summarise(max_impact = max(impact_code), 
                  max_impact_avg_exp = avg_exp[which.max(impact_code)], 
                  max_impact_med_exp = med_exp[which.max(impact_code)] ) %>% 
        mutate(subtissue = tissue)
    
}

vep_files <- paste0('old_data/data/vep/', subtissues,'/variant_summary.txt')
all(file.exists(vep_files))
vep_by_tissue <- lapply(seq_along(subtissues), function(i)read_VEP(vep_files[i], subtissues[i])) %>% bind_rows
apply(vep_by_tissue[,-1], 2, function(x) sum(is.na(x)))

vep_impact_matrix <- vep_by_tissue %>% 
    select(allele_id, subtissue, max_impact) %>% 
    spread(key = subtissue, value = max_impact) %>% 
    mutate(one_count = rowSums(.[,subtissues] == 1), 
           two_count = rowSums(.[,subtissues] == 2), 
           three_count = rowSums(.[,subtissues] == 3),
           var = rowVars(as.matrix(.[,subtissues]) ) )%>% 
    arrange(desc(var))

vep_impact_matrix[is.na(vep_impact_matrix)] <- 1
eye_tissues <- c("RPE_Fetal.Tissue",  "Retina_Adult.Tissue",  "RPE_Adult.Tissue", "Cornea_Adult.Tissue", "Cornea_Fetal.Tissue", "Retina_Fetal.Tissue")
body_tissues <- filter(sample_table, !subtissue%in% eye_tissues) %>% pull(subtissue) %>% unique
save(vep_impact_matrix, vep_by_tissue, file = '/data/swamyvs/ocular_transcriptomes_paper/clean_data/rdata/vep_results.Rdata')
