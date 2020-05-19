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
clinvar_anno <- fread('ref/clinvar_variant_summary.txt.gz', sep = '\t') %>% 
    as_tibble %>% rename(allele_id = `#AlleleID`)
VUS_ids <- clinvar_anno %>% filter(ClinicalSignificance == "Uncertain significance") %>% pull(allele_id) %>% unique
clinsig_rank <- tibble(IMPACT = c('HIGH', 'MODERATE', 'LOW', "MODIFIER")) %>%  
    mutate(impact_code  = case_when(IMPACT == 'HIGH' ~ 3, 
                                    IMPACT == 'MODERATE' ~2, 
                                    TRUE ~ 1) )
sample_table <- read_tsv(files$sample_table)
subtissues <- unique(sample_table$subtissue)

vep_files <- paste0('data/vep/', subtissues,'/variant_summary.txt')
subtissues <- subtissues[file.exists(vep_files)]
vep_files <- vep_files[file.exists(vep_files)]
vep_by_tissue <- lapply(seq_along(subtissues), function(i)read_VEP(vep_files[i], subtissues[i])) %>% reduce(full_join)

eye_tissues <- c("RPE_Fetal.Tissue",  "Retina_Adult.Tissue",  "RPE_Adult.Tissue", "Cornea_Adult.Tissue", "Cornea_Fetal.Tissue", "Retina_Fetal.Tissue")
vep_exp_all_eye <- vep_by_tissue %>% 
    filter(rowSums(.[,eye_tissues] == 3) == length(eye_tissues) ) #%>% 
    mutate(var= rowVars(as.matrix(.[,-1])) )

vep_no_eye <- vep_by_tissue %>% 
    filter(rowSums(.[,eye_tissues] == 1) == length(eye_tissues) ) %>% 
    mutate(var= rowVars(as.matrix(.[,-1])) )
vep_example  <- bind_rows( vep_exp_all_eye %>% arrange(desc(var)) %>% .[1:ncol(.),] %>% select(allele_id, eye_tissues, everything()),
                vep_no_eye %>% arrange(desc(var)) %>% .[1:ncol(.),] %>% select(allele_id, eye_tissues, everything()) )

save(vep_by_tissue, file = all_variant_results_file)
save(vep_example, file = example_variant_result_file)




