library(tidyverse)
library(data.table)
library(matrixStats)
library(argparse)
library(yaml)
library(RBedtools)
library(parallel)
parser <- ArgumentParser()
parser$add_argument('--dataDir', action = 'store', dest = 'data_dir')
parser$add_argument('--filesYaml', action = 'store', dest = 'files_yaml')
####
# data_dir <- '/data/swamyvs/ocular_transcriptomes_pipeline/'
# files_yaml <- '/data/swamyvs/ocular_transcriptomes_paper/files.yaml'
####
list2env(parser$parse_args(), .GlobalEnv)
files <- read_yaml(files_yaml)

setwd(data_dir)


clinvar_anno <- fread(files$clinvar_summary, sep = '\t') %>% 
    as_tibble %>% rename(allele_id = `#AlleleID`) 
clinvar_vus <- clinvar_anno %>% filter(ClinicalSignificance == "Uncertain significance")
clinvar_vus_eye <- clinvar_vus %>% filter(grepl('macula|retin|leber|cone|cornea|bardet|ocular|optic|ocular|joubert', PhenotypeList))

VUS_ids <-clinvar_vus  %>% pull(allele_id) %>% unique
clinsig_rank <- tibble(IMPACT = c('HIGH', 'MODERATE', 'LOW', "MODIFIER")) %>%  
    mutate(impact_code  = case_when(IMPACT == 'HIGH' ~ 3, 
                                    IMPACT == 'MODERATE' ~2, 
                                    TRUE ~ 1) )

sample_table <- read_tsv(files$sample_table)
subtissues <- unique(sample_table$subtissue)

load(files$all_tissue_quant)
load(files$gencode_quant)
gencode_quant[is.na(gencode_quant)] <- 0
all_quant[is.na(all_quant)] <- 0

sample_table <- sample_table %>% filter(sample %in% colnames(all_quant))

eye_tissues <-  c("RPE_Fetal.Tissue",  "Retina_Adult.Tissue",  "RPE_Adult.Tissue", "Cornea_Adult.Tissue", "Cornea_Fetal.Tissue", "Retina_Fetal.Tissue")
subtissues <- unique(sample_table$subtissue)
body_tissues <- sample_table %>% filter(!subtissue%in% eye_tissues) %>% pull(subtissue) %>% unique

make_tissue_level_quant <- function(quant, sum_type, col_name){
    res <- lapply(subtissues, function(x) filter(sample_table, subtissue == x) %>% 
               pull(sample) %>% 
               {select(quant, transcript_id, .)} %>% 
               filter(!is.na(.[,2])) %>% 
               mutate(!!col_name := get(sum_type)(as.matrix(.[,-1]) )) %>% 
               select(transcript_id, !!col_name)
    )
    names(res) <- subtissues
    return(res)
}


avg_dntx_tissue_tpm <-make_tissue_level_quant(all_quant, 'rowMeans', 'avg_exp')
med_dntx_tissue_tpm <- make_tissue_level_quant(all_quant, 'rowMedians', 'med_exp')
avg_gencode_tissue_tpm <- make_tissue_level_quant(gencode_quant, 'rowMeans', 'avg_exp')
med_gencode_tissue_tpm <- make_tissue_level_quant(gencode_quant, 'rowMedians', 'med_exp')


read_VEP <- function(file, tissue, c_med_t_tpm, c_avg_t_tpm){
    fread(file, sep = '\t', skip = '##') %>% 
        as_tibble %>% rename(allele_id = `#Uploaded_variation`) %>% 
        filter(allele_id %in% VUS_ids) %>% 
        select(allele_id,transcript_id = Feature,Consequence, Extra) %>% 
        mutate(IMPACT = str_split(Extra, '=|;') %>% sapply(function(x) x[2]),
               Consequence= str_split(Consequence, ',')) %>% 
        unnest(Consequence) %>% 
        inner_join(clinsig_rank) %>% 
        left_join(c_med_t_tpm[[tissue]]) %>% 
        left_join(c_avg_t_tpm[[tissue]]) %>%
        select(allele_id,impact_code, avg_exp, med_exp) %>% 
        distinct %>% 
        group_by(allele_id) %>% #select the max variant for each allele:transcript pair 
        summarise(max_impact = max(impact_code), 
                  max_impact_avg_exp = avg_exp[which.max(impact_code)], 
                  max_impact_med_exp = med_exp[which.max(impact_code)] ) %>% 
        mutate(subtissue = tissue)
    
}

vep_files <- paste0(files$VEP_dir, subtissues,'/variant_summary.txt')
all(file.exists(vep_files))
dntx_vep_by_tissue <- lapply(seq_along(subtissues), function(i)read_VEP(vep_files[i], subtissues[i],med_dntx_tissue_tpm, avg_dntx_tissue_tpm )) %>%
    bind_rows
gencode_vep_file <- paste0(files$VEP_dir, 'gencode/variant_summary.txt')
#slow but i lazy
gencode_vep_by_tissue <- lapply(seq_along(subtissues), function(i)read_VEP(gencode_vep_file, subtissues[i],med_gencode_tissue_tpm, avg_gencode_tissue_tpm )) %>%
    bind_rows

save.image('testing/process_vep.Rdata')
dntx_global_impacts <- dntx_vep_by_tissue %>% 
    group_by(allele_id) %>% 
    summarise(g_max_impact = max(max_impact), 
              g_min_impact = min(max_impact))


gencode_global_impacts <- gencode_vep_by_tissue %>% 
    select(allele_id, gencode_impact=max_impact) %>% distinct 
impact_diff <- inner_join(dntx_global_impacts, gencode_global_impacts) %>% 
    mutate(max_diff = g_max_impact - gencode_impact, min_diff = gencode_impact - g_min_impact )


dntx_vep_impact_matrix <- dntx_vep_by_tissue %>% 
    select(allele_id, subtissue, max_impact) %>% 
    spread(key = subtissue, value = max_impact) %>% 
    mutate(one_count = rowSums(.[,subtissues] == 1), 
           two_count = rowSums(.[,subtissues] == 2), 
           three_count = rowSums(.[,subtissues] == 3),
           var = rowVars(as.matrix(.[,subtissues]) ) )%>% 
    arrange(desc(var))
dntx_vep_impact_matrix[is.na(dntx_vep_impact_matrix)] <- 1

n=10
n=10
vep_exp_all_eye <- filter(dntx_vep_impact_matrix, allele_id %in% clinvar_vus_eye$allele_id) %>%
    filter(rowSums(.[,eye_tissues] == 3) >=1, rowSums(.[,body_tissues] !=3) ==length(body_tissues))
vep_no_eye <- filter(dntx_vep_impact_matrix, allele_id %in% clinvar_vus_eye$allele_id) %>% 
    filter(rowSums(.[,eye_tissues] != 3) == length(eye_tissues),three_count >4  ) 
impact_diff_changes <- impact_diff %>% filter(min_diff!=0 | max_diff!=0)
eye_impact_increase <- impact_diff_changes %>% filter(max_diff>0, allele_id%in% clinvar_vus_eye$allele_id)
vep_impact_increase <- filter(dntx_vep_impact_matrix, allele_id %in% eye_impact_increase$allele_id)


vep_example <- bind_rows( vep_exp_all_eye %>% arrange(desc(var))%>% 
                              #.[1:n,] %>% 
                              select(allele_id, eye_tissues, everything()),
                          vep_no_eye %>% arrange(desc(var)) %>% 
                              #.[1:n,] %>% 
                              select(allele_id, eye_tissues, everything()) ,
                          vep_impact_increase
)
clinvar_example <- filter(clinvar_vus_eye, allele_id %in% vep_example$allele_id) %>% select(allele_id, contains('pheno'))
#vep_example <- filter(vep_impact_matrix, allele_id %in% clinvar_vus_eye$allele_id) %>% arrange(desc(var))
gencode_vep_by_tissue_simple <- gencode_vep_by_tissue %>% 
    mutate(subtissue = 'gencode', max_impact_avg_exp=1, max_impact_med_exp=1) %>% 
    distinct()
merge_bt <- bind_rows(dntx_vep_by_tissue %>% mutate(build= 'dntx'),
                      gencode_vep_by_tissue_simple%>% mutate(build= 'gencode'))
vep_plot_df <- merge_bt %>%  
    filter(allele_id %in% vep_example$allele_id) %>%
    mutate(avg_exp = log2(max_impact_avg_exp+1), 
           med_exp = log2(max_impact_med_exp + 1),
           avg_exp = replace(avg_exp , avg_exp >5, 5),
           allele_id = as.character(allele_id), 
           `Variant Impact` = case_when(max_impact == 1 ~'Low',
                                        max_impact == 2 ~ 'Moderate', 
                                        max_impact == 3 ~ 'High') %>% factor(levels=c('Low', 'Moderate', 'High')),
           subtissue = factor(subtissue, levels = c('gencode',eye_tissues, body_tissues))
    ) %>% 
    rename(`log2(TPM + 1)` = avg_exp)
vep_plot_df[is.na(vep_plot_df)] <- 0

save(vep_plot_df, impact_diff, clinvar_vus_eye, file = files$variant_results_rdata)
