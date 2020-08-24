library(tidyverse)
library(data.table)
library(matrixStats)
library(argparse)
library(yaml)
library(RBedtools)
library(parallel)
library(glue)
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
    as_tibble %>% rename(allele_id = `#AlleleID`, 
                         var_id = VariationID) 
clinvar_vus <- clinvar_anno %>% filter(ClinicalSignificance == "Uncertain significance")
clinvar_vus_eye <- clinvar_vus %>% filter(grepl('macula|retin|leber|cone|cornea|bardet|ocular|optic|ocular|joubert', PhenotypeList))

VUS_ids <-clinvar_vus  %>% pull(var_id) %>% unique
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
        select(allele_id,impact_code, avg_exp, med_exp, Consequence) %>% 
        distinct %>% 
        group_by(allele_id) %>% #select the max variant for each allele:transcript pair 
        summarise(max_impact = max(impact_code),
                  max_consequence = Consequence[which.max(impact_code)],
                  max_impact_avg_exp = avg_exp[which.max(impact_code)], 
                  max_impact_med_exp = med_exp[which.max(impact_code)] ) %>% 
        mutate(subtissue = tissue)
    
}

vep_files <- paste0(files$VEP_dir, subtissues,'/variant_summary.txt')
all(file.exists(vep_files))
dntx_vep_by_tissue <- lapply(seq_along(subtissues), function(i)read_VEP(vep_files[i], subtissues[i],med_dntx_tissue_tpm, avg_dntx_tissue_tpm )) %>%
    bind_rows %>% rename(var_id = allele_id)
gencode_vep_file <- paste0(files$VEP_dir, 'gencode/variant_summary.txt')
#slow but i lazy
gencode_vep_by_tissue <- lapply(seq_along(subtissues), function(i)read_VEP(gencode_vep_file, subtissues[i],med_gencode_tissue_tpm, avg_gencode_tissue_tpm )) %>%
    bind_rows %>% rename(var_id = allele_id)

save.image('testing/process_vep.Rdata')
dntx_global_impacts <- dntx_vep_by_tissue %>% 
    group_by(var_id) %>% 
    summarise(g_max_impact = max(max_impact), 
              g_min_impact = min(max_impact))


gencode_global_impacts <- gencode_vep_by_tissue %>% 
    select(var_id, gencode_impact=max_impact) %>% distinct 
impact_diff <- inner_join(dntx_global_impacts, gencode_global_impacts) %>% 
    mutate(max_diff = g_max_impact - gencode_impact, min_diff = gencode_impact - g_min_impact )


dntx_vep_impact_matrix <- dntx_vep_by_tissue %>% 
    select(var_id, subtissue, max_impact) %>% 
    spread(key = subtissue, value = max_impact) %>% 
    mutate(one_count = rowSums(.[,subtissues] == 1), 
           two_count = rowSums(.[,subtissues] == 2), 
           three_count = rowSums(.[,subtissues] == 3),
           var = rowVars(as.matrix(.[,subtissues]) ) )%>% 
    arrange(desc(var))
dntx_vep_impact_matrix[is.na(dntx_vep_impact_matrix)] <- 1

n=10
n=10
vep_exp_all_eye <- filter(dntx_vep_impact_matrix, var_id %in% clinvar_vus_eye$var_id) %>%
    filter(rowSums(.[,eye_tissues] == 3) >=1, rowSums(.[,body_tissues] !=3) ==length(body_tissues))

vep_full_priority_increase <- dntx_vep_impact_matrix %>% 
    #filter( var_id %in% clinvar_vus_eye$var_id) %>%
    inner_join(gencode_global_impacts)%>% 
    filter(three_count>=1, gencode_impact==1  )

    vep_no_eye <- filter(dntx_vep_impact_matrix, var_id %in% clinvar_vus_eye$var_id) %>%
    inner_join(gencode_global_impacts) %>% 
    filter(rowSums(.[,eye_tissues] == 1) == length(eye_tissues), gencode_impact == 3, three_count >0) %>% select(-gencode_impact)
impact_diff_changes <- impact_diff %>% filter(min_diff!=0 | max_diff!=0)
eye_impact_increase <- impact_diff_changes %>% filter(max_diff>0, var_id%in% clinvar_vus_eye$var_id)
vep_impact_increase <- filter(dntx_vep_impact_matrix, var_id %in% eye_impact_increase$var_id)


vep_example <- bind_rows( vep_exp_all_eye %>% arrange(desc(var))%>% 
                              .[1:min(n,nrow(.)),] %>% 
                              select(var_id, eye_tissues, everything()),
                          vep_no_eye %>% arrange(desc(var)) %>% 
                              .[1:min(n,nrow(.)),] %>% 
                              select(var_id, eye_tissues, everything()) ,
                          vep_impact_increase
)
clinvar_example <- filter(clinvar_vus_eye, var_id %in% vep_example$var_id) %>% select(var_id, contains('pheno'))
#vep_example <- filter(vep_impact_matrix, var_id %in% clinvar_vus_eye$var_id) %>% arrange(desc(var))
gencode_vep_by_tissue_simple <- gencode_vep_by_tissue %>% 
    mutate(subtissue = 'gencode', max_impact_avg_exp=1, max_impact_med_exp=1) %>% 
    distinct()
merge_bt <- bind_rows(dntx_vep_by_tissue %>% mutate(build= 'dntx'),
                      gencode_vep_by_tissue_simple%>% mutate(build= 'gencode'))
vep_plot_df <- merge_bt %>%  
    filter(var_id %in% vep_example$var_id) %>%
    mutate(avg_exp = log2(max_impact_avg_exp+1), 
           med_exp = log2(max_impact_med_exp + 1),
           avg_exp = replace(avg_exp , avg_exp >5, 5),
           var_id = as.character(var_id), 
           `Variant Impact` = case_when(max_impact == 1 ~'Low',
                                        max_impact == 2 ~ 'Moderate', 
                                        max_impact == 3 ~ 'High') %>% factor(levels=c('Low', 'Moderate', 'High')),
           subtissue = factor(subtissue, levels = c('gencode',eye_tissues, body_tissues))
    ) %>% 
    rename(`log2(TPM + 1)` = avg_exp)
vep_plot_df[is.na(vep_plot_df)] <- 0
#####################################################################################################################################  
#make 

example_var_ids <- c(631672, 631715)
gtf <- rtracklayer::readGFF(files$anno_gtf)

write(example_var_ids %>% paste0('^',.), '/tmp/varids.txt', sep = '\b')
read_VEP_ttx <- function(file, tissue){
    cmd <- glue('grep -Ef /tmp/varids.txt {file}')
    fread(cmd =  cmd, sep = '\t') %>% 
        as_tibble %>% 
        select(var_id = V1, transcript_id = V5, cons = V7, impact = V14 ) %>% 
        left_join(med_dntx_tissue_tpm[[tissue]]) %>% 
        # filter(allele_id %in% VUS_ids) %>% 
        # select(allele_id,transcript_id = Feature,Consequence, Extra) %>% 
        mutate(subtissue = tissue)
    
}

dntx_ttx_bytissue <- lapply(seq_along(subtissues), function(i)read_VEP_ttx(vep_files[i], subtissues[i] )) %>%
    bind_rows
res <- dntx_ttx_bytissue %>% mutate(med_exp = replace_na(med_exp, 0)) %>% 
    group_by(var_id, subtissue) %>% 
    summarise(max_tpm = max(med_exp), 
              transcript_id = transcript_id[which.max(med_exp)],
              cons = cons[which.max(med_exp)], 
              impact = impact[which.max(med_exp)])

k <- filter()

var_bed <- clinvar_vus_eye %>% 
    filter(var_id%in% example_var_ids, Assembly == 'GRCh38') %>% 
    select(Chromosome, Start, Stop, var_id) %>%
    mutate(Chromosome = paste0('chr', Chromosome)) %>% 
    from_data_frame
exon2var <- gtf %>% 
    filter(type == 'exon') %>% 
    select(seqid, start, end, transcript_id) %>% 
    from_data_frame %>% 
    RBedtools('sort', output = 'stdout', i=.) %>% 
    RBedtools('intersect', options = '-wa -wb', a=.,b=var_bed) %>% 
    to_data_frame()
    
conv_tab <- fread(files$tcons2mstrg)  
refid2dntx <- conv_tab %>% mutate(new_id = replace(transcript_id, class_code == '=',refid[class_code == '='])) %>% 
    select(transcript_id, new_id)


t_var <- 631672
t_gene <- "CACNA2D4"
other_tx <- gtf%>% filter(gene_name ==t_gene) %>% select(transcript_id) %>% inner_join(conv_tab) %>% 
    select(-contains('Retina'), -RPE_Fetal.Tissue) %>% filter(rowSums(.[,-(1:4)]!='') >0) %>% distinct %>%
    filter(class_code == '=') %>% pull(transcript_id)
hit_tx <- filter(res, var_id == t_var, transcript_id != '-') %>% pull(transcript_id)
all_tx <- c(hit_tx, other_tx)

all_tx_cleanids <- filter(refid2dntx, transcript_id %in% all_tx) %>% 
    select(transcript_id, new_id) %>% inner_join(tibble(transcript_id = all_tx),.) %>% pull(new_id)
ctab_plotting <- filter(conv_tab, transcript_id %in% other_tx) %>% 
    select(-refid, -class_code, -gene_id, -contains('Retina')) %>% 
    gather(key = 'tissue', value = 'oId', -transcript_id) %>% filter(oId != '') %>% 
    group_by(transcript_id) %>% summarise(subtissue = first(tissue)) %>%
    bind_rows( filter(res %>% ungroup, var_id == t_var, transcript_id != '-') %>% select(transcript_id, subtissue)) %>% 
    inner_join(refid2dntx) %>% select(-transcript_id) %>% rename(transcript_id = new_id)

library(DBI)
con <-  dbConnect(RSQLite::SQLite(), 'data/shiny_data/app_data/DNTX_db.sql')
pg <- con %>% tbl('plotting_gtf') %>% filter(gene_name == t_gene) %>% collect  #%>% filter(transcript_id%in%k$transcript_id)
pg_ctx <- pg %>% filter(transcript_id %in% all_tx_cleanids)


var_plot <- filter(exon2var, X4 == all_tx[1], X8 ==t_var) %>% 
    select(seqid = X1, start = X2, end = X3, transcript_id = X4, var_start = X6, var_end = X7, var_id = X8) %>% 
    inner_join(refid2dntx) %>% select(-transcript_id) %>% rename(transcript_id = new_id) %>% 
    inner_join(pg_ctx) %>% mutate(Xmin = Xmin + (var_start - start), Xmax = Xmin +1) %>% 
    select(seqid, start, end, transcript_id,  var_x = Xmin,var_y = Ymax)

save(vep_plot_df, pg_ctx, all_tx_cleanids,ctab_plotting,var_plot,
     example_var_ids, impact_diff, clinvar_vus_eye, file = files$variant_results_rdata)


dbDisconnect(con)

