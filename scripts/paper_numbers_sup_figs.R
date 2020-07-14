library(tidyverse)
library(RBedtools)
library(argparse)
library(yaml)
library(glue)
parser <- ArgumentParser()
parser$add_argument('--workingDir', action = 'store', dest = 'working_dir')
parser$add_argument('--dataDir', action = 'store', dest = 'data_dir')
parser$add_argument('--filesYaml', action = 'store', dest = 'files_yaml')
list2env(parser$parse_args(), .GlobalEnv)

########################
# working_dir <- '/data/swamyvs/ocular_transcriptomes_paper/'
# data_dir <- '/data/swamyvs/ocular_transcriptomes_pipeline/'
# files_yaml <- '/data/swamyvs/ocular_transcriptomes_paper/files.yaml'
########################


setwd(working_dir)
files <- read_yaml(files_yaml)
######### base numbers 
sample_table <- read_tsv(files$sample_table)
NUM_TOTAL_SAMPLES <- nrow(sample_table)
NUM_EYE_SAMPLES <- sample_table %>% filter(!body_location %in% c('Brain', 'Body')) %>% nrow 
NUM_BODY_SAMPLES <- sample_table %>% filter(body_location %in% c('Brain', 'Body')) %>% nrow 
NUM_BODY_TISSUES <- sample_table %>% filter(body_location %in% c('Brain', 'Body') | subtissue%in%c('Lens_Stem.Cell.Line', 'ESC_Stem.Cell.Line') ) %>% 
    pull(subtissue) %>% unique %>% length 
NUM_EYE_TISSUES <- sample_table %>% filter(!(body_location %in% c('Brain', 'Body') | subtissue%in%c('Lens_Stem.Cell.Line', 'ESC_Stem.Cell.Line')) ) %>% 
    pull(subtissue) %>% unique %>% length 
load(files$core_tight_rdata)
core_tight <- core_tight %>% as_tibble %>%  filter(sample_accession %in% sample_table$sample)
NUM_STUDIES = length(unique(core_tight$study_accession))
NUM_INITIAL_ANOO_TX <- scan(files$union_enst_ids,character(), sep = '\n') %>% n_distinct
NUM_TOTAL_GENCODE_TX <- system2('awk', glue(' \'$3 == "transcript"\' {files$ref_gtf} | wc -l '), stdout = T) %>% as.numeric
pan_body_gtf <- rtracklayer::readGFF(files$anno_gtf)
pan_eye_gtf <- rtracklayer::readGFF(files$pan_eye_gtf)
NUM_REF_TX_BODY <- pan_body_gtf %>% filter( type == 'transcript', class_code=='=')%>% nrow 
NUM_REF_TX_EYE  <- pan_eye_gtf  %>% filter( type == 'transcript', class_code=='=')%>% nrow 

NUM_NOVEL_TX_BODY <- pan_body_gtf %>% filter( type == 'transcript', class_code!='=')%>% nrow 
NUM_NOVEL_TX_EYE  <- pan_eye_gtf  %>% filter( type == 'transcript', class_code!='=')%>% nrow 
NUM_TOTAL_TX_BODY <- pan_body_gtf %>% filter(type == 'transcript') %>% nrow()
NUM_TOTAL_TX_EYE <- pan_eye_gtf %>% filter(type == 'transcript') %>% nrow()
load(files$ref_tx_exon_rdata)
all_exon_bed <- all_exons %>% mutate(score= 888) %>% select(seqid, start, end, origin, score, strand) %>% 
    from_data_frame %>% 
    RBedtools('sort', i=.)
dntx_exon_bed <- pan_body_gtf %>% filter(type == 'exon') %>% mutate(score = 999) %>% 
    select(seqid, start, end, transcript_id, score, strand) %>% from_data_frame %>% 
    RBedtools('sort', i=.)
novel_sequence <- RBedtools('subtract', options = '-s', output = 'stdout', a=dntx_exon_bed, b=all_exon_bed) %>% 
    RBedtools('sort', output = 'stdout', i=.) %>% 
    RBedtools('merge', options = '-s -c 6 -o distinct',i=. )%>% 
    to_data_frame
NUM_TOTAL_NOVEL_SEQ <- novel_sequence %>% mutate(length=X3-X2) %>% pull(length) %>% sum
NUM_REF_TX_BASE_TXOME <- system2('grep', "-c -e 'class_code \"=\"' clean_data/all_base_tx.combined.gtf",stdout = T) %>% 
    as.numeric
NUM_NOVEL_TX_BASE_TXOME <- system2('awk', " '$3 == \"transcript\"' clean_data/all_base_tx.combined.gtf | grep -c -v -e 'class_code \"=\"' - ",
                                   stdout = T) %>% 
    as.numeric
NUM_FRAC_REF_RECOVERED_BASE_TXOME <- {NUM_REF_TX_BASE_TXOME / NUM_TOTAL_GENCODE_TX *100} %>% round(3)
#######ORFs
cds_origin_df <- fread(files$CDS_origin_df)
NUM_REF_ORF <-  cds_origin_df %>% select(-transcript_id) %>% distinct %>% pull(is_annotaed_cds) %>% {sum(.)}
NUM_NOVEL_ORF <- cds_origin_df %>% select(-transcript_id) %>% distinct %>% pull(is_annotaed_cds) %>% {sum(!.)}
########

######## Txome size
load(files$txome_stats_rdata)
NUM_AVG_BASE_TXOME_SIZE = mean(gencode_tx_size$all_exp) %>% round(digits = 2)
size_red <- inner_join(gencode_tx_size, tx_counts) %>% mutate(frac_decrease = (all_exp - filtered)/all_exp)
NUM_AVG_TXOME_SIZE_DECREASE = mean(size_red$frac_decrease) %>% round(digits = 5)
NUM_AVG_MAP_RATE_DIF = mean(all_sample_mapping_rate_difference$mapping_rate_diff) %>% round(digits = 5)
########
######## Tx counts 
conv_tab <- fread(files$tcons2mstrg)
cds_origin <- fread(files$CDS_origin_df)
long_det_df <- apply(conv_tab[,-(1:4)],2, function(col) col!='') %>% as_tibble %>% bind_cols(conv_tab[,1:4], .) %>% 
    select(-refid, -gene_id) %>% gather(key = 'subtissue', value = 'det', -transcript_id,-class_code ) %>% 
    left_join(cds_origin)
tx_count_by_tissue <- long_det_df %>%
    mutate(class_code = replace(class_code,!class_code %in% c('=', 'u'), 'N')) %>% 
    group_by(subtissue, class_code) %>% 
    summarise(count = sum(det))
NUM_AVG_REF_TX_PER_TISSUE = tx_count_by_tissue %>% filter(class_code == '=') %>% pull(count) %>% mean %>% round(digits = 2)
NUM_AVG_NOVEL_ISOFORM_PER_TISSUE = tx_count_by_tissue %>% filter(class_code == 'N') %>% pull(count) %>% mean %>% round(digits = 2)
NUM_AVG_NOVEL_LOCI_PER_TISSUE = tx_count_by_tissue %>% filter(class_code == 'u') %>% pull(count) %>% mean %>% round(digits = 2)

orf_count_per_tissue <- long_det_df %>% filter(det, !is.na(cds_id)) %>% 
    select(subtissue, cds_id, is_annotaed_cds) %>% 
    distinct %>% 
    group_by(subtissue) %>% 
    summarise(count = sum(!is_annotaed_cds))
NUM_AVG_NOVEL_ORF_PER_TISSUE <- mean(orf_count_per_tissue$count) %>% round(digits = 2)

########
########long read data numbers 
load(files$long_read_results_rdata)
NUM_BASE_TXOME_CONST_ACC <- {sum(length_df_plotting$intersection_case == 'stringtie-pacbio') / sum(length_df_plotting$in_stringtie)} %>% round(3)
NUM_TXOME_CONST_ACC_2000_ABOVE <- acc_df_by_length %>% filter(as.numeric(qbin) >2000) %>% pull(acc) %>% mean %>% round(digits = 3)
NUM_TXOME_CONST_ACC_2000_BELOW <- acc_df_by_length %>% filter(as.numeric(qbin) <=2000) %>% pull(acc) %>% mean %>% round(digits = 3)
NUM_FILT_TXOME_ACC_2000_ABOVE <- tpm_df_plotting %>% filter(filtering_type == 'mean',co == 1 ) %>% 
    filter(as.numeric(qbin) >2000) %>% pull(acc) %>% mean %>% round(digits = 3)
NUM_FILT_TXOME_ACC_2000_BELOW <- tpm_df_plotting %>% filter(filtering_type == 'mean',co == 1 ) %>% 
    filter(as.numeric(qbin) <=2000) %>% pull(acc) %>% mean %>% round(digits = 3)
NUM_TOTAL_LR_TX <- sum(length_df_plotting$in_pacbio)
NUM_TOTAL_SR_TX <- sum(length_df_plotting$in_stringtie)
NUM_TXOME_MERGE_ACC_2000_ABOVE <- tpm_df_plotting %>% 
    filter(co == 1,  filtering_type == 'merge', !qbin%in%c('1000', '2000')) %>% pull(acc) %>% mean %>% round(digits = 3)
NUM_TXOME_MERGE_ACC_2000_BELOW <- tpm_df_plotting %>% 
    filter(co == 1,  filtering_type == 'merge', qbin%in%c('1000', '2000')) %>% pull(acc) %>% mean %>% round(digits = 3)
NUM_TXOME_FILTER_SIZE <- tpm_df_plotting %>% filter(co == 1, filtering_type == 'mean') %>% pull(tx_in_build) %>% sum

#######
####### novel isoform numbers 
load(files$novel_isoform_analysis_rdata)

appris_totals <- primary_isoforms_per_tissue %>% 
    spread(key = primary_isoform_origin, value = count) %>% 
    ungroup %>% 
    mutate(total =  rowSums(.[,-1]),
           percent_appris = appris / total,
           percent_gencode = gencode / total, 
           percent_novel = total )
NUM_AVG_APPRIS <- mean(appris_totals$percent_appris) %>% round(digits = 3)
NUM_AVG_NONAPPRIS_GENCODE = mean(appris_totals$percent_gencode) %>% round(digits = 3)
NUM_AVG_NOVEL = mean(appris_totals$percent_novel) %>% round(digits = 3)
plot_list <- novel_eye_tx_by_tissue[!names(novel_eye_tx_by_tissue) %in% c('Lens_Stem.Cell.Line', 'ESC_Stem.Cell.Line') ]
all_tx <- tibble(transcript_id = reduce(plot_list, union))
noveliso_ddf <- bind_cols(all_tx, lapply(plot_list, function(x) all_tx$transcript_id%in% x)) %>% 
    mutate(total = rowSums(.[,-1])) 
NUM_TSPEC_ISO = {sum(noveliso_ddf$total == 1) / nrow(noveliso_ddf)*100} %>% round(2)
NUM_AVG_PIU <- {mean(piu_df$piu)*100} %>% round(2)

#######
####### VEP numbers
load(files$variant_results_rdata)
impact_diff_changes <- impact_diff %>% filter(min_diff!=0 | max_diff!=0)
NUM_TOTAL_VUS <- nrow(impact_diff)
NUM_IMPACT_INCREASE <- impact_diff_changes %>% filter(max_diff>0) %>% nrow
NUM_IMPACT_DECREASE <- impact_diff_changes %>% filter(min_diff > 0) %>% nrow 
##########EXON count numbers 
load(files$build_results_rdata)
exon_count_df <- exon_type_by_transcript_type %>% group_by(label= nv_type_rc) %>% summarise(count= n())
NUM_TES_TSS <- filter(exon_count_df ,grepl('T\\wS', label) ) %>% 
    {sum(.$count) / sum(exon_count_df$count) *100} %>% round(digit =3)

####### Correlation between expression and transcript length 
load(files$gencode_quant)
subtissues <- unique(sample_table$subtissue)
ref_gtf <- rtracklayer::readGFF(files$ref_gtf)
ref_tx_lengths <- filter(ref_gtf, type == 'exon') %>% select(seqid, strand, start, end, transcript_id) %>% 
    mutate(length = end-start) %>% 
    group_by(transcript_id) %>% 
    summarise(length = sum(length) + n())

calc_txome_size <- function(t_tissue){
    samples <- filter(sample_table, subtissue == t_tissue) %>% pull(sample) 
    keep <- rowSums(gencode_quant[,samples]) != 0
    tissue_quant <- gencode_quant %>% filter(keep) %>% 
        mutate(avg_quant  = rowMeans(.[,samples])) %>% 
        inner_join(ref_tx_lengths) %>% 
        select(transcript_id, avg_quant, length)
    return(cor(tissue_quant$avg_quant, tissue_quant$length, method = 'spearman') )
    
}
res <- sapply(subtissues, calc_txome_size)
NUM_AVG_EXP_TX_LENGTH_COR <- mean(res)
############
#fetal retina dev diffexp
load(files$fetal_retina_diffexp_results)
NUM_DEVRET_SAMPLES = nrow(dev_retina_core_tight)
NUM_DEVRET_TIMESPOINTS = length(unique(dev_retina_core_tight$age_str))
NUM_DTU_TX = length(all_dtu_tx)
NUM_DTU_GENE = length(all_dtu_genes)
############
#validation_pvals 
load(files$CAGE_polyA_rdata)
gc_cage <- filter(all_CAGE_phase12,build == 'gencode') %>% pull(abs_dist)
dntx_cage <- filter(all_CAGE_phase12, build == 'dntx') %>% pull(abs_dist)
NUM_PVAL_CAGE <- wilcox.test(dntx_cage, gc_cage, alternative = 'less') %>% .[['p.value']] %>% replace(., . == 0, 2.2e-16)
rcompanion::wilcoxonR(dntx_cage, gc_cage)

gc_polya <- filter(all_polya_closest,label == 'gencode') %>% pull(abs_dist)
dntx_polya <- filter(all_polya_closest, label == 'dntx') %>% pull(abs_dist)
NUM_PVAL_POLYA <- wilcox.test(dntx_polya, gc_polya, alternative = 'less') %>% .[['p.value']] %>% round(digits=3) %>%  replace(., . == 0, 2.2e-16)
rcompanion::wilcoxonR(dntx_polya, gc_polya)

gc_phylop <- all_phylop %>% filter(build == 'gencode') %>% pull(mean_phylop_score)
dntx_phylop <- all_phylop %>% filter(build == 'dntx') %>% pull(mean_phylop_score)
NUM_PVAL_PHYLOP <- wilcox.test(dntx_phylop, gc_phylop, alternative = 'greater') %>% .[['p.value']] %>% replace(., . == 0, 2.2e-16)

rcompanion::wilcoxonR(dntx_phylop, gc_phylop)
#########
genes_with_novel_iso <- conv_tab %>% filter(Retina_Fetal.Tissue != '' | Retina_Adult.Tissue != '') %>% select(transcript_id) %>% 
    inner_join(pan_body_gtf) %>%
    filter(type == 'transcript', !class_code %in%c('=', 'u')) %>% 
    pull(gene_name) %>% unique
all_genes <- conv_tab %>% filter(Retina_Fetal.Tissue != '' | Retina_Adult.Tissue != '') %>% select(transcript_id) %>% 
    inner_join(pan_body_gtf) %>%
    pull(gene_name) %>% unique
ret_net_genes <- scan('/data/swamyvs/ocular_transcriptomes_pipeline/ref/retnet_hgncIDs_2017-03-28.txt',
                      what = character(), sep = '\n')

x=sum(ret_net_genes%in% genes_with_novel_iso)
m=length(ret_net_genes)
n = length(all_genes) - m
k = length(genes_with_novel_iso)
NUM_RET_NET_GENES <- m
NUM_RETNET_IN_DNTX <- x
NUM_RETNET_ENRICH_PVAL= dhyper(x, m, n, k) %>% formatC(format = 'e', digits = 1)

save(list= ls()[grepl('NUM_', ls())], file = files$paper_numbers_rdata)


