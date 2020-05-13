library(tidyverse)
library(data.table)
library(argparse)

parser <- ArgumentParser()
parser$add_argument ('--longReadDir', action = 'store', dest = 'long_read_dir')
parser$add_argument('--cleanedData', action = 'store', dest = 'all_lf_clean_data_file')

list2env(parser$parse_args(), .GlobalEnv)
setwd(long_read_dir)

source('~/scripts/read_salmon.R')

load_and_preprocess <- function(gtf_file, track_file, conv_file){
    gtf <- rtracklayer::readGFF(gtf_file)
    track <- fread(track_file) %>% as_tibble
    conv_tab <- fread(conv_file) %>% as_tibble 
    
    txid2tlid <- conv_tab %>% select(transcript_id,code,  pacbio) %>% filter(!is.na(pacbio)) %>% rename(annot_transcript_id=pacbio)
    talon_abundance <- fread('data/talon_results/RPE_Fetal.Tissue/RPE_Fetal.Tissue_talon_abundance.tsv') %>% 
        inner_join(txid2tlid) %>% mutate(total_ab= RPE_D42_small+ RPE_D42_large)
    no_ism_ab <- talon_abundance %>% 
        filter(transcript_novelty != 'ISM') 
    ism_tx <- talon_abundance %>% filter(transcript_novelty == 'ISM') %>% pull(transcript_id)
    
    
    track <- track %>% filter(!transcript_id %in% ism_tx)
    stringtie_transcripts <- track %>% select(contains('st_')) %>% rowSums(.) %>% {. > 0} %>% {filter(track, .) } %>%
        
        pull(transcript_id)
    filt_stringtie_tx <- track %>% select(contains('stf_')) %>% rowSums(.) %>% {. > 0} %>% {filter(track, .) } %>%
        pull(transcript_id)
    sr_stringtie_stringtie_tx <- track %>% select(contains('stsr_')) %>% rowSums(.) %>% {. > 0} %>% {filter(track, .) } %>%
        pull(transcript_id)
    
    scallop_transcripts <- track %>% select(contains('sp_')) %>% rowSums(.) %>% {. > 0} %>% {filter(track, .) } %>%
        pull(transcript_id)
    pacbio_transcripts <- track %>% filter(.$pacbio) %>% pull(transcript_id)
    
    
    length_df <- gtf %>% filter(type == 'exon', !transcript_id %in% ism_tx) %>%
        mutate(ex_length=end-start) %>%
        group_by(transcript_id) %>%
        summarise(length=sum(ex_length) + n() ) %>%
        filter(transcript_id %in% c(scallop_transcripts, pacbio_transcripts, stringtie_transcripts)) %>% 
        mutate(in_stringtie=transcript_id %in% stringtie_transcripts ,
               in_scallop=transcript_id %in% scallop_transcripts,
               in_pacbio=transcript_id %in% pacbio_transcripts,
               in_stf = transcript_id %in% filt_stringtie_tx,
               in_stsr = transcript_id %in% sr_stringtie_stringtie_tx,
               build_type= case_when( in_stringtie  ~ 'stringtie',
                                      in_scallop   ~ 'scallop' ,
                                      in_pacbio   ~ 'pacbio'),
               intersection_case = case_when( in_stringtie & in_scallop & in_pacbio ~ 'all',
                                              in_stringtie & in_scallop & !in_pacbio ~ 'stringtie-scallop',
                                              in_stringtie & !in_scallop & in_pacbio ~ 'stringtie-pacbio',
                                              !in_stringtie & in_scallop & in_pacbio ~ 'scallop-pacbio',
                                              in_stringtie & !in_scallop & !in_pacbio ~ 'stringtie',
                                              !in_stringtie & in_scallop & !in_pacbio ~ 'scallop',
                                              !in_stringtie & !in_scallop & in_pacbio ~ 'pacbio')
        )
    rm(stringtie_transcripts, scallop_transcripts, pacbio_transcripts)
    res= list(gtf=gtf, length_df=length_df, conv_tab = conv_tab, track = track, no_ism_ab = no_ism_ab)
    return(res)
}

loose_res_list <- load_and_preprocess('data/combined_gtfs/all_RPE_loose.combined.gtf',
                                      'data/gtf_info/all_RPE_loose_detdf.tsv.gz',
                                      'data/gtf_info/all_RPE_loose_convtab.tsv.gz')
length_df <- loose_res_list$length_df

stringtie_tx <- length_df %>% filter(in_stringtie) %>% pull(transcript_id)
scallop_tx <- length_df %>% filter(in_scallop) %>% pull(transcript_id)
pacbio_tx <- length_df %>% filter(in_pacbio) %>% pull(transcript_id)
tx_below_min_length <- length_df %>% filter(length< 300) %>% pull(transcript_id)

st_tx_below_min_length <- round(sum(tx_below_min_length %in% stringtie_tx)/length(stringtie_tx) *100 , 3)
sp_tx_below_min_length <- round(sum(tx_below_min_length %in% scallop_tx)/length(scallop_tx) *100, 3)
gtf <- loose_res_list$gtf
ref_tx <- filter(gtf, class_code == '=', type == 'transcript') %>% pull(transcript_id)
inter_genic_tx <-  filter(gtf, class_code == 'u', type == 'transcript') %>% pull(transcript_id)
filter(length_df, transcript_id %in% ref_tx, in_pacbio) %>% nrow 


sample_table <- read_tsv('sampleTableRPE_V2.tsv')
true_rpe_samples <- sample_table %>% filter(origin == 'true_RPE') %>% pull(sample)
gencode_quant <-read_salmon('data/salmon_quant/gencode_comp_ano/', which_return = 'transcript', quant_type = 'abundance')
gencode_quant <- gencode_quant %>% select(transcript_id, true_rpe_samples) %>% 
    filter(rowSums(.[,-1] > 0) >0) %>% 
    mutate(avg_quant = rowMeans(.[,true_rpe_samples]))
sum(rowSums(gencode_quant[,true_rpe_samples] >=10) >= 1)

no_ism_ab <- loose_res_list$no_ism_ab
total_flnc_hq_no_ism_reads <- sum(no_ism_ab$total_ab)
total_ref_annot_tx <- filter(no_ism_ab, code == '=') %>% nrow 
total_intergenic_tx <- filter(no_ism_ab, code == 'u') %>% nrow 
total_unannotated_tx <- nrow(no_ism_ab) - total_ref_annot_tx  - total_intergenic_tx
long_read_results_table <- tibble(label = c('FLNC, high quality, Non-ISM reads', 'Gencode Annotated Transcripts', 
                                            'Novel Gene Isoforms', 'Novel Transcribed Loci' ) ,
                                  value = c(total_flnc_hq_no_ism_reads, total_ref_annot_tx, total_unannotated_tx, total_intergenic_tx))


length_df_all_valid_tx <- length_df %>% filter(!transcript_id %in% tx_below_min_length)

stringtie_accuracy <- sum(length_df_all_valid_tx$in_stringtie & length_df_all_valid_tx$in_pacbio ) / sum(length_df_all_valid_tx$in_stringtie)
scallop_accuracy <- sum(length_df_all_valid_tx$in_scallop & length_df_all_valid_tx$in_pacbio ) / sum(length_df_all_valid_tx$in_scallop)
acc_df_alltx <- tibble(build = c('stringtie', 'scallop'), 
                 accuracy = c(stringtie_accuracy, scallop_accuracy))

max_tx_length <- length_df_all_valid_tx %>% arrange(length) %>% tail(100) %>% pull(length) %>% .[1]
max_tx_length <- 10000
length_df_plotting <- length_df_all_valid_tx %>% filter(length <= max_tx_length,!(!in_stringtie & !in_pacbio) )
#save(length_df_plotting, file = '/data/swamyvs/ocular_transcriptomes_paper/CSHL_poster/data/pacbio_vs_st_upset.Rdata')

ref_length_df <- length_df_plotting %>% filter(transcript_id %in% ref_tx )
novel_length_df <- length_df_plotting %>% filter(!transcript_id %in% ref_tx, length <= max_tx_length)


acc_by_length <-  length_df_all_valid_tx %>% 
    mutate(qbin = cut(length, breaks=seq(0,7000,1000) ,seq(1000,7000,1000), include.lowest = T,  ), 
           qbin  = as.character(qbin),
           qbin = replace_na(qbin, '7000+')) #%>% 

acc_df_by_length <- acc_by_length %>% group_by(qbin) %>% 
    summarise(stringtie_acc = sum(in_stringtie & in_pacbio)/sum(in_stringtie), 
              scallop_acc =sum(in_scallop & in_pacbio)/sum(in_scallop), 
              number_of_tx_in_bin = n())  %>%  
    mutate(qbin_pretty= c(paste0(seq(0,6000, 1000), '-', qbin)[1:7],'7000+' ))


conv_tab <- loose_res_list$conv_tab 
st_files <- paste0( 'data/stringtie/',true_rpe_samples,'.gtf' )
names <- paste0('st_',true_rpe_samples)

read_stringtie <- function(file, name ){
    df <- data.table::fread(file, sep='\t', header = F) %>% as_tibble %>% filter(V3 == "transcript")
    tibble(
        !!name := str_extract(df$V9, 'transcript_id \".+\";') %>% str_split('\"') %>% sapply(function(x)x[2]),
        TPM=str_extract(df$V9, 'TPM \".+\";') %>% str_split('\"') %>% sapply(function(x)x[2]) %>% as.numeric 
    )
    
}


stringtie_quant <- lapply(seq_along(st_files), function(i) read_stringtie(st_files[i], names[i]) %>% 
                              inner_join(conv_tab[c('transcript_id', names[i])]) %>% 
                              select(-(!!names[i]))  %>% rename(!!names[i]:=TPM)) %>% 
    reduce(full_join) %>% 
    select(transcript_id, everything()) %>% 
    filter(transcript_id %in% length_df_all_valid_tx$transcript_id)
stringtie_quant[is.na(stringtie_quant)] <- 0
pacbio_tx_valid <- length_df %>% filter(in_pacbio) %>% pull(transcript_id)

filter_and_calc_acc <- function(lvl, quant, ftype ){
    if(ftype == 'mean'){
        build_tx <- quant %>% filter(rowMeans(.[,-1]) >=lvl) %>% pull(transcript_id)
    }
    
    else {
        build_tx <- quant %>% filter(rowSums(.[,-1] >= lvl) >0) %>% pull(transcript_id)
    }
    
    acc <- sum(build_tx %in% pacbio_tx_valid )/length(build_tx) 
    return(tibble(co=lvl, acc=acc, tx_in_build = length(build_tx)))
}




filter_and_calc_acc_by_length <- function(q, ftype){
    tx <- filter(acc_by_length, qbin == q) %>% pull(transcript_id)
    c_quant <- filter(stringtie_quant, transcript_id %in% tx )
    lapply(seq(0,25,1),  function(i) filter_and_calc_acc(i, c_quant, ftype)) %>% 
        bind_rows() %>% 
        bind_cols(., tibble(qbin=rep(q, nrow(.))))
    
}

tpm_thresholding_by_length_stmerge <- lapply(unique(acc_by_length$qbin),function(i) filter_and_calc_acc_by_length(i, 'merge') ) %>%
    bind_rows() %>% mutate(filtering_type = 'merge')
tpm_thresholding_by_length_mean <- lapply(unique(acc_by_length$qbin),function(i) filter_and_calc_acc_by_length(i, 'mean') ) %>%
    bind_rows() %>% mutate(filtering_type = 'mean')

tpm_df <- bind_rows(tpm_thresholding_by_length_stmerge, tpm_thresholding_by_length_mean) %>% 
    inner_join(acc_df_by_length %>% select(qbin, qbin_pretty)) %>% rename(`Length Interval` = qbin_pretty)
labdf <- tpm_df %>% group_by(co, filtering_type) %>% summarise(total_tx = sum(tx_in_build), y=max(acc)) %>% 
    filter(co %in%c(0,1,5,10,15,20,25))


below_1tpm <-  stringtie_quant %>% filter(rowMeans(.[,names])<1) %>% .[, names] %>%  colSums()
avg_acc <- tpm_df %>% filter(co ==1, filtering_type == 'mean', !qbin%in%c(1000,2000)) %>% pull(acc) %>% mean
total_exp <- stringtie_quant %>% .[, names] %>%  colSums()
avg_below_exp <- round(mean(below_1tpm/total_exp), 3)*100
ref_above_1tpm <- stringtie_quant %>% filter(rowMeans(.[,names])>=1, transcript_id %in% ref_tx)  %>% nrow  
novel_above_1tpm <- stringtie_quant %>% filter(rowMeans(.[,names])>=1, !transcript_id %in% ref_tx)  %>% nrow   
tpm_df_plotting <- tpm_df
labdf_plotting <- labdf 
#save(tpm_df_plotting, labdf_plotting, file = '/data/swamyvs/ocular_transcriptomes_paper/CSHL_poster/data/TPM_co.Rdata')

save(tpm_df_plotting, labdf_plotting, acc_df_by_length, acc_df_alltx, length_df_plotting, ref_tx, max_tx_length, 
     file =all_lf_clean_data_file )




