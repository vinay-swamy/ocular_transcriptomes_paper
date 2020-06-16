library(tidyverse)
library(data.table)
library(argparse)
library(parallel)
library(yaml)
library(glue)

parser <- ArgumentParser()
parser$add_argument ('--longReadDir', action = 'store', dest = 'long_read_dir')
parser$add_argument('--filesYaml', action = 'store', dest = 'files_yaml')
list2env(parser$parse_args(), .GlobalEnv)
#############
# long_read_dir <- '/data/swamyvs/ocular_transcriptome_longread_analysis/'
# files_yaml <- '/data/swamyvs/ocular_transcriptomes_paper/files.yaml'
#############


setwd(long_read_dir)
files <- read_yaml(files_yaml)
nm_col_clean <- function(col){
    raw_name <- col %>% .[.!='-'] %>% .[!grepl('ENST', .)] %>% head
    name=raw_name %>% str_replace_all('XLOC_\\d+\\|','') %>% 
        str_split('\\.\\d|:|\\|') %>% .[[1]] %>% .[2] %>% 
        str_split('_MSTRG|_\\d+$') %>% .[[1]] %>% .[1]
    return(name)
}


proc_trackfile <- function(track_tab){
    det_df <- apply(track_tab[,-(1:4)],2, function(x) x!='-') %>% as_tibble 
    det_df <- det_df %>% 
        mutate(num_det=rowSums(det_df) ) %>%   
        select(num_det, everything()) 
    
    bind_cols(track_tab[,1:4],det_df) %>% 
        mutate(refid=str_replace_all(refid, '^-$', '.|.'),
               ref_gene_name = str_split(refid, '\\|') %>% sapply(function(x)x[1]),
               ref_transcript_id = str_split(refid, '\\|') %>% sapply(function(x)x[2])  ) %>% 
        select('transcript_id', 'gene_id', 'code', 'ref_gene_name', 'ref_transcript_id', 'num_det', everything()) ->k
    
}
process_columns <- function(tab,col_name){
    tab <- tab %>% filter((!! rlang::sym(col_name)) != '-')
    col <- tab %>% pull(!!col_name)
    name_col <- tab %>% pull(transcript_id)
    det <- suppressWarnings(str_split(col, ':|\\|'))
    z <- '
    most of the oIds have a 1:1 mapping to tx ids, so split into 2 chunks for faster run time, like a lot faster
    _s denotes simple case, _c dentotes complex case
    '
    
    d_lengths <- sapply(det, length)
    base <- min(d_lengths)
    simple <- d_lengths == base
    det_s <- det[simple]
    name_col_s <- name_col[simple]
    tx_simple <-  lapply(1:length(det_s), function(i)  det_s[[i]][3] %>%
                             c(name_col_s[i], .)) %>%
        do.call(rbind, .) %>% as.data.frame(stringsAsFactors=F)
    colnames(tx_simple) <- c('transcript_id', col_name)
    #%>% rename(!!col_name:=oId) %>% distinct
    tx_comp=data.frame()
    if(sum(!simple) > 0){
        det_c <- det[!simple]
        name_col_c <- name_col[!simple]
        
        tx_comp <- lapply(1:length(det_c), function(i)  det_c[[i]][-1] %>%
                              .[grepl('MSTRG\\.\\d+\\.\\d+', .) | grepl('ENST', .)] %>%
                              {tibble(transcript_id=rep(name_col_c[i], length(.)),oId= . )}) %>%
            bind_rows() %>% rename(!!col_name:=oId) %>% distinct
    }
    return(list(simple=tx_simple, comp=tx_comp))
}

process_long_read_results <- function(track_file, gtf_file){
    raw_track_tab <- data.table::fread(track_file, header = F, sep='\t') %>% as_tibble
    colnames(raw_track_tab) <- c('transcript_id', 'gene_id','refid','code', apply(raw_track_tab[,-(1:4)], 2, nm_col_clean))
    colnames(raw_track_tab)[ncol(raw_track_tab)] <- 'pacbio'
    det_df <- proc_trackfile(raw_track_tab)
    raw_track_tab <- raw_track_tab %>% select(transcript_id, gene_id, refid, code, contains('st_'), pacbio)
    cn <- colnames(raw_track_tab)[-(1:4)]
    tcons2mstrg <- mclapply(cn, function(col) process_columns(raw_track_tab,col), mc.cores = min(length(cn), 36))
    convtab <- lapply(tcons2mstrg, function(x) x[['simple']]) %>% reduce(full_join)
    gtf <- rtracklayer::readGFF(gtf_file)
    det_df <- det_df %>% mutate(in_stringtie = rowSums( select(., contains('st_'))) >0) %>% 
        rename(in_pacbio = pacbio)
    stringtie_transcripts <- det_df %>% filter(in_stringtie) %>% pull(transcript_id)
    pacbio_transcripts <- det_df %>% filter(in_pacbio) %>% pull(transcript_id)
    length_df <- gtf %>% filter(type == 'exon') %>%
        mutate(ex_length=end-start) %>%
        group_by(transcript_id) %>%
        summarise(length=sum(ex_length) + n() ) %>%
        mutate(in_stringtie=transcript_id %in% stringtie_transcripts ,
               in_pacbio=transcript_id %in% pacbio_transcripts,
               intersection_case = case_when( 
                                              in_stringtie &  in_pacbio ~ 'stringtie-pacbio',
                                              in_stringtie &  !in_pacbio ~ 'stringtie',
                                              !in_stringtie & in_pacbio ~ 'pacbio')
        )
    acc_by_length <-  length_df%>% 
        mutate(qbin = cut(length, breaks=seq(0,7000,1000) ,seq(1000,7000,1000), include.lowest = T,  ), 
               qbin  = as.character(qbin),
               qbin = replace_na(qbin, '7000+')) #%>% 
    
    acc_df <- acc_by_length %>% group_by(qbin) %>% 
        summarise(stringtie_acc = sum(in_stringtie & in_pacbio)/sum(in_stringtie), 
                  number_of_tx_in_bin = n())  %>%  
        mutate(qbin_pretty= c(paste0(seq(0,6000, 1000), '-', qbin)[1:7],'7000+' ))
    sample_table <- read_tsv('sampleTableRPE_V2.tsv')
    true_rpe_samples <- sample_table %>% filter(origin == 'true_RPE') %>% pull(sample)
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
                                  inner_join(convtab[c('transcript_id', names[i])]) %>% 
                                  select(-(!!names[i]))  %>% rename(!!names[i]:=TPM)) %>% 
        reduce(full_join) %>% 
        select(transcript_id, everything()) 
    stringtie_quant[is.na(stringtie_quant)] <- 0
    
    
    filter_and_calc_acc <- function(lvl, quant, ftype ){
        if(ftype == 'mean'){
            build_tx <- quant %>% filter(rowMeans(.[,-1]) >=lvl) %>% pull(transcript_id)
        }
        
        else {
            build_tx <- quant %>% filter(rowSums(.[,-1] >= lvl) >0) %>% pull(transcript_id)
        }
        
        acc <- sum(build_tx %in% pacbio_transcripts )/length(build_tx) 
        return(tibble(co=lvl, acc=acc, tx_in_build = length(build_tx)))
    }
    
    
    
    
    filter_and_calc_acc_by_length <- function(q, quant, ftype){
        tx <- filter(acc_by_length, qbin == q) %>% pull(transcript_id)
        c_quant <- filter(quant, transcript_id %in% tx )
        lapply(seq(0,25,1),  function(i) filter_and_calc_acc(i, c_quant, ftype)) %>% 
            bind_rows() %>% 
            bind_cols(., tibble(qbin=rep(q, nrow(.))))
        
    }
    
    tpm_thresholding_by_length_stmerge <- lapply(unique(acc_by_length$qbin),function(i) filter_and_calc_acc_by_length(i, stringtie_quant, 'merge') ) %>%
        bind_rows() %>% mutate(filtering_type = 'merge')
    tpm_thresholding_by_length_mean <- lapply(unique(acc_by_length$qbin),function(i) filter_and_calc_acc_by_length(i,stringtie_quant, 'mean') ) %>%
        bind_rows() %>% mutate(filtering_type = 'mean')
    
    tpm_df <- bind_rows(tpm_thresholding_by_length_stmerge, tpm_thresholding_by_length_mean) %>% 
        inner_join(acc_df %>% select(qbin, qbin_pretty)) %>% rename(`Length Interval` = qbin_pretty)
    labdf <- tpm_df %>% group_by(co, filtering_type) %>% summarise(total_tx = sum(tx_in_build), y=max(acc)) %>% 
        filter(co %in%c(0,1,5,10,15,20,25))
    return(list(gtf=gtf, convtab = convtab, det_df=det_df, tpm_df=tpm_df, labdf=labdf, length_df=length_df, 
                quant = stringtie_quant))

}
lr_data <- process_long_read_results(track_file = 'data/combined_gtfs/all_merge_lr_RPE_loose.tracking', 
                                              gtf_file = 'data/combined_gtfs/all_merge_lr_RPE_loose.combined.gtf')

length_df_plotting <- inner_join(lr_data$length_df, lr_data$det_df %>% select(transcript_id, code)) %>% 
    mutate(`Transcript Origin` = ifelse(code == '=', 'Gencode', 'Novel'))
tpm_df_plotting <- lr_data$tpm_df
labdf_plotting <- lr_data$labdf
flnc_lengths <- bind_rows(fread('data/fasta_lengths/RPE_D42_large.tsv', header=F, select = 'V2') %>% 
                              as_tibble %>% mutate(lib='large'),
                          fread('data/fasta_lengths/RPE_D42_small.tsv',header=F, select = 'V2') %>%
                              as_tibble %>% mutate(lib='small'))


acc_df_by_length <- filter(tpm_df_plotting , co ==0, filtering_type == 'mean') %>% 
    distinct %>% 
    select(-co, -filtering_type)
save(acc_df_by_length, length_df_plotting, tpm_df_plotting,flnc_lengths, labdf_plotting, file = files$long_read_results_rdata)

