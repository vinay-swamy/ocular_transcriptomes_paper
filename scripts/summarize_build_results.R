library(yaml)
library(tidyverse)
library(data.table)
library(RBedtools)
library(RColorBrewer)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--workingDir', action='store', dest='working_dir')
parser$add_argument('--dataDir', action = 'store', dest = 'data_dir')
parser$add_argument('--filesYaml', action = 'store', dest = 'files_yaml')
parser$add_argument('--distinctCDStrack', action = 'store', dest = 'distinct_CDS_track_file')
parser$add_argument('--compCDStrack', action = 'store', dest = 'comp_CDS_track_file')
list2env(parser$parse_args(), .GlobalEnv)

####
working_dir <- '/data/swamyvs/ocular_transcriptomes_paper/'
data_dir <- '/data/swamyvs/ocular_transcriptomes_pipeline/'
files_yaml <- '/data/swamyvs/ocular_transcriptomes_paper/files.yaml'
CDS_track_file <- '/data/swamyvs/ocular_transcriptomes_paper/clean_data/CDS_gtf/CDS_comp.tracking'
####
files <- read_yaml(files_yaml)
setwd(data_dir)

load(files$exon_class_rdata) 
gtf <- rtracklayer::readGFF(files$anno_gtf)
tcons2mstrg <- fread(files$tcons2mstrg, sep = '\t') %>% as_tibble 
sample_table <- read_tsv(files$sample_table) %>% filter(subtissue != 'synth')
subtissues <- unique(sample_table$subtissue)
#Count number of novel transcripts per tissue
#---- 
tissue_color_mapping_df <- tibble(body_location=c('Cornea_Adult.Tissue', 'Cornea_Fetal.Tissue', 
                                 'Retina_Adult.Tissue', 'Retina_Fetal.Tissue', 
                                 'RPE_Adult.Tissue', 'RPE_Fetal.Tissue',
                                 "ESC_Stem.Cell.Line","Lens_Stem.Cell.Line", 
                                 'Brain', 'Body'),
                 color=c(brewer.pal(8,'Greens')[c(6,8)], brewer.pal(8,'Blues')[c(6,8)], brewer.pal(8,'Reds')[c(6,8)],
                         'pink','purple',  'orange', 'yellow')
)

tcons2mstrg <- gtf %>% 
    filter(type == "transcript") %>% 
    select(transcript_id, transcript_type) %>% 
    distinct %>% 
    inner_join(tcons2mstrg) %>% 
    mutate(is.PC = transcript_type == 'protein_coding')
novel_transcripts_per_tissue <- lapply(subtissues, function(x) 
                                                tcons2mstrg %>% filter(!class_code %in% c('=', 'u')) %>%
                                                select(is.PC, !!x) %>% filter( .[[x]] !='') %>% 
                                                                    {tibble(subtissue = x,  
                                                                            novel_pc_count=sum(.$`is.PC`),
                                                                            novel_nc_count = sum(!.$`is.PC`)) } ) %>%  bind_rows %>% 
    left_join(sample_table %>% select(subtissue, body_location) %>% distinct ) %>% group_by(body_location) %>% 
    summarise(novel_pc_count=mean(novel_pc_count), novel_nc_count = mean(novel_nc_count)) %>% left_join(tissue_color_mapping_df) %>% 
 mutate(body_location=replace(body_location, body_location == 'Body', 'Body(avg)'),
        body_location=replace(body_location, body_location == 'Brain', 'Brain(avg)')
        )

#----
#count the number of novel ORFs 
#----
## NOTE: I'm considering a  CDS  as annotated if it exists already exists in the gencode annotation, OR 
## if its a unannotated CDS associated with annotated transcript

dntx_cds_track_tab <- fread(comp_CDS_track_file, sep = '\t', header = F) %>% as_tibble %>% 
    mutate(cds_id = str_split(V3, '\\|') %>% sapply(function(x) x[2]), 
           transcript_id = str_split(V5, '\\|') %>% sapply(function(x)x[2])) %>% 
    select(cds_id, transcript_id)

ref_CDS_ids_enst <- fread(distinct_CDS_track_file, sep = '\t', header = F) %>% as_tibble %>% 
    filter(V6!='-') %>% pull(V1) %>% unique()


ref_cds_ids_dntx <- tcons2mstrg %>% 
    inner_join(dntx_cds_track_tab) %>% 
    filter(!is.na(cds_id), class_code == '=') %>%
    pull(cds_id) %>% 
    unique()
all_annotaed_cdsid = c(ref_cds_ids_dntx, ref_CDS_ids_enst) %>% unique
cds_origin_df <- dntx_cds_track_tab %>% mutate(is_annotaed_cds = cds_id%in% all_annotaed_cdsid)
fwrite(cds_origin_df, file = files$CDS_origin_df, sep = '\t')

transcripts_with_novel_cds <- filter(dntx_cds_track_tab, !cds_id %in% all_annotaed_cdsid) %>% pull(transcript_id)
tcons2mstrg <- tcons2mstrg %>% mutate(has_novel_orf = (transcript_id %in% transcripts_with_novel_cds ) & transcript_type =='protein_coding' )
novel_orfs_per_tissue <- lapply(subtissues, function(x) 
    tcons2mstrg %>% filter(!class_code %in% c('=', 'u')) %>%
        select(has_novel_orf,is.PC, !!x) %>% filter(.[[x]] !='') %>% 
        {tibble(subtissue = x,  
                novel_orf_count=sum(.$has_novel_orf & .$`is.PC`),
                ref_orf_count = sum(!(.$has_novel_orf) & .$`is.PC`) )}) %>%  
    bind_rows %>% 
    left_join(sample_table %>% select(subtissue, body_location) %>% distinct ) %>% 
    group_by(body_location) %>% 
    summarise(novel_orf_count=mean(novel_orf_count), ref_orf_count = mean(ref_orf_count)) %>% 
    left_join(tissue_color_mapping_df) %>% 
    mutate(body_location=replace(body_location, body_location == 'Body', 'Body(avg)'),
           body_location=replace(body_location, body_location == 'Brain', 'Brain(avg)')
    )




#----
#count the number of novel loci per tissue
#----
tc2m_novel_loci <- filter(tcons2mstrg, transcript_id %in% novel_loci_distinct$transcript_id)
novel_loci_per_tissue <- lapply(subtissues, function(i) tc2m_novel_loci %>% 
                                filter( .[[i]]!='')  %>% 
                                {tibble(subtissue=i, num_pc=sum(.$`is.PC`), num_nc=nrow(.) - sum(.$`is.PC`) )}  ) %>%
    bind_rows() %>% 
    left_join( sample_table %>% select(subtissue, body_location) %>% distinct ) %>% group_by(body_location) %>% 
    summarise(noncoding_tx=mean(num_nc), pc_tx=mean(num_pc)) %>% left_join(tissue_color_mapping_df) %>% 
    mutate(body_location=replace(body_location, body_location == 'Body', 'Body(avg)'),
           body_location=replace(body_location, body_location == 'Brain', 'Brain(avg)')
    ) 
#---- 
all_novel_tx <- filter(gtf, !class_code %in% c('=', 'u') ) %>% pull(transcript_id) %>% unique
exons_per_nvtx <- gtf %>% filter(type == 'exon', transcript_id %in% all_novel_tx ) %>%
    group_by(transcript_id) %>% summarise(n_exons = n())
novel_transcripts_no_novel_exon <- gtf %>% filter(type == 'exon', transcript_id %in% all_novel_tx) %>%
    anti_join(novel_exons_TSES, by = c("seqid", "start", "end", "strand")) %>% 
    group_by(transcript_id) %>% 
    summarise(n_annot_exons = n()) %>% 
    left_join(exons_per_nvtx) %>% 
    mutate(n_annot_exons = replace_na(n_annot_exons, 0), diff = n_exons - n_annot_exons) %>% #if the number exons did not change when remove novel exons, all exons must be annot 
    filter(diff == 0) %>% pull(transcript_id) %>% unique 
exon_count_df <- tibble(group = 'Novel Isoform', 
                        label = c('Novel Exon', 'Exon Omission'), 
                        count = c(sum(!all_novel_tx %in% novel_transcripts_no_novel_exon), sum(all_novel_tx %in% novel_transcripts_no_novel_exon)  ) 
                        )
exon_count_df <- novel_exons_TSES %>% 
    select(id, nv_type_rc) %>% 
    distinct %>% group_by(nv_type_rc) %>%
    summarise(count = n()) %>% 
    mutate(group = 'Novel Exon Type') %>% 
    select(group, label = nv_type_rc, count) %>% bind_rows(exon_count_df)


#----
#make make pc/nc -> pc-nvorf/pc-annorf table 
#----

pcnc_df <- tcons2mstrg %>% filter( !class_code %in%c('=', 'u')) %>% 
    {tibble(group = 'Transcipt Type', 
            label = c('Protein Coding', 'Noncoding'), 
            count = c(sum(.$transcript_type == 'protein_coding'), sum(.$transcript_type != 'protein_coding')  )
            )}
orfdf <- tcons2mstrg %>% filter(!class_code %in%c('=', 'u')) %>% 
    {tibble(group = 'ORF Type', 
            label = c('Novel ORF', 'Annotated ORF') ,
            count = c(sum(.$has_novel_orf), sum(!.$has_novel_orf & .$transcript_type == 'protein_coding') )
            )}
novel_isoform_anno_df <- bind_rows(pcnc_df, orfdf)

#---- 
#novel exon location analysis
### is the novel exon in the cds part 3 -  only check if the novel regions of novel exons are actually in the CDS

cds_bed <- gtf %>% filter(type == 'CDS') %>%
    group_by(transcript_id) %>%
    summarise(seqid=first(seqid), strand=first(strand), start=min(start), end=max(end)) %>%
    mutate(score=999) %>%
    select(seqid, start, end, transcript_id, score, strand) %>%
    from_data_frame %>%
    RBedtools('sort', i=.)
exon_bed <- all_exons %>% mutate(score=888) %>%
    select(seqid, start, end, origin, score, strand) %>%
    from_data_frame %>%
    RBedtools('sort', i=.)

novel_exons_tx <- novel_exons_TSES %>%
    inner_join(gtf %>% filter(type == 'exon') %>% select( seqid,strand, start, end, transcript_id ))
novel_seq_bed <- novel_exons_tx  %>%
    mutate(score= 777, cp_id=paste(id, transcript_id, sep = ';')) %>% # join novel exon id and txid
    select(seqid, start, end, cp_id, score, strand) %>%
    from_data_frame %>%
    RBedtools('sort',output = 'stdout', i=.) %>%
    RBedtools('subtract', options = '-s',output = 'stdout',  a=., b=exon_bed) %>%
    RBedtools('intersect', options = '-s -wo', a=., b=cds_bed) %>%
    to_data_frame

protein_coding_tx <- gtf %>% 
    filter(type == 'transcript', transcript_type == 'protein_coding') %>% pull(transcript_id)
exon_locations <- novel_seq_bed %>%
    mutate(transcript_id= str_split(X4, ';') %>% sapply(function(x) x[2]),
           id=str_split(X4, ';') %>% sapply(function(x) x[1])) %>%
    filter(transcript_id == X10) %>%
    select(id, transcript_id) %>%
    mutate(exon_location='CDS') %>%
    left_join(novel_exons_tx,.) %>%
    mutate(exon_location= case_when(is.na(exon_location) & (!transcript_id %in% protein_coding_tx) ~ 'NC',
                                    is.na(exon_location) & (transcript_id %in% protein_coding_tx)  ~ 'UTR',
                                    TRUE ~ exon_location)
    ) %>% distinct






#----

novel_exon_annotation <- gtf %>% 
    filter(type =='exon') %>%
    select(seqid, strand, start,end, transcript_id) %>% 
    distinct %>% 
    inner_join(tcons2mstrg %>% select(transcript_id,has_novel_orf, transcript_type)) %>% 
    left_join(dntx_cds_track_tab) %>% 
    inner_join(novel_exons_TSES) %>% 
    left_join(exon_locations %>% select(seqid, strand, start, end, transcript_id, exon_location)) %>% 
    select(transcript_id, nv_exon_id = id, cds_id, has_novel_orf,transcript_type, exon_location, nv_type_rc ) %>% 
    distinct()
novel_exon_annotation <-tcons2mstrg %>% filter(!class_code %in% c('=', 'u'), transcript_id%in% novel_transcripts_no_novel_exon) %>% 
    select(transcript_id, has_novel_orf, transcript_type) %>% 
    mutate(nv_type_rc = 'Exon Omission') %>% bind_rows(novel_exon_annotation,.)
    
exon_type_by_transcript_type <-novel_exon_annotation %>% 
    select(transcript_id, transcript_type, has_novel_orf, nv_type_rc) %>% distinct %>% 
    mutate(comp_transcript_type = case_when(transcript_type == 'protein_coding' & has_novel_orf ~ 'Novel Isoform Novel Protein', 
                                            transcript_type == 'protein_coding' & !has_novel_orf ~ 'Novel Isoform Known ORF',
                                            transcript_type != 'protein_coding' ~ 'Noncoding')
           ) %>% 
    select(comp_transcript_type, nv_type_rc)
    



#Determine whether a novel exon is in the coding region of a transcript
# novel exons protein coding, noncoding, or UTR
#----
# library(UpSetR)
# novel_exons_tx_added <- novel_exons_TSES %>% 
#     inner_join(gtf %>% select(seqid, strand, start,end ,transcript_id, gene_name)) %>% 
#     mutate(is.PC= transcript_id %in% gff3$ID)
# cds_bed <- gff3 %>% filter(type == 'CDS') %>% 
#     select(seqid, strand, start, end, transcript_id=ID) %>% 
#     group_by(transcript_id) %>% summarise(seqid=first(seqid), strand=first(strand), start=min(start), end=max(end)) %>% 
#     mutate(score=1000) %>% 
#     select(seqid, start, end, transcript_id, score, strand)
# novel_exon_bed <- novel_exons_tx_added %>% mutate(score = 999) %>% 
#     select(seqid, start, end, transcript_id, score, strand) %>% distinct 
# intersect <-  cds_bed %>% from_data_frame %>% 
#     RBedtools('sort', output = 'stdout', i=.) %>% 
#     RBedtools('intersect', options = '-s -wao', a= novel_exon_bed %>% from_data_frame %>% RBedtools('sort', i=. ), b=.) %>% 
#     to_data_frame()
# 
# 
# filter(intersect, X4==X10, X13 > 1) %>% select(X1, X6, X2, X3) %>% distinct %>% nrow 
# novel_exons_ano <- filter(intersect, X4==X10) %>% select(transcript_id=X4, seqid=X1, strand=X6, start=X2, end=X3) %>% 
#     mutate(is.CDS=T) %>% distinct %>% left_join(novel_exons_tx_added,.) %>%
#     mutate(is.CDS=replace_na(is.CDS, F))# make sure exons match the transcript they are found in 
# 
# novel_exon_location_analysis <- list(protein_coding=filter(novel_exons_ano, is.PC, is.CDS) %>% pull(id) %>% unique, 
#             non_coding=filter(novel_exons_ano, !is.PC) %>% pull(id) %>% unique,
#             UTR_change=filter(novel_exons_ano, is.PC, !is.CDS) %>% pull(id) %>% unique)


save(tissue_color_mapping_df, file = files$color_mapping_rdata)
save(novel_transcripts_per_tissue, 
     novel_orfs_per_tissue,
     novel_loci_per_tissue, 
     #novel_exon_location_analysis,
     novel_isoform_anno_df,
     exon_type_by_transcript_type,
     exon_count_df, file=files$build_results_rdata)






