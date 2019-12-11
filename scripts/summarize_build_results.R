library(tidyverse)
library(RBedtools)
library(RColorBrewer)
args <- c('/Volumes/data/occular_transcriptomes_paper/',
          '/Volumes/data/eyeintegration_splicing/',
          'rdata/novel_exon_classification.rdata',
          'data/gtfs/all_tissues.combined_NovelAno.gtf', 
          'sampleTableFull.tsv',
          'data/seqs/transdecoder_results/all_tissues.combined_transdecoderCDS.gff3',
          'data/misc/TCONS2MSTRG.tsv', 
          '/Volumes/data/occular_transcriptomes_paper/clean_data/tissue_to_colors.Rdata',
          '/Volumes/data/occular_transcriptomes_paper/clean_data/buildResultsSummary.Rdata')
args <- commandArgs(trailingOnly = T)
working_dir <- args[1]
data_dir <- args[2]
exon_classifcation_file <- args[3]
gtf_file <- args[4]
sample_table_file <- args[5]
gff3_file <- args[6]
tcons2mstrg_file <- args[7]
color_mapping_df <- args[8]
data_to_plot <- args[9]
setwd(data_dir)

load(exon_classifcation_file) 
gtf <- rtracklayer::readGFF(gtf_file)
gff3 <- rtracklayer::readGFF(gff3_file) %>% as_tibble  %>% mutate(ID=str_extract(ID,'DNTX_[0-9]+|ENSG[0-9]+'))
tcons2mstrg <- read_tsv(tcons2mstrg_file)
sample_table <- read_tsv(sample_table_file) %>% filter(subtissue != 'synth')
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
save(tissue_color_mapping_df, file = color_mapping_df )

novel_transcripts_per_tissue <- tibble(subtissue=subtissues, 
                                       novel_tx_count=sapply(subtissues, 
                                                             function(x) tcons2mstrg %>% 
                                                                 filter(!class_code %in% c('=', 'u')) %>%  .[[x]] %>% 
                                                                 .[!is.na(.)] %>% length)) %>% 
    left_join(sample_table %>% select(subtissue, body_location) %>% distinct ) %>% group_by(body_location) %>% 
    summarise(novel_transcript_count=mean(novel_tx_count)) %>% left_join(tissue_color_mapping_df) %>% 
 mutate(body_location=replace(body_location, body_location == 'Body', 'Body(avg)'),
        body_location=replace(body_location, body_location == 'Brain', 'Brain(avg)')
        )



# color_list<- novel_transcripts_per_tissue$color
# names(color_list) <- novel_transcripts_per_tissue$body_location
# ggplot(data = novel_transcripts_per_tissue) + 
#     geom_col(aes(x=body_location, y=novel_transcript_count, fill=body_location)) +
#     scale_fill_manual(values = color_list)+
#     theme(axis.text.x=element_text(angle=45, hjust = 1))

#----
#count the number of novel loci per tissue
#----
tc2m_novel_loci <- filter(tcons2mstrg, transcript_id %in% novel_loci_distinct$transcript_id) %>% 
    mutate(is.PC=transcript_id %in% gff3$ID)
novel_loci_per_tissue <- lapply(subtissues, function(i) tc2m_novel_loci %>% 
                                filter( !is.na(.[,i]))  %>% 
                                {tibble(subtissue=i, num_pc=sum(.$`is.PC`), num_nc=nrow(.) - sum(.$`is.PC`) )}  ) %>%
    bind_rows() %>% 
    left_join( sample_table %>% select(subtissue, body_location) %>% distinct ) %>% group_by(body_location) %>% 
    summarise(noncoding_tx=mean(num_nc), pc_tx=mean(num_pc)) %>% left_join(tissue_color_mapping_df) %>% 
    mutate(body_location=replace(body_location, body_location == 'Body', 'Body(avg)'),
           body_location=replace(body_location, body_location == 'Brain', 'Brain(avg)')
    ) %>% 
    gather(transcript_type, counts, -body_location)
    
# color_list<- novel_loci_per_tissue$color
# names(color_list) <- novel_loci_per_tissue$body_location
# ggplot(data = novel_loci_per_tissue) +
#     geom_col(aes(x=body_location, y=counts, fill=body_location, alpha=transcript_type), position = 'dodge') +
#     scale_fill_manual(values = color_list)+
#     scale_alpha_discrete(range=c(.5,1))+
#     theme_minimal()






#Determine whether a novel exon is in the coding region of a transcript
# novel exons protein coding, noncoding, or UTR
#----
library(UpSetR)
novel_exons_tx_added <- novel_exons_TSES %>% 
    inner_join(gtf %>% select(seqid, strand, start,end ,transcript_id, gene_name)) %>% 
    mutate(is.PC= transcript_id %in% gff3$ID)
cds_bed <- gff3 %>% filter(type == 'CDS') %>% 
    select(seqid, strand, start, end, transcript_id=ID) %>% 
    group_by(transcript_id) %>% summarise(seqid=first(seqid), strand=first(strand), start=min(start), end=max(end)) %>% 
    mutate(score=1000) %>% 
    select(seqid, start, end, transcript_id, score, strand)
novel_exon_bed <- novel_exons_tx_added %>% mutate(score = 999) %>% 
    select(seqid, start, end, transcript_id, score, strand) %>% distinct 
intersect <-  cds_bed %>% from_data_frame %>% 
    RBedtools('sort', output = 'stdout', i=.) %>% 
    RBedtools('intersect', options = '-s -wao', a= novel_exon_bed %>% from_data_frame %>% RBedtools('sort', i=. ), b=.) %>% 
    to_data_frame()


filter(intersect, X4==X10, X13 > 1) %>% select(X1, X6, X2, X3) %>% distinct %>% nrow 
novel_exons_ano <- filter(intersect, X4==X10) %>% select(transcript_id=X4, seqid=X1, strand=X6, start=X2, end=X3) %>% 
    mutate(is.CDS=T) %>% distinct %>% left_join(novel_exons_tx_added,.) %>%
    mutate(is.CDS=replace_na(is.CDS, F))# make sure exons match the transcript they are found in 

novel_exon_location_analysis <- list(protein_coding=filter(novel_exons_ano, is.PC, is.CDS) %>% pull(id) %>% unique, 
            non_coding=filter(novel_exons_ano, !is.PC) %>% pull(id) %>% unique,
            UTR_change=filter(novel_exons_ano, is.PC, !is.CDS) %>% pull(id) %>% unique)



save(novel_transcripts_per_tissue, novel_loci_per_tissue, novel_exon_location_analysis, file=data_to_plot)

# 
# 
# 
# pc_dist <- reduce(novel_exon_location_analysis, setdiff)
# nc_dist <- reduce(novel_exon_location_analysis[c('non_coding', 'UTR_change', 'protein_coding')], setdiff)
# utr_dist <- reduce(novel_exon_location_analysis[c('UTR_change','non_coding',  'protein_coding')], setdiff)
# 
# 
# changing_exons <- setdiff(unlist(novel_exon_location_analysis), c(pc_dist, nc_dist, utr_dist))
# eye_tissues <- c('Cornea_Adult.Tissue', 'Cornea_Fetal.Tissue', 
#                  'Retina_Adult.Tissue', 'Retina_Fetal.Tissue', 
#                  'RPE_Adult.Tissue', 'RPE_Fetal.Tissue')
# eye_tx <- tcons2mstrg %>% select(transcript_id, eye_tissues) %>% filter( rowSums(apply(.,2,is.na))<6) %>% pull(transcript_id)
# 
# filter(novel_exons_tx_added,transcript_id %in% eye_tx, id %in% changing_exons) %>% group_by(gene_name) %>% 
#     summarise(count=n()) %>% filter(count>1) %>% 
#     pull(gene_name ) %>% unique() %>% clipr::write_clip()






# novel exons, Promoter, terminal, or splicing











