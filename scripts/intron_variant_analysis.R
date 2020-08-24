library(tidyverse)
library(yaml)
library(RBedtools)
library(glue)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--workingDir', action = 'store',dest = 'working_dir')
parser$add_argument('--dataDir', action = 'store', dest = 'data_dir')
parser$add_argument('--fileYaml', action = 'store', dest = 'files_yaml')

####
# working_dir <- '/data/swamyvs/ocular_transcriptomes_paper/'
# data_dir <- '/data/swamyvs/ocular_transcriptomes_pipeline/'
# files_yaml <- '/data/swamyvs/ocular_transcriptomes_paper/files.yaml'
#####
list2env(parser$parse_args(), .GlobalEnv)

files <- read_yaml(files_yaml)
panel_variants <- read_tsv(files$intron_variant_panel, 
                           col_names = c('gene_name', 'raw_loc', 'variant','NM_accession', 'rs_number')) %>% 
    mutate(seqid = str_split(raw_loc,':') %>% sapply(function(x)x[1]) %>% tolower,
           start = str_split(raw_loc,':') %>% sapply(function(x)x[2]) %>% as.numeric, 
           end = start+1, 
           var_id = paste0('pnl',0:(nrow(.)-1) ))
complete_gtf <- rtracklayer::readGFF(files$anno_gtf)

conv_tab <- fread(files$tcons2mstrg)
eye_tissues <- c('Retina_Adult.Tissue', 'Retina_Fetal.Tissue', 'RPE_Adult.Tissue', 'RPE_Fetal.Tissue', 'Cornea_Adult.Tissue', 'Cornea_Fetal.Tissue')
det_df <- conv_tab[,-(1:4)] %>% apply(2, function(x) x!='') %>% 
    as.data.frame %>% 
    mutate(num_det_eye =rowSums(select(., all_of(eye_tissues)) ), 
           num_det_body=rowSums(select(., -all_of(eye_tissues)) ))  %>% 
    select(num_det_eye,num_det_body, all_of(eye_tissues), everything()) %>% 
    bind_cols(conv_tab[,1:4] ,.)

low_priority_vars <- c('intron_variant', "intron_variant,non_coding_transcript_variant",
                       "splice_region_variant,intron_variant,non_coding_transcript_variant",  "downstream_gene_variant",
                       "upstream_gene_variant", "intergenic_variant", 'splice_region_variant,intron_variant')
t2g <- complete_gtf %>% filter(type == 'transcript') %>% select(transcript_id, gene_name)
t2c <- conv_tab %>% select(transcript_id, class_code)


adult_retina_panel_var_preds <-read_tsv( glue('{files$VEP_dir}/Retina_Adult.Tissue/variant_summary.txt'),  comment = '##') %>% 
    filter( !Consequence%in% low_priority_vars) %>% 
    rename(transcript_id = Feature) %>% 
    inner_join(t2g, .) %>% inner_join(t2c) %>% filter(class_code != '=') %>% 
    mutate(Adult = T) %>% inner_join(det_df) 





fetal_retina_panel_var_preds <-read_tsv(glue('{files$VEP_dir}/Retina_Fetal.Tissue/variant_summary.txt'),  comment = '##') %>% 
    filter( !Consequence%in% low_priority_vars) %>% 
    rename(transcript_id = Feature) %>% 
    inner_join(t2g, .) %>% inner_join(t2c) %>% filter(class_code != '=') %>% 
    mutate(Fetal = T) 


gencode_panel <- read_tsv(glue('{files$VEP_dir}/gencode/variant_summary.txt'),  comment = '##') %>% 
    rename(transcript_id = Feature, gencode_consequence = Consequence)

all_vars <- full_join(adult_retina_panel_var_preds %>% select(gene_name,`#Uploaded_variation`, Location, Consequence,Adult ) %>% distinct,
                      fetal_retina_panel_var_preds %>% select(gene_name,`#Uploaded_variation`, Location, Consequence,Fetal ) %>% distinct) %>% 
    left_join(gencode_panel %>% select( Location, gencode_consequence) %>% distinct)

keep_locs <- c('chr16:1526594', 'chr1:94002463', 'chr4:15988237', 'chr14:21321429')

clean_variant_results <-  filter(all_vars, Location %in% keep_locs | gene_name == 'ABCA4') %>% 
    group_by(Location, `#Uploaded_variation`, gene_name, Consequence) %>% 
    summarise(`Gencode Predicted Consequence` = str_split(gencode_consequence, ',') %>% 
                                                unlist %>% unique %>% paste(collapse = ',\n')) %>%
    ungroup %>% 
    rename(var_id = `#Uploaded_variation`) %>%
    inner_join(panel_variants) %>% 
    select( `Gene Name` = gene_name, `Location(hg19)` = raw_loc, variant, NM_accession,
            `DNTX Predicted Consequence` = Consequence,`Gencode Predicted Consequence` ) %>% 
    distinct 


var_to_study <-list(c('Chr1:94468019','Bauwens et al.', "ABCA4-associated maculopathy", "Chr1:94468019 G>T"),
                    c('Chr1:94481967', 'Bauwens et al.', "ABCA4-associated maculopathy",'Chr1:94481967 C>T'),
                    c('Chr1:94546814','Bauwens et al.',"ABCA4-associated maculopathy",'Chr1:94546814 G>C' ),
                    c('Chr1:94484001', 'Braun et al.\nZernant et al.', "Stargardt disease","Chr1:94484001 C>T" ),
                    #c('Chr1:94484001', 'Braun et al.\nZernant et al.', "Stargardt disease","Chr1:94484001 C>A" ),
                    c('Chr1:94484082', 'Braun et al.\nZernant et al.',"Stargardt disease",'Chr1:94484082 T>G' ),
                    c('Chr1:94526934', 'Zernant et al.', "Stargardt disease",'Chr1:94526934 T>G' ),
                    c('Chr1:94527698', 'Sangermano et al.',"Stargardt disease", 'Chr1:94527698 G>C'),
                    c('Chr1:94546780', 'Sangermano et al.',"Stargardt disease",'Chr1:94546780 C>G'),
                    c('Chr14:21789588', 'Jamshidi et al.',"RPGRIP1-mediated inherited retinal degeneration","Chr14:21789588 G>A"),
                    c('Chr4:15989860', 'Mayer et al.', "Cone–rod dystrophy",  "Chr4:15989860 T>G"),
                    c('Chr16:1576595','Geoffroy et al.',"Ciliopathy", "Chr16:1576595 C>A" )
) %>% do.call(rbind, .) %>% as_tibble %>%  rename( `Location(hg19)`=V1, `Published Study` =V2, `Associated Disease` = V3, 
                                                   loc_with_var = V4)
clean_variant_results <- clean_variant_results %>% 
    mutate(`DNTX Predicted Consequence` = str_replace_all(`DNTX Predicted Consequence`, '_', ' '),
           `Gencode Predicted Consequence` = str_replace_all(`Gencode Predicted Consequence`, '_', ' ')) %>% 
    inner_join(var_to_study)

save(clean_variant_results, file = files$intron_variant_analysis_rdata)

