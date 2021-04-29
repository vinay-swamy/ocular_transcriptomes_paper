## this script is supposed to run as part of the Snakefile, so no working directory will be specified
library(tidyverse)
library(data.table)
library(tximport)
library(edgeR)
library(limma)
library(glue)
library(patchwork)
library(ggplotify)
library(ComplexHeatmap)
library(yaml)
library(scattermore)


source('scripts/read_salmon.R')

args <- commandArgs(trailingOnly = T)
# files_yaml <- '/data/swamyvs/ocular_transcriptomes_paper/files.yaml'
files_yaml <- args[1]
files <- read_yaml(files_yaml)

replace_nan <- function(df) {
    df_fixed <- lapply(colnames(df),function(col) pull(df, col) %>%  
                           {replace(., is.nan(.), 0)}) %>% bind_cols %>% as_tibble
    colnames(df_fixed) <- colnames(df)
    return(df_fixed)
    
}

sample_table <- read_tsv(files$sample_table)
conv_tab <- fread(files$tcons2mstrg, sep ='\t') %>% select(all_of(colnames(.)[1:4]), Retina_Fetal.Tissue) %>% 
    filter(Retina_Fetal.Tissue!= '')
gtf <- rtracklayer::readGFF(files$anno_gtf) %>% as_tibble %>% filter(transcript_id %in% conv_tab$transcript_id)
excluded_quant_files <- list.files(files$excluded_retina_quant, pattern = 'quant.sf', recursive = T,
                                   full.names = T)



txi_raw <- read_salmon(which_return = 'txi', normalize_counts = F, qfiles = excluded_quant_files)
txi_normed <- read_salmon(which_return = 'txi', qfiles = excluded_quant_files)

load(files$core_tight_rdata)

retina_dev_samples <- filter(sample_table, origin %in% c('Retina_FetalTissue', 'Retina_Organoid')) 
fetal_retina_core_tight <- core_tight %>% 
    ungroup %>%  
    filter(Tissue == 'Retina', grepl('fetal|organoid',Sub_Tissue,ignore.case = T)) %>% 
    mutate(subtissue = 'Retina_Fetal.Tissue', 
           origin = ifelse(grepl('Fetal',Sub_Tissue, ignore.case = T), 'Retina_FetalTissue', 'Retina_Organoid'), 
           in_dntx = sample_accession%in% retina_dev_samples$sample) %>% 
    select(sample = sample_accession, tissue = Tissue, subtissue, origin, age_str=Age_Days, study=study_accession, in_dntx) %>% 
    distinct %>% 
    mutate(age_num = as.numeric(age_str)) %>% 
    arrange(age_num) %>% 
    mutate(age_fac = factor(age_str, levels = unique(age_str) ))

dev_retina_core_tight <- fetal_retina_core_tight %>% filter(origin =='Retina_FetalTissue') %>% 
    filter(study == 'SRP105756') #%>% 
#filter(study == 'SRP119766') %>% 

tx_counts_raw <- txi_raw$counts %>% as.data.frame %>% 
    mutate(transcript_id = rownames(.)) %>% 
    select(transcript_id, all_of(dev_retina_core_tight$sample)) %>% 
    filter(rowSums(.[,-1] >10) >0)

tx_abundance_normed <- txi_normed$abundance %>% 
    as.data.frame %>% 
    mutate(transcript_id = rownames(.)) %>% 
    select(transcript_id, all_of(dev_retina_core_tight$sample)) %>% 
    filter(transcript_id %in% tx_counts_raw$transcript_id)


edist <- function(a,b){
    return( sqrt((a[1] -b[1])^2 + (a[2] - b[2])^2 ) )
}
pca <- tx_abundance_normed %>% select(-transcript_id) %>% as.matrix %>% t %>% prcomp
pca_df <- pca$x %>% 
    as.data.frame %>% 
    mutate(sample = rownames(.)) %>% 
    select(sample, PC1, PC2) %>% 
    inner_join(dev_retina_core_tight)

center <- c(mean(pca_df$PC1), mean(pca_df$PC2) )
pca_df$c_dist <- apply(pca_df[,2:3], 1, function(x)edist(center, x) )
pca_df <- pca_df %>% arrange(desc(c_dist)) %>%  mutate(outlier = c(rep('OL', 5), rep('NOL', nrow(pca_df) -5) ))

ggplot(pca_df)+
    geom_point(aes(x= PC1, y=PC2, color = outlier), size = 7)+
    geom_text(aes(x = PC1, y=PC2,label =sample))+
    theme_minimal()
ol_samps <- filter(pca_df, outlier == 'OL') %>% pull(sample)
dev_retina_core_tight <- dev_retina_core_tight %>% filter(!sample %in% ol_samps) %>% 
    arrange(age_num) %>% 
    mutate(pseudo_age = factor(age_str, levels = unique(age_str)))
tx_counts_raw <- tx_counts_raw  %>% select(-all_of(ol_samps))
tx_abundance_normed <- tx_abundance_normed %>% select(-all_of(ol_samps))

t2g <- gtf %>% filter(type == "transcript") %>% select(transcript_id, gene_name) %>% distinct
gene_abundance_normed <- tximport::summarizeToGene(txi_normed, tx2gene = t2g) %>% 
    .[['abundance']] %>% 
    as.data.frame %>% 
    mutate(gene_name = rownames(.)) %>% 
    select(gene_name,all_of(dev_retina_core_tight$sample)) %>% 
    inner_join(t2g, .) %>%
    inner_join(tx_abundance_normed%>% select(transcript_id),.)
stopifnot(all(gene_abundance_normed$transcript_id == tx_abundance_normed$transcript_id))# critical 
samples <- dev_retina_core_tight$sample
tx_piu <- {tx_abundance_normed[,samples] / gene_abundance_normed[,samples] } %>%
    replace_nan %>% 
    {bind_cols(tx_abundance_normed[,'transcript_id'] ,.) }
colnames(tx_piu)[1] <- 'transcript_id'

counts_mat <- tx_counts_raw %>% select(-transcript_id, all_of(dev_retina_core_tight$sample))

rownames(counts_mat) <- tx_counts_raw$transcript_id
design_mat <- model.matrix(~dev_retina_core_tight$pseudo_age + 0)
colnames(design_mat) <- levels(dev_retina_core_tight$pseudo_age)
rownames(design_mat) <- dev_retina_core_tight$sample
dge <- calcNormFactors(DGEList(counts_mat))
dge <- voom(dge, design_mat)
fit <- lmFit(dge, design = design_mat)
efit <- eBayes(fit)



tpm_co <-10
piu_diff_co <- .25
min_qval <- .01
tps <- levels(dev_retina_core_tight$pseudo_age)
tx_ab_pseudot <- tx_abundance_normed %>%   
    gather(key = 'sample', value = 'ab', -transcript_id) %>% 
    inner_join(dev_retina_core_tight%>% select(sample, pseudo_age)) %>% 
    group_by(transcript_id, pseudo_age) %>% 
    summarise(avg_ab = mean(ab)) %>% 
    ungroup %>% 
    spread(key =  pseudo_age, value = avg_ab) %>% 
    filter(rowSums(.[,tps] >=tpm_co) >0)




multi_tx_genes <-t2g  %>% filter(transcript_id %in% tx_counts_raw$transcript_id) %>% 
    distinct %>% 
    select(transcript_id, gene_name) %>% group_by(gene_name) %>% summarise(n=n()) %>% pull(gene_name)


tx_piu_pseudot <- tx_piu %>% 
    gather(key = 'sample', value = 'ab', -transcript_id) %>% 
    inner_join(dev_retina_core_tight%>% select(sample, pseudo_age)) %>% 
    group_by(transcript_id, pseudo_age) %>% 
    summarise(avg_ab = mean(ab)) %>% 
    ungroup %>% 
    spread(key =  pseudo_age, value = avg_ab) %>% 
    filter(transcript_id %in% tx_ab_pseudot$transcript_id)

tx_labdf <- t2g %>% filter(gene_name %in% multi_tx_genes, transcript_id %in% tx_ab_pseudot$transcript_id)

novel_single_exons <- gtf %>% filter(is.singleExon == 'TRUE', class_code!= '=') %>% pull(transcript_id) %>% unique
#save.image('testing/retfet0616.rdata')
analyze_DTU<- function(coef){
    #message(glue('processing {coef}'))
    
    top_tx_table <- topTable(efit, coef =coef,number = 100000000) %>%
        as.data.frame( ) %>% 
        mutate(transcript_id = rownames(.)) %>% 
        filter(transcript_id %in%tx_ab_pseudot$transcript_id)
    #message( glue('{nrow(top_tx_table)} initially signficant transcripts') )
    
    other_coefs <- colnames(tx_piu_pseudot)[!colnames(tx_piu_pseudot)%in%c(coef, 'transcript_id')]
    
    above_tpm_co <- tx_ab_pseudot %>% dplyr::filter( !!as.symbol(coef) >=tpm_co) %>% pull(transcript_id)
    #message(glue('{length(above_tpm_co)} transcripts above {tpm_co} TPM'))
    piu_coef <- tx_piu_pseudot %>% select(transcript_id, c_coef_piu :=!!coef)
    
    other_piu <- tx_piu_pseudot %>%  dplyr::select(transcript_id, all_of(other_coefs)) %>% gather(key='sample', value = 'piu', -transcript_id)
    
    piu_diff_tx <- inner_join(other_piu, piu_coef) %>% mutate(piu_diff = c_coef_piu - piu)
    above_piu_diff_tx <- piu_diff_tx %>% 
        dplyr::filter(abs(piu_diff) >= piu_diff_co) %>% 
        pull(transcript_id) %>% unique
    # message(glue('{length(above_piu_diff_tx)} transcripts with over {piu_diff_co} piu diff from mean'))
    
    
    best_results <- top_tx_table %>% 
        dplyr::filter(adj.P.Val < .01,
                      transcript_id%in% tx_labdf$transcript_id,#only genes that have multiple tx
                      transcript_id %in% above_piu_diff_tx, 
                      !transcript_id %in% novel_single_exons) # large change in piu 
    #message(glue('{nrow(best_results)} transcripts that have significant DTU'))
    if(nrow( best_results ) >0){
        res <- best_results %>% 
            #inner_join(t2g, .) #%>% 
            inner_join(tx_labdf)
        
        all_results <- inner_join(top_tx_table, piu_diff_tx) %>% select(transcript_id, adj.P.Val, piu_diff)
        return(list(filtered =res, all=all_results)  )
        
    }   
}

all_dtu_list <- lapply(levels(dev_retina_core_tight$pseudo_age), analyze_DTU) 
#all_dtu_genes <- lapply(all_dtu, function(x) x$filtered$gene_name) %>% reduce(union)
all_dtu <- lapply(all_dtu_list, function(x) x$filtered %>% select(gene_name, transcript_id)) %>% bind_rows() %>% distinct()
dtu_genes_multi_tx <- table(all_dtu$gene_name) %>% .[.>1] %>% names
all_dtu <- all_dtu %>% filter(gene_name %in% dtu_genes_multi_tx)
### for suppdata
names(all_dtu_list) <- levels(dev_retina_core_tight$pseudo_age)
all_dtu_tab <- lapply(1:length(all_dtu_list), function(i) {all_dtu_list[[i]][['filtered']] } %>% mutate(tp = names(all_dtu_list)[i]) ) %>% 
    bind_rows() %>% distinct()
all_dtu_tab <- all_dtu_tab %>% filter(gene_name %in% dtu_genes_multi_tx)
save(all_dtu_tab, file = 'clean_data/rdata/all_dtu.Rdata')


all_dtu_genes <- all_dtu$gene_name %>% unique
all_dtu_tx <- all_dtu$transcript_id


names(all_dtu_list) <- levels(dev_retina_core_tight$pseudo_age)

tx_above_tpm_co <- tx_ab_pseudot %>% filter(rowSums(.[,-1] >=tpm_co)>0 ) %>% pull(transcript_id)
vp_data <- lapply(all_dtu_list, function(x) x[['all']] ) %>% bind_rows() %>% 
    mutate(abdiff = abs(piu_diff))
# 
# vp_data_sig <-vp_data %>%  filter(adj.P.Val < min_qval, abdiff >= piu_diff_co, transcript_id %in% tx_above_tpm_co) %>% 
#     filter(!duplicated(transcript_id), transcript_id %in% all_dtu_tx) %>% 
#     mutate(sig = 'DTU')
# vp_data_not_sig <- vp_data %>% filter(!transcript_id%in% vp_data_sig$transcript_id) %>% 
#     filter(!duplicated(transcript_id)) %>% 
#     mutate(sig = 'NOT DTU')
# 
# vp_data <-bind_rows(vp_data_sig, vp_data_not_sig) %>% mutate(qval= -1* log2(adj.P.Val), 
#                                                              sig = factor(sig, levels = c('NOT DTU', 'DTU')) ) %>% 
#     arrange(sig)
vp_data <- vp_data %>% mutate(sig = ifelse(adj.P.Val < min_qval & abdiff >= piu_diff_co & transcript_id %in% tx_above_tpm_co,
                                           'DTU',
                                           'NOT DTU'),
                              qval= -1* log2(adj.P.Val),
                              sig = factor(sig, levels = c('NOT DTU', 'DTU'))) %>% 
    arrange(sig)

vp <- ggplot(vp_data) +
    geom_scattermore(aes(x=piu_diff, y=qval, color = sig), pointsize = 0)+
    scale_color_manual(values = c( 'DTU'= 'red', 'NOT DTU'='lightgrey'),)+
    guides(color = guide_legend(title = '', override.aes = list(pointsize = 3)))+
    scale_alpha_manual(values =  c( 'DTU'= 1, 'NOT DTU'=1))+
    ylab('-log2(adjusted p value)') + 
    xlab('change in fraction isoform usage')+
    #ggtitle('Differentially expressed transcripts \nwith differential transcript usage(DTU)')+
    cowplot::theme_cowplot() +
    theme(legend.position = 'top')

#table(vp_data$sig)


write(all_dtu_genes, '/tmp/sdfsdfg3', sep = '\n')
system2('/usr/local/bin/Rscript', args = 'scripts/enrichGO.R /tmp/sdfsdfg3  /tmp/gse3') # in repo/scripts/
gse <- readRDS('/tmp/gse3')

sz <- 15

eye_interseting <- gse@result %>% 
    filter(p.adjust < .05, grepl('photo|eye|visual|retin', Description)) %>% 
    mutate(tag = 'eye')
neuro_interseting <- gse@result %>% filter(p.adjust <.05, grepl('neur', Description)) %>% 
    mutate(tag = 'neuro') %>% .[1:5,]
k <- sz-nrow(eye_interseting)-nrow(neuro_interseting) 
other_interseting <- gse@result %>% filter(!ID %in% eye_interseting$ID, !ID%in%neuro_interseting$ID ) %>% .[1:k ,] %>% 
    mutate(tag = 'other')
gse_interesting_only <- gse
gse_interesting_only@result <- bind_rows(other_interseting, eye_interseting, neuro_interseting)

#View(gse@result %>% filter(p.adjust < .05 ))
#View(gse@result %>% filter(p.adjust < .05, grepl('photo|eye|visual|neur|retina', Description) ))

add_newline <- function(str){
    str_list <- str_split(str, ' ') %>% .[[1]]
    split_val <- {length(str_list)/2 } %>% trunc
    paste0(c( paste(str_list[1:split_val], collapse = ' ') , 
              '\n', 
              paste(str_list[(split_val+1) : length(str_list)], collapse = ' ') 
    ), collapse = '')
}
## pulled from clusterProfiler package
## Note 03/2021: running on a fresh R installation, output of cluster profiler script is now a
## an S4 object,  changed fortify command below. May cause a problem if running in biowulf
mod_dotplot_internal <- function(object, x = "geneRatio", color = "p.adjust",
                                 showCategory=10, size=NULL, split = NULL,
                                 font.size=12, title = "", orderBy="x", decreasing=TRUE) {
    
    colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
    if (x == "geneRatio" || x == "GeneRatio") {
        x <- "GeneRatio"
        if (is.null(size))
            size <- "Count"
    } else if (x == "count" || x == "Count") {
        x <- "Count"
        if (is.null(size))
            size <- "GeneRatio"
    } else if (is(x, "formula")) {
        x <- as.character(x)[2]
        if (is.null(size))
            size <- "Count"
    } else {
        ## message("invalid x, setting to 'GeneRatio' by default")
        ## x <- "GeneRatio"
        ## size <- "Count"
        if (is.null(size))
            size  <- "Count"
    }
    
    df <- fortify(object@result, showCategory = showCategory, split=split)
    #View(df)
    ## already parsed in fortify
    ## df$GeneRatio <- parse_ratio(df$GeneRatio)
    
    if (orderBy !=  'x' && !orderBy %in% colnames(df)) {
        message('wrong orderBy parameter; set to default `orderBy = "x"`')
        orderBy <- "x"
    }
    
    if (orderBy == "x") {
        df <- dplyr::mutate(df, x = object@result$GeneRatio %>% str_split('/') %>% do.call(rbind,.) %>% apply(2, as.numeric) %>% {.[,1]/.[,2]} )
    }
    
    idx <- order(df[[orderBy]], decreasing = decreasing)
    df$Description <- sapply(df$Description , add_newline)
    df$Description <- factor(df$Description, levels=rev(unique(df$Description[idx])))
    ggplot(df, aes_string(x="x", y="Description", size=size, color=colorBy)) +
        geom_point() +
        scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
        ## scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
        xlab('GeneRatio') +
        ylab(NULL) + ggtitle(title) + theme_minimal(font.size) + scale_size(range=c(3, 8)) +
        theme(axis.text.y = element_text(size=8))
    
}

dp <- mod_dotplot_internal(gse_interesting_only, showCategory = sz)
pr_gene_table <- gse@result %>% 
    filter(p.adjust <.05, grepl('photo|eye|visual|retin', Description)) %>% 
    mutate(gene_name = str_split(geneID, '/')) %>% 
    unnest(gene_name)
pr_genes <- gse@result %>% filter(p.adjust <.05, grepl('photo|eye|visual|retin', Description)) %>% pull(geneID) %>% str_split('/') %>% unlist %>% unique
lab_map_text<- "visual perception,visual perception
photoreceptor cell differentiation,photoreceptor development
eye morphogenesis,eye development
eye development,eye development
visual system development,eye development
photoreceptor cell development,photoreceptor development
retina development in camera-type eye,eye development
eye photoreceptor cell differentiation,photoreceptor development
eye photoreceptor cell development,photoreceptor development"
lab_map<- read.csv(text = lab_map_text, col.names = c('Description', 'lab'), header = F, stringsAsFactors = F) %>% 
    mutate(lab = str_replace(lab,' ', '\n'))

neuro_genes <- gse@result %>% filter(p.adjust <.05, grepl('neur', Description)) %>% pull(geneID) %>% str_split('/') %>%
    unlist%>% unique


hm_genes <- c(pr_genes, neuro_genes)
dtu_df<- filter(tx_piu, transcript_id %in% all_dtu_tx) %>% inner_join(t2g, .) %>% filter(gene_name %in% pr_genes) %>%
    arrange(gene_name) #%>% select(-transcript_id, -gene_name)
dtu_mat <- dtu_df[,samples] %>% as.matrix
rownames(dtu_mat) <- dtu_df$transcript_id
colnames(dtu_mat) <- dev_retina_core_tight$age_str
cf <- circlize::colorRamp2(c(0,.25,.5,.75, 1), viridisLite::viridis(5))


novel_transcripts <- filter(gtf, !class_code %in%c('=', 'u')) %>% pull(transcript_id) %>% unique
pr_gene_table_long <- pr_gene_table %>% 
    left_join(lab_map) %>% 
    select(gene_name,lab ) %>% 
    distinct %>% 
    mutate(det = T) %>% 
    spread(key = lab, value = lab) %>% 
    inner_join(dtu_df %>% select(gene_name, transcript_id),.) %>% 
    mutate(transcipt_origin = ifelse(transcript_id %in% novel_transcripts, 'Novel\nIsoform', 'Annotated\nIsoform'))
all(dtu_df$transcript_id == pr_gene_table_long$transcript_id)

col = list(`Photoreceptor\nDevelopment` = c('photoreceptor\ndevelopment' = 'yellowgreen', 'white'),
           `Eye\nDevelopment` = c('eye\ndevelopment' = 'sienna', 'white') , 
           `Visual\nPerception` =c('visual\nperception' = 'orange', 'white'),
           `Isoform\nOrigin` = c('Novel\nIsoform' = 'orchid', 'Annotated\nIsoform' = 'springgreen' ))
ha <- HeatmapAnnotation(`Photoreceptor\nDevelopment` = pr_gene_table_long$`photoreceptor\ndevelopment`,
                        `Eye\nDevelopment` = pr_gene_table_long$`eye\ndevelopment` , 
                        `Visual\nPerception` = pr_gene_table_long$`visual\nperception`,
                        `Isoform\nOrigin` = pr_gene_table_long$transcipt_origin,
                        col=col,na_col = 'white',
                        simple_anno_size = unit(.20, "cm"),
                        which = 'row',
                        annotation_legend_param = list(title = '', grid_height = unit(8, "mm"),labels_gp = gpar(fontsize = 10)), 
                        show_annotation_name = F)

ha_col <- HeatmapAnnotation(`Photoreceptor\nDevelopment` = pr_gene_table_long$`photoreceptor\ndevelopment`,
                        `Eye\nDevelopment` = pr_gene_table_long$`eye\ndevelopment` , 
                        `Visual\nPerception` = pr_gene_table_long$`visual\nperception`,
                        `Isoform\nOrigin` = pr_gene_table_long$transcipt_origin,
                        col=col,na_col = 'white',
                        simple_anno_size = unit(.20, "cm"),
                        which = 'column',
                        annotation_legend_param = list(title = '', grid_height = unit(8, "mm"),labels_gp = gpar(fontsize = 10)), 
                        show_annotation_name = F)

#save.image('testing/frdd.Rdata')


heatmap_fe <- Heatmap(dtu_mat, col=cf, show_row_dend = F, show_row_names = F,cluster_columns = T , cluster_rows = T,
              name = 'FIU', right_annotation = ha, column_names_gp = gpar(fontsize=10) )

heatmap_rot <- Heatmap(t(dtu_mat), col=cf, show_column_dend =  F, show_column_names =  F,cluster_columns = T , cluster_rows = T,
                      name = 'FIU', top_annotation = ha_col, row_names_gp =  gpar(fontsize=10),  )

gghm <- as.ggplot(heatmap_fe)
gghm_transpose <- as.ggplot(heatmap_rot)

custom_filter_df <- function(df){
    df <- df %>% select(-gene_name) %>% 
        gather(key = 'sample', value = 'piu', -transcript_id) %>% group_by(transcript_id) %>% 
        summarise(max_piu = max(piu)) %>% 
        inner_join(refid2dntx) %>% 
        arrange(desc(max_piu)) 
    
    res <- df %>% head(3) %>% pull(transcript_id) %>% {tibble(transcript_id = .)}
    if(all(grepl('DNTX', df$pretty_txid[1:3]))){
        res[3,'transcript_id'] <- df[grep('ENST', df$pretty_txid)[1],'transcript_id']
    }
    return(res)
}

refid2dntx <-conv_tab %>% 
    select(transcript_id, class_code,  pretty_txid =refid) %>% 
    mutate(pretty_txid = replace(pretty_txid, pretty_txid == '',transcript_id[pretty_txid == '']), 
           pretty_txid = replace(pretty_txid, class_code !='=',transcript_id[ class_code !='='] )) %>% 
    select(-class_code)




plot_piu_bar <- function(gene,white_list, black_list=''){
    
    #tx_piu_filter <- tx_piu %>% inner_join(t2g,.) %>%  filter(gene_name == gene) %>% 
    #group_by(gene_name) %>% do(custom_filter_df (.))
    piu <- tx_piu %>% 
        inner_join(t2g,.) %>%  
        filter(gene_name == gene, 
               transcript_id %in% white_list,
        ) %>% 
        select(-gene_name) %>% 
        gather(key = 'sample', value = 'piu', -transcript_id) %>% 
        inner_join(dev_retina_core_tight %>% select(sample, age_str, age_num)) %>% 
        group_by(transcript_id, age_str, age_num) %>% 
        summarise(avg_piu = mean(piu)) %>% 
        ungroup %>% 
        inner_join(tx_labdf) %>% 
        left_join(refid2dntx)  %>% 
        filter(!transcript_id %in% black_list) %>% 
        arrange(age_num) %>% 
        mutate(age_fac = factor(age_str, levels = unique(age_str)))
    #%>% 
    #left_join(dev_retina_core_tight %>% select(-sample) %>% distinct)
    
    
    #tx_labdf %>% left_join(dntx2ref) %>% mutate(new=replace(transcript_id))
    tx_ab <- tx_abundance_normed %>% 
        inner_join(t2g,.) %>%  filter(gene_name == gene) %>% select(-gene_name) %>% 
        filter(transcript_id %in% piu$transcript_id) %>%
        gather(key = 'sample', value = 'ab', -transcript_id) %>% 
        inner_join(dev_retina_core_tight %>% select(sample, age_str, age_num)) %>% 
        group_by(transcript_id, age_str, age_num) %>% summarise(avg_ab = mean(ab)) %>% 
        arrange(age_num) %>% 
        left_join(refid2dntx) 
    ab <- gene_abundance_normed %>% filter(gene_name == gene) %>% select(-transcript_id) %>% 
        distinct %>% 
        gather(key = 'sample', value = 'ab', -gene_name) %>% 
        inner_join(dev_retina_core_tight %>% select(sample, age_str, age_num)) %>% 
        group_by(gene_name, age_str, age_num) %>% 
        summarise(avg_ab = mean(ab)) %>% 
        ungroup %>%
        rename(pretty_txid = gene_name) %>% 
        bind_rows(tx_ab) %>% 
        filter(!transcript_id %in% black_list) %>% 
        mutate(ab=log2(avg_ab+1)) %>% 
        arrange(age_num) %>% 
        mutate(age_fac = factor(age_str, levels = unique(age_str)))
    color_list <- c('red','blue', 'green', 'purple')
    names(color_list) <- unique(ab$pretty_txid) %>% {c(.[. == gene], .[grepl('ENST',.)], .[grepl('DNTX',.)] )}
    
    piu_plot <- ggplot(piu %>% rename(`Transcript ID` = pretty_txid) )+
        geom_col(aes(x=age_fac, y=avg_piu, fill = `Transcript ID`), position = 'dodge')+
        scale_fill_manual(values = color_list)+
        xlab('Days Post Fertilization')+
        ylab('fraction of total\ngene expression')+
        cowplot::theme_cowplot() +
        theme(legend.position = c(1,1))
    ab_plot <- ggplot(ab %>% rename(`Transcript ID` = pretty_txid) )+
        geom_line(aes(x = age_fac, y=ab, group = `Transcript ID`, color =`Transcript ID`)) +
        scale_color_manual(values = color_list)+
        xlab('Days Post Fertilization')+
        ylab('log2(TPM+1)') +
        cowplot::theme_cowplot() +
        theme(legend.title = element_blank())
    p <- piu_plot / ab_plot +plot_annotation(tag_level = 'A') #+plot_layout(guides = 'collect')
    p
    return(p)
    
}

plot_piu_bar_listout <- function(gene,white_list, black_list=''){
    
    #tx_piu_filter <- tx_piu %>% inner_join(t2g,.) %>%  filter(gene_name == gene) %>% 
    #group_by(gene_name) %>% do(custom_filter_df (.))
    piu <- tx_piu %>% 
        inner_join(t2g,.) %>%  
        filter(gene_name == gene, 
               transcript_id %in% white_list,
        ) %>% 
        select(-gene_name) %>% 
        gather(key = 'sample', value = 'piu', -transcript_id) %>% 
        inner_join(dev_retina_core_tight %>% select(sample, age_str, age_num)) %>% 
        group_by(transcript_id, age_str, age_num) %>% 
        summarise(avg_piu = mean(piu)) %>% 
        ungroup %>% 
        inner_join(tx_labdf) %>% 
        left_join(refid2dntx)  %>% 
        filter(!transcript_id %in% black_list) %>% 
        arrange(age_num) %>% 
        mutate(age_fac = factor(age_str, levels = unique(age_str)))
    #%>% 
    #left_join(dev_retina_core_tight %>% select(-sample) %>% distinct)
    
    
    #tx_labdf %>% left_join(dntx2ref) %>% mutate(new=replace(transcript_id))
    tx_ab <- tx_abundance_normed %>% 
        inner_join(t2g,.) %>%  filter(gene_name == gene) %>% select(-gene_name) %>% 
        filter(transcript_id %in% piu$transcript_id) %>%
        gather(key = 'sample', value = 'ab', -transcript_id) %>% 
        inner_join(dev_retina_core_tight %>% select(sample, age_str, age_num)) %>% 
        group_by(transcript_id, age_str, age_num) %>% summarise(avg_ab = mean(ab)) %>% 
        arrange(age_num) %>% 
        left_join(refid2dntx) 
    ab <- gene_abundance_normed %>% filter(gene_name == gene) %>% select(-transcript_id) %>% 
        distinct %>% 
        gather(key = 'sample', value = 'ab', -gene_name) %>% 
        inner_join(dev_retina_core_tight %>% select(sample, age_str, age_num)) %>% 
        group_by(gene_name, age_str, age_num) %>% 
        summarise(avg_ab = mean(ab)) %>% 
        ungroup %>%
        rename(pretty_txid = gene_name) %>% 
        bind_rows(tx_ab) %>% 
        filter(!transcript_id %in% black_list) %>% 
        mutate(ab=log2(avg_ab+1)) %>% 
        arrange(age_num) %>% 
        mutate(age_fac = factor(age_str, levels = unique(age_str)))
    color_list <- c('red','blue', 'green', 'purple')
    names(color_list) <- unique(ab$pretty_txid) %>% {c(.[. == gene], .[grepl('ENST',.)], .[grepl('DNTX',.)] )}
    
    piu_plot <- ggplot(piu %>% rename(`Transcript ID` = pretty_txid) )+
        geom_col(aes(x=age_fac, y=avg_piu, fill = `Transcript ID`), position = 'dodge')+
        scale_fill_manual(values = color_list)+
        xlab('Days Post Fertilization')+
        ylab('fraction of total\ngene expression')+
        cowplot::theme_cowplot() +
        theme(legend.position = c(1,1))
    ab_plot <- ggplot(ab %>% rename(`Transcript ID` = pretty_txid) )+
        geom_line(aes(x = age_fac, y=ab, group = `Transcript ID`, color =`Transcript ID`)) +
        scale_color_manual(values = color_list)+
        xlab('Days Post Fertilization')+
        ylab('log2(TPM+1)') +
        cowplot::theme_cowplot() +
        theme(legend.title = element_blank())
    return(list(ab=ab_plot, piu = piu_plot))
    
}

draw_all_transcripts_static <- function(gene, gtf,keep_tx, black_list = ''){
    main_cols <- c("seqid", "type",  "start", "end", "strand", "transcript_id", "gene_name", "exon_number" )
    gtf_gene <- filter(gtf, gene_name == gene, transcript_id %in% keep_tx)
    unique_tx <- length(unique(gtf_gene$transcript_id))
    gtf_exons <- filter(gtf_gene, type == 'exon') %>% 
        select(seqid, strand, start, end) %>% 
        distinct %>% 
        arrange(start) %>% 
        mutate(length=end-start, length=sqrt(length), Xmin=0, Xmax=0, Ymin=-1, Ymax=1) %>% 
        mutate(lab1=paste0('sdfs',1:nrow(.)), 
               lab2=paste0('wef',1:nrow(.)), 
               lab3=paste0('gfr',1:nrow(.)) )
    gap=mean(gtf_exons$length)
    gtf_exons[1,'Xmax'] <- gtf_exons[1,'Xmin'] + gtf_exons[1,'length']
    if(nrow(gtf_exons) > 1){
        for(i in 2:nrow(gtf_exons)){
            gtf_exons[i,'Xmin'] <- gtf_exons[(i-1),'Xmax'] + gap 
            gtf_exons[i,'Xmax'] <- gtf_exons[i,'Xmin'] + gtf_exons[i,'length']
        }
    }
    
    plot_data <- filter(gtf_gene, type == 'exon') %>%  
        select(main_cols,novel_exon_id) %>% 
        inner_join(gtf_exons) %>% 
        mutate(`exon type`=ifelse(is.na(novel_exon_id), 'ref', 'novel')) %>% 
        inner_join(refid2dntx)%>% 
        filter(!transcript_id %in% black_list)
    # mutate(pretty_txid = replace(pretty_txid, transcript_id%in%all_dtu_tx, 
    #                            paste0(pretty_txid[transcript_id%in%all_dtu_tx], '*'  )) )
    color_list <- c( 'blue', 'green', 'purple')
    names(color_list) <- unique(plot_data$pretty_txid) %>% {c(gene, .[grepl('ENST',.)], .[grepl('DNTX',.)] )}
    #print(plot_data)
    #print(color_list)
    plot <- ggplot(data = plot_data %>% rename(`Transcript ID` = pretty_txid) ) +
        geom_rect( aes(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax,fill=`Transcript ID`))+
        scale_fill_manual(values = color_list) +
        facet_wrap(~`Transcript ID`, ncol=1)+
        #ggtitle(gene) +
        theme_void() +
        theme( strip.background = element_blank(),strip.text.x = element_blank())
    #print(nchar(plot$data$transcript_id))
    return(plot)
    #return(girafe(ggobj = plot, width_svg = wsvg, height_svg = hsvg))
    
}

df <- conv_tab %>% inner_join(t2g) %>% inner_join(all_dtu) %>% inner_join(refid2dntx) %>% 
    filter(class_code != '=', 
           #gene_name %in% pr_genes
           )
#t_gene <- unique(df$gene_name)[8]#MYO9A 
#1 DNTX_00140848 LRP2      ENST00000263816.7
#2 DNTX_00142728 LRP2      DNTX_00142728    
#3 DNTX_00140849 LRP2      DNTX_00140849   
t_gene <-  "LRP2"
#gm_fetret <- draw_all_transcripts_static(t_gene , gtf, all_dtu_tx, black_list = 'DNTX_00080651')
piu_fetret <- plot_piu_bar(t_gene,white_list =c('DNTX_00142728','DNTX_00140848' ),  black_list = 'DNTX_00080651')

###new code 
library(DBI)
db_name <- '/data/swamyvs/ocular_transcriptomes_pipeline/data/shiny_data/app_data/DNTX_db.sql'
con <- dbConnect(RSQLite::SQLite(),db_name)
keep_tx <- c('DNTX_00142728','DNTX_00140848')
pretty_keep_tx <- filter(refid2dntx, transcript_id %in% keep_tx) %>% pull(pretty_txid) %>% unique
s_gtf <- con %>% tbl('plotting_gtf') %>% filter(gene_name == t_gene) %>% collect %>% filter(transcript_id %in% pretty_keep_tx)

c_merge <- function(mat){
    #mat must bed sorted by start and interval length 
    res <- list()
    cur <- mat[1, ]
    j <- 1
    for( i in 2:nrow(mat)){
        #1=start, 2=end#
        if(mat[i, 1]<= cur[1,  2]){
            #merge the interval 
            if(mat[i, 2] > cur[1, 2]){
                cur[1, 2] <- mat[i, 2]
                cur[1, 3] <- paste(cur[1, 3],mat[i, 3], sep='-')
            }else{
                cur[1, 3] <- paste(cur[1, 3],mat[i, 3], sep='-')
            }
            
        }else{
            # interval is good 
            if(any(cur %in% res)) {cur <- mat[i, ] }
            res[[j]] <- cur 
            cur <- mat[i, ]
            j <- j+1
        }
        
    }
    if(!any(cur %in% res)) {
        #cur[1, 3] <- paste(cur[1, 3],mat[i, 3], sep='-')
        res[[j+1]] <- cur
    }
    return(do.call(rbind, res))
    
} 



simple_scale_introns <- function(df, sf){
    df <- df %>% mutate(xlength = Xmax-Xmin)
    exon_length <- mean(df$xlength)
    intron_length <- exon_length*sf
    new_xmin = rep(-1, nrow(df))
    new_xmax  = rep(-1, nrow(df))
    
    new_xmin[1] <- 0
    new_xmax[1] <- exon_length
    for (i in 2:nrow(df)){
        if (df$start[i] == df$end[(i-1)] + 1) {
            new_xmin[i] = new_xmax[(i-1)]
            new_xmax[i] = new_xmin[i] + exon_length
        }else{
            new_xmin[i] = new_xmax[(i-1)] + intron_length
            new_xmax[i] = new_xmin[i] + exon_length
        }
    }
    df <- df %>% mutate(Xmin = new_xmin, Xmax = new_xmax)
    return(df)
}


draw_all_transcripts_csm <- function(gene, gtf, keep_tx, g, black_list, sf){
    #manually scale intron and exons, because it doesnt work well precopmuting it for genes with a lot of exons
    # original code dervied from webapp
    gtf_gene <- filter(gtf, gene_name == gene, transcript_id %in% keep_tx, !transcript_id %in% black_list)
    distinct_exons <- gtf_gene %>% filter(type == 'exon') %>% 
        select(seqid, strand, start, end, Xmin, Xmax) %>% 
        distinct %>% 
        mutate(length = sqrt(end-start),) %>% 
        arrange(start, length) %>% 
        mutate(id= as.character(1:nrow(.))) %>% 
        select(-length)
    
    all_correct <- simple_scale_introns(distinct_exons, sf=sf)
    
    
    all_plot_data <- distinct_exons %>% 
        select(-Xmax, -Xmin) %>%  
        inner_join(all_correct) %>% 
        inner_join(gtf_gene %>% select(-Xmin, -Xmax, -length)) %>% 
        mutate(`exon type`=ifelse(is.na(novel_exon_id), 'ref', 'novel')) 
    plot_data <- all_plot_data %>% filter(type == 'exon') %>% mutate(exon_number = as.numeric(exon_number))

    color_list <- c('black', 'red')
    names(color_list) <- c('ref', 'novel')
    color_list <- c('red', 'blue', 'darkred', 'purple')
    names(color_list) <- unique(plot_data$transcript_id) %>% {c(gene, .[grepl('ENST',.)], .[grepl('DNTX',.)] )}
    plot_data$duptxid <- plot_data$transcript_id
    m <- replace_na(plot_data$novel_exon_id == 'nvl_exon120703', F)
    plot_data[m, 'duptxid'] = 'nvl_exon120703'
    arrow_data <- plot_data %>% filter(duptxid == 'nvl_exon120703') %>% select(Xmin,Xmax, duptxid,transcript_id) %>% 
        mutate(x=mean(c(Xmin, Xmax)), y=1.5)
    plot <- ggplot(data = plot_data) +
        geom_rect( aes(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax,fill=duptxid))+
        scale_fill_manual(values = color_list) +
        geom_segment(data = arrow_data, aes(x=x, xend=x),y=2, yend=1,size=.5, arrow=arrow(length = unit(0.1, "npc")))+
        geom_text(data = arrow_data,aes(x=x+1), y=1.8, label='Novel Exon', size =4, hjust = 'left')+
        ylim(c(-1, 2))+
        facet_wrap(~transcript_id, ncol=1)+
        ggtitle(gene)+
        theme_void() +
        theme( strip.background = element_blank(),strip.text.x = element_blank(), legend.position = 'none', 
               plot.title = element_text(hjust = 0.5))
    return(plot)
    
}
gm_fetret <- draw_all_transcripts_csm(t_gene,s_gtf, pretty_keep_tx, 1, black_list = 'DNTX_00080651', sf=1)
####
bdes <- '
AB
AC
AC
AC
AC
AC
AC'
bottom <- gghm  +gm_fetret + piu_fetret +plot_layout(design = bdes)

top <- vp |dp
full_des <-'
AAA
BBB
BBB'
fp <- top/bottom + 
    plot_layout(design = full_des) +
    plot_annotation(tag_levels = 'A') &  
    theme(legend.margin = margin(0,0,0,0), legend.box.margin=margin(-10,-10,-10,-10))
ggsave('/data/swamyvs/ocular_transcriptomes_paper/testing/fig5.png', fp, height = 15, width = 13)


br_des <-'
AA
AA
AA
BB
CD
CD
CD'
piu_fetret_list <- plot_piu_bar_listout(t_gene,white_list =c('DNTX_00142728','DNTX_00140848' ),  black_list = 'DNTX_00080651')
ab_line <-  piu_fetret_list$ab
piu_bar <- piu_fetret_list$piu


bottom_rework <- gghm_transpose + gm_fetret + piu_bar + ab_line +plot_layout(design = br_des, guides = 'collect')
fp_rw <- top/bottom_rework + 
    plot_layout(design = full_des) +
    plot_annotation(tag_levels = 'A') &  
    theme(legend.margin = margin(0,0,0,0), legend.box.margin=margin(-10,-10,-10,-10))

ggsave('/data/swamyvs/ocular_transcriptomes_paper/testing/fig5_rq.png', fp_rw, height = 15, width = 13)

save(vp, dp, ab_line, piu_bar, gm_fetret, gghm_transpose, top, bottom_rework,fp_rw, full_des, br_des, t_gene, dev_retina_core_tight, all_dtu_genes, all_dtu_tx,
    file =files$fetal_retina_diffexp_results )



