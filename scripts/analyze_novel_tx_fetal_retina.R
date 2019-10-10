library(tidyverse)
library(limma)
library(edgeR)
library(ComplexHeatmap)
library(viridis)
args <- c('~/NIH/occular_transcriptomes_paper/',
          '~/NIH/eyeintegration_splicing/sampleTableV6.tsv', 
          '~/NIH/occular_transcriptomes_paper/clean_data/EiaD_quant.Rdata',
          '~/NIH/eyeintegration_splicing/dl_data/gfc_TCONS_to_st_MSTRG.tsv.gz',
          '~/NIH/eyeintegration_splicing/all_tissues.combined.gtf',
          '~/NIH/eyeintegration_splicing/dl_data/all_tissues.combined_transdecoderCDS.gff3.gz',
          '~/NIH/occular_transcriptomes_paper/clean_data/fetal_novel_tx_diffexp_results.Rdata',
          '~/NIH/occular_transcriptomes_paper/clean_data/fetal_novel_tx_diffexp.txt',
          '~/NIH/occular_transcriptomes_paper/clean_data/fetal_novel_tx_diffexp_hm.Rdata'
          )

args <- commandArgs(trailingOnly = T)

wd <- args[1]
sample_table_file <- args[2]
eiad_quant_data <- args[3]
t2m_file <- args[4]
gtf_file <- args[5]
gff3_file <- args[6]
diffexp_data <- args[7]
diff_exp_tx <- args[8]
heatmap_file <- args[9]

setwd(wd)

sample_table <- read_tsv(sample_table_file)
load(eiad_quant_data)
tcons2mstrg <- read_tsv(t2m_file)
gtf <- rtracklayer::readGFF(gtf_file)
gff3 <- rtracklayer::readGFF(gff3_file) %>% 
    as_tibble  %>% mutate(ID=str_extract(ID,'TCONS_[0-9]+|ENSG[0-9]+'))
retina_fetal_exp <- tcons2mstrg %>% select(transcript_id, Retina_Fetal.Tissue) %>% 
    filter(!is.na(Retina_Fetal.Tissue)) %>%
    pull(transcript_id)


samples <- gene_quant$sample_accession
fetal_gene_quant <- gene_quant[,-1] %>% t() %>% as.data.frame %>% mutate(gene_name=rownames(.)) %>% 
    select(gene_name, everything())
colnames(fetal_gene_quant) <- c('gene_name', samples)

samples <- tx_quant$sample_accession
fetal_tx_quant <- tx_quant[,-1] %>% t() %>% as.data.frame %>%
    mutate(transcript_id=str_extract(rownames(.), '\\((.*?)\\)') %>% str_remove('\\(') %>% str_remove('\\)'), 
           gene_name=str_remove(rownames(.), '\\((.*?)\\)')) %>% select(transcript_id,gene_name, everything())
colnames(fetal_tx_quant) <- c('transcript_id','gene_name', samples)


core_tight_fetal <- filter(core_tight_2019, sample_accession %in% colnames(fetal_gene_quant)) %>% 
    select(colnames(.)[1:6], study_accession) %>%
    mutate(Sub_Tissue= ifelse(Sub_Tissue =="Retina - Fetal Tissue", 'Retina_Fetal.Tissue', 'Retina_Organoid'),
           Age_Days=as.numeric(Age_Days))
colnames(core_tight_fetal) <- c('sample', 'tissue', 'subtissue', 'origin', 'age', 'kept', 'study_accession')
fetal_tissue_samples <- core_tight_fetal %>% filter(subtissue == 'Retina_Fetal.Tissue', kept=='Kept') %>% arrange(age)
organoid_samples <- core_tight_fetal %>% filter(subtissue == 'Retina_Organoid', kept=='Kept') %>% arrange(age)


library(limma)
library(edgeR)
#gdata::keep(fetal_tissue_samples, fetal_gene_quant)
sample_design <- fetal_tissue_samples %>% select(sample, age, study_accession) %>%
    mutate(stage= case_when(age <=84 ~ 'early',
                            age >84 & age<119 ~ 'mid', 
                            age >=119 ~ 'late'))
table(sample_design$stage)
sample_design %>% group_by(stage) %>% summarise(range=max(age) - min(age))

gene_names <- fetal_gene_quant$gene_name
exp_mat <- fetal_gene_quant %>% .[,sample_design$sample] 
rownames(exp_mat) <- gene_names
keep <- rowSums(exp_mat) > 5*ncol(exp_mat)
exp_mat <- exp_mat[keep,]
nrow(exp_mat)
table(sample_design$study_accession)
dge <- calcNormFactors(DGEList(exp_mat))
stage <- factor(sample_design$stage)
study <- factor(sample_design$study_accession)
design_mat <- model.matrix(~0 + study + stage)
colnames(design_mat) <-colnames(design_mat) %>% str_remove('study|stage')
voom_dge <- voom(dge, design = design_mat)
design.pairs <-function(levels) {
    n <- length(levels)
    design <- matrix(0,n,choose(n,2))
    rownames(design) <- levels
    colnames(design) <- 1:choose(n,2)
    k <- 0
    for (i in 1:(n-1))
        for (j in (i+1):n) {
            k <- k+1
            design[i,k] <- 1
            design[j,k] <- -1
            colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
        }
    design
}

cont_mat <- design.pairs(c(levels(stage),levels(study))) %>% .[,grepl('early|mid|late', colnames(.))]
model_fitted <-lmFit(voom_dge, design = design_mat)
cont_mat <- cont_mat[colnames(model_fitted),]
model_results <- contrasts.fit(model_fitted, cont_mat) %>% eBayes
colnames(model_results)
targets <- c("early-late","early-mid", 'late-mid')
upregulated <- list(early=c(), mid=c(), late=c())
limma_de <- list()

for( i in c("early-late","early-mid", 'late-mid')){
    conds <- str_split(i, '-')[[1]]
    top <- conds[1]
    bottom <- conds[2]
    df <-  topTable(model_results, adjust.method = 'fdr', coef = i, number = 30000000, p.value = .05) %>% 
        as.data.frame %>% mutate(gene_name=rownames(.)) %>% 
        filter(grepl('TCONS', gene_name), gene_name %in% retina_fetal_exp)
    upregulated[[top]] <- c(upregulated[[top]],
                            df %>% filter(logFC >0) %>% arrange(desc(logFC)) %>%  pull(gene_name))
    upregulated[[bottom]] <- c(upregulated[[bottom]], 
                               df %>% filter(logFC <0) %>% arrange(desc(abs(logFC))) %>% pull(gene_name) )
    limma_de[[i]] <- df
}

shared_upregulated <- reduce(upregulated, intersect)  
distinct_upregulated <- combine(upregulated) %>% unique 
des <- model.matrix(~0 + stage)
batch_cor_exp <- removeBatchEffect(exp_mat, study, design = des)
batch_cor_exp[batch_cor_exp<0] <- 0



mat <- log2(batch_cor_exp[distinct_upregulated,] + 1 )
sum(mat>10)/length(mat)
#mat[mat>10] <- 10
topAno <- HeatmapAnnotation(age=sample_design$age, 
                            stage=sample_design$stage, 
                            study=sample_design$study_accession,
                            col=list(stage=c('early'='yellow', 'mid'='orange', 'late'='red'),
                                     study=c("SRP119766"='green', "SRP105756"='purple', "SRP090040"='blue')),
                            which = 'col')
rightAno <- HeatmapAnnotation(protein_coding=ifelse(distinct_upregulated %in% gff3$ID, 'pc', 'nc'),
                              col=list(protein_coding=c('pc'='pink', 'nc'='white')),
                              which = 'row')

hm <- Heatmap(mat, cluster_rows = F, cluster_columns = T, 
        show_row_names = F,
        name = 'log2(TPM+1)',
        top_annotation = topAno, 
        right_annotation = rightAno,
        heatmap_width = unit(7, 'in'),
        heatmap_height = unit(12, 'in'),
        col = viridis(100), column_labels = sample_design$age)
fetal_retina_novel_tx_heatmap <- draw(hm)

save(limma_de, upregulated, file = diffexp_data)
write(distinct_upregulated, file = diff_exp_tx, sep = '\n')
save(fetal_retina_novel_tx_heatmap, file = heatmap_file)


#fetal_gene_quant[,c('gene_name', sample_design$sample)] %>% View












