library(tidyverse)
library(limma)
library(edgeR)
library(ComplexHeatmap)
library(viridis)
args <- c('/Volumes/data/occular_transcriptomes_paper/',
          '/Volumes/data/eyeintegration_splicing/sampleTableFull.tsv', 
          '/Volumes/data/eyeintegration_splicing/data/all_tissue_quant.Rdata',
          '/Volumes/data/eyeintegration_splicing/data/misc/TCONS2MSTRG.tsv',
          '/Volumes/data/eyeintegration_splicing/data/gtfs/all_tissues.combined.gtf',
          '/Volumes/data/eyeintegration_splicing/data/seqs/transdecoder_results/all_tissues.combined_transdecoderCDS.gff3',
          '/Volumes/data/eyeintegration_splicing/ref/core_tight.Rdata',
          '~/NIH/occular_transcriptomes_paper/clean_data/fetal_novel_tx_diffexp_results.Rdata',
          '~/NIH/occular_transcriptomes_paper/clean_data/fetal_novel_tx_diffexp.txt',
          '~/NIH/occular_transcriptomes_paper/clean_data/fetal_novel_tx_diffexp_hm.Rdata'
          )

args <- commandArgs(trailingOnly = T)
wd <- args[1]
sample_table_file <- args[2]
salmon_quant_data <- args[3]
t2m_file <- args[4]
gtf_file <- args[5]
gff3_file <- args[6]
eiad_core_tight <- args[7]
diffexp_data <- args[8]
diff_exp_tx <- args[9]
heatmap_file <- args[10]

setwd(wd)

sample_table <- read_tsv(sample_table_file)
load(salmon_quant_data)
all_quant[is.na(all_quant)] <- 0
tcons2mstrg <- read_tsv(t2m_file)
gtf <- rtracklayer::readGFF(gtf_file)
gff3 <- rtracklayer::readGFF(gff3_file) %>% 
    as_tibble  %>% mutate(ID=str_extract(ID,'DNTX_[0-9]+|ENSG[0-9]+'))
retina_fetal_exp <- tcons2mstrg %>% select(transcript_id, Retina_Fetal.Tissue) %>% 
    filter(!is.na(Retina_Fetal.Tissue)) %>%
    pull(transcript_id)
retina_fetal_samples <- sample_table %>% filter(subtissue == 'Retina_Fetal.Tissue') %>% pull(sample)
fetal_tx_quant <- all_quant %>% filter(transcript_id %in% retina_fetal_exp) %>% select(transcript_id,retina_fetal_samples)
load(eiad_core_tight)
core_tight_fetal <- filter(core_tight, sample_accession %in% colnames(fetal_tx_quant), Tissue == 'Retina', Kept=='Kept') %>% 
    select(colnames(.)[1:6], study_accession) %>%
    mutate(Sub_Tissue= ifelse(Sub_Tissue =="Retina - Fetal Tissue", 'Retina_Fetal.Tissue', 'Retina_Organoid'),
           Age_Days=as.numeric(Age_Days))# there's a rogue RPE sample marked as retina organoid
fetal_tx_quant <- fetal_tx_quant %>% select(transcript_id, core_tight_fetal$sample_accession)
colnames(core_tight_fetal) <- c('sample', 'tissue', 'subtissue', 'origin', 'age', 'kept', 'study_accession')
fetal_tissue_samples <- core_tight_fetal %>% filter(subtissue == 'Retina_Fetal.Tissue') %>% arrange(age)
organoid_samples <- core_tight_fetal %>% filter(subtissue == 'Retina_Organoid') %>% arrange(age)


library(limma)
library(edgeR)
#gdata::keep(fetal_tissue_samples, fetal_gene_quant)
sample_design <- fetal_tissue_samples %>% select(sample, age, study_accession) %>%
    mutate(stage= case_when(age <=107 ~ 'early',
                            age >107 ~ 'late'))
table(sample_design$stage)
sample_design %>% group_by(stage) %>% summarise(range=max(age) - min(age))

tx_names <- fetal_tx_quant$transcript_id
exp_mat <- fetal_tx_quant %>% .[,sample_design$sample] 
rownames(exp_mat) <- tx_names
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

cont_mat <- design.pairs(c(levels(stage),levels(study))) %>% .[,grepl('early|late', colnames(.))]
model_fitted <-lmFit(voom_dge, design = design_mat)
cont_mat <- cont_mat[colnames(model_fitted),]
model_results <- contrasts.fit(model_fitted, cont_mat) %>% eBayes
colnames(model_results)
targets_cont <- topTable(model_results, adjust.method = 'fdr', coef = "early-late", number = 30000000, p.value = .05)
sum(rownames(targets_cont) %in% novel_loci_distinct$transcript_id)
upregulated <-  target_cont %>% filter(logFC >=1)
downregulated <- target_cout %>% filter(logfc<=-1)

des <- model.matrix(~0 + stage)
batch_cor_exp <- removeBatchEffect(exp_mat, study, design = des)
batch_cor_exp[batch_cor_exp<0] <- 0

novel_loci <- filter(gtf, grepl('DNTX', gene_name), type == "transcript") 
sum(rownames(targets_cont) %in% novel_loci_distinct$transcript_id)
keep <- rownames(batch_cor_exp) %in% novel_loci_distinct$transcript_id & rownames(batch_cor_exp) %in% rownames(targets_cont)
mat <- log2(batch_cor_exp[keep,]+1)
# exp_by_stage <- lapply(c('early', 'mid', 'late'), function(x) filter(sample_design, stage == x) %>% pull(sample) %>% {batch_cor_exp[,.]} %>%
#            rowMeans) %>% bind_cols %>% mutate()
# 
# 
# sum(mat>500)/length(mat)
# mat[mat>500] <- 500
topAno <- HeatmapAnnotation(age=sample_design$age, 
                            stage=sample_design$stage, 
                            study=sample_design$study_accession,
                            col=list(stage=c('early'='yellow', 'late'='red'),
                                     study=c("SRP119766"='green',  "SRP090040"='blue')),
                            which = 'col')
rightAno <- HeatmapAnnotation(protein_coding=ifelse(rownames(mat)%in% gff3$ID, 'pc', 'nc'),
                              col=list(protein_coding=c('pc'='pink', 'nc'='white')),
                              which = 'row')

hm <- Heatmap(mat, cluster_rows = F, cluster_columns = T, 
        show_row_names = F,
        name = 'log2(TPM+1)',
        top_annotation = topAno, 
        right_annotation = rightAno,
        col = viridis(100), column_labels = sample_design$age)
fetal_retina_novel_tx_heatmap <- draw(hm)

save(limma_de, upregulated, file = diffexp_data)
write(distinct_upregulated, file = diff_exp_tx, sep = '\n')
save(fetal_retina_novel_tx_heatmap, file = heatmap_file)


#fetal_gene_quant[,c('gene_name', sample_design$sample)] %>% View












