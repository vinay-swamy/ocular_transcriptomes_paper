library(tidyverse)
library(ComplexHeatmap)
library(viridis)
library(matrixStats)
 setwd('~/NIH/')
 psi_file <- 'eyeintegration_splicing/dl_data/all_tissues_psi.tsv.gz'
 sample_file <- 'eyeintegration_splicing/sampleTableV6.tsv'
args <- commandArgs(trailingOnly = T)
wd <- args[1]
psi_file <- args[2]
sample_file <- args[3]
outfile <- args[4]
#save(args, file='/tmp/sphm_args.rdata')
setwd(wd)
psi_tab <- read_tsv(psi_file)
colnames(psi_tab) <- str_remove(colnames(psi_tab), '_psi')
sample_table <- read_tsv(sample_file) %>% filter(sample %in% colnames(psi_))
set.seed(34543)
psi_tab[is.na(psi_tab)] <- 0
#mat <- psi_tab[,-(1:4)] %>% mutate(vars=rowVars(as.matrix(.)) ) %>% arrange(vars) %>% tail(10000) %>% .[,sample_table$sample]
mat <- psi_tab[,-(1:4)] %>% filter(rowSums(.)< (1*ncol(.)), rowSums(.)>(.05*ncol(.))) %>% sample_n(10000) %>% .[,sample_table$sample] %>% t()


ht_opt(fast_hclust = TRUE)
ha <- HeatmapAnnotation(tissue=sample_table$body_location, which = 'row')
hm <- Heatmap(mat,col = viridis(100)  ,name = 'PSI', right_annotation = ha, show_row_dend = F, show_column_dend =F,
              heatmap_height = unit(20,'cm'), heatmap_width = unit(40,'cm'), show_column_names = F, show_row_names = F)
splicing_heatmap <- draw(hm)
splicing_heatmap
save(splicing_heatmap,file = outfile)
