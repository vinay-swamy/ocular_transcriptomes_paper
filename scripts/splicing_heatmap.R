library(tidyverse)
library(ComplexHeatmap)
library(viridis)
library(matrixStats)
 # setwd('~/NIH/dev_eyeintegration_splicing/')
 # psi_file <- 'data/rmats/all_tissues_psi.tsv'
 # sample_file <- 'sampleTableDev.tsv'
args <- commandArgs(trailingOnly = T)
wd <- args[1]
psi_file <- args[2]
sample_file <- args[3]
outfile <- args[4]
#save(args, file='/tmp/sphm_args.rdata')
setwd(wd)
psi_tab <- read_tsv(psi_file)
colnames(psi_tab) <- str_split(colnames(psi_tab), '_psi') %>% sapply(function(x) x[[1]])
sample_table <- read_tsv(sample_file, col_names = c('sample', 'run', 'paired', 'tissue', 'subtissue', 'origin')) %>%
  filter( subtissue!='synth', paired=='y')
set.seed(34543)
psi_tab[is.na(psi_tab)] <- 0
#mat <- psi_tab[,-(1:4)] %>% mutate(vars=rowVars(as.matrix(.)) ) %>% arrange(vars) %>% tail(10000) %>% .[,sample_table$sample]
mat <- psi_tab[,-(1:4)] %>% filter(rowSums(.)< (1*ncol(.)), rowSums(.)>(.05*ncol(.))) %>% sample_n(10000) %>% .[,sample_table$sample] %>% t()

#mat[is.na(mat)] <- 0
ht_opt$fast_hclust = TRUE
ha <- HeatmapAnnotation(tissue=sample_table$subtissue, which = 'row')
hm <- Heatmap(mat,col = viridis(100)  ,name = 'PSI', right_annotation = ha, show_row_dend = F, show_column_dend =F,
              heatmap_height = unit(20,'cm'), heatmap_width = unit(40,'cm'), show_column_names = F, show_row_names = F)
splicing_heatmap <- draw(hm)
splicing_heatmap
save(splicing_heatmap,file = outfile)
