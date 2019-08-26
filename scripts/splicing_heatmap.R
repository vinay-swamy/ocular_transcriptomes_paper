library(tidyverse)
library(ComplexHeatmap)
library(viridis)
setwd('~/NIH/dev_eyeintegration_splicing/')
psi_file <- 'data/rmats/all_tissues_psi.tsv'
sample_file <- 'sampleTableDev.tsv'
psi_tab <- read_tsv(psi_file)
colnames(psi_tab) <- str_split(colnamFes(psi_tab), '_psi') %>% sapply(function(x) x[[1]])
sample_table <- read_tsv(sample_file, col_names = c('sample', 'run', 'paired', 'tissue', 'subtissue', 'origin')) %>% 
  filter( subtissue!='synth')
set.seed(34543)
mat <- psi_tab[,-(1:4)] %>% sample_n(10000) %>% .[,sample_table$sample]
mat[is.na(mat)] <- 0
ht_opt$fast_hclust = TRUE
ha <- HeatmapAnnotation(tissue=sample_table$subtissue, which = 'col')
hm <- Heatmap(mat,col = viridis(100)  ,name = 'PSI', top_annotation = ha, show_row_dend = F, show_column_dend =F,
              heatmap_height = unit(30,'cm'), heatmap_width = unit(20,'cm'))
splicing_heatmap <- draw(hm)
splicing_heatmap
save(splicing_heatmap,file = '~/NIH/occular_transcriptomes_paper/data/splicing_heatmap.Rdata')











