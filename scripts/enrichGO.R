library(clusterProfiler)
args=commandArgs(trailingOnly = T)

genes <- scan(args[1], character(), sep ='\n')
gse <- enrichGO(genes,  OrgDb = 'org.Hs.eg.db', ont = 'BP', keyType = 'SYMBOL', qvalueCutoff = .05)
saveRDS(gse, file = args[2])