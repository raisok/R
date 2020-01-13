
library(pheatmap)

args <- commandArgs(TRUE)
input1<-args[1]
output<-args[2]

gene = read.table(input1,stringsAsFactors=F,header = T,row.names = 1)
gene1 = log10(gene)
gene1[gene1 == "-Inf"] =0
gene1[gene1 >= 20] =20
pdf(output)
pheatmap(gene1,cluster_cols = F,show_rownames = F)
dev.off()
